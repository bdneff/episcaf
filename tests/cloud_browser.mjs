// Drive both private cloud paths against the real local API with in-memory Google services.
import assert from 'node:assert/strict';
import {spawn} from 'node:child_process';
import {randomBytes} from 'node:crypto';
import {once} from 'node:events';
import {setTimeout as delay} from 'node:timers/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE || 'playwright-core');
const manual=process.env.EPISCAF_PILOT_MODE==='manual-gemini';
const token='alice.'+randomBytes(30).toString('hex');
const setup=`
import csv, io, itertools, os, sys
from types import SimpleNamespace
from zipfile import ZipFile
sys.path.insert(0,"tests")
from test_cloud_backbones import Bucket, CloudPilot
from episcaf_service import cloud, pilot, disco
from episcaf_service.app import app
import uvicorn
os.environ.pop("K_SERVICE",None)
requests,results=Bucket(),Bucket()
cloud.store=lambda:requests
cloud.outputs=lambda:results
fake=CloudPilot()
fake.jobs={}
# Advance the test clock past the status cache between requests.
clock=itertools.count(1800000000,31)
cloud.time=SimpleNamespace(time=lambda:next(clock))
def batch(method, suffix="", **kwargs):
    response=fake.batch(method,suffix,**kwargs)
    if method=="POST":
        rid=kwargs["params"]["jobId"][3:]
        run=cloud.read(cloud.record_blob("alice",rid))
        table=io.StringIO()
        csv_rows=csv.writer(table)
        csv_rows.writerow(["source","design_id","epitope","route","structure_file"])
        archive=io.BytesIO()
        with ZipFile(archive,"w") as z:
            for index,region in enumerate(run["regions"]):
                for i in range(8):
                    name=f"fixture_{index}_{i}"
                    path=f"backbones/{name}.cif.gz"
                    csv_rows.writerow(["cloud-rfd3",name,region,run["route"],path])
                    z.writestr(path,b"browser test fixture")
            z.writestr("backbones.csv",table.getvalue())
        results.blob(rid+"/backbones.csv").upload_from_string(table.getvalue())
        results.blob(rid+"/designs.zip").upload_from_string(archive.getvalue())
        cloud.write(results.blob(rid+"/done.json"),dict(workflow="rfd3-backbones-v1",generated=8*len(run["regions"]),completed=next(clock)))
    return response
cloud.batch_request=batch
disco._from_biolib=disco._from_local=lambda path,chain,kind: {"%d"%i:i/129 for i in range(1,130)}
cloud.write(requests.blob("testers/alice"),{"digest":pilot.digest(os.environ.pop("BROWSER_TEST_TOKEN"))})
uvicorn.run(app,host="127.0.0.1",port=0)
`;
const server=spawn(process.env.PYTHON || 'python3',['-c',setup],{env:{...process.env,
  EPISCAF_PILOT_MODE:manual?'manual-gemini':'cloud-lysozyme', EPISCAF_DISCOTOPE:manual?'biolib':'local',EPISCAF_GPU_IMAGE:'fixture',
  EPISCAF_RESULT_BUCKET:'results',EPISCAF_WORKER_SA:'worker',BROWSER_TEST_TOKEN:token}});
let logs='',browser;
server.stderr.on('data',s=>logs+=s);
try{
  for(let i=0;i<150 && !/Uvicorn running on/.test(logs);i++){
    if(server.exitCode!==null) throw Error(logs);
    await delay(100);
  }
  const origin=logs.match(/Uvicorn running on (http:\/\/127\.0\.0\.1:\d+)/)[1];
  browser=await chromium.launch({channel:'chrome',headless:true});
  const page=await browser.newPage();
  for(const scenario of (manual?['antibody','growpatch']:['antibody','growpatch','review'])){
    const route=scenario==='review'?'growpatch':scenario;
    await page.goto(origin+'/home.html');
    if(route==='antibody') await page.getByRole('link',{name:'Find antibody contacts'}).click();
    else await page.getByRole('link',{name:'Find candidate regions',exact:true}).click();
    await page.locator('#live-note').waitFor();
    const ref=await page.locator('#live-note a').getAttribute('href');
    let structure=await (await page.request.get(origin+ref)).body();
    if(scenario==='review')structure=Buffer.from(structure.toString().replaceAll(' Y ',' A '));
    await page.locator('#file').setInputFiles({name:'test.pdb',mimeType:'chemical/x-pdb',buffer:structure});
    await page.getByRole('radio',{name:scenario==='review'?/Chain A/:/Chain Y/}).waitFor();
    if(route==='growpatch'){
      if(manual){
        assert.match(await page.locator('#provider').innerText(),/selected protein chain to BioLib/);
        assert.equal(await page.locator('#disco').isChecked(),false);
        await page.locator('#disco').check();
      }else{
        assert.match(await page.locator('#provider').innerText(),/episcaf service/);
        assert.ok(await page.locator('#disco').isChecked());
        assert.ok(await page.locator('#disco-type').isVisible());
      }
    }
    await page.locator('#run').click();
    await page.locator('.pick').first().waitFor();
    if(route==='growpatch') await page.locator('.disco:not(.none)').first().waitFor();
    if(route==='antibody') assert.ok(await page.locator('.pick').first().isChecked());
    else {
      const picks=page.locator('.pick');
      await picks.first().check();
      if(!manual && scenario!=='review'){
        await picks.nth(1).check();
        await picks.nth(2).check();
        assert.ok(await picks.nth(3).isDisabled());
        await picks.nth(1).uncheck();
        assert.ok(await picks.nth(3).isEnabled());
        await picks.nth(3).check();
      }
    }
    if(await page.locator('#signin').isVisible()){
      await page.locator('#token').fill(token);
      await page.locator('#signin button').click();
      await page.locator('#signout').waitFor();
    }
    await page.locator('#stage').click();
    const launch=page.getByRole('button',{name:manual?'Request a Gemini run':scenario==='review'?'Request lab review':'Start cloud run'});
    await launch.waitFor();
    const total=manual?16:scenario==='growpatch'?24:8;
    assert.match(await page.locator('#reviewtext').innerText(),new RegExp("Total: "+total+" backbone attempts"));
    if(manual)assert.match(await page.locator('#reviewtext').innerText(),/up to 8 backbones/);
    await launch.click();
    if(manual){
      await page.locator('#runlist').getByText(/awaiting Brandon’s submission/).first().waitFor();
      const link=page.getByRole('link',{name:'Download submitted request (ZIP)'});
      const [download]=await Promise.all([page.waitForEvent('download'),link.click()]);
      assert.match(download.suggestedFilename(),/^request-[a-f0-9]{32}\.zip$/);
      assert.equal(await download.failure(),null);
      assert.match(await page.locator('#request-regions').innerText(),/Y:/);
    }else if(scenario==='review'){
      await page.locator('#runlist tr[aria-current="true"]').getByText('awaiting lab review').waitFor();
      assert.ok(await page.locator('#models').isHidden());
    }else{
      await page.getByRole('link',{name:'Download backbone list (CSV)'}).waitFor();
      assert.equal(await page.locator('#tbl tbody tr').count(),total);
      assert.equal(await page.locator('#stats').innerText(),total+' RFdiffusion3 backbones');
      assert.ok(await page.locator('#screenhelp').isHidden());
      assert.ok(await page.locator('#passonly').isHidden());
      const [download]=await Promise.all([page.waitForEvent('download'),
        page.getByRole('link',{name:'Download backbone coordinates (ZIP)'}).click()]);
      assert.equal(await download.failure(),null);
      if(scenario==='growpatch'){
        await page.locator('#file').setInputFiles({name:'historical.csv',mimeType:'text/csv',
          buffer:Buffer.from('design_id,source,s_p1_free\nexample,cloud-backbone,2.0\n')});
        await page.getByText('1 retained backbone',{exact:true}).waitFor();
        await page.locator('#file').setInputFiles({name:'scores.csv',mimeType:'text/csv',
          buffer:Buffer.from('design_id,overall_rmsd\nexample,1.5\n')});
        await page.locator('#screenhelp').waitFor();
        assert.ok(await page.locator('#passonly').isVisible());
      }
    }
  }
  assert.equal((await (await page.request.get(origin+'/pilot/runs')).json()).length,manual?2:3);
  assert.ok(!logs.includes(token));
  console.log(manual?'Both manual Gemini reference routes passed.':'Cloud demo passed antibody contacts, ranked region selection, 8- and 24-backbone result downloads, and review of other uploads.');
}finally{
  await browser?.close();
  if(server.exitCode===null){const done=once(server,'exit');server.kill();await done;}
}
