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
import os, sys
sys.path.insert(0,"tests")
from test_cloud_backbones import Bucket, CloudPilot
from episcaf_service import cloud, pilot
from episcaf_service.app import app
import uvicorn
os.environ.pop("K_SERVICE",None)
requests,results=Bucket(),Bucket()
cloud.store=lambda:requests
cloud.outputs=lambda:results
fake=CloudPilot()
fake.jobs={}
cloud.batch_request=fake.batch
cloud.write(requests.blob("testers/alice"),{"digest":pilot.digest(os.environ.pop("BROWSER_TEST_TOKEN"))})
uvicorn.run(app,host="127.0.0.1",port=0)
`;
const server=spawn(process.env.PYTHON || 'python3',['-c',setup],{env:{...process.env,
  EPISCAF_PILOT_MODE:manual?'manual-gemini':'cloud-lysozyme', EPISCAF_DISCOTOPE:'off',EPISCAF_GPU_IMAGE:'fixture',
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
  for(const route of ['antibody','growpatch']){
    await page.goto(origin+'/home.html');
    if(route==='antibody') await page.getByRole('link',{name:'Find antibody contacts'}).click();
    else await page.getByRole('link',{name:'Find candidate regions',exact:true}).click();
    await page.locator('#live-note').waitFor();
    const ref=await page.locator('#live-note a').getAttribute('href');
    const structure=await (await page.request.get(origin+ref)).body();
    await page.locator('#file').setInputFiles({name:'test.pdb',mimeType:'chemical/x-pdb',buffer:structure});
    await page.getByRole('radio',{name:/Chain Y/}).waitFor();
    if(route==='growpatch') await page.locator('#disco').uncheck();
    await page.locator('#run').click();
    await page.locator('.pick').first().waitFor();
    if(route==='antibody') assert.ok(await page.locator('.pick').first().isChecked());
    else await page.locator('.pick').first().check();
    if(await page.locator('#signin').isVisible()){
      await page.locator('#token').fill(token);
      await page.locator('#signin button').click();
      await page.locator('#signout').waitFor();
    }
    await page.locator('#stage').click();
    await page.getByRole('button',{name:manual?'Request a Gemini run':'Start cloud run'}).waitFor();
    assert.match(await page.locator('#reviewtext').innerText(),/Total: 16 backbone attempts/);
    assert.match(await page.locator('#reviewtext').innerText(),/up to 8 backbones/);
    await page.getByRole('button',{name:manual?'Request a Gemini run':'Start cloud run'}).click();
    await page.locator('#runlist').getByText(manual?/awaiting Brandon’s submission/:/Google Batch job/).first().waitFor();
    if(manual){
      const link=page.getByRole('link',{name:'Download submitted request (ZIP)'});
      const [download]=await Promise.all([page.waitForEvent('download'),link.click()]);
      assert.match(download.suggestedFilename(),/^request-[a-f0-9]{32}\.zip$/);
      assert.equal(await download.failure(),null);
      assert.match(await page.locator('#request-regions').innerText(),/Y:/);
    }
  }
  assert.equal((await (await page.request.get(origin+'/pilot/runs')).json()).length,2);
  assert.ok(!logs.includes(token));
  console.log('Both epitope routes passed the 16-trial review and submission in '+(manual?'manual Gemini':'Google Batch')+' mode.');
}finally{
  await browser?.close();
  if(server.exitCode===null){const done=once(server,'exit');server.kill();await done;}
}
