// Run from the repo root. Requires playwright-core, Chrome, and the web service Python environment.
import assert from 'node:assert/strict';
import {mkdtemp, readFile, rm} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {spawn, execFileSync} from 'node:child_process';
import {once} from 'node:events';
import {setTimeout as delay} from 'node:timers/promises';
const {chromium} = await import(process.env.PLAYWRIGHT_MODULE || 'playwright-core');
const dir = await mkdtemp(join(tmpdir(), 'episcaf-pilot-'));
const python = process.env.PYTHON || 'python3';
const env = {...process.env, EPISCAF_PILOT_DB:join(dir,'runs.sqlite'), EPISCAF_WS:dir,
  EPISCAF_PILOT_MODE:'simulated', EPISCAF_DISCOTOPE:'', K_SERVICE:'', EPISCAF_SIM_FAIL:'0'};
const admin=(action,owner)=>execFileSync(python,['-m','episcaf_service.pilot',action,owner],
  {env,encoding:'utf8'}).trim();
let service,browser,origin,logs='';
async function start(port=0) {
  service=spawn(python,['-m','uvicorn','episcaf_service.app:app','--host','127.0.0.1','--port',String(port)],
    {env,stdio:['ignore','pipe','pipe']});
  let startup='';
  service.stderr.on('data',chunk=>{startup+=chunk;logs+=chunk;});
  service.stdout.on('data',chunk=>{logs+=chunk;});
  for(let i=0;i<150;i++){
    const match=startup.match(/Uvicorn running on (http:\/\/127\.0\.0\.1:\d+)/);
    if(match){origin=match[1];return;}
    if(service.exitCode!==null)throw Error(startup);
    await delay(100);
  }
  throw Error('Service startup timed out.');
}
async function stop() {
  if(service && service.exitCode===null){
    const stopped=once(service,'exit');service.kill();await stopped;
  }
}
async function launchBrowser(){
  browser=await chromium.launch({channel:'chrome',headless:true});
}
async function login(page,token){
  await page.locator('#token').fill(token);
  await page.getByRole('button',{name:'Sign in',exact:true}).click();
  await page.locator('#signout').waitFor();
}
try{
  const alice=admin('issue','alice'), bob=admin('issue','bob');
  await start();
  await launchBrowser();
  let page=await browser.newPage();
  await page.goto(origin+'/growpatch_tool.html?backend=https://unconfigured.invalid');
  assert.equal(await page.locator('#file').inputValue(),'');
  // Strip partners from the checked-in fixture to exercise the lone-antigen path.
  const fixture=await readFile('docs/pdb/8PWW.pdb','utf8');
  assert.equal(await page.locator('#structure').isVisible(),false);
  let enumerations=0;
  page.on('request',req=>{if(new URL(req.url()).pathname==='/enumerate')enumerations++;});
  let finishInspect;
  const inspectionGate=new Promise(resolve=>{finishInspect=resolve;});
  await page.route('**/inspect?*',async route=>{await inspectionGate;await route.continue();},{times:1});
  await page.locator('#file').setInputFiles({name:'complex.pdb',mimeType:'text/plain',buffer:Buffer.from(fixture)});
  try{
    await page.locator('#viewer canvas').first().waitFor({timeout:10000});
    const preview=await page.locator('#viewer').boundingBox();
    const form=await page.locator('#formcard').boundingBox();
    assert.ok(preview.x>=form.x+form.width && preview.y<form.y+form.height);
    assert.ok(preview.width===340 && preview.height===255);
    await page.setViewportSize({width:390,height:844});
    const mobile=await page.locator('#viewer').boundingBox();
    assert.ok(mobile.x>=0 && mobile.x+mobile.width<=390 && mobile.height<=255);
    await page.setViewportSize({width:1280,height:720});
    assert.equal(await page.locator('#results').isVisible(),false);
    assert.ok(await page.evaluate(()=>viewer.getModel().selectedAtoms({chain:'B'}).length>0));
  }finally{finishInspect();}
  await page.waitForFunction(()=>!document.getElementById('run').disabled);
  await page.locator('#chain').selectOption('B');
  assert.ok(await page.evaluate(()=>viewer.getModel().selectedAtoms({chain:'B'}).some(a=>a.style.cartoon)));
  assert.equal(await page.evaluate(()=>viewer.getModel().selectedAtoms({chain:'A'}).some(a=>a.style.cartoon)),false);
  assert.equal(enumerations,0);
  const pdb=fixture.split('\n')
    .filter(l=>l.startsWith('ATOM  ') && l[21]==='A').join('\n')+'\nEND\n';
  await page.locator('#file').setInputFiles({name:'8PWW-A.pdb',mimeType:'text/plain',buffer:Buffer.from(pdb)});
  await page.locator('#disco').uncheck();
  await page.locator('#run').click();
  for(let i=0;i<3;i++)await page.locator('.pick').nth(i).check();
  const expanded=await page.locator('#viewer').boundingBox();
  assert.ok(expanded.width>340 && expanded.height>=480);
  await login(page,alice);
  const reviewResponse=page.waitForResponse(r=>r.url().endsWith('/pilot/reviews'));
  await page.locator('#stage').click();
  const review=await (await reviewResponse).json();
  assert.match(await page.locator('#reviewtext').innerText(),/Total: 24 backbone attempts/);
  assert.equal(review.rounds,1); assert.equal(review.partial_diffusion,false);
  assert.equal(review.regions.length,3);
  for(const region of review.regions){
    const rows=review.contigs.filter(r=>r.growpatch_epitope===region);
    assert.ok(rows.length>1);
    assert.equal(rows.reduce((n,r)=>n+r.backbone_trials,0),8);
  }
  const path='/pilot/runs/'+review.id;
  const result=page.waitForResponse(r=>r.url().endsWith(path));
  const [duplicate]=await Promise.all([page.request.post(origin+path),page.locator('#launch').click()]);
  assert.equal((await duplicate.json()).id,(await (await result).json()).id);
  assert.equal((await (await page.request.get(origin+'/pilot/runs')).json()).length,1);
  await page.waitForURL(origin+'/results.html?run='+review.id+'#result-title');
  await page.locator('#runlist tr[aria-current="true"]').waitFor();
  assert.match(await page.locator('#runlist').innerText(),/8PWW-A.pdb/);
  await browser.close();browser=null;
  // Read the durable state directly, without polling the web API or leaving a browser open.
  for(let i=0;i<20;i++){
    const stage=execFileSync(python,['-c',
      'import sqlite3,sys;print(sqlite3.connect(sys.argv[1]).execute("SELECT stage FROM runs").fetchone()[0])',
      env.EPISCAF_PILOT_DB],{encoding:'utf8'}).trim();
    if(stage==='completed')break;
    assert.ok(i<19,'Monitor did not complete while the browser was closed');
    await delay(500);
  }
  const port=new URL(origin).port;
  await stop();env.EPISCAF_SIM_FAIL='1';await start(port);
  await launchBrowser();page=await browser.newPage();
  await page.goto(origin+'/results.html?run='+review.id);
  await login(page,alice);
  await page.getByRole('link',{name:'My runs',exact:true}).click();
  await page.getByRole('link',{name:'View results',exact:true}).click();
  assert.equal(await page.locator('#tbl tbody tr').first().innerText().then(t=>t.includes('simulated')),true);
  assert.equal(await page.locator('#tbl tbody tr.na').count(),3);
  assert.equal((await page.request.get(origin+path+'/results?download=true')).status(),200);
  await page.locator('#file').setInputFiles({name:'legacy.csv',mimeType:'text/csv',
    buffer:Buffer.from('epitope_chunk_rmsd,overall_rmsd,mean_pae,aomclash_pred\n0.1,1,2,0.1\n')});
  await page.getByText('legacy.csv',{exact:true}).waitFor();
  assert.equal(await page.locator('#tbl tbody tr.na').count(),1);
  await page.getByRole('link',{name:'View results',exact:true}).click();
  await page.waitForURL(origin+'/results.html?run='+review.id+'#result-title');
  await page.locator('#original').waitFor();
  assert.equal(await page.locator('#tbl tbody tr').count(),3);
  const other=await browser.newPage();
  await other.goto(origin+'/results.html#myruns');
  await login(other,bob);
  assert.deepEqual(await (await other.request.get(origin+'/pilot/runs')).json(),[]);
  for(const url of [path,path+'/results',path+'/results?download=true','/results/'+review.id,'/status/'+review.id])
    assert.equal((await other.request.get(origin+url)).status(),404);
  assert.equal((await other.request.post(origin+path)).status(),404);
  const failedReview=await (await other.request.post(origin+'/pilot/reviews',{
    data:{filename:review.filename,chain:'A',regions:review.regions,structure:pdb}})).json();
  await other.request.post(origin+'/pilot/runs/'+failedReview.id);
  await other.getByText(/simulated worker stopped/).waitFor();
  admin('revoke','bob');
  assert.equal((await other.request.get(origin+'/pilot/runs')).status(),401);
  const rotated=admin('issue','alice');
  assert.equal((await page.request.get(origin+path+'/results')).status(),401);
  assert.equal((await page.request.post(origin+'/pilot/session',{
    headers:{Authorization:'Bearer '+rotated,Origin:'https://unconfigured.invalid'}})).status(),403);
  await page.request.post(origin+'/pilot/session',{headers:{Authorization:'Bearer '+rotated}});
  assert.equal((await (await page.request.get(origin+path)).json()).id,review.id);
  // Exercise the deployed mode guard and the non-persistent preparation endpoint.
  await stop();env.K_SERVICE='cloud-run-staging-test';await start(port);
  assert.equal((await (await page.request.get(origin+'/health')).json()).pilot_mode,'staging');
  assert.equal((await page.request.get(origin+'/pilot/runs')).status(),404);
  const selection={filename:review.filename,chain:'A',regions:review.regions};
  const staged=await (await page.request.post(origin+'/prepare',{data:selection})).json();
  assert.deepEqual(staged.contigs,review.contigs);
  assert.equal(staged.total,24);
  const tight=await (await page.request.post(origin+'/prepare',{data:{...selection,regions:['A:1-101']}})).json();
  assert.equal(tight.contigs.length,1);
  assert.equal(tight.contigs[0].backbone_trials,8);
  assert.equal((await page.request.post(origin+'/prepare',{data:{...selection,regions:['A:1-103']}})).status(),422);
  assert.ok(!logs.includes(alice) && !logs.includes(bob) && !logs.includes(rotated));
  console.log('Verified 3 regions × 8 trials, concurrent duplicate launch, background completion, restart, fresh browser sign-in, Results, manual CSV, ownership, failure and revocation.');
}finally{
  await browser?.close();
  await stop();
  await rm(dir,{recursive:true,force:true});
}
