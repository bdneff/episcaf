// Exercise the sequence picker with real pointer events and a local PDB fixture.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
const {chromium}=await import(process.env.PLAYWRIGHT_MODULE || 'playwright-core');
const browser=await chromium.launch({channel:'chrome',headless:true});
try{
  const page=await browser.newPage({viewport:{width:900,height:850}});
  const html=await readFile('docs/scaffold.html','utf8');
  const pdb=await readFile('episcaf_v3/energetics/md/3hfm_apo/structures/3hfm.pdb');
  let complex=false, staging=false;
  let stagedRequest;
  const errors=[];
  page.on('pageerror',error=>errors.push(error.message));
  await page.route('**/*',route=>{
    const path=new URL(route.request().url()).pathname;
    if(path==='/scaffold.html')return route.fulfill({contentType:'text/html',body:html});
    if(path==='/health')return route.fulfill({json:{service:'episcaf',pilot_mode:staging?'staging':'manual-gemini'}});
    if(path==='/submit'){
      stagedRequest=new URL(route.request().url());
      return route.fulfill({json:{run_id:'example',n_contigs:1,n_epitopes:1,run_dir:'example',submit_sh:'example'}});
    }
    if(path==='/inspect')return route.fulfill({json:{antigen_chain:'Y',is_complex:complex,
      partner_chains:complex?['H','L']:[],epitope:'',n_epitope:0}});
    return route.abort();
  });
  async function upload(){
    await page.goto('http://episcaf.test/scaffold.html');
    await page.locator('#pdbfile').setInputFiles({name:'test.pdb',mimeType:'text/plain',buffer:pdb});
    await page.locator('#seqbody .seqrow').first().waitFor();
  }
  const letters=chain=>page.locator('#seqbody span[data-c="'+chain+'"]');
  async function point(chain,index){
    const bounds=await letters(chain).nth(index).boundingBox();
    return [bounds.x+bounds.width/2,bounds.y+bounds.height/2];
  }
  async function drag(chain,start,end,via){
    await letters(chain).nth(start).scrollIntoViewIfNeeded();
    await page.mouse.move(...await point(chain,start));
    await page.mouse.down();
    if(via!==undefined)await page.mouse.move(...await point(chain,via));
    await page.mouse.move(...await point(chain,end));
    await page.mouse.up();
  }
  await upload();
  assert.equal(await page.locator('input[type=number],input[type=range],#backend,#strategy,input[type=checkbox]').count(),0);
  assert.equal(await page.locator('#simulation').isVisible(),false);
  await page.getByRole('button',{name:'Simulate a run'}).click();
  await page.locator('#rbody tr').first().waitFor();
  assert.equal(stagedRequest,undefined);
  staging=true;
  await upload();
  await drag('Y',4,8,11);
  const picked=()=>page.locator('#seqbody .epi').count();
  assert.equal(await picked(),5);
  await drag('Y',7,5);
  assert.equal(await picked(),2);
  await letters('Y').nth(8).click();
  assert.equal(await picked(),1);
  await letters('Y').nth(11).click({modifiers:['Shift']});
  assert.equal(await picked(),5);
  assert.equal(await page.evaluate(()=>MANUAL.epi.length),5);
  await page.getByRole('button',{name:'Prepare design inputs'}).click();
  await page.locator('#fetchres').waitFor();
  assert.equal(stagedRequest.searchParams.get('chain'),'Y');
  assert.ok(stagedRequest.searchParams.get('epitopes').includes('Y:'));
  for(const [key,value] of Object.entries({rfd3:'8',mpnn:'8',ipd:'false',prefilter:'true',keepfrac:'0.5'})){
    assert.equal(stagedRequest.searchParams.get(key),value);
  }
  assert.equal(await page.locator('#simulation').isVisible(),false);
  complex=true;
  await upload();
  assert.equal(await page.locator('#residue-help').isVisible(),false);
  await page.locator('#moderes').click();
  await drag('H',2,5);
  assert.equal(await page.locator('#seqbody .absel').count(),4);
  await drag('H',5,2);
  const wrap=await letters('H').evaluateAll(nodes=>nodes.findIndex((el,i)=>
    i>0 && el.offsetTop!==nodes[i-1].offsetTop));
  assert.ok(wrap>2,'Fixture sequence should wrap');
  await drag('H',wrap-2,wrap+2);
  assert.equal(await page.locator('#seqbody .absel').count(),5);
  await drag('H',wrap+2,wrap-2);
  assert.deepEqual(await page.evaluate(()=>CFG.__upload.epi),[]);
  assert.equal(await page.locator('#seqbody .absel').count(),0);
  assert.deepEqual(errors,[]);
  console.log('Passed fixed-default staging, simulation, drag selection, reversal, removal, click and Shift-click in the sequence picker.');
}finally{await browser.close();}
