/* Per-design structure viewer — a docked 3Dmol panel for the results page. Self-contained: injects
   its own styles + dock markup, owns the viewer and the layer toggles. Driven by the viewer/<id>.pdb
   (aligned scaffold) + viewer/native.pdb (shared antigen+antibody) + viewer/<id>.json that
   rfd4/viewer.py emits. Needs 3Dmol loaded globally ($3Dmol).

   API:  DesignView.open({native, scaffold, cfg, title, metrics})  DesignView.close()
   metrics: {epitope_rmsd, overall_rmsd, overall_pae, backbone_clashes, composite} (any may be null). */
(function(){
  "use strict";
  var G="#2f9e6e", NG="#707a86", RED="#c0392b", BLUE="#2b47d6", AM="#d98b1f";
  var S={scaffold:true, native:true, epi:true, ab:true, clash:true, aom:false, spin:false};
  var built=false, viewer=null, scaffModel=null, nativeText=null, aom=[], CFG=null, els={};

  function css(){
    var s=document.createElement("style");
    s.textContent=
    ".dv-dock{position:fixed;top:0;right:0;width:440px;max-width:100vw;height:100vh;z-index:1000;"+
      "background:#fff;border-left:1px solid #d5deea;box-shadow:-6px 0 22px rgba(20,40,80,.10);"+
      "transform:translateX(100%);transition:transform .22s ease;display:flex;flex-direction:column;overflow:auto}"+
    ".dv-dock.open{transform:none}@media(prefers-reduced-motion:reduce){.dv-dock{transition:none}}"+
    "@media(min-width:1000px){.dv-open .wrap{margin-right:460px}}"+
    "@media(max-width:620px){.dv-dock{width:100vw}}"+
    ".dv-h{font-size:12px;font-weight:bold;letter-spacing:.03em;text-transform:uppercase;color:#3b6ea5;"+
      "background:#eef3f9;border-bottom:1px solid #d5deea;padding:9px 13px;display:flex;align-items:center;justify-content:space-between;gap:8px;position:sticky;top:0}"+
    ".dv-x{border:0;background:none;color:#555;cursor:pointer;font-size:18px;line-height:1;padding:2px 5px}.dv-x:hover{color:#1a1a1a}"+
    ".dv-stats{display:flex;gap:5px;flex-wrap:wrap;padding:10px 13px 2px}"+
    ".dv-chip{display:flex;flex-direction:column;align-items:flex-end;border:1px solid #d5deea;border-radius:5px;padding:3px 7px;min-width:52px}"+
    ".dv-chip .k{font-size:9px;text-transform:uppercase;letter-spacing:.04em;color:#555}"+
    ".dv-chip .v{font-variant-numeric:tabular-nums;font-size:12.5px;font-weight:bold}"+
    ".dv-body{padding:10px 13px 14px}"+
    ".dv-stage{position:relative;width:100%;height:340px;border:1px solid #ccc;background:#fff;border-radius:6px;overflow:hidden}"+
    ".dv-stage>div.v{position:absolute;inset:0}"+
    ".dv-hint{position:absolute;left:9px;bottom:8px;font-size:10.5px;color:#555;background:rgba(255,255,255,.82);border:1px solid #d5deea;border-radius:5px;padding:2px 6px;z-index:2}"+
    ".dv-legend{position:absolute;right:9px;bottom:8px;background:rgba(255,255,255,.92);border:1px solid #d5deea;border-radius:6px;padding:6px 8px;font-size:11px;z-index:2;display:flex;flex-direction:column;gap:3px}"+
    ".dv-legend .r{display:flex;align-items:center;gap:6px}.dv-legend .r[hidden]{display:none}"+
    ".dv-sw{width:10px;height:10px;border-radius:2px;flex:none}.dv-sw.mesh{border:1.5px solid "+BLUE+";background:rgba(43,71,214,.14);border-radius:50%}"+
    ".dv-fallback{display:none;padding:16px;color:#555;text-align:center;font-size:12px}"+
    ".dv-ctrls{display:flex;gap:6px;flex-wrap:wrap;margin-top:10px}"+
    ".dv-pill{font:inherit;font-size:11.5px;padding:4px 10px;border:1.5px solid #cfd8e6;background:#fff;color:#1a1a1a;border-radius:13px;cursor:pointer;display:inline-flex;align-items:center;gap:5px;user-select:none}"+
    ".dv-pill .dot{width:8px;height:8px;border-radius:2px;flex:none;opacity:.35}"+
    ".dv-pill[aria-pressed=true]{border-color:#1b4a86;background:#eef3f9;font-weight:bold;box-shadow:inset 0 0 0 1px #296DC0}"+
    ".dv-pill[aria-pressed=true] .dot{opacity:1}.dv-pill:focus-visible{outline:2px solid #296DC0;outline-offset:2px}.dv-pill.plain{border-radius:6px}"+
    ".dv-foot{font-size:10.5px;color:#555;margin-top:8px}"+
    ".dv-seq{margin-top:9px;font-family:'Courier New',monospace;font-size:11px;line-height:1.7;word-break:break-all;border:1px solid #d5deea;border-radius:4px;padding:7px 9px;background:#fcfdfe;max-height:66px;overflow:auto}"+
    ".dv-seq .e{background:#f7d9d4;outline:1px solid "+RED+";border-radius:2px}.dv-seq .c{background:#fdefcf;outline:1px solid "+AM+";border-radius:2px}"+
    ".dv-pan{position:absolute;top:8px;left:8px;z-index:3;display:grid;grid-template-columns:repeat(3,22px);grid-template-rows:repeat(3,22px);gap:2px}"+
    ".dv-pan button{border:1px solid #d5deea;background:rgba(255,255,255,.9);border-radius:4px;cursor:pointer;font-size:10px;line-height:1;color:#3b6ea5;padding:0;display:flex;align-items:center;justify-content:center}"+
    ".dv-pan button:hover{background:#eef3f9;color:#1b4a86}.dv-pan button:focus-visible{outline:2px solid #296DC0;outline-offset:1px}"+
    ".dv-pan .u{grid-area:1/2}.dv-pan .l{grid-area:2/1}.dv-pan .r{grid-area:2/3}.dv-pan .d{grid-area:3/2}";
    document.head.appendChild(s);
  }

  function pill(id,label,color,pressed){
    return "<button class='dv-pill' id='"+id+"' aria-pressed='"+(pressed?"true":"false")+"'>"+
      (color?"<span class='dot' style='background:"+color+"'></span>":"")+label+"</button>";
  }
  function build(){
    css();
    var d=document.createElement("div"); d.className="dv-dock"; d.setAttribute("role","dialog");
    d.setAttribute("aria-label","design structure viewer");
    d.innerHTML=
    "<div class='dv-h'><span id='dv-title'>Design</span><button class='dv-x' id='dv-close' aria-label='Close viewer'>&times;</button></div>"+
    "<div class='dv-stats'>"+
      "<div class='dv-chip'><span class='k'>Epitope RMSD</span><span class='v' id='dv-epi'>&ndash;</span></div>"+
      "<div class='dv-chip'><span class='k'>Overall RMSD</span><span class='v' id='dv-overall'>&ndash;</span></div>"+
      "<div class='dv-chip'><span class='k'>Mean PAE</span><span class='v' id='dv-pae'>&ndash;</span></div>"+
      "<div class='dv-chip'><span class='k'>Bb clashes</span><span class='v' id='dv-bb'>&ndash;</span></div>"+
      "<div class='dv-chip'><span class='k'>Composite</span><span class='v' id='dv-comp'>&ndash;</span></div></div>"+
    "<div class='dv-body'><div class='dv-stage'><div class='v' id='dv-viewer'></div>"+
      "<div class='dv-pan' id='dv-pan' aria-label='pan view'>"+
        "<button class='u' data-dx='0' data-dy='-1' title='Pan up' aria-label='Pan up'>&#9650;</button>"+
        "<button class='l' data-dx='-1' data-dy='0' title='Pan left' aria-label='Pan left'>&#9664;</button>"+
        "<button class='r' data-dx='1' data-dy='0' title='Pan right' aria-label='Pan right'>&#9654;</button>"+
        "<button class='d' data-dx='0' data-dy='1' title='Pan down' aria-label='Pan down'>&#9660;</button></div>"+
      "<div class='dv-hint'>drag to rotate &middot; scroll to zoom</div>"+
      "<div class='dv-legend' id='dv-legend'>"+
        "<div class='r' id='dv-lg-scaf'><span class='dv-sw' style='background:"+G+"'></span><span>scaffold (design)</span></div>"+
        "<div class='r' id='dv-lg-epi'><span class='dv-sw' style='background:"+RED+"'></span><span id='dv-lg-epi-t'>epitope on scaffold</span></div>"+
        "<div class='r' id='dv-lg-clash'><span class='dv-sw' style='background:"+AM+"'></span><span>steric clashes</span></div>"+
        "<div class='r' id='dv-lg-nat'><span class='dv-sw' style='background:"+NG+"'></span><span>native antigen</span></div>"+
        "<div class='r' id='dv-lg-ab'><span class='dv-sw' style='background:"+BLUE+"'></span><span>antibody</span></div>"+
        "<div class='r' id='dv-lg-aom' hidden><span class='dv-sw mesh'></span><span>antibody occupancy</span></div></div>"+
      "<div class='dv-fallback' id='dv-fallback'>WebGL isn't available here, so the 3D viewer can't load.</div></div>"+
      "<div class='dv-ctrls' role='group' aria-label='viewer layers'>"+
        pill("dv-t-scaffold","Scaffold design",G,true)+pill("dv-t-native","Native structure",NG,true)+
        pill("dv-t-epi","Epitope",RED,true)+pill("dv-t-ab","Antibody",BLUE,true)+
        pill("dv-t-clash","Steric clashes",AM,true)+pill("dv-t-aom","AOM occupancy mesh",BLUE,false)+
        "<button class='dv-pill plain' id='dv-t-spin' aria-pressed='false'>Spin</button>"+
        "<button class='dv-pill plain' id='dv-t-reset'>Reset view</button></div>"+
      "<div class='dv-foot' id='dv-foot'></div><div class='dv-seq' id='dv-seq'></div></div>";
    document.body.appendChild(d);
    ["title","close","epi","overall","pae","bb","comp","viewer","legend","lg-scaf","lg-epi","lg-epi-t",
     "lg-clash","lg-nat","lg-ab","lg-aom","fallback","foot","seq"].forEach(function(k){ els[k]=d.querySelector("#dv-"+k); });
    els.dock=d;
    d.querySelector("#dv-close").onclick=close;
    tog("dv-t-scaffold","scaffold"); tog("dv-t-native","native"); tog("dv-t-epi","epi");
    tog("dv-t-ab","ab"); tog("dv-t-clash","clash"); tog("dv-t-aom","aom");
    var sp=d.querySelector("#dv-t-spin");
    sp.onclick=function(){ S.spin=!S.spin; sp.setAttribute("aria-pressed",String(S.spin)); if(viewer) viewer.spin(S.spin?"y":false); };
    d.querySelector("#dv-t-reset").onclick=fit;
    d.querySelectorAll("#dv-pan button").forEach(function(b){ b.onclick=function(){
      if(viewer) viewer.translate(+b.dataset.dx*30, +b.dataset.dy*30, 120); }; });  // pan the view ~30px
    document.addEventListener("keydown",function(e){ if(e.key==="Escape" && d.classList.contains("open")) close(); });
    built=true;
  }
  function tog(id,key){ var b=els.dock.querySelector("#"+id); if(b) b.onclick=function(){
    S[key]=!S[key]; b.setAttribute("aria-pressed",String(S[key])); paint(); }; }
  function webgl(){ try{var c=document.createElement("canvas");
    return !!(window.WebGLRenderingContext&&(c.getContext("webgl")||c.getContext("experimental-webgl")));}catch(e){return false;} }
  function nativeByChain(){ var m={}; (CFG.epitopeNative||[]).forEach(function(p){ (m[p[0]]=m[p[0]]||[]).push(p[1]); }); return m; }

  function paint(){
    if(!viewer||!CFG) return;
    viewer.setStyle({},{});
    if(S.scaffold){
      viewer.setStyle({chain:CFG.scaffoldChain},{cartoon:{color:G,thickness:0.8,arrows:true}});
      if(S.epi) viewer.setStyle({chain:CFG.scaffoldChain,resi:CFG.epitopeScaffold},{cartoon:{color:RED,thickness:1.0,arrows:true}});
      if(S.clash) viewer.setStyle({chain:CFG.scaffoldChain,resi:CFG.clash},{cartoon:{color:AM,thickness:1.0,arrows:true},stick:{color:AM,radius:0.2}});
    }
    if(S.native){
      (CFG.antigenChains||[]).forEach(function(ch){ viewer.setStyle({chain:ch},{cartoon:{color:NG,opacity:0.85,thickness:0.6,arrows:true}}); });
      if(S.epi && !S.scaffold){ var m=nativeByChain(); for(var ch in m) viewer.setStyle({chain:ch,resi:m[ch]},{cartoon:{color:RED,opacity:0.95,thickness:0.9,arrows:true}}); }
    }
    if(S.ab) (CFG.antibodyChains||[]).forEach(function(ch){ viewer.setStyle({chain:ch},{cartoon:{color:BLUE,thickness:0.75,arrows:true}}); });
    while(aom.length) viewer.removeShape(aom.pop());
    if(S.aom && CFG.aomCenter){ var c={x:CFG.aomCenter[0],y:CFG.aomCenter[1],z:CFG.aomCenter[2]};
      aom.push(viewer.addSphere({center:c,radius:CFG.aomRadius,color:BLUE,opacity:0.22}));
      aom.push(viewer.addSphere({center:c,radius:CFG.aomRadius+1.4,color:BLUE,opacity:0.6,wireframe:true})); }
    els["lg-scaf"].hidden=!S.scaffold;
    els["lg-epi"].hidden=!(S.epi&&(S.scaffold||S.native));
    els["lg-epi-t"].textContent=S.scaffold?"epitope on scaffold":"epitope (native)";
    els["lg-clash"].hidden=!(S.scaffold&&S.clash);
    els["lg-nat"].hidden=!S.native; els["lg-ab"].hidden=!S.ab; els["lg-aom"].hidden=!S.aom;
    viewer.render();
  }
  function fit(){ if(!viewer) return; viewer.zoomTo(); viewer.zoom(0.9); viewer.render(); }  // frame the whole complex, not just the scaffold
  function seqHtml(){
    var q=(CFG.sequence||""), e={}, c={}, h="", i;
    (CFG.epitopeScaffold||[]).forEach(function(r){e[r]=1;}); (CFG.clash||[]).forEach(function(r){c[r]=1;});
    for(i=0;i<q.length;i++){ var r=i+1, cl=c[r]?" class=c":e[r]?" class=e":""; h+="<span"+cl+">"+q[i]+"</span>"; }
    els.seq.innerHTML=h;
  }
  function num(v,d){ return (v===null||v===undefined||v===""||isNaN(+v))?"–":(+v).toFixed(d); }

  function open(o){
    if(typeof $3Dmol==="undefined"){ return; }
    if(!built) build();
    CFG=o.cfg||{};
    var m=o.metrics||{};
    els.title.textContent=o.title||"Design";
    els.epi.textContent=num(m.epitope_rmsd,2)+(m.epitope_rmsd!=null?" Å":"");
    els.overall.textContent=num(m.overall_rmsd,2)+(m.overall_rmsd!=null?" Å":"");
    els.pae.textContent=num(m.overall_pae,2);
    els.bb.textContent=(m.backbone_clashes==null?"–":m.backbone_clashes);
    els.comp.textContent=num(m.composite,2);
    var cl=(CFG.clash||[]).length;
    els.foot.textContent=(CFG.epitopeRMSD!=null?"Epitope aligns at "+(+CFG.epitopeRMSD).toFixed(2)+" Å. ":"")+
      (CFG.antibodyChains&&CFG.antibodyChains.length?cl+" scaffold residue"+(cl===1?"":"s")+" overlap the antibody.":"No antibody in this structure.");
    els.dock.classList.add("open"); document.documentElement.classList.add("dv-open");
    if(!viewer){
      if(!webgl()){ els.fallback.style.display="block"; return; }
      viewer=$3Dmol.createViewer(els.viewer,{backgroundColor:"white",ambientOcclusion:{strength:0.9}});
    }
    if(nativeText!==o.native){ viewer.removeAllModels(); scaffModel=null;
      if(o.native) viewer.addModel(o.native,"pdb"); nativeText=o.native; }
    else if(scaffModel){ viewer.removeModel(scaffModel); scaffModel=null; }
    scaffModel=viewer.addModel(o.scaffold,"pdb");
    seqHtml(); paint(); fit();
  }
  function close(){ if(els.dock){ els.dock.classList.remove("open"); } document.documentElement.classList.remove("dv-open");
    if(viewer && S.spin){ S.spin=false; viewer.spin(false); var sp=els.dock&&els.dock.querySelector("#dv-t-spin"); if(sp) sp.setAttribute("aria-pressed","false"); } }

  window.DesignView={ open:open, close:close, isOpen:function(){ return built && els.dock.classList.contains("open"); } };
})();
