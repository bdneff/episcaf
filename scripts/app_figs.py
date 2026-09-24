#!/usr/bin/env python3
"""Application-concept illustrations for the "Why this matters" slide (diagnostics + vaccine).

Clean science-illustration SVGs in the deck palette, then rasterized by scripts/render_app_figs.sh
(headless Chrome). Programmatic so the well grid and nanoparticle display come from loops.

  Diagnostics: a microplate whose columns are scaffolded variant epitopes; serum lights up only its
               match (variant ID without sequencing). A magnifier zooms one well to the binding event.
  Vaccine:     syringe -> magnifier -> a nanoparticle displaying many related scaffolded epitopes, so
               the immune system learns the whole family at once.

Run:  /usr/bin/python3 scripts/app_figs.py     # writes build/appfig_{diagnostics,vaccine}.svg
"""
import math, pathlib, base64

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT/"build"; OUT.mkdir(exist_ok=True)
def datauri(p):
    return "data:image/png;base64,"+base64.b64encode(pathlib.Path(p).read_bytes()).decode()
CHIP_BIND=OUT/"chips/chip_binding.png"   # VMD scaffold+epitope+scFv (transparent)
CHIP_UNIT=OUT/"chips/chip_unit.png"      # VMD scaffold+epitope (transparent)
BLUE="#4285F4"; BLUE_D="#2456E6"; RED="#C0392B"; ORANGE="#E8820E"; PINK="#D6336C"
SILVER="#C7CDD6"; SILVER_L="#E9EDF2"; SILVER_D="#9AA3AE"
GREY_BG="#EEF2F7"; INK="#212121"; MUTE="#5b6570"; RIM="#8A94A3"; WELL="#B9C2CE"

DEFS = f"""
<defs>
  <filter id="soft" x="-30%" y="-30%" width="160%" height="160%">
    <feDropShadow dx="0" dy="5" stdDeviation="7" flood-color="#1b2a4a" flood-opacity="0.16"/>
  </filter>
  <filter id="soft2" x="-40%" y="-40%" width="180%" height="180%">
    <feDropShadow dx="0" dy="3" stdDeviation="4" flood-color="#1b2a4a" flood-opacity="0.22"/>
  </filter>
  <radialGradient id="wellPos" cx="38%" cy="34%" r="70%">
    <stop offset="0%" stop-color="#8fb4ff"/><stop offset="55%" stop-color="{BLUE}"/>
    <stop offset="100%" stop-color="{BLUE_D}"/>
  </radialGradient>
  <radialGradient id="core" cx="36%" cy="30%" r="75%">
    <stop offset="0%" stop-color="#F2F5F9"/><stop offset="45%" stop-color="{SILVER}"/>
    <stop offset="100%" stop-color="#727C89"/>
  </radialGradient>
  <radialGradient id="glass" cx="34%" cy="28%" r="80%">
    <stop offset="0%" stop-color="#ffffff" stop-opacity="0.85"/>
    <stop offset="60%" stop-color="#eaf0f8" stop-opacity="0.35"/>
    <stop offset="100%" stop-color="#dbe3ee" stop-opacity="0.30"/>
  </radialGradient>
  <linearGradient id="scaf" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0%" stop-color="{SILVER_L}"/><stop offset="100%" stop-color="{SILVER_D}"/>
  </linearGradient>
  <linearGradient id="barrel" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0%" stop-color="#ffffff"/><stop offset="45%" stop-color="#eef3f9"/>
    <stop offset="100%" stop-color="#d7dee8"/>
  </linearGradient>
  <linearGradient id="fluid" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0%" stop-color="#6fa0ff"/><stop offset="100%" stop-color="{BLUE_D}"/>
  </linearGradient>
</defs>
"""

def txt(x,y,s,size=22,fill=INK,anchor="middle",weight="normal",style="normal",spacing=None):
    sp=f' letter-spacing="{spacing}"' if spacing else ""
    return (f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
            f'fill="{fill}" text-anchor="{anchor}" font-weight="{weight}" font-style="{style}"{sp}>{s}</text>')

def binding_unit(cx, cy, s=1.0):
    """scaffold mound (silver) + red epitope + blue antibody clamped on top. Local, scaled."""
    g=[f'<g transform="translate({cx},{cy}) scale({s})">']
    g.append(f'<ellipse cx="0" cy="72" rx="82" ry="34" fill="url(#scaf)" stroke="{SILVER_D}" stroke-width="2"/>')
    g.append(f'<rect x="-48" y="26" width="96" height="22" rx="11" fill="{RED}" stroke="#8f2a1f" stroke-width="1.5"/>')
    # antibody Y (blue): stem + two Fab arms, clamped onto the red epitope
    g.append(f'<path d="M 0 24 L 0 -8 M 0 -8 L -34 -48 M 0 -8 L 34 -48" fill="none" stroke="{BLUE}" '
             f'stroke-width="17" stroke-linecap="round" stroke-linejoin="round"/>')
    g.append(f'<circle cx="-34" cy="-48" r="9" fill="{BLUE_D}"/><circle cx="34" cy="-48" r="9" fill="{BLUE_D}"/>')
    # signal spark
    for a in (0,60,120):
        rad=math.radians(a); x2=64+22*math.cos(rad); y2=-2+22*math.sin(rad)
        g.append(f'<line x1="60" y1="0" x2="{x2:.0f}" y2="{y2:.0f}" stroke="#F4B400" stroke-width="6" stroke-linecap="round"/>')
    g.append('</g>')
    return "".join(g)

def diagnostics():
    W,H=1040,860
    cols=8; rows=6; x0=145; dx=60; y0=270; dy=56
    labels=["Wu","α","β","γ","δ","κ","λ","ο"]; pos=4  # detected variant column
    posrows=set(range(rows))  # the whole detected column lights up
    extra=set()               # no background: keep the specificity message clean
    px,py,pw,ph=78,182,510,445
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',DEFS]
    s.append(f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>')
    # plate
    s.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" rx="30" fill="{GREY_BG}" stroke="{"#CDD6E0"}" stroke-width="3" filter="url(#soft)"/>')
    s.append(f'<path d="M {px+18} {py} q -18 0 -18 18 l 0 22 z" fill="#DFE6EE"/>')  # A1 corner notch
    # column headers (variants)
    for c in range(cols):
        cx=x0+c*dx; hot=(c==pos)
        s.append(txt(cx,py+38,labels[c],size=23 if hot else 20,fill=BLUE if hot else MUTE,weight="bold" if hot else "normal"))
    s.append(txt(x0+pos*dx,py+60,"detected",size=13,fill=BLUE,weight="bold"))
    # wells
    for r in range(rows):
        for c in range(cols):
            cx=x0+c*dx; cy=y0+r*dy
            on=(c==pos and r in posrows) or (r,c) in extra
            if on:
                s.append(f'<circle cx="{cx}" cy="{cy}" r="20" fill="url(#wellPos)" stroke="{BLUE_D}" stroke-width="2.5"/>')
                s.append(f'<circle cx="{cx-6}" cy="{cy-7}" r="5" fill="#ffffff" opacity="0.55"/>')
            else:
                s.append(f'<circle cx="{cx}" cy="{cy}" r="20" fill="#ffffff" stroke="{WELL}" stroke-width="2"/>')
    s.append(txt(px+pw/2,py+ph+34,"each column = one variant's scaffolded epitope",size=17,fill=MUTE,style="italic"))
    # zoom cone from the detected column to the lens
    Lx,Ly,Lr=812,330,158
    s.append(f'<polygon points="{x0+pos*dx+16},{y0} {Lx-Lr*0.72},{Ly-Lr*0.72} {Lx-Lr*0.72},{Ly+Lr*0.72} {x0+pos*dx+16},{y0+2*dy}" fill="{BLUE}" opacity="0.06"/>')
    # magnifier
    s.append(f'<line x1="{Lx+Lr*0.66:.0f}" y1="{Ly+Lr*0.66:.0f}" x2="{Lx+Lr*0.66+96:.0f}" y2="{Ly+Lr*0.66+96:.0f}" stroke="{RIM}" stroke-width="26" stroke-linecap="round"/>')
    s.append(f'<clipPath id="lensD"><circle cx="{Lx}" cy="{Ly}" r="{Lr-8}"/></clipPath>')
    s.append(f'<circle cx="{Lx}" cy="{Ly}" r="{Lr-8}" fill="url(#glass)"/>')
    bw=232; bh=int(bw*1219/805)
    s.append(f'<g clip-path="url(#lensD)"><image href="{datauri(CHIP_BIND)}" x="{Lx-bw/2:.0f}" y="{Ly-bh*0.56:.0f}" width="{bw}" height="{bh}"/></g>')
    s.append(f'<circle cx="{Lx}" cy="{Ly}" r="{Lr}" fill="none" stroke="{RIM}" stroke-width="12"/>')
    s.append(f'<circle cx="{Lx}" cy="{Ly}" r="{Lr-8}" fill="none" stroke="#ffffff" stroke-width="3" opacity="0.7"/>')
    s.append(txt(Lx,Ly+Lr+30,"serum antibody binds its epitope, and reports",size=15,fill=MUTE,style="italic"))
    s.append('</svg>')
    (OUT/"appfig_diagnostics.svg").write_text("\n".join(s))
    return W,H

def vaccine():
    W,H=1240,860
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',DEFS]
    s.append(f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>')
    # ---- syringe (points right) ----
    by,bh=352,120; bx,bw=96,352
    s.append(f'<g filter="url(#soft2)">')
    s.append(f'<rect x="{bx-56}" y="{by+24}" width="56" height="{bh-48}" rx="8" fill="#cfd7e2"/>')  # thumb flange
    s.append(f'<rect x="{bx-40}" y="{by+bh/2-9}" width="150" height="18" rx="9" fill="#c3ccd8"/>')   # plunger rod
    s.append(f'<rect x="{bx}" y="{by}" width="{bw}" height="{bh}" rx="16" fill="url(#barrel)" stroke="{"#AEB8C4"}" stroke-width="2.5"/>')
    s.append(f'<rect x="{bx+96}" y="{by+14}" width="212" height="{bh-28}" rx="10" fill="url(#fluid)"/>')  # fluid
    for i in range(6):  # graduation ticks
        tx=bx+120+i*34; s.append(f'<line x1="{tx}" y1="{by+8}" x2="{tx}" y2="{by+26}" stroke="{"#8FA0B4"}" stroke-width="2"/>')
    s.append(f'<polygon points="{bx+bw},{by+18} {bx+bw+40},{by+bh/2} {bx+bw},{by+bh-18}" fill="url(#barrel)" stroke="{"#AEB8C4"}" stroke-width="2"/>')  # tip cone
    s.append(f'<rect x="{bx+bw+38}" y="{by+bh/2-3}" width="118" height="6" rx="3" fill="{"#9AA6B4"}"/>')  # needle
    s.append('</g>')
    # ---- injection arrow to the lens ----
    tipx=bx+bw+156; cy=by+bh/2
    Lx,Ly,Lr=880,410,206
    ax2=Lx-Lr-14
    s.append(f'<path d="M {tipx} {cy} C {tipx+70} {cy-10}, {ax2-70} {Ly+30}, {ax2} {Ly}" fill="none" stroke="{BLUE}" stroke-width="4" stroke-dasharray="2 10" stroke-linecap="round"/>')
    s.append(f'<polygon points="{ax2},{Ly-9} {ax2+18},{Ly} {ax2},{Ly+9}" fill="{BLUE}"/>')
    # ---- magnifier with nanoparticle ----
    s.append(f'<line x1="{Lx+Lr*0.68:.0f}" y1="{Ly+Lr*0.68:.0f}" x2="{Lx+Lr*0.68+110:.0f}" y2="{Ly+Lr*0.68+110:.0f}" stroke="{RIM}" stroke-width="30" stroke-linecap="round"/>')
    s.append(f'<circle cx="{Lx}" cy="{Ly}" r="{Lr-8}" fill="url(#glass)"/>')
    s.append(f'<clipPath id="lensV"><circle cx="{Lx}" cy="{Ly}" r="{Lr-8}"/></clipPath>')
    s.append(f'<g clip-path="url(#lensV)">')
    core_r=84; Rc=128; n=8
    uw=96; uh=int(uw*1220/1005)  # chip aspect (~116); kept upright so the display reads cleanly
    unit=datauri(CHIP_UNIT)
    for i in range(n):
        a=math.radians(i*360/n - 90)
        cxu=Lx+Rc*math.cos(a); cyu=Ly+Rc*math.sin(a)
        ex=Lx+core_r*math.cos(a); ey=Ly+core_r*math.sin(a)
        s.append(f'<line x1="{ex:.0f}" y1="{ey:.0f}" x2="{cxu:.0f}" y2="{cyu:.0f}" stroke="{SILVER_D}" stroke-width="7" stroke-linecap="round"/>')
        s.append(f'<image href="{unit}" x="{cxu-uw/2:.0f}" y="{cyu-uh/2:.0f}" width="{uw}" height="{uh}"/>')
    s.append(f'<circle cx="{Lx}" cy="{Ly}" r="{core_r}" fill="url(#core)" stroke="#5f6875" stroke-width="2"/>')
    s.append(f'<ellipse cx="{Lx-24}" cy="{Ly-26}" rx="24" ry="15" fill="#ffffff" opacity="0.35"/>')
    s.append('</g>')
    s.append(f'<circle cx="{Lx}" cy="{Ly}" r="{Lr}" fill="none" stroke="{RIM}" stroke-width="13"/>')
    s.append(f'<circle cx="{Lx}" cy="{Ly}" r="{Lr-8}" fill="none" stroke="#ffffff" stroke-width="3" opacity="0.7"/>')
    # labels
    s.append(txt(Lx,Ly+Lr+40,"one nanoparticle, many related epitopes",size=19,fill=INK,weight="bold"))
    s.append(f'<g transform="translate({bx+40},{by-56})">')
    s.append(txt(0,0,"vaccine",size=19,fill=MUTE,anchor="start",style="italic"))
    s.append('</g>')
    s.append(txt(bx+40,by+bh+96,"many related scaffolded epitopes on one particle",size=17,fill=MUTE,anchor="start",style="italic"))
    s.append('</svg>')
    (OUT/"appfig_vaccine.svg").write_text("\n".join(s))
    return W,H

if __name__=="__main__":
    d=diagnostics(); v=vaccine()
    print(f"diagnostics {d}  vaccine {v}  -> {OUT}")
