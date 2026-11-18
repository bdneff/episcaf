#!/usr/bin/env python3
"""Trim the white border of a render and make near-white pixels transparent.

Used by render_webtool_figs.sh to turn Tachyon's white-background PNGs into trimmed,
transparent-background figures that sit cleanly on any page color.

    python trim_transparent.py fig1.png fig2.png ...

Overwrites each file in place. White (and near-white) becomes transparent; everything
else is kept. A small padding is added so marks are not flush to the edge.
"""
import sys
from PIL import Image

PAD = 6
THRESH = 244   # pixels with min(r,g,b) >= THRESH are treated as background

def process(path):
    im = Image.open(path).convert("RGBA")
    px = im.load()
    w, h = im.size
    minx, miny, maxx, maxy = w, h, 0, 0
    # make near-white transparent, track content bbox
    for y in range(h):
        for x in range(w):
            r, g, b, a = px[x, y]
            if min(r, g, b) >= THRESH:
                px[x, y] = (r, g, b, 0)
            else:
                if x < minx: minx = x
                if y < miny: miny = y
                if x > maxx: maxx = x
                if y > maxy: maxy = y
    if maxx < minx:            # nothing found; leave as-is
        im.save(path)
        return
    minx = max(0, minx - PAD); miny = max(0, miny - PAD)
    maxx = min(w, maxx + PAD); maxy = min(h, maxy + PAD)
    im.crop((minx, miny, maxx, maxy)).save(path)
    print("  trimmed %s -> %dx%d" % (path, maxx - minx, maxy - miny))

if __name__ == "__main__":
    for p in sys.argv[1:]:
        process(p)
