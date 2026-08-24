# -*- coding: utf-8 -*-
"""Debug 2 : ticks = segments s'etendant a gauche de la spine (x<116) ou a l'interieur."""
import numpy as np
from PIL import Image

SCRATCH = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
img = np.asarray(Image.open(SCRATCH + r"\p12_img0.png").convert("RGB")).astype(int)

def nonwhite(a):
    return a.sum(axis=2) < 3 * 225

# Panel (a) pression : crop rows 648..946. Spine verte vers x=116-118.
# Ticks sortants -> x ~104..115 ; entrants -> x ~121..132.
for label, xs in (("sortant", (100, 115)), ("entrant", (121, 136))):
    band = nonwhite(img[648:946, xs[0]:xs[1]])
    rows = np.where(band.any(axis=1))[0]
    # cluster
    if rows.size:
        groups, start, prev = [], rows[0], rows[0]
        for rr in rows[1:]:
            if rr - prev <= 2:
                prev = rr
            else:
                groups.append((start, prev)); start = prev = rr
        groups.append((start, prev))
        info = [(a, bb, int(band[a:bb+1].sum())) for a, bb in groups]
        print(f"pression(a) {label}:", info)
    else:
        print(f"pression(a) {label}: rien")

# Meme chose panel (a) torque (rows 351..647)
for label, xs in (("sortant", (100, 115)), ("entrant", (121, 136))):
    band = nonwhite(img[351:647, xs[0]:xs[1]])
    rows = np.where(band.any(axis=1))[0]
    if rows.size:
        groups, start, prev = [], rows[0], rows[0]
        for rr in rows[1:]:
            if rr - prev <= 2:
                prev = rr
            else:
                groups.append((start, prev)); start = prev = rr
        groups.append((start, prev))
        info = [(a, bb, int(band[a:bb+1].sum())) for a, bb in groups]
        print(f"torque(a) {label}:", info)
    else:
        print(f"torque(a) {label}: rien")
