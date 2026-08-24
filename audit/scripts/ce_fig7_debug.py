# -*- coding: utf-8 -*-
"""Debug: ou sont les ticks ? inspection de bandes autour du cadre gauche."""
import numpy as np
from PIL import Image

SCRATCH = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
img = np.asarray(Image.open(SCRATCH + r"\p12_img0.png").convert("RGB")).astype(int)

# colonne du cadre gauche panel (a) : cherchons la colonne noire verticale vers x=120
col_dark = []
for x in range(100, 140):
    col = img[651:940, x]
    dark = ((col.max(axis=1) < 120)).sum()
    col_dark.append((x, dark))
print("colonnes sombres (pression a, rows 651-940):", [c for c in col_dark if c[1] > 50])

# echantillon de couleurs dans la bande x=104..122 aux rows 651..945 ou il y a des pixels non blancs
band = img[648:946, 100:124]
nonwhite = np.where(band.sum(axis=2) < 3 * 235)
rows = sorted(set(nonwhite[0].tolist()))
print("rows non-blancs dans bande gauche pression(a):", rows[:80])
# couleurs presentes
for rr in rows[:10]:
    cols = nonwhite[1][nonwhite[0] == rr]
    print(" row", rr, [(cc, tuple(band[rr, cc])) for cc in cols[:6]])
