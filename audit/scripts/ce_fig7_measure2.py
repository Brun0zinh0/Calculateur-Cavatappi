# -*- coding: utf-8 -*-
"""Contre-expertise figure 7 : calibration par etiquettes de ticks + mesure des pics.

Sortie :
 - pour chaque panneau : rows des etiquettes, fit lineaire value(row), valeurs vraies
   aux bords du crop, et biais de la calibration PANEL_SPECS a divers niveaux ;
 - pics par cycle (courbe noire = theorie, marqueurs couleur = experience,
   vert = pression) en unites VRAIES.
"""
import numpy as np
from PIL import Image

SCRATCH = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
img = np.asarray(Image.open(SCRATCH + r"\p12_img0.png").convert("RGB")).astype(int)
r, g, b = img[..., 0], img[..., 1], img[..., 2]

def mask_color(name):
    if name == "blue":
        return (b > 135) & (g > 80) & (r < 135) & ((b - r) > 55)
    if name == "orange":
        return (r > 180) & (g > 60) & (g < 185) & (b < 125) & ((r - b) > 80)
    if name == "green":
        return (g > 135) & (r < 120) & (b < 150) & ((g - r) > 45)
    if name == "black":
        return (r < 105) & (g < 105) & (b < 105)

MASKS = {c: mask_color(c) for c in ("blue", "orange", "green", "black")}
NONWHITE = (img.sum(axis=2) < 3 * 225)

PANELS = {
    ("a", "force"):   dict(crop=(120, 51, 704, 350),  ticks=[2200, 2000, 1800, 1600, 1400, 1200], ps=(1200., 2200.)),
    ("a", "torque"):  dict(crop=(120, 351, 704, 647), ticks=[1650, 1100, 550, 0],                 ps=(-250., 1650.)),
    ("a", "pressure"):dict(crop=(120, 648, 704, 946), ticks=[1.5, 1.0, 0.5, 0.0],                ps=(0., 1.5)),
    ("b", "force"):   dict(crop=(899, 51, 1484, 350), ticks=[2800, 2600, 2400, 2200, 2000],      ps=(2000., 2800.)),
    ("b", "torque"):  dict(crop=(899, 351, 1484, 647),ticks=[1440, 960, 480, 0],                 ps=(-250., 1800.)),
    ("b", "pressure"):dict(crop=(899, 648, 1484, 946),ticks=[1.5, 1.0, 0.5, 0.0],                ps=(0., 1.5)),
}

def label_rows(crop, n_expected):
    """Clusters de texte dans la bande des etiquettes (a gauche de la spine)."""
    x0, y0, x1, y1 = crop
    band = NONWHITE[y0:y1, x0 - 48:x0 - 6]
    rows = np.where(band.any(axis=1))[0]
    if rows.size == 0:
        return None
    groups, start, prev = [], rows[0], rows[0]
    for rr in rows[1:]:
        if rr - prev <= 3:
            prev = rr
        else:
            groups.append((start, prev)); start = prev = rr
    groups.append((start, prev))
    # texte : hauteur 12-30 px
    groups = [gg for gg in groups if 8 <= (gg[1] - gg[0]) <= 34]
    if len(groups) != n_expected:
        return ("MISMATCH", groups)
    return [0.5 * (a + bb) for a, bb in groups]

print("=" * 90)
for key, spec in PANELS.items():
    crop = spec["crop"]
    ticks = spec["ticks"]
    res = label_rows(crop, len(ticks))
    print(f"\n=== panneau {key} crop={crop} ===")
    if res is None or (isinstance(res, tuple) and res[0] == "MISMATCH"):
        print("  clusters d'etiquettes inattendus:", res)
        continue
    rows = np.asarray(res, float)
    vals = np.asarray(ticks, float)  # ordonnes du haut vers le bas
    A, Bc = np.polyfit(rows, vals, 1)
    resid = vals - (A * rows + Bc)
    hgt = crop[3] - crop[1]
    spec["cal"] = (A, Bc)
    print(f"  rows etiquettes: {np.round(rows,1)} -> {vals}; resid max {np.abs(resid).max():.2f}")
    print(f"  VRAI  : bord haut={Bc:.2f}  bord bas={A*hgt+Bc:.2f}  (pente {A:.4f}/px)")
    lo, hi = spec["ps"]
    print(f"  SCRIPT: bord haut={hi}  bord bas={lo}  (pente {-(hi-lo)/(hgt-1):.4f}/px)")
    # biais script - vrai a plusieurs niveaux vrais
    levels = np.linspace(lo if lo < 0 else 0, vals.max(), 6)
    Aps = -(hi - lo) / (hgt - 1)
    for lv in [vals.max(), 0.75 * vals.max(), 0.5 * vals.max()]:
        row_true = (lv - Bc) / A
        script_val = hi + Aps * row_true
        print(f"    valeur vraie {lv:9.2f} -> le script lit {script_val:9.2f}  (biais {script_val - lv:+.2f})")

def cycle_peaks(mask, crop, cal, cut_true=None, col_margin=10, min_gap=18):
    x0, y0, x1, y1 = crop
    A, Bc = cal
    sub = mask[y0:y1, x0:x1].copy()
    hgt, wid = sub.shape
    sub[:4, :] = False; sub[-4:, :] = False
    if cut_true is not None:
        row_cut = int((cut_true - Bc) / A)  # rows < row_cut sont au-dessus de cut_true
        sub[:max(0, row_cut), :] = False
    top = np.full(wid, np.nan)
    for cc in range(col_margin, wid - col_margin):
        rr = np.where(sub[:, cc])[0]
        if rr.size:
            top[cc] = A * rr.min() + Bc
    peaks = []
    finite = np.isfinite(top)
    med = np.nanmedian(top[finite]) if finite.any() else np.nan
    cc = col_margin
    while cc < wid - col_margin:
        if np.isfinite(top[cc]):
            lo_ = max(col_margin, cc - min_gap); hi_ = min(wid - col_margin, cc + min_gap + 1)
            seg = top[lo_:hi_]
            if top[cc] == np.nanmax(seg) and top[cc] > med:
                peaks.append((cc, top[cc])); cc += min_gap; continue
        cc += 1
    return top, peaks

CUTS = {("a", "force"): 1870., ("a", "torque"): 1020., ("b", "force"): 2500., ("b", "torque"): 1080.,
        ("a", "pressure"): None, ("b", "pressure"): None}
EXPC = {"force": "blue", "torque": "orange", "pressure": "green"}

print("\n" + "=" * 90)
for key, spec in PANELS.items():
    if "cal" not in spec:
        continue
    panel, kind = key
    unit = {"force": "mN", "torque": "uNm", "pressure": "MPa"}[kind]
    fmt = "%.3f" if kind == "pressure" else "%.0f"
    print(f"\n### {key} (unites vraies) ###")
    _, pk_exp = cycle_peaks(MASKS[EXPC[kind]], spec["crop"], spec["cal"], CUTS[key])
    vexp = np.asarray([p[1] for p in pk_exp])
    if vexp.size:
        print(f"  exp/pression ({EXPC[kind]}): n={vexp.size}  max={fmt % vexp.max()}  moy={fmt % vexp.mean()}  {unit}")
        print("   pics:", [fmt % v for _, v in pk_exp])
    if kind != "pressure":
        _, pk_th = cycle_peaks(MASKS["black"], spec["crop"], spec["cal"], CUTS[key])
        vth = np.asarray([p[1] for p in pk_th])
        if vth.size:
            print(f"  theorie (noir): n={vth.size}  max={fmt % vth.max()}  moy={fmt % vth.mean()}  {unit}")
            print("   pics:", [fmt % v for _, v in pk_th])
