# -*- coding: utf-8 -*-
"""Contre-expertise : mesure directe des pics de la figure 7 (image embarquee p12_img0.png).

1. Detection des ticks majeurs sur l'axe gauche de chaque panneau (segments colores
   a gauche du cadre) -> calibration lineaire valeur(row).
2. Mesure des pics par cycle : courbe noire (theorie), marqueurs colores (exp),
   courbe verte (pression).
3. Comparaison avec la calibration PANEL_SPECS de validation_figure7.py.
"""
import numpy as np
from PIL import Image

SCRATCH = r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
img = np.asarray(Image.open(SCRATCH + r"\p12_img0.png").convert("RGB")).astype(float)
H, W, _ = img.shape
print("image", img.shape)

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
    raise ValueError(name)

MASKS = {c: mask_color(c) for c in ("blue", "orange", "green", "black")}

# Panneaux : crop du cadre trace (x0, y0, x1, y1) comme PANEL_SPECS + couleur axe + valeurs ticks majeurs
PANELS = {
    ("a", "force"):   dict(crop=(120, 51, 704, 350),  axcolor="blue",   ticks=[1200, 1400, 1600, 1800, 2000, 2200]),
    ("a", "torque"):  dict(crop=(120, 351, 704, 647), axcolor="orange", ticks=[0, 550, 1100, 1650]),
    ("a", "pressure"):dict(crop=(120, 648, 704, 946), axcolor="green",  ticks=[0.0, 0.5, 1.0, 1.5]),
    ("b", "force"):   dict(crop=(899, 51, 1484, 350), axcolor="blue",   ticks=[2000, 2200, 2400, 2600, 2800]),
    ("b", "torque"):  dict(crop=(899, 351, 1484, 647),axcolor="orange", ticks=[0, 480, 960, 1440]),
    ("b", "pressure"):dict(crop=(899, 648, 1484, 946),axcolor="green",  ticks=[0.0, 0.5, 1.0, 1.5]),
}

def detect_ticks(crop, axcolor, n_expected):
    """Ticks majeurs = segments colores dans une bande a gauche du cadre."""
    x0, y0, x1, y1 = crop
    best = None
    # cherche la bande la plus informative entre 4 et 14 px a gauche du cadre
    for dx in (14, 12, 10, 8, 6):
        strip = MASKS[axcolor][y0:y1, max(0, x0 - dx):x0 - 2]
        rows = np.where(strip.any(axis=1))[0]
        if rows.size == 0:
            continue
        # clustering des rangees contigues
        groups = []
        start = prev = rows[0]
        for rr in rows[1:]:
            if rr - prev <= 2:
                prev = rr
            else:
                groups.append((start, prev))
                start = prev = rr
        groups.append((start, prev))
        centers = [0.5 * (a + bb) for a, bb in groups]
        lengths = []
        for a, bb in groups:
            seg = MASKS[axcolor][y0 + a:y0 + bb + 1, max(0, x0 - 20):x0 - 2]
            lengths.append(seg.sum())
        if len(centers) >= n_expected:
            # garde les n_expected plus "longs" (ticks majeurs vs mineurs)
            idx = np.argsort(lengths)[::-1][:n_expected]
            sel = sorted(np.asarray(centers)[sorted(idx)])
            if best is None:
                best = sel
        if best is not None:
            break
    return best  # rows relatifs au crop

for key, spec in PANELS.items():
    crop = spec["crop"]
    ticks_vals = spec["ticks"]
    centers = detect_ticks(crop, spec["axcolor"], len(ticks_vals))
    print("\n=== panel", key, "===")
    if centers is None:
        print("  ECHEC detection ticks")
        continue
    # valeurs decroissantes avec row : la plus haute valeur = row le plus petit
    rows = np.asarray(centers, float)
    vals = np.asarray(sorted(ticks_vals, reverse=True), float)
    # regression lineaire value = A*row + B
    A, B = np.polyfit(rows, vals, 1)
    resid = vals - (A * rows + B)
    print("  tick rows (crop):", np.round(rows, 1), "-> vals", vals, " resid max", np.abs(resid).max())
    spec["cal"] = (A, B)
    # valeur au bord haut (row 0) et bas (row = hauteur du crop)
    hgt = crop[3] - crop[1]
    print(f"  valeur au bord haut du crop = {B:.1f} ; au bord bas = {A * hgt + B:.1f}")
    # comparaison avec PANEL_SPECS (mapping lineaire bord haut / bord bas)
    ps_ylim = {("a","force"):(1200.,2200.),("a","torque"):(-250.,1650.),("a","pressure"):(0.,1.5),
               ("b","force"):(2000.,2800.),("b","torque"):(-250.,1800.),("b","pressure"):(0.,1.5)}[key]
    print(f"  PANEL_SPECS suppose bord haut = {ps_ylim[1]}, bord bas = {ps_ylim[0]}")

def cycle_peaks(mask, crop, cal, row_min_value=None, col_margin=10, min_gap=15):
    """Pic (valeur max) par colonne puis extraction des maxima locaux (cycles)."""
    x0, y0, x1, y1 = crop
    A, B = cal
    sub = mask[y0:y1, x0:x1].copy()
    hgt, wid = sub.shape
    if row_min_value is not None:
        # exclut les rows dont la valeur depasse row_min_value (zone legende)
        row_vals = A * np.arange(hgt) + B
        sub[row_vals > row_min_value, :] = False
    top = np.full(wid, np.nan)
    for cc in range(col_margin, wid - col_margin):
        rr = np.where(sub[:, cc])[0]
        if rr.size:
            top[cc] = A * rr.min() + B
    # maxima locaux espaces de min_gap colonnes
    peaks = []
    cc = col_margin
    while cc < wid - col_margin:
        if np.isfinite(top[cc]):
            lo = max(col_margin, cc - min_gap)
            hi = min(wid - col_margin, cc + min_gap + 1)
            seg = top[lo:hi]
            if np.isfinite(seg).any() and top[cc] == np.nanmax(seg) and top[cc] > np.nanmedian(top[np.isfinite(top)]):
                peaks.append((cc, top[cc]))
                cc += min_gap
                continue
        cc += 1
    return top, peaks

def summarize(name, peaks, unit):
    if not peaks:
        print(f"  {name}: aucun pic")
        return
    v = np.asarray([p[1] for p in peaks])
    print(f"  {name}: n={v.size} pics ; max={v.max():.1f} ; moyenne={v.mean():.1f} ; min={v.min():.1f} {unit}")

# seuils anti-legende (valeurs au-dessus desquelles on ignore les pixels)
LEGEND_CUT = {("a","force"): 1900., ("a","torque"): 1050., ("b","force"): 2550., ("b","torque"): 1100.,
              ("a","pressure"): None, ("b","pressure"): None}
EXP_COLOR = {"force": "blue", "torque": "orange", "pressure": "green"}

for key, spec in PANELS.items():
    if "cal" not in spec:
        continue
    panel, kind = key
    crop, cal = spec["crop"], spec["cal"]
    unit = {"force": "mN", "torque": "uNm", "pressure": "MPa"}[kind]
    print(f"\n### {key} ###")
    cut = LEGEND_CUT[key]
    # courbe experimentale (couleur)
    topc, pk = cycle_peaks(MASKS[EXP_COLOR[kind]], crop, cal, row_min_value=cut)
    summarize("exp (couleur)", pk, unit)
    if kind != "pressure":
        topb, pkb = cycle_peaks(MASKS["black"], crop, cal, row_min_value=cut)
        summarize("theorie (noir)", pkb, unit)
        # liste des pics noirs
        print("   pics noirs:", [f"{c}:{v:.0f}" for c, v in pkb])
        print("   pics exp :", [f"{c}:{v:.0f}" for c, v in pk])
    else:
        print("   pics pression:", [f"{c}:{v:.3f}" for c, v in pk])
