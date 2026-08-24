# -*- coding: utf-8 -*-
"""Synthese de la matrice figure 7 : tableau CSV/markdown + figures.

- Lit matrix_results.jsonl (runs de ce_matrix_runner.py).
- Cibles : theorie article et ancrages texte exp (unites vraies).
- Figures : overlay meilleure config vs config conforme article ; barres des ecarts.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

SCRATCH = Path(
    r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude"
    r"\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-"
    r"artifici-le-mod-le-mod-le-chinois-Espace-de-travail"
    r"\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
)

# Cibles (mesure pixel avec calibration vraie, cf. ce_fig7_measure2)
THEORY = {
    (0.8, "force"): (1758.0, 1806.0, 1783.0),
    (1.0, "force"): (2389.0, 2466.0, 2427.0),
    (0.8, "torque"): (591.0, 733.0, 664.0),
    (1.0, "torque"): (702.0, 876.0, 789.0),
}
EXP_TEXT = {
    (0.8, "force"): 1824.60,
    (1.0, "force"): 2424.67,
    (0.8, "torque"): 877.01,
    (1.0, "torque"): 873.41,
}

recs = [json.loads(line) for line in open(SCRATCH / "matrix_results.jsonl", encoding="utf-8")]

rows = []
for r in recs:
    eps = r["eps"]
    out = dict(
        tag=r["tag"], bias=r["bias"], section=r["section"], eps=eps, pmode=r["pmode"],
        disc=f"L{r['n_layers']}/p{r['n_phi']}/dt{r['dt']}",
        p_peak=r["pressure"]["peak_mean"],
        F_peak=r["force"]["peak_mean"], F_ampl=r["force"]["amplitude"],
        T_peak=r["torque"]["peak_mean"], T_ampl=r["torque"]["amplitude"],
    )
    for kind, key in (("force", "F"), ("torque", "T")):
        lo, hi, mid = THEORY[(eps, kind)]
        out[f"{key}_vs_theorie_pct"] = 100.0 * (out[f"{key}_peak"] - mid) / mid
        out[f"{key}_vs_exp_pct"] = 100.0 * (out[f"{key}_peak"] - EXP_TEXT[(eps, kind)]) / EXP_TEXT[(eps, kind)]
    rows.append(out)

# ---- CSV ----
cols = ["tag", "bias", "section", "eps", "pmode", "disc", "p_peak",
        "F_peak", "F_vs_theorie_pct", "F_vs_exp_pct", "F_ampl",
        "T_peak", "T_vs_theorie_pct", "T_vs_exp_pct", "T_ampl"]
with open(SCRATCH / "matrix_summary.csv", "w", encoding="utf-8") as fh:
    fh.write(";".join(cols) + "\n")
    for o in rows:
        fh.write(";".join(
            f"{o[c]:.2f}" if isinstance(o[c], float) else str(o[c]) for c in cols
        ) + "\n")

# ---- markdown ----
md = ["# Matrice figure 7 - pics etablis (t >= 185 s)", "",
      "Cibles theorie article (pics, unites vraies) : force 1758-1806 mN (e0.8), "
      "2389-2466 mN (e1.0) ; couple 591-733 uNm (e0.8), 702-876 uNm (e1.0).",
      "Ancrages texte exp : 1824.60 / 2424.67 mN ; 877.01 / 873.41 uNm.", "",
      "| tag | bias | section | eps | pression | disc | P pic (MPa) | F pic (mN) | vs th. | vs exp | ampl F | T pic (uNm) | vs th. | vs exp | ampl T |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for o in rows:
    md.append(
        f"| {o['tag']} | {o['bias']} | {o['section']} | {o['eps']} | {o['pmode']} | {o['disc']} "
        f"| {o['p_peak']:.3f} | {o['F_peak']:.0f} | {o['F_vs_theorie_pct']:+.1f}% | {o['F_vs_exp_pct']:+.1f}% | {o['F_ampl']:.0f} "
        f"| {o['T_peak']:.0f} | {o['T_vs_theorie_pct']:+.1f}% | {o['T_vs_exp_pct']:+.1f}% | {o['T_ampl']:.0f} |"
    )
(SCRATCH / "matrix_summary.md").write_text("\n".join(md), encoding="utf-8")
print("\n".join(md))

# ---- calibration vraie des panneaux force/couple (etiquettes de ticks) ----
img = np.asarray(Image.open(SCRATCH / "p12_img0.png").convert("RGB")).astype(int)
NONWHITE = img.sum(axis=2) < 3 * 225

PANELS = {
    (0.8, "force"):  dict(crop=(120, 51, 704, 350),  ticks=[2200, 2000, 1800, 1600, 1400, 1200], ps=(1200., 2200.)),
    (0.8, "torque"): dict(crop=(120, 351, 704, 647), ticks=[1650, 1100, 550, 0],                 ps=(-250., 1650.)),
    (1.0, "force"):  dict(crop=(899, 51, 1484, 350), ticks=[2800, 2600, 2400, 2200, 2000],      ps=(2000., 2800.)),
    (1.0, "torque"): dict(crop=(899, 351, 1484, 647), ticks=[1440, 960, 480, 0],                ps=(-250., 1800.)),
}


def label_rows(crop, n_expected):
    x0, y0, x1, y1 = crop
    band = NONWHITE[y0:y1, x0 - 48:x0 - 6]
    rr = np.where(band.any(axis=1))[0]
    groups, start, prev = [], rr[0], rr[0]
    for v in rr[1:]:
        if v - prev <= 3:
            prev = v
        else:
            groups.append((start, prev)); start = prev = v
    groups.append((start, prev))
    groups = [g for g in groups if 8 <= (g[1] - g[0]) <= 34]
    assert len(groups) == n_expected, (crop, groups)
    return [0.5 * (a + b) for a, b in groups]


def script_to_true(eps, kind, y_script):
    spec = PANELS[(eps, kind)]
    rows_l = np.asarray(label_rows(spec["crop"], len(spec["ticks"])), float)
    A, B = np.polyfit(rows_l, np.asarray(spec["ticks"], float), 1)
    lo, hi = spec["ps"]
    h = spec["crop"][3] - spec["crop"][1]
    height = max(1.0, float(h - 1))
    row = (hi - np.asarray(y_script)) * height / (hi - lo)
    return A * row + B


# ---- overlay : meilleure config vs conforme article ----
import sys
WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
sys.path.insert(0, str(WORKSPACE))
import validation_figure7 as f7  # noqa: E402
REAL_PDF = WORKSPACE / "Article support" / f7.PDF_FILENAME
_orig = f7.extract_figure7_image
f7.extract_figure7_image = lambda pdf_path=REAL_PDF: _orig(REAL_PDF)
exp_targets = f7.digitize_figure7()
th_targets = f7.digitize_figure7_theory()

BEST = {0.8: "lf_e08_rs_L6", 1.0: "lf_e10_rs_L6"}
CONF = {0.8: "au_e08_rs_L6", 1.0: "arc_upd_e10_rescale_L3p16d10"}

fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.5), sharex=True, constrained_layout=True)
for col, eps in enumerate((0.8, 1.0)):
    for rowi, kind in enumerate(("force", "torque")):
        ax = axes[rowi, col]
        tx, ty = exp_targets[eps][kind]
        ax.scatter(tx, script_to_true(eps, kind, ty), s=4, alpha=0.3,
                   color="#1f77b4" if kind == "force" else "#ff7f0e",
                   label="Exp. numerisee (unites vraies)")
        thx, thy = th_targets[eps][kind]
        ax.plot(thx, script_to_true(eps, kind, thy), "k-", lw=1.2, alpha=0.8,
                label="Theorie article (unites vraies)")
        d = np.load(SCRATCH / f"mx_{BEST[eps]}.npz")
        key = "force" if kind == "force" else "torque"
        ax.plot(d["time"], d[key], color="#d62728", lw=1.6,
                label="Meilleure config (lin/fixed, P recalee, L6)")
        d2 = np.load(SCRATCH / f"mx_{CONF[eps]}.npz")
        ax.plot(d2["time"], d2[key], color="#9467bd", lw=1.4, ls="--",
                label="Conforme article (arctan/updated, P recalee)")
        ax.set_ylabel("Force bloquee (mN)" if kind == "force" else "Couple actionnement (uNm)")
        if rowi == 0:
            ax.set_title(f"eps = {eps}")
        if rowi == 1:
            ax.set_xlabel("Temps (s)")
        ax.grid(alpha=0.3)
        if rowi == 0 and col == 0:
            ax.legend(fontsize=7.5, loc="lower right")
fig.suptitle("Figure 7 : matrice de runs - meilleure config vs configuration conforme article")
fig.savefig(SCRATCH / "fig_matrix_overlay.png", dpi=150)
print("\nfig_matrix_overlay.png ecrit")

# ---- barres des ecarts (runs corr, L3p16, les 4 combos) ----
combos = [("paper_linear", "fixed"), ("uniform_twist", "fixed"),
          ("paper_linear", "updated"), ("uniform_twist", "updated")]
labels = ["lin/fixed", "arctan/fixed", "lin/updated", "arctan/updated"]
fig2, axes2 = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
width = 0.2
for axi, kind in enumerate(("F", "T")):
    ax = axes2[axi]
    for ci, (bias, sec) in enumerate(combos):
        vals = []
        for eps in (0.8, 1.0):
            sel = [o for o in rows
                   if o["bias"] == bias and o["section"] == sec and o["eps"] == eps
                   and o["pmode"] == "corr" and o["disc"] == "L3/p16/dt1.0"]
            vals.append(sel[0][f"{kind}_vs_theorie_pct"] if sel else np.nan)
        ax.bar(np.arange(2) + (ci - 1.5) * width, vals, width, label=labels[ci])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["eps=0.8", "eps=1.0"])
    ax.set_ylabel("Ecart pics vs theorie article (%)")
    ax.set_title("Force" if kind == "F" else "Couple")
    ax.grid(axis="y", alpha=0.3)
axes2[0].legend(fontsize=8)
fig2.suptitle("Ecarts aux pics theorie article - pression corrigee, L3/p16/dt1")
fig2.savefig(SCRATCH / "fig_matrix_deficits.png", dpi=150)
print("fig_matrix_deficits.png ecrit")
