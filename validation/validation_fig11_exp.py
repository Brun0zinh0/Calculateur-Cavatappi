# -*- coding: utf-8 -*-
"""Validation quantitative du mode suspendu contre la figure 11 de l'article EXP.

Article de référence : Liu, Zhang & Liu, « Twisted and coiled tube actuators
driven by hydraulic pressure: Experiment and theory », Thin-Walled Structures
209 (2025) 112957. La figure 11 (p. 9) est la validation TCPA de l'article :

* panneau (a) — fluage sous charge morte de 1 N, P = 0, 600 s ;
* panneau (b) — 10 cycles d'injection 0-3 mL à 5 mL/min (72 s/cycle, 720 s)
  sous 1 N ; l'entrée du modèle de l'article est la pression MESURÉE, tracée
  dans la figure 5(a) (pics 1,21 -> 1,18 MPa), numérisée ici.

Protocole du moteur (item 4.1 du plan de correction, audit 2026-08) :
Base.run_suspended_actuation, load_N = 1.0, eps = 0 (aucune précontrainte dans
le protocole EXP), equilibrate_load_before_pressure = False (le fluage sous le
poids fait partie de l'essai), et les DEUX modes de précontrainte
('elastic_tk_reference' et 'viscoelastic_history') — avec eps = 0 la phase de
précontrainte est vide et les deux modes coïncident exactement ; le script le
vérifie numériquement à chaque exécution.

Géométrie et matériau : ceux de l'article EXP (Tables 1-3) — TCPA 2R_TO0 =
2,00 mm, 2R_TI0 = 0,80 mm, alpha_f0 = 35,69°, R_h0 = 2,66 mm, beta_h0 = 6,87°,
L_T0 = 25,07 mm, nylon d = 0,77 mm ; Maxwell m = 3 (E0 = 6,36, E1 = 20,67,
E2 = 5,98, E3 = 4,75 MPa ; eta1 = 154,57, eta2 = 977,79, eta3 = 11044,83
MPa·s) ; E_radial = 9,1, G12 = 7,2, E_nylon = 3694, G_nylon = 790 MPa.
Seuls E_radial (9,1 vs 8,82), G12 (7,2 vs 7,24) et E_nylon (3694 vs 3690)
diffèrent des défauts alpha V3.

CONVENTION DE DÉFORMATION (point crucial) : la figure 11 démarre à 0 alors que
le poids est déjà accroché — sa référence est la longueur CHARGÉE à t = 0+
(saut élastique inclus dans la référence, invisible sur la figure), et la
normalisation est L_T0 = 25,07 mm. Vérification : beta_h0 = 6,87° et
R_h0 = 2,66 mm donnent un pas d'hélice 2*pi*R_h0*tan(beta_h0) = 2,01 mm =
2R_TO0, c'est-à-dire des spires jointives au repos ; le moteur prédit un
allongement élastique de 25,07 -> ~57,3 mm sous 1 N, cohérent avec les photos
de la figure 5(b) (~55 mm d'hélice chargée). Le script calcule donc
eps_fig11(t) = -(L(t) - L(0+)) / L_T0 à partir de axial_length_mm — convention
indépendante de la version du moteur (free_actuation_percent a changé de
normalisation en phase4-14).

Résultats attendus (moteur 2026.08.21-audit-phase4-14, dt = 1 s, 4 couches,
16 secteurs) :

* fluage : modèle -0,297 @ 72 s / -0,375 @ 600 s, contre -0,283 / -0,345 pour
  la théorie de l'article et -0,260 / -0,340 pour l'expérience (RMSE ~0,023 vs
  théorie, ~0,040 vs expérience) — amplitude et échelle de temps reproduites,
  légère surestimation du fluage (même signe d'écart que la théorie de
  l'article vs son expérience : anisotropie visqueuse supposée égale) ;
* cycles : 1er cycle reproduit (pic ~0,00, amplitude ~0,31) puis le moteur
  sous-prédit la récupération sous pression (amplitude ~0,15 aux cycles 2-10
  contre 0,28-0,287 dans le texte de l'article) et dérive vers le bas — le
  moteur n'a pas de butée d'auto-contact, que l'article invoque comme état de
  contraction maximale.

Usage :

    python validation_fig11_exp.py                 # les deux panneaux, 2 modes
    python validation_fig11_exp.py --quick         # discrétisation réduite
    python validation_fig11_exp.py --panel a       # fluage seulement
    python validation_fig11_exp.py --show          # affiche les figures

Sorties : PNG et métriques JSON dans validation/out/. La numérisation reste
transparente et approximative (± ~0,005 de déformation, ± ~0,02 MPa) : les
RMSE sont des ordres de grandeur pour diagnostiquer les écarts, pas un
substitut aux données brutes. Ancrage de qualité : l'expérience numérisée du
panneau (a) donne -0,260 à 72 s, contre « elongation 0.262 at 72 s » dans le
texte de l'article (+0,8 %).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from io import BytesIO
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
APP_DIR = SCRIPT_DIR.parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import Base  # noqa: E402
from parametres import DEFAULT_SETTINGS, build_config  # noqa: E402


PDF_FILENAME = "Twisted and coiled tube actuators driven by hydraulic pressure Experiment .pdf"
PDF_CANDIDATES = (
    APP_DIR.parent / "Article support" / PDF_FILENAME,
    APP_DIR / PDF_FILENAME,
    SCRIPT_DIR / PDF_FILENAME,
)

OUT_DIR = SCRIPT_DIR / "out"

# ---------------------------------------------------------------------------
# Paramètres de l'article EXP (Tables 1-3)
# ---------------------------------------------------------------------------

L_T0_MM = 25.07          # longueur initiale du TCPA (p. 3)
LOAD_N = 1.0             # charge morte (p. 2)
CREEP_DURATION_S = 600.0
CYCLES_DURATION_S = 720.0

# Ancrages chiffrés du texte de l'article.
TEXT_CREEP_ELONGATION_72S = 0.262        # p. 10 : élongation de fluage à 72 s
TEXT_CYCLE_ACTUATION = {2: 0.28, 5: 0.287, 10: 0.284}   # p. 10 : ε_AF à 3 mL

EXP_SETTINGS = {
    "eps": 0.0,                  # aucune précontrainte dans le protocole EXP
    "rout_mm": 1.00,             # 2R_TO0 = 2,00 mm
    "rin_mm": 0.40,              # 2R_TI0 = 0,80 mm
    "nylon_diameter_mm": 0.77,
    "rho0_mm": 2.66,             # R_h0
    "alpha0_deg": 6.87,          # beta_h0 (angle d'hélice)
    "theta_f_deg": 35.69,        # alpha_f0 (angle de biais en surface)
    "initial_length_mm": L_T0_MM,
    "uncoiled_length_mm": 0.0,
    "E_radius_mpa": 9.1,
    "G12_mpa": 7.2,
    "E_nylon_mpa": 3694.0,
    "G_nylon_mpa": 790.0,
    # Maxwell : les défauts alpha V3 sont déjà la Table 2 de l'article EXP.
}

PRESTRAIN_MODES = ("elastic_tk_reference", "viscoelastic_history")

# ---------------------------------------------------------------------------
# Calibrations des images embarquées (mesurées sur les cadres et graduations)
# ---------------------------------------------------------------------------

# Figure 11 = image unique de la page 9 (1900 x 1653 px).
FIG11_PAGE_INDEX = 8
FIG11A = {
    "frame": (131.0, 907.5, 14.0, 639.0),   # x0, x1, y0(haut), y1(bas) en px
    "xlim": (0.0, 600.0),
    "ylim": (0.0, -0.4),                     # haut, bas
    "xticks": 7,
    "yticks": 5,
    "legend": (420, 900, 30, 180),           # boîte de légende à exclure
}
FIG11B = {
    "frame": (1089.0, 1866.0, 14.0, 639.0),
    "xlim": (0.0, 720.0),
    "ylim": (0.1, -0.4),
    "xticks": 7,
    "yticks": 6,
    "legend": (1350, 1860, 25, 185),
}
# Figure 5 = image unique de la page 5 (1300 x 1988 px), panneau pression.
FIG5_PAGE_INDEX = 4
FIG5P = {
    "xticks_px": [199.0, 377.5, 556.0, 735.0, 913.0, 1091.0, 1270.0],
    "xticks_val": [0.0, 120.0, 240.0, 360.0, 480.0, 600.0, 720.0],
    "yticks_px": [430.5, 526.5, 621.5, 714.0],
    "yticks_val": [1.2, 0.8, 0.4, 0.0],
    "frame_x": (199.0, 1270.0),
    "crop_rows": (390, 764),
    "tick_len": 16,
}

EDGE_MARGIN_PX = 12   # exclusion des graduations mineures le long des cadres


def resolve_pdf_path() -> Path:
    override = os.environ.get("TCPA_EXP_PDF", "")
    candidates = ((Path(override),) if override else ()) + PDF_CANDIDATES
    for candidate in candidates:
        if candidate.exists():
            return candidate
    tried = "\n  - ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        "PDF de l'article EXP introuvable (la variable d'environnement TCPA_EXP_PDF "
        "permet d'indiquer un chemin explicite). Chemins essayés :\n  - " + tried
    )


# ---------------------------------------------------------------------------
# Numérisation
# ---------------------------------------------------------------------------


def page_image(page_index: int, pdf_path: Path | None = None) -> np.ndarray:
    """Image raster unique embarquée dans la page (pypdf page.images)."""
    from PIL import Image
    from pypdf import PdfReader

    reader = PdfReader(str(pdf_path or resolve_pdf_path()))
    page = reader.pages[page_index]
    if len(page.images) < 1:
        raise RuntimeError(f"Aucune image embarquée sur la page {page_index + 1} du PDF.")
    return np.asarray(Image.open(BytesIO(page.images[0].data)).convert("RGB"))


def color_mask(rgb: np.ndarray, name: str) -> np.ndarray:
    r, g, b = [rgb[..., i].astype(float) for i in range(3)]
    if name == "red":
        return (r > 130) & ((r - g) > 60) & ((r - b) > 60)
    if name == "black":
        return (r < 100) & (g < 100) & (b < 100)
    if name == "blue":
        return (b > 130) & ((b - r) > 60) & ((b - g) > 60)
    raise ValueError(f"Masque de couleur inconnu : {name}")


def panel_interior_mask(shape: Tuple[int, ...], spec: dict) -> np.ndarray:
    """Intérieur du cadre, graduations (majeures et mineures) et légende exclues."""
    fx0, fx1, fy0, fy1 = spec["frame"]
    inner = np.zeros(shape[:2], dtype=bool)
    inner[int(fy0) + EDGE_MARGIN_PX:int(fy1) - EDGE_MARGIN_PX + 1,
          int(fx0) + EDGE_MARGIN_PX:int(fx1) - EDGE_MARGIN_PX + 1] = True
    lx0, lx1, ly0, ly1 = spec["legend"]
    inner[ly0:ly1, lx0:lx1] = False
    # bandes des graduations majeures (longueur ~16 px vers l'intérieur)
    for c in np.linspace(fx0, fx1, spec["xticks"]):
        c0, c1 = int(c - 4), int(c + 5)
        inner[int(fy1) - 16:int(fy1), c0:c1] = False
        inner[int(fy0):int(fy0) + 16, c0:c1] = False
    for rr in np.linspace(fy0, fy1, spec["yticks"]):
        r0, r1 = int(rr - 4), int(rr + 5)
        inner[r0:r1, int(fx0):int(fx0) + 16] = False
        inner[r0:r1, int(fx1) - 16:int(fx1)] = False
    return inner


def column_band_curve(mask: np.ndarray, peak_restore: bool = True):
    """Courbe par centre de bande verticale, pics restitués par l'enveloppe.

    Pour chaque colonne de pixels : centre = (rangée min + rangée max)/2.
    Sur un segment monotone raide, le centre correspond à la valeur au milieu
    de la colonne ; au sommet d'un pic, il sous-estime — les maxima locaux de
    l'enveloppe supérieure (rangée min) remplacent alors le centre.
    """
    cols, ctr, top, bot = [], [], [], []
    for c in range(mask.shape[1]):
        rows = np.where(mask[:, c])[0]
        if len(rows) == 0:
            continue
        cols.append(c)
        top.append(rows.min())
        bot.append(rows.max())
        ctr.append(0.5 * (rows.min() + rows.max()))
    cols = np.asarray(cols, dtype=float)
    ctr = np.asarray(ctr)
    top = np.asarray(top, dtype=float)
    bot = np.asarray(bot, dtype=float)
    if peak_restore and len(cols) > 10:
        for i in range(2, len(cols) - 2):
            if top[i] <= top[i - 1] and top[i] <= top[i + 1] and (bot[i] - top[i]) > 8:
                lo, hi = max(0, i - 2), min(len(cols), i + 3)
                if top[i] == top[lo:hi].min():
                    ctr[i] = top[i]
    return cols, ctr


def px_to_data(xs: np.ndarray, ys: np.ndarray, spec: dict):
    fx0, fx1, fy0, fy1 = spec["frame"]
    x = spec["xlim"][0] + (xs - fx0) / (fx1 - fx0) * (spec["xlim"][1] - spec["xlim"][0])
    y = spec["ylim"][0] + (ys - fy0) / (fy1 - fy0) * (spec["ylim"][1] - spec["ylim"][0])
    return x, y


def digitize_fig11(arr9: np.ndarray) -> Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """Numérise les panneaux (a) et (b) : marqueurs rouges (exp), trait noir (théorie)."""
    out: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    for key, spec in (("a", FIG11A), ("b", FIG11B)):
        inner = panel_interior_mask(arr9.shape, spec)
        out[key] = {}
        for curve, color in (("exp", "red"), ("theory", "black")):
            cols, ctr = column_band_curve(color_mask(arr9, color) & inner)
            t, e = px_to_data(cols, ctr, spec)
            out[key][curve] = (t, e)
    return out


def digitize_fig5_pressure(arr5: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Pression mesurée de la figure 5(a) (courbe bleue), en MPa vs s."""
    r0, r1 = FIG5P["crop_rows"]
    x0 = int(FIG5P["frame_x"][0]) + 8
    x1 = int(FIG5P["frame_x"][1]) - 3
    blue = color_mask(arr5, "blue")
    keep = np.zeros_like(blue)
    keep[r0:r1, x0:x1] = True
    # graduations bleues vers l'intérieur (celle de 0.0 MPa est confondue avec
    # la courbe au départ : on la garde)
    for rr in FIG5P["yticks_px"][:-1]:
        keep[int(rr - 5):int(rr + 6), x0:x0 + FIG5P["tick_len"]] = False
        keep[int(rr - 5):int(rr + 6), x1 - FIG5P["tick_len"]:x1] = False
    blue &= keep
    cols, ctr = column_band_curve(blue)
    Ax, Bx = np.polyfit(FIG5P["xticks_px"], FIG5P["xticks_val"], 1)
    Ay, By = np.polyfit(FIG5P["yticks_px"], FIG5P["yticks_val"], 1)
    t = Ax * cols + Bx
    p = np.clip(Ay * ctr + By, 0.0, None)
    # ancrage protocole : la pression part de 0 à t = 0 (t < ~6 s perdu par la
    # marge d'exclusion des graduations)
    if t[0] > 0.0:
        t = np.concatenate(([0.0], t))
        p = np.concatenate(([0.0], p))
    return t, p


# ---------------------------------------------------------------------------
# Modèle
# ---------------------------------------------------------------------------


def make_config(mode: str, dt: float, n_layers: int, n_phi: int):
    settings = dict(DEFAULT_SETTINGS)
    settings.update(EXP_SETTINGS)
    settings.update(
        {
            "dt": dt,
            "n_layers": n_layers,
            "n_phi": n_phi,
            "prestrain_reference_mode": mode,
        }
    )
    return build_config(settings)


def run_protocol(mode: str, t: np.ndarray, p: np.ndarray, dt: float, n_layers: int, n_phi: int):
    """Exécute le mode suspendu et renvoie (t, eps_fig11, L_mm).

    eps_fig11(t) = -(L(t) - L(0+)) / L_T0 : référence = état chargé à t = 0+
    (convention de la figure 11), normalisation par la longueur naturelle L_T0.
    """
    cfg = make_config(mode, dt=dt, n_layers=n_layers, n_phi=n_phi)
    _, arr = Base.run_suspended_actuation(
        cfg,
        load_N=LOAD_N,
        pressure_time=t,
        pressure_MPa=p,
        equilibrate_load_before_pressure=False,
    )
    L = np.asarray(arr["axial_length_mm"], dtype=float)
    tt = np.asarray(arr["time"], dtype=float)
    eps_fig = -(L - L[0]) / L_T0_MM
    return tt, eps_fig, L


def rmse_to_digitized(target_x, target_y, model_x, model_y) -> float:
    m = (np.asarray(target_x) >= np.nanmin(model_x)) & (np.asarray(target_x) <= np.nanmax(model_x))
    if np.count_nonzero(m) < 2:
        return float("nan")
    err = np.interp(np.asarray(target_x)[m], model_x, model_y) - np.asarray(target_y)[m]
    return float(np.sqrt(np.nanmean(err * err)))


def cycle_metrics(t: np.ndarray, e: np.ndarray, period: float = 72.0) -> Dict[str, dict]:
    """Pic, vallée et amplitude (vallée -> pic) par cycle de 72 s."""
    out = {}
    n = int(np.floor(t.max() / period + 1e-9))
    for cyc in range(1, n + 1):
        w = (t >= (cyc - 1) * period) & (t <= cyc * period)
        if np.count_nonzero(w) < 4:
            continue
        out[str(cyc)] = {
            "peak": float(e[w].max()),
            "valley": float(e[w].min()),
            "amplitude": float(e[w].max() - e[w].min()),
        }
    return out


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def check_prestrain_modes(dt: float, n_layers: int, n_phi: int) -> float:
    """Vérifie que les deux modes de précontrainte coïncident (eps = 0).

    Fluage court de 60 s : renvoie max|dL| entre les deux modes (attendu : 0).
    """
    t = np.arange(0.0, 60.0 + 0.5 * dt, dt)
    p = np.zeros_like(t)
    lengths = []
    for mode in PRESTRAIN_MODES:
        _, _, L = run_protocol(mode, t, p, dt=dt, n_layers=n_layers, n_phi=n_phi)
        lengths.append(L)
    return float(np.max(np.abs(lengths[0] - lengths[1])))


def run_validation(panel: str, dt_creep: float, dt_cycles: float, n_layers: int,
                   n_phi: int, modes, show: bool = False) -> dict:
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = resolve_pdf_path()
    print(f"PDF : {pdf_path}")
    print(f"Moteur : Base {Base.MODEL_VERSION} — dt fluage {dt_creep:g} s, "
          f"dt cycles {dt_cycles:g} s, {n_layers} couches, {n_phi} secteurs")

    dig = digitize_fig11(page_image(FIG11_PAGE_INDEX, pdf_path))
    p_t, p_mpa = digitize_fig5_pressure(page_image(FIG5_PAGE_INDEX, pdf_path))

    # qualité de la numérisation : ancrage sur le texte de l'article
    a_exp_t, a_exp_e = dig["a"]["exp"]
    exp72 = float(np.interp(72.0, a_exp_t, a_exp_e))
    print(f"Numérisation (a) : exp @ 72 s = {exp72:+.4f} "
          f"(texte article : -{TEXT_CREEP_ELONGATION_72S}) ; "
          f"pression fig 5 : pics {p_mpa.max():.2f} MPa, {len(p_t)} points")

    metrics: dict = {
        "model_version": Base.MODEL_VERSION,
        "discretization": {"dt_creep": dt_creep, "dt_cycles": dt_cycles,
                           "n_layers": n_layers, "n_phi": n_phi},
        "digitization": {
            "a_exp_at_72s": exp72,
            "text_creep_elongation_72s": TEXT_CREEP_ELONGATION_72S,
            "fig5_pressure_peak_MPa": float(p_mpa.max()),
        },
        "panels": {},
    }

    # les deux modes de précontrainte doivent coïncider avec eps = 0
    dmax = check_prestrain_modes(dt=2.0, n_layers=max(2, n_layers // 2), n_phi=max(8, n_phi // 2))
    metrics["prestrain_modes_max_dL_mm"] = dmax
    print(f"Modes de précontrainte (eps = 0) : max|dL| = {dmax:.3e} mm "
          f"({'identiques' if dmax < 1e-9 else 'ATTENTION : différents'})")

    runs: dict = {}
    if panel in ("a", "both"):
        t = np.arange(0.0, CREEP_DURATION_S + 0.5 * dt_creep, dt_creep)
        p = np.zeros_like(t)
        runs["a"] = {mode: run_protocol(mode, t, p, dt_creep, n_layers, n_phi) for mode in modes}
    if panel in ("b", "both"):
        tc = np.arange(0.0, CYCLES_DURATION_S + 0.5 * dt_cycles, dt_cycles)
        pc = np.clip(np.interp(tc, p_t, p_mpa), 0.0, None)
        runs["b"] = {mode: run_protocol(mode, tc, pc, dt_cycles, n_layers, n_phi) for mode in modes}
        runs["b_pressure"] = (tc, pc)

    n_rows = len([k for k in ("a", "b") if k in runs])
    fig, axes = plt.subplots(n_rows, 1, figsize=(11.5, 4.6 * n_rows), squeeze=False)
    axes = axes[:, 0]
    fig.suptitle(
        "Validation figure 11 (EXP) : mode suspendu alpha V3, charge 1 N, eps = 0\n"
        f"Base {Base.MODEL_VERSION} — référence = état chargé t = 0+, normalisation L_T0 = {L_T0_MM} mm"
    )
    row = 0

    if "a" in runs:
        ax = axes[row]
        row += 1
        ax.plot(*dig["a"]["theory"], "k-", lw=1.1, label="Théorie article (numérisée)")
        ax.plot(*dig["a"]["exp"], "r.", ms=3.2, label="Expérience (numérisée)")
        panel_metrics: dict = {"modes": {}}
        for mode, style in zip(modes, ("-", "--")):
            tt, ee, LL = runs["a"][mode]
            ax.plot(tt, ee, style, color="#1f4fd8", lw=1.6, label=f"Alpha V3 ({mode})")
            key_points = {str(int(tq)): float(np.interp(tq, tt, ee))
                          for tq in (18.0, 72.0, 300.0, 600.0)}
            panel_metrics["modes"][mode] = {
                "L_loaded_0plus_mm": float(LL[0]),
                "elastic_elongation_strain": float((LL[0] - L_T0_MM) / L_T0_MM),
                "key_points": key_points,
                "rmse_vs_experiment": rmse_to_digitized(*dig["a"]["exp"], tt, ee),
                "rmse_vs_paper_theory": rmse_to_digitized(*dig["a"]["theory"], tt, ee),
            }
        panel_metrics["digitized"] = {
            "exp": {str(int(tq)): float(np.interp(tq, *dig["a"]["exp"]))
                    for tq in (72.0, 300.0)},
            "theory": {str(int(tq)): float(np.interp(tq, *dig["a"]["theory"]))
                       for tq in (72.0, 300.0)},
            "exp_end": [float(dig["a"]["exp"][0][-1]), float(dig["a"]["exp"][1][-1])],
            "theory_end": [float(dig["a"]["theory"][0][-1]), float(dig["a"]["theory"][1][-1])],
        }
        metrics["panels"]["a"] = panel_metrics
        ax.set_title("(a) Fluage sous 1 N, P = 0")
        ax.set_xlabel("Temps (s)")
        ax.set_ylabel("Déformation d'actionnement (-)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    if "b" in runs:
        ax = axes[row]
        ax.plot(*dig["b"]["theory"], "k-", lw=1.0, label="Théorie article (numérisée)")
        ax.plot(*dig["b"]["exp"], "r.", ms=2.5, label="Expérience (numérisée)")
        panel_metrics = {"modes": {}}
        for mode, style in zip(modes, ("-", "--")):
            tb, eb, Lb = runs["b"][mode]
            ax.plot(tb, eb, style, color="#1f4fd8", lw=1.3, label=f"Alpha V3 ({mode})")
            cyc = cycle_metrics(tb, eb)
            panel_metrics["modes"][mode] = {
                "L_loaded_0plus_mm": float(Lb[0]),
                "rmse_vs_experiment": rmse_to_digitized(*dig["b"]["exp"], tb, eb),
                "rmse_vs_paper_theory": rmse_to_digitized(*dig["b"]["theory"], tb, eb),
                "cycles": cyc,
                "text_cycle_actuation": TEXT_CYCLE_ACTUATION,
            }
        te, ee_ = dig["b"]["exp"]
        panel_metrics["digitized_exp_cycles"] = cycle_metrics(te, ee_)
        metrics["panels"]["b"] = panel_metrics
        tp, pp = runs["b_pressure"]
        ax2 = ax.twinx()
        ax2.plot(tp, pp, color="#2ca02c", lw=0.7, alpha=0.6)
        ax2.set_ylabel("Pression (MPa)", color="#2ca02c")
        ax.set_title("(b) Cycles 0-3 mL à 5 mL/min sous 1 N (pression mesurée, fig. 5)")
        ax.set_xlabel("Temps (s)")
        ax.set_ylabel("Déformation d'actionnement (-)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower left")

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig_path = OUT_DIR / f"validation_fig11_exp_{panel}_L{n_layers}p{n_phi}.png"
    fig.savefig(fig_path, dpi=160)

    # ---- synthèse console
    print()
    if "a" in metrics["panels"]:
        pm = metrics["panels"]["a"]
        print("Panneau (a) — fluage 1 N :")
        for mode in modes:
            mm = pm["modes"][mode]
            kp = mm["key_points"]
            print(f"  [{mode}] L(0+) = {mm['L_loaded_0plus_mm']:.2f} mm "
                  f"(saut élastique {mm['elastic_elongation_strain']:+.3f}, hors figure)")
            print(f"    modèle  : {kp['72']:+.3f} @ 72 s ; {kp['300']:+.3f} @ 300 s ; "
                  f"{kp['600']:+.3f} @ 600 s")
            print(f"    RMSE    : {mm['rmse_vs_paper_theory']:.4f} vs théorie article ; "
                  f"{mm['rmse_vs_experiment']:.4f} vs expérience")
        dd = pm["digitized"]
        print(f"    article : théorie {dd['theory']['72']:+.3f} @ 72 s, "
              f"{dd['theory_end'][1]:+.3f} @ {dd['theory_end'][0]:.0f} s ; "
              f"exp {dd['exp']['72']:+.3f} @ 72 s, "
              f"{dd['exp_end'][1]:+.3f} @ {dd['exp_end'][0]:.0f} s")
    if "b" in metrics["panels"]:
        pm = metrics["panels"]["b"]
        print("Panneau (b) — cycles 0-3 mL :")
        for mode in modes:
            mm = pm["modes"][mode]
            print(f"  [{mode}] RMSE : {mm['rmse_vs_paper_theory']:.4f} vs théorie ; "
                  f"{mm['rmse_vs_experiment']:.4f} vs expérience")
            for cyc, target in TEXT_CYCLE_ACTUATION.items():
                cm = mm["cycles"].get(str(cyc))
                if cm:
                    print(f"    cycle {cyc:2d} : amplitude modèle {cm['amplitude']:.3f} "
                          f"vs texte article {target:.3f} "
                          f"({100.0 * (cm['amplitude'] / target - 1.0):+.0f} %)")
    print()
    print("Rappel : moteur sans butée d'auto-contact — l'article invoque l'auto-contact")
    print("comme état de contraction maximale (spires jointives au repos, pas = 2R_TO0).")

    metrics_path = OUT_DIR / f"validation_fig11_exp_{panel}_L{n_layers}p{n_phi}.json"
    with open(metrics_path, "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, ensure_ascii=False, indent=2)
    print(f"Figure : {fig_path}\nMétriques : {metrics_path}")

    if show:
        plt.show()
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validation du mode suspendu contre la figure 11 de l'article EXP."
    )
    parser.add_argument("--panel", choices=("a", "b", "both"), default="both",
                        help="a = fluage, b = cycles, both = les deux (défaut)")
    parser.add_argument("--quick", action="store_true",
                        help="discrétisation réduite (2 couches, 8 secteurs, dt 2/1 s)")
    parser.add_argument("--layers", type=int, default=None,
                        help="nombre de couches (défaut : 4, ou 2 avec --quick)")
    parser.add_argument("--mode", choices=PRESTRAIN_MODES + ("both",), default="both",
                        help="mode(s) de précontrainte simulé(s) (défaut : les deux)")
    parser.add_argument("--show", action="store_true", help="affiche les figures")
    args = parser.parse_args()

    n_layers = args.layers if args.layers is not None else (2 if args.quick else 4)
    n_phi = 8 if args.quick else 16
    dt_creep = 2.0 if args.quick else 1.0
    dt_cycles = 1.0 if args.quick else 0.5
    modes = PRESTRAIN_MODES if args.mode == "both" else (args.mode,)
    run_validation(args.panel, dt_creep, dt_cycles, n_layers, n_phi, modes, show=args.show)


if __name__ == "__main__":
    main()
