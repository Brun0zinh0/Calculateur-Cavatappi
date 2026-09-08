# -*- coding: utf-8 -*-
"""Validation quantitative contre la figure 7 de l'article « Blocked actuation ».

Version alpha V3 — portée depuis le script historique de la racine avec les
corrections issues de l'audit du 17-18/08/2026 (plan de correction, item 0.1) :

1. moteur : importe le Base et le parametres d'alpha V3 (le script historique
   importait le moteur du dossier « livrable », plus ancien) ;
2. chemin du PDF : résolu vers « Article support/ » (le chemin historique,
   cassé, supposait le PDF à la racine de l'espace de travail) ;
3. calibration : la conversion pixel->pression est corrigée par un ajustement
   sur les positions réelles des étiquettes d'axe (les ticks sont en retrait
   des bords du cadrage ; l'ancienne conversion sous-lisait les pics :
   1,32-1,37 MPa au lieu de 1,41-1,43) ;
4. recalage : en mode « rescale » (défaut), les pics de la pression numérisée
   sont recalés sur la valeur du texte de l'article (1,41-1,43 MPa, cible
   1,42) — sans ce recalage les pics restent limités par la résolution de la
   numérisation.

Usage :

    python validation_figure7.py                # run complet (6 couches)
    python validation_figure7.py --quick        # discrétisation réduite
    python validation_figure7.py --baseline     # (ré)génère la baseline JSON
    python validation_figure7.py --show         # affiche les figures
    python validation_figure7.py --pmode raw    # pression numérisée d'origine

Sorties : figures PNG et métriques JSON dans validation/out/. La baseline de
non-régression est lue/écrite dans validation/baseline_figure7.json et
consommée par test_scientifique.py (test_validation_figure7_non_regression).

La numérisation reste volontairement transparente et approximative : les RMSE
sont des ordres de grandeur pour diagnostiquer les écarts, pas un substitut
aux données expérimentales brutes.
"""

from __future__ import annotations

import argparse
import json
from io import BytesIO
from pathlib import Path
import sys
from typing import Dict, Tuple

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
APP_DIR = SCRIPT_DIR.parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import Base  # noqa: E402
from parametres import DEFAULT_SETTINGS, build_config, cycle_period_seconds  # noqa: E402


PDF_FILENAME = "Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure.pdf"
PDF_CANDIDATES = (
    APP_DIR.parent / "Article support" / PDF_FILENAME,
    APP_DIR / PDF_FILENAME,
    SCRIPT_DIR / PDF_FILENAME,
)

OUT_DIR = SCRIPT_DIR / "out"
BASELINE_PATH = SCRIPT_DIR / "baseline_figure7.json"

PAPER_T_MIN = 180.0
PAPER_T_MAX = 380.0
PAPER_PMAX_MPA = 1.50
PAPER_FLOW_RATE_ML_MIN = 10.0
PAPER_VOLUME_ML = 1.50

# Pics de pression du texte de l'article (p. 10 : « peak pressure 1.41-1.43 MPa »).
TEXT_PRESSURE_PEAK_MPA = 1.42
# Ancrages chiffrés du texte (valeurs expérimentales, 11e cycle).
TEXT_ANCHORS = {
    0.8: {"force_peak_mN": 1824.60, "torque_peak_microNm": 877.01},
    1.0: {"force_peak_mN": 2424.67, "torque_peak_microNm": 873.41},
}

# Fenêtre « établie » pour la détection des pics (écarte le transitoire de
# jonction préfixe généré -> pression numérisée à t = 180 s).
ESTABLISHED_T0 = 185.0


def resolve_pdf_path() -> Path:
    import os

    override = os.environ.get("TCPA_FIGURE7_PDF", "")
    candidates = ((Path(override),) if override else ()) + PDF_CANDIDATES
    for candidate in candidates:
        if candidate.exists():
            return candidate
    tried = "\n  - ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        "PDF de l'article « Blocked actuation » introuvable (la variable d'environnement "
        "TCPA_FIGURE7_PDF permet d'indiquer un chemin explicite). Chemins essayés :\n  - " + tried
    )


# --------------------------------------------------------------------------
# Numérisation de la figure 7 (image embarquée page 12 du PDF)
# --------------------------------------------------------------------------

PANEL_SPECS = {
    0.8: {
        "xlim": (180.0, 380.0),
        "force": {
            "crop": (120, 51, 704, 350),
            "ylim": (1200.0, 2200.0),
            "valid_ylim": (1200.0, 1900.0),
            "color": "blue",
        },
        "torque": {
            "crop": (120, 351, 704, 647),
            "ylim": (-250.0, 1650.0),
            "valid_ylim": (-250.0, 1050.0),
            "color": "orange",
        },
        "pressure": {"crop": (120, 648, 704, 946), "ylim": (0.0, 1.5), "color": "green"},
    },
    1.0: {
        "xlim": (180.0, 380.0),
        "force": {
            "crop": (899, 51, 1484, 350),
            "ylim": (2000.0, 2800.0),
            "valid_ylim": (1950.0, 2550.0),
            "color": "blue",
        },
        "torque": {
            "crop": (899, 351, 1484, 647),
            "ylim": (-250.0, 1800.0),
            "valid_ylim": (-250.0, 1100.0),
            "color": "orange",
        },
        "pressure": {"crop": (899, 648, 1484, 946), "ylim": (0.0, 1.5), "color": "green"},
    },
}

# Calibration vraie des axes de pression : positions (en lignes de pixels du
# cadrage) des étiquettes 1.5 / 1.0 / 0.5 / 0.0, mesurées sur l'image embarquée
# (audit, ce_fig7_measure2.py). L'ancienne conversion supposait 1.5 et 0.0 aux
# bords du cadrage, d'où des pics sous-lus d'environ 0,05-0,09 MPa.
PRESSURE_TICK_FITS = {
    0.8: ([21.0, 107.5, 195.0, 282.0], [1.5, 1.0, 0.5, 0.0], (120, 648, 704, 946)),
    1.0: ([15.0, 104.0, 192.5, 282.0], [1.5, 1.0, 0.5, 0.0], (899, 648, 1484, 946)),
}


def extract_figure7_image(pdf_path: Path | None = None):
    """Retourne l'image raster de la figure 7 embarquée dans le PDF."""
    from PIL import Image
    from pypdf import PdfReader

    reader = PdfReader(str(pdf_path or resolve_pdf_path()))
    page = reader.pages[11]
    if len(page.images) < 1:
        raise RuntimeError("Image de la figure 7 introuvable sur la page 12 du PDF.")
    image = Image.open(BytesIO(page.images[0].data))
    return image.convert("RGB")


def color_mask(rgb: np.ndarray, color: str) -> np.ndarray:
    r = rgb[..., 0].astype(float)
    g = rgb[..., 1].astype(float)
    b = rgb[..., 2].astype(float)
    if color == "blue":
        return (b > 135) & (g > 80) & (r < 135) & ((b - r) > 55)
    if color == "orange":
        return (r > 180) & (g > 60) & (g < 185) & (b < 125) & ((r - b) > 80)
    if color == "green":
        return (g > 135) & (r < 120) & (b < 150) & ((g - r) > 45)
    if color == "black":
        return (r < 105) & (g < 105) & (b < 105)
    raise ValueError(f"Masque de couleur inconnu : {color}")


def pixels_to_data(
    xs: np.ndarray,
    ys: np.ndarray,
    crop: Tuple[int, int, int, int],
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
) -> Tuple[np.ndarray, np.ndarray]:
    x0, y0, x1, y1 = crop
    width = max(1.0, float(x1 - x0 - 1))
    height = max(1.0, float(y1 - y0 - 1))
    x = xlim[0] + xs / width * (xlim[1] - xlim[0])
    y = ylim[1] - ys / height * (ylim[1] - ylim[0])
    return x, y


def bin_curve(
    x: np.ndarray,
    y: np.ndarray,
    bin_width: float = 0.35,
    statistic: str = "median",
) -> Tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if x.size == 0:
        return x, y
    bins = np.arange(np.nanmin(x), np.nanmax(x) + bin_width, bin_width)
    centers = []
    values = []
    for left, right in zip(bins[:-1], bins[1:]):
        mask = (x >= left) & (x < right)
        if np.count_nonzero(mask) >= 2:
            centers.append(0.5 * (left + right))
            if statistic == "median":
                values.append(float(np.nanmedian(y[mask])))
            elif statistic == "p90":
                values.append(float(np.nanpercentile(y[mask], 90.0)))
            elif statistic == "p99":
                values.append(float(np.nanpercentile(y[mask], 99.0)))
            else:
                raise ValueError("statistic doit être 'median', 'p90' ou 'p99'.")
    return np.asarray(centers), np.asarray(values)


def digitize_figure7(image=None) -> Dict[float, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """Numérise les marqueurs expérimentaux colorés (force, couple, pression)."""
    if image is None:
        image = extract_figure7_image()
    arr = np.asarray(image)
    digitized: Dict[float, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    for eps, panel in PANEL_SPECS.items():
        digitized[eps] = {}
        xlim = panel["xlim"]
        for key in ("force", "torque", "pressure"):
            spec = panel[key]
            x0, y0, x1, y1 = spec["crop"]
            crop_rgb = arr[y0:y1, x0:x1]
            mask = color_mask(crop_rgb, spec["color"])
            ys, xs = np.where(mask)
            width = x1 - x0
            height = y1 - y0
            inside = (xs > 8) & (xs < width - 8) & (ys > 4) & (ys < height - 4)
            xs = xs[inside]
            ys = ys[inside]
            x, y = pixels_to_data(xs, ys, spec["crop"], xlim, spec["ylim"])
            if "valid_ylim" in spec:
                ylo, yhi = spec["valid_ylim"]
                valid_value = (y >= ylo) & (y <= yhi)
                x = x[valid_value]
                y = y[valid_value]
            if key == "pressure":
                bx, by = bin_curve(x, y, bin_width=0.20, statistic="p99")
            else:
                bx, by = bin_curve(x, y, bin_width=0.45)
            digitized[eps][key] = (bx, by)
    return digitized


def digitize_figure7_theory(image=None) -> Dict[float, Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """Numérise les courbes théoriques noires, séparées des marqueurs.

    Attention (audit) : aux sommets des pics de couple, les marqueurs
    expérimentaux recouvrent le trait noir — la cible « théorie » est donc
    tronquée aux pics et son maximum numérisé sous-estime le vrai sommet.
    """
    if image is None:
        image = extract_figure7_image()
    arr = np.asarray(image)
    digitized: Dict[float, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}
    for eps, panel in PANEL_SPECS.items():
        digitized[eps] = {}
        for key in ("force", "torque"):
            spec = panel[key]
            x0, y0, x1, y1 = spec["crop"]
            crop_rgb = arr[y0:y1, x0:x1]
            ys, xs = np.where(color_mask(crop_rgb, "black"))
            inside = (xs > 8) & (xs < (x1 - x0) - 8) & (ys > 4) & (ys < (y1 - y0) - 4)
            x, y = pixels_to_data(xs[inside], ys[inside], spec["crop"], panel["xlim"], spec["ylim"])
            if key == "force":
                ylo, yhi = ((1380.0, 1900.0) if eps == 0.8 else (1950.0, 2500.0))
            else:
                ylo, yhi = (-180.0, 900.0)
            valid = (y >= ylo) & (y <= yhi)
            digitized[eps][key] = bin_curve(x[valid], y[valid], bin_width=0.45)
    return digitized


def script_to_true_pressure(eps: float, y_script: np.ndarray) -> np.ndarray:
    """Convertit la pression issue du mappage naïf vers la calibration vraie."""
    rows_t, vals_t, crop = PRESSURE_TICK_FITS[eps]
    A, B = np.polyfit(rows_t, vals_t, 1)
    height = max(1.0, float(crop[3] - crop[1] - 1))
    row = (PANEL_SPECS[eps]["pressure"]["ylim"][1] - np.asarray(y_script)) * height / 1.5
    return A * row + B


def build_targets(pmode: str = "rescale"):
    """Cibles numérisées, avec pression corrigée puis recalée selon pmode.

    pmode : "raw" (mappage naïf historique), "corr" (calibration des ticks),
    "rescale" (corr + pics recalés sur la valeur du texte, 1,42 MPa).
    """
    if pmode not in ("raw", "corr", "rescale"):
        raise ValueError("pmode doit être 'raw', 'corr' ou 'rescale'.")
    targets = digitize_figure7()
    info: Dict[float, dict] = {}
    if pmode == "raw":
        return targets, info
    out = {e: dict(targets[e]) for e in targets}
    for eps in (0.8, 1.0):
        tp, pp = targets[eps]["pressure"]
        pp_true = np.clip(script_to_true_pressure(eps, pp), 0.0, None)
        info[eps] = {}
        if pmode == "rescale":
            ext = cycle_extrema(tp, pp_true, t0=ESTABLISHED_T0)
            scale = TEXT_PRESSURE_PEAK_MPA / ext["peak_mean"]
            pp_true = np.clip(pp_true * scale, 0.0, None)
            info[eps]["scale"] = scale
        out[eps]["pressure"] = (tp, pp_true)
        ext2 = cycle_extrema(tp, pp_true, t0=ESTABLISHED_T0)
        info[eps]["p_peak_mean"] = ext2["peak_mean"]
        info[eps]["p_peak_max"] = ext2["peak_max"]
    return out, info


# --------------------------------------------------------------------------
# Modèle
# --------------------------------------------------------------------------


def make_paper_config(eps: float, dt: float = 1.0, n_layers: int = 6, n_phi: int = 16):
    n_cycles = int(np.ceil(PAPER_T_MAX / (2.0 * 60.0 * PAPER_VOLUME_ML / PAPER_FLOW_RATE_ML_MIN))) + 1
    settings = dict(DEFAULT_SETTINGS)
    settings.update(
        {
            "eps": eps,
            "n_cycles": n_cycles,
            "use_fixed_duration": False,
            "p_max_mpa": PAPER_PMAX_MPA,
            "dt": dt,
            "n_layers": n_layers,
            "n_phi": n_phi,
            "pre_steps": 24,
            "flow_rate_mL_min": PAPER_FLOW_RATE_ML_MIN,
            "volume_mL": PAPER_VOLUME_ML,
            "nonlinear_pressure": True,
        }
    )
    return build_config(settings)


def digitized_pressure_history(
    eps: float,
    targets,
    prefix_dt: float = 1.0,
    figure_dt: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray]:
    """Préfixe d'entraînement généré (0-180 s) puis pression numérisée (180-380 s).

    Le dt de la config est ignoré quand une histoire de pression est fournie :
    le pas de temps effectif du solveur est la grille construite ici.
    """
    period = 2.0 * 60.0 * PAPER_VOLUME_ML / PAPER_FLOW_RATE_ML_MIN
    prefix_cycles = int(np.ceil(PAPER_T_MIN / period)) + 1
    prefix_t, prefix_p = Base.cyclic_pressure_history(
        n_cycles=prefix_cycles,
        Pmax=PAPER_PMAX_MPA,
        flow_rate_mL_min=PAPER_FLOW_RATE_ML_MIN,
        volume_mL=PAPER_VOLUME_ML,
        dt=prefix_dt,
        nonlinear=True,
    )
    keep = prefix_t < PAPER_T_MIN
    prefix_t = prefix_t[keep]
    prefix_p = prefix_p[keep]

    target_t, target_p = targets[eps]["pressure"]
    figure_t = np.arange(PAPER_T_MIN, PAPER_T_MAX + 0.5 * figure_dt, figure_dt)
    figure_p = np.interp(figure_t, target_t, target_p)
    figure_p = np.clip(figure_p, 0.0, PAPER_PMAX_MPA)
    return np.concatenate([prefix_t, figure_t]), np.concatenate([prefix_p, figure_p])


def window_model_data(data: dict) -> dict:
    time = np.asarray(data["time"], dtype=float)
    mask = (time >= PAPER_T_MIN) & (time <= PAPER_T_MAX)
    return {key: np.asarray(value)[mask] for key, value in data.items() if np.asarray(value).shape == time.shape}


def cycle_extrema(t, v, t0: float = ESTABLISHED_T0, min_gap: float = 8.0) -> dict:
    """Pics et vallées par cycle sur la fenêtre établie (> t0)."""
    t = np.asarray(t, float)
    v = np.asarray(v, float)
    m = t >= t0
    t, v = t[m], v[m]
    med = np.median(v)

    def _find(sig, ref):
        out = []
        for i in range(1, len(sig) - 1):
            if sig[i] >= sig[i - 1] and sig[i] > sig[i + 1] and sig[i] > ref:
                if not out or t[i] - out[-1][0] > min_gap:
                    out.append((t[i], sig[i]))
                elif sig[i] > out[-1][1]:
                    out[-1] = (t[i], sig[i])
        return out

    peaks = _find(v, med)
    valleys = [(tt, -vv) for tt, vv in _find(-v, -med)]
    pv = np.asarray([p[1] for p in peaks]) if peaks else np.asarray([np.nan])
    vv = np.asarray([p[1] for p in valleys]) if valleys else np.asarray([np.nan])
    return dict(
        n_peaks=int(len(peaks)),
        peak_mean=float(np.nanmean(pv)),
        peak_min=float(np.nanmin(pv)),
        peak_max=float(np.nanmax(pv)),
        valley_mean=float(np.nanmean(vv)),
        amplitude=float(np.nanmean(pv) - np.nanmean(vv)),
        sig_min=float(v.min()),
        sig_max=float(v.max()),
    )


def run_case(eps: float, targets, dt: float = 1.0, n_layers: int = 6, n_phi: int = 16):
    """Exécute le protocole figure 7 et retourne (config, fenêtre, extrema)."""
    config = make_paper_config(eps=eps, dt=dt, n_layers=n_layers, n_phi=n_phi)
    pressure_time, pressure_mpa = digitized_pressure_history(
        eps, targets, prefix_dt=dt, figure_dt=0.5 * dt
    )
    _, data = Base.run_blocked_actuation(config, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    win = window_model_data(data)
    extrema = {
        "force": cycle_extrema(win["time"], win["force_total_mN"]),
        "torque": cycle_extrema(win["time"], win["torque_act_microNm"]),
        "pressure": cycle_extrema(win["time"], win["pressure_MPa"]),
    }
    return config, win, extrema


def rmse_to_digitized(target_x, target_y, model_x, model_y) -> float:
    if len(target_x) < 2 or len(model_x) < 2:
        return float("nan")
    mask = (target_x >= np.nanmin(model_x)) & (target_x <= np.nanmax(model_x))
    if np.count_nonzero(mask) < 2:
        return float("nan")
    interp = np.interp(target_x[mask], model_x, model_y)
    err = interp - target_y[mask]
    return float(np.sqrt(np.nanmean(err * err)))


# --------------------------------------------------------------------------
# Baseline de non-régression
# --------------------------------------------------------------------------

BASELINE_CONFIGS = {
    "reduced": {"dt": 1.0, "n_layers": 2, "n_phi": 8},
    "production": {"dt": 1.0, "n_layers": 6, "n_phi": 16},
}


def compute_baseline_entry(targets, disc: dict) -> dict:
    entry = {}
    for eps in (0.8, 1.0):
        _, _, extrema = run_case(eps, targets, **disc)
        entry[str(eps)] = {
            "force_peak_mean_mN": round(extrema["force"]["peak_mean"], 2),
            "force_valley_mean_mN": round(extrema["force"]["valley_mean"], 2),
            "torque_peak_mean_microNm": round(extrema["torque"]["peak_mean"], 2),
            "torque_valley_mean_microNm": round(extrema["torque"]["valley_mean"], 2),
            "pressure_peak_mean_MPa": round(extrema["pressure"]["peak_mean"], 4),
        }
    return entry


def write_baseline(pmode: str = "rescale") -> dict:
    targets, _ = build_targets(pmode)
    payload = {
        "model_version": Base.MODEL_VERSION,
        "pmode": pmode,
        "established_t0_s": ESTABLISHED_T0,
        "text_anchors": {str(k): v for k, v in TEXT_ANCHORS.items()},
        "configs": {},
    }
    for name, disc in BASELINE_CONFIGS.items():
        payload["configs"][name] = {"discretization": disc, **compute_baseline_entry(targets, disc)}
    BASELINE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(BASELINE_PATH, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print(f"Baseline écrite : {BASELINE_PATH}")
    return payload


def load_baseline() -> dict:
    with open(BASELINE_PATH, encoding="utf-8") as fh:
        return json.load(fh)


# --------------------------------------------------------------------------
# Validation graphique complète
# --------------------------------------------------------------------------


def plot_validation(pmode: str = "rescale", dt: float = 1.0, n_layers: int = 6, n_phi: int = 16, show: bool = False):
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    experimental_targets, pinfo = build_targets(pmode)
    theory_targets = digitize_figure7_theory()
    runs = {eps: run_case(eps, experimental_targets, dt=dt, n_layers=n_layers, n_phi=n_phi) for eps in (0.8, 1.0)}

    fig, axes = plt.subplots(3, 2, figsize=(12.0, 8.6), sharex="col", constrained_layout=True)
    fig.suptitle(
        "Validation figure 7 : cibles numérisées vs moteur alpha V3\n"
        f"volume {PAPER_VOLUME_ML:.2f} mL, débit {PAPER_FLOW_RATE_ML_MIN:.1f} mL/min, "
        f"Pmax {PAPER_PMAX_MPA:.2f} MPa — pression {pmode}, "
        f"dt {dt:g} s, {n_layers} couches, {n_phi} secteurs\n"
        f"Base {Base.MODEL_VERSION}"
    )

    row_defs = [
        ("force", "force_total_mN", "Force bloquée (mN)", "#1f77b4"),
        ("torque", "torque_act_microNm", "Couple d'actionnement (µN·m)", "#ff7f0e"),
        ("pressure", "pressure_MPa", "Pression (MPa)", "#2ca02c"),
    ]

    metrics = {}
    for col, eps in enumerate((0.8, 1.0)):
        _, win, extrema = runs[eps]
        t = win["time"]
        metrics[eps] = {"extrema": extrema, "pressure_info": pinfo.get(eps, {})}
        for row, (target_key, model_key, ylabel, color) in enumerate(row_defs):
            ax = axes[row, col]
            tx, ty = experimental_targets[eps][target_key]
            ax.scatter(tx, ty, s=5.0, alpha=0.35, color=color, label="Expérience numérisée")
            if target_key in theory_targets[eps]:
                theory_t, theory_y = theory_targets[eps][target_key]
                ax.plot(theory_t, theory_y, color="black", lw=1.3, alpha=0.75, label="Théorie article numérisée")
            else:
                theory_t, theory_y = np.array([]), np.array([])
            ax.plot(t, win[model_key], color="#d62728", lw=1.7, label="Alpha V3")
            if row == 0:
                ax.set_title(f"eps = {eps:.1f}")
                ax.legend(loc="best", fontsize=8)
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)
            if row == 2:
                ax.set_xlabel("Temps depuis le début de la pression (s)")
            metrics[eps][target_key] = {
                "rmse_experiment": rmse_to_digitized(tx, ty, t, win[model_key]),
                "rmse_theory": rmse_to_digitized(theory_t, theory_y, t, win[model_key]),
            }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = OUT_DIR / f"validation_figure7_{pmode}_L{n_layers}p{n_phi}.png"
    fig.savefig(fig_path, dpi=160)

    print("Validation figure 7 (moteur alpha V3) terminée")
    for eps in (0.8, 1.0):
        config, win, extrema = runs[eps]
        anchors = TEXT_ANCHORS[eps]
        f_peak = extrema["force"]["peak_mean"]
        t_peak = extrema["torque"]["peak_mean"]
        print(f"eps={eps:.1f}  (période de cycle {cycle_period_seconds(config):.1f} s)")
        print(f"  pression modèle : pics {extrema['pressure']['peak_mean']:.3f} MPa "
              f"(texte article : 1.41-1.43)")
        print(f"  force : pics établis {f_peak:.1f} mN "
              f"[{extrema['force']['peak_min']:.1f}..{extrema['force']['peak_max']:.1f}], "
              f"vallées {extrema['force']['valley_mean']:.1f} — "
              f"ancrage texte {anchors['force_peak_mN']:.1f} mN "
              f"({100.0 * (f_peak / anchors['force_peak_mN'] - 1.0):+.1f} %)")
        print(f"  couple : pics établis {t_peak:.1f} µN·m "
              f"[{extrema['torque']['peak_min']:.1f}..{extrema['torque']['peak_max']:.1f}], "
              f"vallées {extrema['torque']['valley_mean']:.1f} — "
              f"ancrage texte {anchors['torque_peak_microNm']:.1f} µN·m "
              f"({100.0 * (t_peak / anchors['torque_peak_microNm'] - 1.0):+.1f} %)")
        print(f"  RMSE force théorie/exp : {metrics[eps]['force']['rmse_theory']:.1f} / "
              f"{metrics[eps]['force']['rmse_experiment']:.1f} mN ; "
              f"couple théorie/exp : {metrics[eps]['torque']['rmse_theory']:.1f} / "
              f"{metrics[eps]['torque']['rmse_experiment']:.1f} µN·m")
    print("Rappel (audit) : la cible « théorie » numérisée est tronquée aux pics de couple "
          "(occlusion par les marqueurs) ; RMSE = ordres de grandeur.")

    metrics_path = OUT_DIR / f"validation_figure7_{pmode}_L{n_layers}p{n_phi}.json"
    with open(metrics_path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "model_version": Base.MODEL_VERSION,
                "pmode": pmode,
                "discretization": {"dt": dt, "n_layers": n_layers, "n_phi": n_phi},
                "metrics": {
                    str(eps): {
                        "extrema": metrics[eps]["extrema"],
                        "pressure_info": metrics[eps]["pressure_info"],
                        "force": metrics[eps]["force"],
                        "torque": metrics[eps]["torque"],
                    }
                    for eps in (0.8, 1.0)
                },
            },
            fh,
            ensure_ascii=False,
            indent=2,
        )
    print(f"Figure : {fig_path}\nMétriques : {metrics_path}")

    if show:
        plt.show()
    return fig, metrics


def main() -> None:
    parser = argparse.ArgumentParser(description="Validation figure 7 (moteur alpha V3).")
    parser.add_argument("--pmode", choices=("raw", "corr", "rescale"), default="rescale")
    parser.add_argument("--quick", action="store_true", help="discrétisation réduite (2 couches, 8 secteurs)")
    parser.add_argument("--layers", type=int, default=None, help="nombre de couches (défaut : 6, ou 2 avec --quick)")
    parser.add_argument("--baseline", action="store_true", help="(ré)génère la baseline de non-régression")
    parser.add_argument("--show", action="store_true", help="affiche les figures")
    args = parser.parse_args()

    if args.baseline:
        write_baseline(args.pmode)
        return

    n_layers = args.layers if args.layers is not None else (2 if args.quick else 6)
    n_phi = 8 if args.quick else 16
    plot_validation(pmode=args.pmode, n_layers=n_layers, n_phi=n_phi, show=args.show)


if __name__ == "__main__":
    main()
