# -*- coding: utf-8 -*-
"""Matrice de runs figure 7 : {bias_profile} x {section_update_mode} x {eps}.

Usage:
    python ce_matrix_runner.py <bias> <section> <eps> <pmode> [dt n_layers n_phi tag]

  bias    : paper_linear | uniform_twist
  section : fixed | updated
  eps     : 0.8 | 1.0
  pmode   : corr (calibration vraie des ticks) | rescale (corr + pics recales a 1.42 MPa)
            | raw (pression numerisee biaisee, comme la validation d'origine)
  dt, n_layers, n_phi : discretisation (defaut 1.0, 3, 16)
  tag     : etiquette du run (defaut auto)

Chaque run ajoute une ligne JSON a matrix_results.jsonl et sauvegarde la trace
fenetree en npz. NE MODIFIE AUCUN FICHIER DU PROJET.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

SCRATCH = Path(
    r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude"
    r"\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-"
    r"artifici-le-mod-le-mod-le-chinois-Espace-de-travail"
    r"\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad"
)
WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
sys.path.insert(0, str(WORKSPACE))

import validation_figure7 as f7  # noqa: E402
from livrable_parametres import DEFAULT_SETTINGS, build_config  # noqa: E402
import Base  # noqa: E402  (moteur livrable, verifie identique a alpha V2 sur ce protocole)

REAL_PDF = WORKSPACE / "Article support" / f7.PDF_FILENAME
_orig_extract = f7.extract_figure7_image
_cache = {}


def patched(pdf_path=REAL_PDF):
    if "img" not in _cache:
        _cache["img"] = _orig_extract(REAL_PDF)
    return _cache["img"]


f7.extract_figure7_image = patched

# Calibration vraie (fit sur rows des etiquettes mesurees dans ce_fig7_measure2)
TICK_FITS = {
    0.8: ([21.0, 107.5, 195.0, 282.0], [1.5, 1.0, 0.5, 0.0], (120, 648, 704, 946)),
    1.0: ([15.0, 104.0, 192.5, 282.0], [1.5, 1.0, 0.5, 0.0], (899, 648, 1484, 946)),
}


def script_to_true_pressure(eps, y_script):
    rows_t, vals_t, crop = TICK_FITS[eps]
    A, B = np.polyfit(rows_t, vals_t, 1)
    h = crop[3] - crop[1]
    height = max(1.0, float(h - 1))
    row = (1.5 - np.asarray(y_script)) * height / 1.5
    return A * row + B


def cycle_extrema(t, v, t0=185.0, min_gap=8.0):
    """Pics et vallees par cycle sur la fenetre etablie (> t0)."""
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
    valls = [(tt, -vv) for tt, vv in _find(-v, -med)]
    pv = np.asarray([p[1] for p in peaks]) if peaks else np.asarray([np.nan])
    vv = np.asarray([p[1] for p in valls]) if valls else np.asarray([np.nan])
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


def make_cfg(eps, dt, n_layers, n_phi, bias, section):
    n_cycles = int(np.ceil(f7.PAPER_T_MAX / (2.0 * 60.0 * f7.PAPER_VOLUME_ML / f7.PAPER_FLOW_RATE_ML_MIN))) + 1
    settings = dict(DEFAULT_SETTINGS)
    settings.update(
        {
            "eps": eps,
            "n_cycles": n_cycles,
            "use_fixed_duration": False,
            "p_max_mpa": f7.PAPER_PMAX_MPA,
            "dt": dt,
            "n_layers": n_layers,
            "n_phi": n_phi,
            "pre_steps": 24,
            "flow_rate_mL_min": f7.PAPER_FLOW_RATE_ML_MIN,
            "volume_mL": f7.PAPER_VOLUME_ML,
            "nonlinear_pressure": True,
            "constitutive_mode": f7.CONSTITUTIVE_MODE,
            "prestrain_reference_mode": f7.PRESTRAIN_REFERENCE_MODE,
            "maxwell_anisotropy_mode": f7.MAXWELL_ANISOTROPY_MODE,
            "nylon_condition_mode": f7.NYLON_CONDITION_MODE,
            "pressure_end_force_mode": f7.PRESSURE_END_FORCE_MODE,
            "pressure_end_force_scale": f7.PRESSURE_END_FORCE_SCALE,
            "bias_angle_profile": bias,
            "section_update_mode": section,
        }
    )
    return build_config(settings)


def build_targets(pmode):
    targets = f7.digitize_figure7()
    info = {}
    if pmode == "raw":
        return targets, info
    out = {e: dict(targets[e]) for e in targets}
    for eps in (0.8, 1.0):
        tp, pp = targets[eps]["pressure"]
        pp_true = np.clip(script_to_true_pressure(eps, pp), 0.0, None)
        if pmode == "rescale":
            ext = cycle_extrema(tp, pp_true, t0=185.0)
            scale = 1.42 / ext["peak_mean"]
            pp_true = np.clip(pp_true * scale, 0.0, None)
            info[eps] = dict(scale=scale, peak_before=ext["peak_mean"])
        out[eps]["pressure"] = (tp, pp_true)
        ext2 = cycle_extrema(tp, pp_true, t0=185.0)
        info.setdefault(eps, {})
        info[eps]["p_peak_mean"] = ext2["peak_mean"]
        info[eps]["p_peak_max"] = ext2["peak_max"]
    return out, info


def main():
    bias, section, eps_s, pmode = sys.argv[1:5]
    eps = float(eps_s)
    dt = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    n_layers = int(sys.argv[6]) if len(sys.argv) > 6 else 3
    n_phi = int(sys.argv[7]) if len(sys.argv) > 7 else 16
    short = {"paper_linear": "lin", "uniform_twist": "arc", "fixed": "fix", "updated": "upd"}
    tag = sys.argv[8] if len(sys.argv) > 8 else (
        f"{short[bias]}_{short[section]}_e{eps_s.replace('.', '')}_{pmode}"
        f"_L{n_layers}p{n_phi}d{str(dt).replace('.', '')}"
    )

    t_start = time.time()
    targets, pinfo = build_targets(pmode)
    cfg = make_cfg(eps, dt, n_layers, n_phi, bias, section)
    # dt effectif = grille de pression fournie (cfg.dt est ignore avec une
    # histoire imposee) : prefix_dt=dt, figure_dt=dt/2 (defaut 1.0/0.5 = f7).
    pressure_time, pressure_mpa = f7.digitized_pressure_history(
        eps, targets, prefix_dt=dt, figure_dt=0.5 * dt
    )
    _, data = Base.run_blocked_actuation(cfg, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    win = f7.window_model_data(data)
    elapsed = time.time() - t_start

    resF = cycle_extrema(win["time"], win["force_total_mN"])
    resT = cycle_extrema(win["time"], win["torque_act_microNm"])
    resP = cycle_extrema(win["time"], win["pressure_MPa"])

    rec = dict(
        tag=tag, bias=bias, section=section, eps=eps, pmode=pmode,
        dt=dt, n_layers=n_layers, n_phi=n_phi, elapsed_s=round(elapsed, 1),
        p_info={str(k): v for k, v in pinfo.items()},
        force=resF, torque=resT, pressure=resP,
    )
    with open(SCRATCH / "matrix_results.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    np.savez_compressed(
        SCRATCH / f"mx_{tag}.npz",
        time=win["time"], force=win["force_total_mN"],
        torque=win["torque_act_microNm"], pressure=win["pressure_MPa"],
    )
    print(f"[{tag}] {elapsed:.1f} s")
    print(f"  P    : pics moy={resP['peak_mean']:.3f} max={resP['peak_max']:.3f} MPa")
    print(f"  F    : pics n={resF['n_peaks']} moy={resF['peak_mean']:.1f} "
          f"[{resF['peak_min']:.1f}..{resF['peak_max']:.1f}] vallee={resF['valley_mean']:.1f} "
          f"ampl={resF['amplitude']:.1f} mN")
    print(f"  T    : pics n={resT['n_peaks']} moy={resT['peak_mean']:.1f} "
          f"[{resT['peak_min']:.1f}..{resT['peak_max']:.1f}] vallee={resT['valley_mean']:.1f} "
          f"ampl={resT['amplitude']:.1f} uNm")


if __name__ == "__main__":
    main()
