"""Compare le moteur alpha V2 (2026.07.30) et le moteur livrable (2026.07.17)
sur le protocole exact de validation_figure7 (pression numerisee figure 7).

Aucun fichier du projet n'est modifie. Les deux Base.py sont charges sous des
noms de module distincts via importlib.
"""
import importlib.util
import sys
import time
from functools import partial
from pathlib import Path

import numpy as np

WORKSPACE = Path(r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle\modèle\modèle chinois\Espace de travail")
LIVRABLE = WORKSPACE / "livrable"
ALPHA = WORKSPACE / "alpha V2"
SCRATCH = Path(r"C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad")
PDF_REAL = WORKSPACE / "Article support" / "Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure.pdf"

for p in (str(LIVRABLE), str(WORKSPACE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import validation_figure7 as f7  # noqa: E402  (Base -> livrable/Base.py)

_orig_extract = f7.extract_figure7_image
f7.PDF_PATH = PDF_REAL
f7.extract_figure7_image = partial(_orig_extract, PDF_REAL)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


alpha_base = load_module("alpha_Base", ALPHA / "Base.py")
alpha_params = load_module("alpha_parametres", ALPHA / "parametres.py")

liv_base = f7.Base
from livrable_parametres import DEFAULT_SETTINGS as LIV_DEFAULTS, build_config as liv_build  # noqa: E402

print(f"livrable Base : {liv_base.__file__} version={liv_base.MODEL_VERSION}")
print(f"alpha   Base : {alpha_base.__file__} version={alpha_base.MODEL_VERSION}")

COMMON = {
    "use_fixed_duration": False,
    "p_max_mpa": f7.PAPER_PMAX_MPA,
    "dt": 1.0,
    "n_layers": 3,
    "n_phi": 16,
    "pre_steps": 24,
    "flow_rate_mL_min": f7.PAPER_FLOW_RATE_ML_MIN,
    "volume_mL": f7.PAPER_VOLUME_ML,
    "nonlinear_pressure": True,
}
N_CYCLES = int(np.ceil(f7.PAPER_T_MAX / (2.0 * 60.0 * f7.PAPER_VOLUME_ML / f7.PAPER_FLOW_RATE_ML_MIN))) + 1


def run_livrable(eps, pressure_time, pressure_mpa):
    settings = dict(LIV_DEFAULTS)
    settings.update(COMMON)
    settings.update(
        {
            "eps": eps,
            "n_cycles": N_CYCLES,
            "constitutive_mode": "generalized_maxwell",
            "prestrain_reference_mode": "elastic_tk_reference",
            "maxwell_anisotropy_mode": "paper_equal",
            "nylon_condition_mode": "bonded_linear",
            "pressure_end_force_mode": "none",
            "pressure_end_force_scale": 0.0,
        }
    )
    config = liv_build(settings)
    _, data = liv_base.run_blocked_actuation(config, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    return data


def run_alpha(eps, pressure_time, pressure_mpa):
    settings = dict(alpha_params.DEFAULT_SETTINGS)
    settings.update(COMMON)
    settings.update({"eps": eps, "n_cycles": N_CYCLES})
    config = alpha_params.build_config(settings)
    _, data = alpha_base.run_blocked_actuation(config, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    return data


if __name__ == "__main__":
    t0 = time.time()
    targets = f7.digitize_figure7()
    theory = f7.digitize_figure7_theory()
    for eps in (0.8, 1.0):
        pt, pp = f7.digitized_pressure_history(eps, targets)
        liv = f7.window_model_data(run_livrable(eps, pt, pp))
        alp = f7.window_model_data(run_alpha(eps, pt, pp))
        print(f"\n=== eps = {eps:.1f} ===")
        for label, win in (("livrable", liv), ("alphaV2 ", alp)):
            tx_f, ty_f = targets[eps]["force"]
            tx_t, ty_t = targets[eps]["torque"]
            th_f = theory[eps]["force"]
            th_t = theory[eps]["torque"]
            rmse_f_exp = f7.rmse_to_digitized(tx_f, ty_f, win["time"], win["force_total_mN"])
            rmse_t_exp = f7.rmse_to_digitized(tx_t, ty_t, win["time"], win["torque_act_microNm"])
            rmse_f_th = f7.rmse_to_digitized(th_f[0], th_f[1], win["time"], win["force_total_mN"])
            rmse_t_th = f7.rmse_to_digitized(th_t[0], th_t[1], win["time"], win["torque_act_microNm"])
            print(
                f"  {label}: F [{np.min(win['force_total_mN']):.1f}, {np.max(win['force_total_mN']):.1f}] mN, "
                f"T [{np.min(win['torque_act_microNm']):.1f}, {np.max(win['torque_act_microNm']):.1f}] uNm"
            )
            print(
                f"            RMSE F theorie/exp = {rmse_f_th:.2f}/{rmse_f_exp:.2f} mN ; "
                f"RMSE T theorie/exp = {rmse_t_th:.2f}/{rmse_t_exp:.2f} uNm"
            )
        # difference directe entre moteurs sur la grille commune
        n = min(len(liv["time"]), len(alp["time"]))
        df = np.abs(liv["force_total_mN"][:n] - alp["force_total_mN"][:n])
        dt_ = np.abs(liv["torque_act_microNm"][:n] - alp["torque_act_microNm"][:n])
        print(
            f"  diff moteurs: max|dF| = {np.max(df):.3f} mN (moy {np.mean(df):.3f}) ; "
            f"max|dT| = {np.max(dt_):.3f} uNm (moy {np.mean(dt_):.3f})"
        )
    print(f"\nDuree totale: {time.time() - t0:.1f} s")
