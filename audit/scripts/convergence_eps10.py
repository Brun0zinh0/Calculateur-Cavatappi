"""Sonde de convergence spatiale pour eps=1.0 : (layers, phi) = (1,8), (3,16), (5,32)."""
import sys
from pathlib import Path

import numpy as np

WORKSPACE = Path(
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents"
    r"\Stage muscle artificièle\modèle\modèle chinois\Espace de travail"
)
sys.path.insert(0, str(WORKSPACE))

import validation_figure7 as figure7  # noqa: E402

REAL_PDF = WORKSPACE / "Article support" / figure7.PDF_FILENAME
_orig = figure7.extract_figure7_image
_cache = {}


def patched(pdf_path=REAL_PDF):
    if "img" not in _cache:
        _cache["img"] = _orig(REAL_PDF)
    return _cache["img"]


figure7.extract_figure7_image = patched

targets = figure7.digitize_figure7()
theory = figure7.digitize_figure7_theory()

for n_layers, n_phi in ((1, 8), (3, 16), (5, 32)):
    _, data = figure7.run_model_for_eps(1.0, targets=targets, dt=1.0, n_layers=n_layers, n_phi=n_phi)
    win = figure7.window_model_data(data)
    f = win["force_total_mN"]
    tq = win["torque_act_microNm"]
    tx, ty = theory[1.0]["force"]
    rmse_th = figure7.rmse_to_digitized(tx, ty, win["time"], f)
    print(
        f"layers={n_layers} phi={n_phi}: force {np.min(f):.1f}..{np.max(f):.1f} mN, "
        f"torque {np.min(tq):.1f}..{np.max(tq):.1f} uNm, RMSE force/theorie {rmse_th:.1f} mN"
    )
