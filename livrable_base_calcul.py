"""Moteur de calcul expose au livrable Streamlit.

Ce module sert de façade claire pour l'interface livrable. Il charge
``Base.py`` par chemin absolu, à partir de l'emplacement de ce fichier.
Ainsi l'interface peut être lancée depuis n'importe quel répertoire.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parent
BASE_CANDIDATES = (
    MODULE_DIR / "Base.py",
    MODULE_DIR.parent / "Base.py",
)


def _load_base_module():
    for candidate in BASE_CANDIDATES:
        if not candidate.exists():
            continue
        spec = importlib.util.spec_from_file_location("_livrable_base_model", candidate)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    searched = "\n".join(f"- {path}" for path in BASE_CANDIDATES)
    raise FileNotFoundError(
        "Impossible de trouver Base.py pour le livrable. Emplacements cherches :\n"
        f"{searched}\n\n"
        "Place Base.py dans le dossier livrable ou dans son dossier parent."
    )


_BASE = _load_base_module()
BASE_PATH = Path(_BASE.__file__).resolve()

TCPAMaxwellBlockedModel = _BASE.TCPAMaxwellBlockedModel
add_corrected_output_conventions = _BASE.add_corrected_output_conventions
cyclic_pressure_history = _BASE.cyclic_pressure_history
default_discretization = _BASE.default_discretization
default_geometry_params = _BASE.default_geometry_params
default_material_params = _BASE.default_material_params
default_maxwell_tensile_params = _BASE.default_maxwell_tensile_params
default_simulation_config = _BASE.default_simulation_config
ramp_hold_pressure_history = _BASE.ramp_hold_pressure_history
run_blocked_actuation = _BASE.run_blocked_actuation
run_suspended_actuation = _BASE.run_suspended_actuation
summary = _BASE.summary
