from __future__ import annotations

import json
import os
import pickle
import tempfile
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np


SettingValue = float | int | str | bool

PSI_TO_MPA = 0.006894757293168361
P_MAX_MPA = 1.50
P_MAX_PSI = P_MAX_MPA / PSI_TO_MPA
DEFAULT_PARALLEL_WORKERS = max(1, min(4, (os.cpu_count() or 2) - 1))

CACHE_DIR = Path(tempfile.gettempdir()) / "tcpa_cavatappi_beta_cache"
SETTINGS_PATH = CACHE_DIR / "cavatappi_beta_settings.json"
TIMING_PROFILE_PATH = CACHE_DIR / "cavatappi_beta_timing_profile.json"
BLOCKED_RESULT_PATH = CACHE_DIR / "cavatappi_beta_blocked.pkl"
RELAXATION_RESULT_PATH = CACHE_DIR / "cavatappi_beta_relaxation.pkl"
PRESTRAIN_RESULT_PATH = CACHE_DIR / "cavatappi_beta_prestrain.pkl"
SUSPENDED_RESULT_PATH = CACHE_DIR / "cavatappi_beta_suspended.pkl"
HYSTERESIS_RESULT_PATH = CACHE_DIR / "cavatappi_beta_hysteresis.pkl"
RESULT_CACHE_PATHS = (
    BLOCKED_RESULT_PATH,
    RELAXATION_RESULT_PATH,
    PRESTRAIN_RESULT_PATH,
    SUSPENDED_RESULT_PATH,
    HYSTERESIS_RESULT_PATH,
)
SETTINGS_SCHEMA_VERSION = 12

INTEGRATION_OPTIONS = ["exponential", "paper_explicit"]
VISUAL_STATE_OPTIONS = ["fabricated", "prestrained"]
SECTION_UPDATE_OPTIONS = ["fixed", "updated"]
BIAS_ANGLE_PROFILE_OPTIONS = ["paper_linear", "uniform_twist"]
CONSTITUTIVE_OPTIONS = ["generalized_maxwell", "instantaneous_elastic"]
AXIAL_MODULUS_OPTIONS = ["paper_table", "maxwell_sum"]
PRESTRAIN_REFERENCE_OPTIONS = ["elastic_tk_reference", "viscoelastic_ramp"]
MAXWELL_ANISOTROPY_OPTIONS = ["axial_test_only", "paper_equal"]
NYLON_CONDITION_OPTIONS = ["bonded_linear", "tension_only", "axially_sliding_confined"]
PRESSURE_INPUT_OPTIONS = ["generated", "measured_csv"]

INTEGRATION_LABELS = {
    "paper_explicit": "Euler explicite",
    "exponential": "Intégration exponentielle stable",
}
VISUAL_STATE_LABELS = {
    "fabricated": "Fabriqué",
    "prestrained": "Précontraint",
}
SECTION_UPDATE_LABELS = {
    "fixed": "Mise à jour hélicoïdale uniquement",
    "updated": "Section radiale et orientation évolutives",
}
BIAS_ANGLE_PROFILE_LABELS = {
    "paper_linear": "Variation linéaire avec le rayon",
    "uniform_twist": "Torsion uniforme (loi en tangente)",
}
CONSTITUTIVE_LABELS = {
    "generalized_maxwell": "Maxwell généralisé",
    "instantaneous_elastic": "Élastique instantané",
}
AXIAL_MODULUS_LABELS = {
    "paper_table": "Module axial saisi manuellement",
    "maxwell_sum": "Somme des modules de Maxwell",
}
PRESTRAIN_REFERENCE_LABELS = {
    "elastic_tk_reference": "Précontrainte élastique conservée",
    "viscoelastic_ramp": "Précontrainte viscoélastique",
}
MAXWELL_ANISOTROPY_LABELS = {
    "axial_test_only": "Relaxation limitée à la direction axiale",
    "paper_equal": "Relaxation proportionnelle dans toutes les directions",
}
NYLON_CONDITION_LABELS = {
    "bonded_linear": "Linéaire bilatéral, lié aux extrémités",
    "tension_only": "Filament non collé, traction seulement",
    "axially_sliding_confined": "Glissant axialement, confiné dans le tube",
}
PRESSURE_INPUT_LABELS = {
    "generated": "Profil généré par le modèle",
    "measured_csv": "Historique pression/temps mesuré (CSV)",
}


DEFAULT_SETTINGS: dict[str, SettingValue] = {
    "_settings_schema_version": SETTINGS_SCHEMA_VERSION,
    "eps": 0.8,
    "eps_study_min": 0.0,
    "eps_study_max": 1.2,
    "eps_study_points": 13,
    "p_max_mpa": 1.5,
    "rout_mm": 1.0,
    "rin_mm": 0.4,
    "nylon_diameter_mm": 0.77,
    "rho0_mm": 2.16,
    "alpha0_deg": 10.53,
    "theta_f_deg": 37.91,
    "initial_length_mm": 32.45,
    "bias_angle_profile": "paper_linear",
    "section_update_mode": "fixed",
    "n_cycles": 11,
    "hysteresis_cycle": 1,
    "hysteresis_cycles": "1",
    "hysteresis_compare_mode": "current",
    "hysteresis_prestrain_values": "0.6, 0.8, 1.0",
    "hysteresis_pressure_rates_mpa_s": "0.05, 0.10, 0.20",
    "use_fixed_duration": False,
    "duration_s": 500.0,
    "flow_rate_mL_min": 10.0,
    "volume_mL": 1.50,
    "nonlinear_pressure": False,
    "pressure_input_mode": "generated",
    "measured_pressure_time_column": "",
    "measured_pressure_column": "",
    "measured_pressure_unit": "MPa",
    "measured_pressure_file_hash": "",
    "measured_pressure_subtract_initial": True,
    "relaxation_ramp_time_s": 9.0,
    "relaxation_hold_time_s": 300.0,
    "suspended_mass_g": 100.0,
    "suspended_duration_s": 120.0,
    "suspended_pressure_rate_mpa_s": 0.10,
    "suspended_hold_pressure": False,
    "suspended_show_geometry_plot": False,
    "dt": 0.5,
    "pre_steps": 24,
    "integration": "exponential",
    "prestrain_reference_mode": "elastic_tk_reference",
    "constitutive_mode": "generalized_maxwell",
    "maxwell_anisotropy_mode": "paper_equal",
    "axial_modulus_mode": "maxwell_sum",
    "E_axial_mpa": 31.24,
    "E_radius_mpa": 8.82,
    "G12_mpa": 7.24,
    "nu12": 0.205,
    "nu23": 0.422,
    "maxwell_E0_mpa": 6.36,
    "maxwell_E1_mpa": 20.67,
    "maxwell_eta1_mpa_s": 154.57,
    "maxwell_E2_mpa": 5.98,
    "maxwell_eta2_mpa_s": 977.79,
    "maxwell_E3_mpa": 4.75,
    "maxwell_eta3_mpa_s": 11044.83,
    "E_nylon_mpa": 3.69e3,
    "G_nylon_mpa": 0.79e3,
    "nylon_axial_prestrain_coupling": 1.0,
    "nylon_axial_actuation_coupling": 1.0,
    "nylon_condition_mode": "bonded_linear",
    "nylon_scale": 1.0,
    "n_layers": 4,
    "n_phi": 16,
    "parallel_workers": DEFAULT_PARALLEL_WORKERS,
    "view_elev_deg": 22.0,
    "view_azim_deg": -58.0,
}


@dataclass
class MaxwellTensileParams:
    E0: float = 6.36
    E1: float = 20.67
    eta1: float = 154.57
    E2: float = 5.98
    eta2: float = 977.79
    E3: float = 4.75
    eta3: float = 11044.83

    @property
    def E(self) -> np.ndarray:
        return np.array([self.E1, self.E2, self.E3], dtype=float)

    @property
    def eta(self) -> np.ndarray:
        return np.array([self.eta1, self.eta2, self.eta3], dtype=float)

    @property
    def E_total(self) -> float:
        return self.E0 + self.E1 + self.E2 + self.E3


@dataclass
class MaterialParams:
    E_axial: float = 31.24
    E_radius: float = 8.82
    G12: float = 7.24
    nu12: float = 0.205
    nu23: float = 0.422
    maxwell: MaxwellTensileParams = field(default_factory=MaxwellTensileParams)
    E_nylon: float = 3.69e3
    G_nylon: float = 0.79e3
    maxwell_anisotropy_mode: str = "paper_equal"
    nylon_condition_mode: str = "bonded_linear"
    nylon_axial_prestrain_coupling: float = 1.0
    nylon_axial_actuation_coupling: float = 1.0


@dataclass
class GeometryParams:
    Rout: float = 1.0
    Rin: float = 0.4
    r_nylon: float = 0.77 / 2.0
    rho0: float = 2.16
    alpha0_deg: float = 10.53
    theta_f_deg: float = 37.91
    initial_length: float = 32.45
    bias_angle_profile: str = "paper_linear"
    section_update_mode: str = "fixed"


@dataclass
class SimulationParams:
    eps: float = 1.0
    n_cycles: int = 11
    duration_s: float | None = 500.0
    Pmax: float = P_MAX_MPA
    dt: float = 0.5
    n_layers: int = 8
    n_phi: int = 24
    pre_steps: int = 24
    integration: str = "exponential"
    prestrain_reference_mode: str = "elastic_tk_reference"
    constitutive_mode: str = "generalized_maxwell"
    axial_modulus_mode: str = "maxwell_sum"
    flow_rate_mL_min: float = 10.0
    volume_mL: float = 1.50
    nonlinear_pressure: bool = True
    nylon_stiffness_scale: float = 1.0
    mat: MaterialParams = field(default_factory=MaterialParams)
    geom: GeometryParams = field(default_factory=GeometryParams)


def option_index(options: list[str], value: object) -> int:
    try:
        return options.index(str(value))
    except ValueError:
        return 0


def load_settings() -> dict[str, SettingValue]:
    try:
        with SETTINGS_PATH.open("r", encoding="utf-8") as file:
            saved = json.load(file)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        saved = {}

    try:
        saved_version = int(saved.get("_settings_schema_version", 1))
    except (TypeError, ValueError, AttributeError):
        saved_version = 1
    if saved_version < SETTINGS_SCHEMA_VERSION:
        saved = dict(saved)
        try:
            if abs(float(saved.get("dt", DEFAULT_SETTINGS["dt"])) - 2.0) < 1e-12:
                saved["dt"] = DEFAULT_SETTINGS["dt"]
        except (TypeError, ValueError):
            saved["dt"] = DEFAULT_SETTINGS["dt"]
        if str(saved.get("integration", "")) in {"paper_incremental", "paper"}:
            saved["integration"] = "exponential"
        saved["nylon_axial_actuation_coupling"] = 1.0
        saved["use_fixed_duration"] = False
        saved["nonlinear_pressure"] = False
        saved["p_max_mpa"] = min(float(saved.get("p_max_mpa", P_MAX_MPA)), P_MAX_MPA)
        saved["section_update_mode"] = DEFAULT_SETTINGS["section_update_mode"]
        saved["bias_angle_profile"] = DEFAULT_SETTINGS["bias_angle_profile"]
        saved.setdefault("constitutive_mode", DEFAULT_SETTINGS["constitutive_mode"])
        saved["axial_modulus_mode"] = DEFAULT_SETTINGS["axial_modulus_mode"]
        saved["prestrain_reference_mode"] = DEFAULT_SETTINGS["prestrain_reference_mode"]
        saved["maxwell_anisotropy_mode"] = DEFAULT_SETTINGS["maxwell_anisotropy_mode"]
        saved["nylon_condition_mode"] = DEFAULT_SETTINGS["nylon_condition_mode"]
        saved["nylon_scale"] = 1.0
        saved["nylon_axial_prestrain_coupling"] = 1.0
        saved["nylon_axial_actuation_coupling"] = 1.0
        saved.setdefault("pressure_input_mode", DEFAULT_SETTINGS["pressure_input_mode"])
        saved["_settings_schema_version"] = SETTINGS_SCHEMA_VERSION

    settings = dict(DEFAULT_SETTINGS)
    settings.update({key: saved[key] for key in settings.keys() & saved.keys()})
    for key in (
        "n_cycles",
        "hysteresis_cycle",
        "eps_study_points",
        "n_layers",
        "n_phi",
        "pre_steps",
        "parallel_workers",
    ):
        settings[key] = int(settings[key])
    for key in (
        "use_fixed_duration",
        "nonlinear_pressure",
        "measured_pressure_subtract_initial",
        "suspended_hold_pressure",
        "suspended_show_geometry_plot",
    ):
        settings[key] = bool(settings[key])
    return settings


def save_settings(settings: dict[str, SettingValue]) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with SETTINGS_PATH.open("w", encoding="utf-8") as file:
        json.dump(settings, file, indent=2, sort_keys=True)


def reset_settings_and_cache() -> dict[str, SettingValue]:
    settings = dict(DEFAULT_SETTINGS)
    save_settings(settings)
    for path in RESULT_CACHE_PATHS:
        try:
            path.unlink(missing_ok=True)
        except OSError:
            pass
    return settings


def load_result_cache(path: Path) -> Any | None:
    try:
        with path.open("rb") as file:
            return pickle.load(file)
    except (FileNotFoundError, OSError, pickle.PickleError, AttributeError, EOFError, ImportError, TypeError, ValueError):
        return None


def save_result_cache(path: Path, payload: Any) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as file:
        pickle.dump(payload, file)


def geometry_error(settings: dict[str, SettingValue]) -> str | None:
    rout = float(settings["rout_mm"])
    rin = float(settings["rin_mm"])
    nylon_diameter = float(settings["nylon_diameter_mm"])
    rho0 = float(settings["rho0_mm"])
    values = np.array(
        [rout, rin, nylon_diameter, rho0, settings["initial_length_mm"], settings["alpha0_deg"]], dtype=float
    )
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        return "Les dimensions et l'angle hélicoïdal doivent être finis et strictement positifs."
    if rin >= rout:
        return "Rin doit être inférieur à Rout."
    if nylon_diameter > 2.0 * rin:
        return "Le diamètre du nylon ne doit pas dépasser le diamètre intérieur du tube."
    if rho0 <= rout:
        return "rho0 doit être supérieur à Rout pour une ligne centrale hélicoïdale."
    return None


def material_error(settings: dict[str, SettingValue]) -> str | None:
    E0 = float(settings["maxwell_E0_mpa"])
    branch_E = np.array(
        [settings["maxwell_E1_mpa"], settings["maxwell_E2_mpa"], settings["maxwell_E3_mpa"]], dtype=float
    )
    branch_eta = np.array(
        [settings["maxwell_eta1_mpa_s"], settings["maxwell_eta2_mpa_s"], settings["maxwell_eta3_mpa_s"]],
        dtype=float,
    )
    if E0 < 0.0 or np.any(branch_E < 0.0) or E0 + float(branch_E.sum()) <= 0.0:
        return "Au moins un module de Maxwell doit être strictement positif et aucun ne peut être négatif."
    if np.any((branch_E > 0.0) & (branch_eta <= 0.0)):
        return "Chaque branche de Maxwell active doit avoir une viscosité strictement positive."

    maxwell_sum = E0 + float(branch_E.sum())
    axial_modulus = (
        maxwell_sum
        if str(settings.get("axial_modulus_mode", "maxwell_sum")) == "maxwell_sum"
        else float(settings["E_axial_mpa"])
    )
    if (
        str(settings.get("maxwell_anisotropy_mode", "paper_equal")) == "axial_test_only"
        and axial_modulus - float(branch_E.sum()) <= 0.0
    ):
        return (
            "La relaxation axiale identifiee exige un module axial instantane superieur "
            "a la somme des branches transitoires."
        )
    positive_moduli = np.array(
        [
            axial_modulus,
            settings["E_radius_mpa"],
            settings["G12_mpa"],
            settings["E_nylon_mpa"],
            settings["G_nylon_mpa"],
        ],
        dtype=float,
    )
    if np.any(~np.isfinite(positive_moduli)) or np.any(positive_moduli <= 0.0):
        return "Tous les modules du tube et du nylon doivent être strictement positifs."

    Ea, Er, G12 = positive_moduli[:3]
    nu12 = float(settings["nu12"])
    nu23 = float(settings["nu23"])
    try:
        nu21 = nu12 * Er / Ea
        system = np.array(
            [
                [1.0, -2.0 * nu12, 0.0, 0.0],
                [0.0, 1.0, -nu12, -nu12],
                [-nu21, 1.0 - nu23, 0.0, 0.0],
                [0.0, -nu21, 1.0, -nu23],
            ],
            dtype=float,
        )
        C11, C12, C22, C23 = np.linalg.solve(system, np.array([Ea, 0.0, 0.0, Er], dtype=float))
        C = np.zeros((6, 6), dtype=float)
        C[0, 0] = C11
        C[0, 1] = C[1, 0] = C12
        C[0, 2] = C[2, 0] = C12
        C[1, 1] = C[2, 2] = C22
        C[1, 2] = C[2, 1] = C23
        C[3, 3] = 0.5 * (C22 - C23)
        C[4, 4] = C[5, 5] = G12
        eig_min = float(np.linalg.eigvalsh(C).min())
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        eig_min = -np.inf
    if not np.isfinite(eig_min) or eig_min <= 1e-10:
        return "Les modules et coefficients de Poisson produisent une matrice de rigidité non physique."
    return None


def numerical_error(settings: dict[str, SettingValue]) -> str | None:
    dt = float(settings["dt"])
    if not np.isfinite(dt) or dt <= 0.0:
        return "Le pas de temps doit être strictement positif."
    pressure = float(settings["p_max_mpa"])
    if not 0.0 <= pressure <= P_MAX_MPA:
        return f"La pression doit rester comprise entre 0 et {P_MAX_MPA:g} MPa pour ce modèle."
    if float(settings["eps"]) < 0.0:
        return "La précontrainte initiale ne peut pas être négative."
    if int(settings["n_cycles"]) < 1 or int(settings["n_layers"]) < 1 or int(settings["n_phi"]) < 4:
        return "Le nombre de cycles et le maillage doivent être strictement positifs."
    if str(settings.get("section_update_mode")) not in SECTION_UPDATE_OPTIONS:
        return "Le mode de mise à jour de section est inconnu."
    if str(settings.get("bias_angle_profile")) not in BIAS_ANGLE_PROFILE_OPTIONS:
        return "Le profil radial de l'angle de biais est inconnu."
    if str(settings.get("constitutive_mode")) not in CONSTITUTIVE_OPTIONS:
        return "Le mode constitutif est inconnu."
    if str(settings.get("axial_modulus_mode")) not in AXIAL_MODULUS_OPTIONS:
        return "La convention du module axial est inconnue."
    if str(settings.get("prestrain_reference_mode")) not in PRESTRAIN_REFERENCE_OPTIONS:
        return "Le mode de traitement de la précontrainte est inconnu."
    if str(settings.get("maxwell_anisotropy_mode")) not in MAXWELL_ANISOTROPY_OPTIONS:
        return "Le mode d'anisotropie viscoelastique est inconnu."
    if str(settings.get("nylon_condition_mode")) not in NYLON_CONDITION_OPTIONS:
        return "La condition physique du nylon est inconnue."
    if str(settings.get("pressure_input_mode", "generated")) not in PRESSURE_INPUT_OPTIONS:
        return "La source de pression est inconnue."
    for key in ("nylon_axial_prestrain_coupling", "nylon_axial_actuation_coupling"):
        if not 0.0 <= float(settings[key]) <= 1.0:
            return "Les coefficients de couplage du nylon doivent rester entre 0 et 1."
    if str(settings["integration"]) == "paper_explicit" and str(settings["constitutive_mode"]) == "generalized_maxwell":
        E = np.array(
            [settings["maxwell_E1_mpa"], settings["maxwell_E2_mpa"], settings["maxwell_E3_mpa"]], dtype=float
        )
        eta = np.array(
            [settings["maxwell_eta1_mpa_s"], settings["maxwell_eta2_mpa_s"], settings["maxwell_eta3_mpa_s"]],
            dtype=float,
        )
        active = E > 0.0
        if np.any(active):
            tau_min = float(np.min(eta[active] / E[active]))
            if dt >= 2.0 * tau_min:
                return f"Euler explicite est instable avec ce pas : utilisez dt < {2.0 * tau_min:.3g} s."
            pre_dt = (
                60.0
                * float(settings["eps"])
                * float(settings["initial_length_mm"])
                / 20.0
                / max(int(settings["pre_steps"]), 1)
            )
            if pre_dt >= 2.0 * tau_min:
                return (
                    "Euler explicite est instable pendant la précontrainte : "
                    f"augmentez le nombre d'étapes pour obtenir un pas inférieur à {2.0 * tau_min:.3g} s."
                )
    if float(settings["suspended_mass_g"]) <= 0.0:
        return "La masse suspendue doit être strictement positive."
    return None


def settings_error(settings: dict[str, SettingValue]) -> str | None:
    return geometry_error(settings) or material_error(settings) or numerical_error(settings)


def derived_geometry(settings: dict[str, SettingValue]) -> dict[str, float]:
    rout = float(settings["rout_mm"])
    rin = float(settings["rin_mm"])
    rho0 = float(settings["rho0_mm"])
    alpha = np.deg2rad(float(settings["alpha0_deg"]))
    eps = float(settings["eps"])
    initial_length = float(settings["initial_length_mm"])

    h0 = rho0 * np.tan(alpha)
    pitch0 = 2.0 * np.pi * h0
    turns = initial_length / pitch0
    centerline_length = initial_length / max(np.sin(alpha), 1e-12)
    wall_area = np.pi * (rout**2 - rin**2)
    inner_area = np.pi * rin**2
    nylon_area = np.pi * (0.5 * float(settings["nylon_diameter_mm"])) ** 2
    tube_internal_volume_ml = inner_area * centerline_length * 1.0e-3
    prestrained_length = (1.0 + eps) * initial_length
    prestrained_pitch = (1.0 + eps) * pitch0

    return {
        "h0_mm_per_rad": h0,
        "pitch0_mm": pitch0,
        "turns": turns,
        "centerline_length_mm": centerline_length,
        "wall_area_mm2": wall_area,
        "inner_area_mm2": inner_area,
        "nylon_area_mm2": nylon_area,
        "nylon_fill_ratio": nylon_area / inner_area,
        "tube_internal_volume_ml": tube_internal_volume_ml,
        "prestrained_length_mm": prestrained_length,
        "prestrained_pitch_mm": prestrained_pitch,
        "spring_index": rho0 / rout,
        "equivalent_mandrel_diameter_mm": max(0.0, 2.0 * (rho0 - rout)),
    }


def with_scaled_nylon(params: SimulationParams) -> SimulationParams:
    mat = params.mat
    scale = params.nylon_stiffness_scale
    return replace(
        params,
        mat=replace(
            mat,
            E_nylon=mat.E_nylon * scale,
            G_nylon=mat.G_nylon * scale,
        ),
    )


def build_config(settings: dict[str, SettingValue]) -> SimulationParams:
    maxwell = MaxwellTensileParams(
        E0=float(settings["maxwell_E0_mpa"]),
        E1=float(settings["maxwell_E1_mpa"]),
        eta1=float(settings["maxwell_eta1_mpa_s"]),
        E2=float(settings["maxwell_E2_mpa"]),
        eta2=float(settings["maxwell_eta2_mpa_s"]),
        E3=float(settings["maxwell_E3_mpa"]),
        eta3=float(settings["maxwell_eta3_mpa_s"]),
    )
    constitutive_mode = str(settings.get("constitutive_mode", "generalized_maxwell"))
    if constitutive_mode == "instantaneous_elastic":
        maxwell = replace(maxwell, E0=maxwell.E_total, E1=0.0, E2=0.0, E3=0.0)
    elif constitutive_mode != "generalized_maxwell":
        raise ValueError("Mode constitutif inconnu.")
    axial_modulus_mode = str(settings.get("axial_modulus_mode", "maxwell_sum"))
    if axial_modulus_mode == "maxwell_sum":
        E_axial = maxwell.E_total
    elif axial_modulus_mode == "paper_table":
        E_axial = float(settings["E_axial_mpa"])
    else:
        raise ValueError("Convention de module axial inconnue.")
    mat = MaterialParams(
        E_axial=E_axial,
        E_radius=float(settings["E_radius_mpa"]),
        G12=float(settings["G12_mpa"]),
        nu12=float(settings["nu12"]),
        nu23=float(settings["nu23"]),
        maxwell=maxwell,
        E_nylon=float(settings["E_nylon_mpa"]),
        G_nylon=float(settings["G_nylon_mpa"]),
        maxwell_anisotropy_mode=str(settings.get("maxwell_anisotropy_mode", "paper_equal")),
        nylon_condition_mode=str(settings.get("nylon_condition_mode", "bonded_linear")),
        nylon_axial_prestrain_coupling=(
            0.0 if str(settings.get("nylon_condition_mode")) == "axially_sliding_confined" else 1.0
        ),
        nylon_axial_actuation_coupling=(
            0.0 if str(settings.get("nylon_condition_mode")) == "axially_sliding_confined" else 1.0
        ),
    )
    geom = GeometryParams(
        Rout=float(settings["rout_mm"]),
        Rin=float(settings["rin_mm"]),
        r_nylon=0.5 * float(settings["nylon_diameter_mm"]),
        rho0=float(settings["rho0_mm"]),
        alpha0_deg=float(settings["alpha0_deg"]),
        theta_f_deg=float(settings["theta_f_deg"]),
        initial_length=float(settings["initial_length_mm"]),
        bias_angle_profile=str(settings.get("bias_angle_profile", "paper_linear")),
        section_update_mode=str(settings.get("section_update_mode", "fixed")),
    )
    duration_s = float(settings["duration_s"]) if bool(settings["use_fixed_duration"]) else None
    return with_scaled_nylon(
        SimulationParams(
            eps=float(settings["eps"]),
            n_cycles=int(settings["n_cycles"]),
            duration_s=duration_s,
            Pmax=float(settings["p_max_mpa"]),
            dt=float(settings["dt"]),
            n_layers=int(settings["n_layers"]),
            n_phi=int(settings["n_phi"]),
            pre_steps=int(settings["pre_steps"]),
            integration=str(settings["integration"]),
            prestrain_reference_mode=str(settings.get("prestrain_reference_mode", "elastic_tk_reference")),
            constitutive_mode=constitutive_mode,
            axial_modulus_mode=axial_modulus_mode,
            flow_rate_mL_min=float(settings["flow_rate_mL_min"]),
            volume_mL=float(settings["volume_mL"]),
            nonlinear_pressure=bool(settings["nonlinear_pressure"]),
            nylon_stiffness_scale=float(settings["nylon_scale"]),
            mat=mat,
            geom=geom,
        )
    )


def cycle_period_seconds(config: SimulationParams) -> float:
    if config.duration_s is not None:
        return config.duration_s / config.n_cycles
    return 2.0 * 60.0 * config.volume_mL / config.flow_rate_mL_min


def effective_flow_rate_mL_min(config: SimulationParams) -> float:
    period = cycle_period_seconds(config)
    return 2.0 * 60.0 * config.volume_mL / period


def make_pressure_history(config: SimulationParams) -> tuple[np.ndarray | None, np.ndarray | None]:
    if config.duration_s is None:
        return None, None
    if config.duration_s <= 0.0:
        raise ValueError("duration_s must be positive or None.")
    if config.n_cycles <= 0:
        raise ValueError("n_cycles must be positive.")

    period = cycle_period_seconds(config)
    half_period = 0.5 * period
    regular = np.arange(0.0, config.duration_s, config.dt, dtype=float)
    transitions = np.arange(0.0, config.duration_s + 0.5 * half_period, half_period, dtype=float)
    t = np.unique(np.concatenate((regular, transitions, np.array([0.0, config.duration_s]))))
    t = t[(t >= 0.0) & (t <= config.duration_s)]
    phase = (t % period) / period
    loading = phase < 0.5
    pressure = np.zeros_like(t)

    if config.nonlinear_pressure:
        gamma_load = 3.5
        gamma_unload = 2.8
        x = phase[loading] / 0.5
        y = (phase[~loading] - 0.5) / 0.5
        pressure[loading] = config.Pmax * x**gamma_load
        pressure[~loading] = config.Pmax * (1.0 - y) ** gamma_unload
    else:
        x = phase[loading] / 0.5
        y = (phase[~loading] - 0.5) / 0.5
        pressure[loading] = config.Pmax * x
        pressure[~loading] = config.Pmax * (1.0 - y)

    pressure[np.isclose(t, config.duration_s, rtol=0.0, atol=1e-12)] = 0.0
    return t, np.clip(pressure, 0.0, config.Pmax)
