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

CACHE_DIR = Path(tempfile.gettempdir()) / "tcpa_livrable_streamlit_cache"
SETTINGS_PATH = CACHE_DIR / "cavatappi_livrable_settings.json"
TIMING_PROFILE_PATH = CACHE_DIR / "cavatappi_livrable_timing_profile.json"
BLOCKED_RESULT_PATH = CACHE_DIR / "cavatappi_livrable_blocked.pkl"
RELAXATION_RESULT_PATH = CACHE_DIR / "cavatappi_livrable_relaxation.pkl"
PRESTRAIN_RESULT_PATH = CACHE_DIR / "cavatappi_livrable_prestrain.pkl"
SUSPENDED_RESULT_PATH = CACHE_DIR / "cavatappi_livrable_suspended.pkl"
RESULT_CACHE_PATHS = (
    BLOCKED_RESULT_PATH,
    RELAXATION_RESULT_PATH,
    PRESTRAIN_RESULT_PATH,
    SUSPENDED_RESULT_PATH,
)
SETTINGS_SCHEMA_VERSION = 5

PRESSURE_END_FORCE_OPTIONS = ["none", "projected_inner_area", "axial_inner_area"]
INTEGRATION_OPTIONS = ["paper_incremental", "paper_explicit", "exponential"]
VISUAL_STATE_OPTIONS = ["fabricated", "prestrained"]

PRESSURE_END_FORCE_LABELS = {
    "none": "Aucune",
    "projected_inner_area": "Poussée projetée sur l'axe",
    "axial_inner_area": "Poussée axiale complète",
}
INTEGRATION_LABELS = {
    "paper_incremental": "Incrémentale article",
    "paper_explicit": "Explicite article",
    "exponential": "Exponentielle stable",
}
VISUAL_STATE_LABELS = {
    "fabricated": "Fabriqué",
    "prestrained": "Précontraint",
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
    "pressure_end_force_mode": "none",
    "pressure_end_force_scale": 0.0,
    "n_cycles": 11,
    "hysteresis_cycle": 1,
    "hysteresis_cycles": "1",
    "hysteresis_compare_mode": "current",
    "hysteresis_prestrain_values": "0.6, 0.8, 1.0",
    "hysteresis_pressure_rates_mpa_s": "0.05, 0.10, 0.20",
    "use_fixed_duration": True,
    "duration_s": 500.0,
    "flow_rate_mL_min": 10.0,
    "volume_mL": 1.50,
    "nonlinear_pressure": True,
    "relaxation_ramp_time_s": 9.0,
    "relaxation_hold_time_s": 300.0,
    "suspended_mass_g": 100.0,
    "suspended_duration_s": 120.0,
    "suspended_pressure_rate_mpa_s": 0.10,
    "suspended_hold_pressure": False,
    "suspended_show_geometry_plot": False,
    "dt": 0.5,
    "pre_steps": 24,
    "integration": "paper_incremental",
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
    "nylon_axial_actuation_coupling": 0.0,
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
    nylon_axial_prestrain_coupling: float = 1.0
    nylon_axial_actuation_coupling: float = 0.0


@dataclass
class GeometryParams:
    Rout: float = 1.0
    Rin: float = 0.4
    r_nylon: float = 0.77 / 2.0
    rho0: float = 2.16
    alpha0_deg: float = 10.53
    theta_f_deg: float = 37.91
    initial_length: float = 32.45
    pressure_end_force_mode: str = "none"
    pressure_end_force_scale: float = 0.0


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
    integration: str = "paper_incremental"
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
    for key in ("use_fixed_duration", "nonlinear_pressure", "suspended_hold_pressure", "suspended_show_geometry_plot"):
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
    except (FileNotFoundError, OSError, pickle.PickleError, AttributeError, EOFError):
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
    if rin >= rout:
        return "Rin doit être inférieur à Rout."
    if nylon_diameter > 2.0 * rin:
        return "Le diamètre du nylon ne doit pas dépasser le diamètre intérieur du tube."
    if rho0 <= rout:
        return "rho0 doit être supérieur à Rout pour une ligne centrale hélicoïdale."
    return None


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


def with_calibrated_material(params: SimulationParams) -> SimulationParams:
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
    mat = MaterialParams(
        E_axial=float(settings["E_axial_mpa"]),
        E_radius=float(settings["E_radius_mpa"]),
        G12=float(settings["G12_mpa"]),
        nu12=float(settings["nu12"]),
        nu23=float(settings["nu23"]),
        maxwell=maxwell,
        E_nylon=float(settings["E_nylon_mpa"]),
        G_nylon=float(settings["G_nylon_mpa"]),
        nylon_axial_prestrain_coupling=float(settings["nylon_axial_prestrain_coupling"]),
        nylon_axial_actuation_coupling=float(settings["nylon_axial_actuation_coupling"]),
    )
    geom = GeometryParams(
        Rout=float(settings["rout_mm"]),
        Rin=float(settings["rin_mm"]),
        r_nylon=0.5 * float(settings["nylon_diameter_mm"]),
        rho0=float(settings["rho0_mm"]),
        alpha0_deg=float(settings["alpha0_deg"]),
        theta_f_deg=float(settings["theta_f_deg"]),
        initial_length=float(settings["initial_length_mm"]),
        pressure_end_force_mode=str(settings["pressure_end_force_mode"]),
        pressure_end_force_scale=float(settings["pressure_end_force_scale"]),
    )
    duration_s = float(settings["duration_s"]) if bool(settings["use_fixed_duration"]) else None
    return with_calibrated_material(
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

    t = np.arange(0.0, config.duration_s, config.dt)
    if t.size == 0 or t[-1] < config.duration_s:
        t = np.append(t, config.duration_s)

    period = cycle_period_seconds(config)
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

    return t, pressure
