from __future__ import annotations

import hashlib
import sys
import re
import time
import os
import inspect
import importlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) in sys.path:
    sys.path.remove(str(APP_DIR))
sys.path.insert(0, str(APP_DIR))

import Base as modele
import affichage as affichage_module
import parametres as parametres_module

# Streamlit conserve les modules importés entre deux exécutions. Ce rechargement
# garantit que les libellés, valeurs par défaut et validations restent synchronisés.
importlib.invalidate_caches()
parametres_module = importlib.reload(parametres_module)
affichage_module = importlib.reload(affichage_module)

from affichage import (
    format_seconds,
    hysteresis_to_csv_bytes,
    make_cavatappi_figure,
    make_cross_section_figure,
    mapping_to_csv_bytes,
    plot_hysteresis_overlay,
    plot_prestrain_study,
    plot_relaxation_response,
    plot_suspended_response,
    plot_time_response_fr,
)
from parametres import (
    BLOCKED_RESULT_PATH,
    AXIAL_MODULUS_LABELS,
    AXIAL_MODULUS_OPTIONS,
    BIAS_ANGLE_PROFILE_LABELS,
    BIAS_ANGLE_PROFILE_OPTIONS,
    CONSTITUTIVE_LABELS,
    CONSTITUTIVE_OPTIONS,
    DEFAULT_SETTINGS,
    HYSTERESIS_RESULT_PATH,
    INTEGRATION_LABELS,
    INTEGRATION_OPTIONS,
    MAXWELL_ANISOTROPY_LABELS,
    MAXWELL_ANISOTROPY_OPTIONS,
    NYLON_CONDITION_LABELS,
    NYLON_CONDITION_OPTIONS,
    PRESTRAIN_REFERENCE_LABELS,
    PRESTRAIN_REFERENCE_OPTIONS,
    PRESSURE_INPUT_LABELS,
    PRESSURE_INPUT_OPTIONS,
    PRESSURE_END_FORCE_LABELS,
    PRESSURE_END_FORCE_OPTIONS,
    PRESTRAIN_RESULT_PATH,
    RELAXATION_RESULT_PATH,
    SECTION_UPDATE_LABELS,
    SECTION_UPDATE_OPTIONS,
    SettingValue,
    SETTINGS_SCHEMA_VERSION,
    SUSPENDED_RESULT_PATH,
    TIMING_PROFILE_PATH,
    VISUAL_STATE_LABELS,
    VISUAL_STATE_OPTIONS,
    build_config,
    cycle_period_seconds,
    derived_geometry,
    settings_error,
    load_result_cache,
    load_settings,
    make_pressure_history,
    option_index,
    reset_settings_and_cache,
    save_result_cache,
    save_settings,
)
from parallel import (
    run_hysteresis_pressure_rate_case,
    run_hysteresis_prestrain_case,
    run_prestrain_case,
)
from pression import (
    estimate_measured_period,
    infer_column,
    infer_pressure_unit,
    measured_pressure_payload,
    parse_uploaded_numeric_csv,
)
from timing import (
    estimate_compute_seconds,
    load_timing_profile,
    model_cost_index,
    record_timing_sample,
)


HYSTERESIS_COMPARE_OPTIONS = ["current", "prestrain", "pressure_rate"]
INTERFACE_DEFAULT_PARALLEL_WORKERS = max(1, min(4, (os.cpu_count() or 2) - 1))
INTERFACE_DEFAULT_SUSPENDED_SHOW_GEOMETRY_PLOT = False
HYSTERESIS_COMPARE_LABELS = {
    "current": "Cycles du calcul courant",
    "prestrain": "Plusieurs précontraintes initiales",
    "pressure_rate": "Plusieurs vitesses de pression injectée",
}


def render_csv_download(payload: bytes | None, filename: str, key: str) -> None:
    st.download_button(
        "Exporter les résultats en CSV",
        data=payload if payload is not None else b"",
        file_name=filename,
        mime="text/csv; charset=utf-8",
        key=key,
        on_click="ignore",
        icon=":material/download:",
        disabled=payload is None,
    )

BLOCKED_RESULT_IGNORE_KEYS = {
    "eps_study_min",
    "eps_study_max",
    "eps_study_points",
    "hysteresis_cycle",
    "hysteresis_cycles",
    "hysteresis_compare_mode",
    "hysteresis_prestrain_values",
    "hysteresis_pressure_rates_mpa_s",
    "relaxation_ramp_time_s",
    "relaxation_hold_time_s",
    "suspended_mass_g",
    "suspended_duration_s",
    "suspended_pressure_rate_mpa_s",
    "suspended_hold_pressure",
    "suspended_show_geometry_plot",
    "parallel_workers",
    "view_elev_deg",
    "view_azim_deg",
}
RELAXATION_RESULT_IGNORE_KEYS = {
    "eps_study_min",
    "eps_study_max",
    "eps_study_points",
    "hysteresis_cycle",
    "hysteresis_cycles",
    "hysteresis_compare_mode",
    "hysteresis_prestrain_values",
    "hysteresis_pressure_rates_mpa_s",
    "n_cycles",
    "use_fixed_duration",
    "duration_s",
    "flow_rate_mL_min",
    "volume_mL",
    "suspended_mass_g",
    "suspended_duration_s",
    "suspended_pressure_rate_mpa_s",
    "suspended_hold_pressure",
    "suspended_show_geometry_plot",
    "parallel_workers",
    "view_elev_deg",
    "view_azim_deg",
}
PRESTRAIN_RESULT_IGNORE_KEYS = {
    "eps",
    "hysteresis_cycle",
    "hysteresis_cycles",
    "hysteresis_compare_mode",
    "hysteresis_prestrain_values",
    "hysteresis_pressure_rates_mpa_s",
    "relaxation_ramp_time_s",
    "relaxation_hold_time_s",
    "suspended_mass_g",
    "suspended_duration_s",
    "suspended_pressure_rate_mpa_s",
    "suspended_hold_pressure",
    "suspended_show_geometry_plot",
    "parallel_workers",
    "view_elev_deg",
    "view_azim_deg",
}
SUSPENDED_RESULT_IGNORE_KEYS = {
    "eps_study_min",
    "eps_study_max",
    "eps_study_points",
    "hysteresis_cycle",
    "hysteresis_cycles",
    "hysteresis_compare_mode",
    "hysteresis_prestrain_values",
    "hysteresis_pressure_rates_mpa_s",
    "relaxation_ramp_time_s",
    "relaxation_hold_time_s",
    "n_cycles",
    "use_fixed_duration",
    "duration_s",
    "flow_rate_mL_min",
    "volume_mL",
    "nonlinear_pressure",
    "suspended_show_geometry_plot",
    "parallel_workers",
    "view_elev_deg",
    "view_azim_deg",
}
HYSTERESIS_COMPARE_IGNORE_KEYS = {
    "eps_study_min",
    "eps_study_max",
    "eps_study_points",
    "hysteresis_cycle",
    "hysteresis_cycles",
    "hysteresis_compare_mode",
    "hysteresis_prestrain_values",
    "hysteresis_pressure_rates_mpa_s",
    "relaxation_ramp_time_s",
    "relaxation_hold_time_s",
    "suspended_mass_g",
    "suspended_duration_s",
    "suspended_pressure_rate_mpa_s",
    "suspended_hold_pressure",
    "suspended_show_geometry_plot",
    "parallel_workers",
    "view_elev_deg",
    "view_azim_deg",
}


def _numbers_from_text(text: object) -> list[float]:
    values = []
    for item in re.findall(r"-?\d+(?:[\.,]\d+)?", str(text)):
        try:
            values.append(float(item.replace(",", ".")))
        except ValueError:
            continue
    return values


def parse_cycle_list(text: object, max_cycle: int) -> list[int]:
    cycles = []
    for value in _numbers_from_text(text):
        cycle = int(round(value))
        if 1 <= cycle <= max_cycle and cycle not in cycles:
            cycles.append(cycle)
    return cycles or [1]


def parse_positive_float_list(text: object, max_count: int = 8) -> list[float]:
    values = []
    for value in _numbers_from_text(text):
        if value > 0.0 and value not in values:
            values.append(value)
        if len(values) >= max_count:
            break
    return values


def settings_signature(settings: dict[str, SettingValue], ignore_keys: set[str] | None = None) -> tuple[tuple[str, str], ...]:
    ignored = ignore_keys or set()
    settings_items = tuple(
        sorted(
            (key, str(value))
            for key, value in settings.items()
            if key not in ignored and not key.startswith("_")
        )
    )
    return (("_model_version", str(modele.MODEL_VERSION)),) + settings_items


def result_matches_settings(
    result: dict | None,
    signature: tuple[tuple[str, str], ...],
    ignore_keys: set[str] | None = None,
) -> bool:
    if result is None:
        return False
    if result.get("signature") == signature:
        return True
    saved_settings = result.get("settings")
    if isinstance(saved_settings, dict):
        return settings_signature(saved_settings, ignore_keys) == signature
    return False


def make_pressure_rate_history(config, rate_mpa_s: float, n_cycles_for_history: int):
    if rate_mpa_s <= 0.0:
        raise ValueError("La vitesse de pression doit être positive.")
    if config.Pmax <= 0.0:
        raise ValueError("La pression maximale doit être positive.")

    half_period = config.Pmax / rate_mpa_s
    period = 2.0 * half_period
    total_time = max(1, int(n_cycles_for_history)) * period
    regular_time = np.arange(0.0, total_time + config.dt, config.dt)
    transition_time = np.arange(0.0, total_time + half_period, half_period)
    time_values = np.unique(np.concatenate((regular_time, transition_time, [total_time])))
    time_values = time_values[(time_values >= 0.0) & (time_values <= total_time)]

    phase = (time_values % period) / period
    pressure = np.empty_like(time_values)
    loading = phase <= 0.5
    pressure[loading] = config.Pmax * (phase[loading] / 0.5)
    pressure[~loading] = config.Pmax * (1.0 - (phase[~loading] - 0.5) / 0.5)
    pressure[np.isclose(time_values, total_time)] = 0.0
    return time_values, np.clip(pressure, 0.0, config.Pmax), period


def apply_view_query_params(settings: dict[str, SettingValue]) -> dict[str, SettingValue]:
    for key, lower, upper in (
        ("view_elev_deg", 0.0, 90.0),
        ("view_azim_deg", -180.0, 180.0),
    ):
        try:
            raw = st.query_params.get(key)
        except Exception:
            raw = None
        if isinstance(raw, list):
            raw = raw[0] if raw else None
        if raw is None:
            continue
        try:
            settings[key] = float(np.clip(float(raw), lower, upper))
        except (TypeError, ValueError):
            continue
    return settings


def inject_keyboard_view_controls(settings: dict[str, SettingValue]) -> None:
    elev = float(settings.get("view_elev_deg", 22.0))
    azim = float(settings.get("view_azim_deg", -58.0))
    renderer = st.components.v2.component(
        "cavatappi_keyboard_view",
        html="<span aria-hidden='true'></span>",
        js="""
        export default function(component) {
            const state = {
                elev: Number(component.data.elev),
                azim: Number(component.data.azim)
            };
            const clamp = (value, minValue, maxValue) => Math.min(maxValue, Math.max(minValue, value));
            const wrapAzim = (value) => {
                const wrapped = ((value + 180) % 360 + 360) % 360 - 180;
                return wrapped === -180 ? 180 : wrapped;
            };
            const handler = (event) => {
                const tagName = event.target?.tagName?.toUpperCase() || "";
                if (["INPUT", "TEXTAREA", "SELECT"].includes(tagName) || event.ctrlKey || event.metaKey || event.altKey) return;
                let changed = true;
                if (event.key === "ArrowLeft") state.azim = wrapAzim(state.azim - 5);
                else if (event.key === "ArrowRight") state.azim = wrapAzim(state.azim + 5);
                else if (event.key === "ArrowUp") state.elev = clamp(state.elev + 3, 0, 90);
                else if (event.key === "ArrowDown") state.elev = clamp(state.elev - 3, 0, 90);
                else changed = false;
                if (!changed) return;
                event.preventDefault();
                const url = new URL(window.location.href);
                url.searchParams.set("view_elev_deg", state.elev.toFixed(1));
                url.searchParams.set("view_azim_deg", state.azim.toFixed(1));
                window.location.href = url.toString();
            };
            document.addEventListener("keydown", handler, true);
            return () => document.removeEventListener("keydown", handler, true);
        }
        """,
        isolate_styles=False,
    )
    renderer(key="keyboard_view", data={"elev": elev, "azim": azim}, width=1, height=1)


def make_suspended_response_figure(data: dict[str, np.ndarray], show_geometry: bool):
    fresh_affichage = importlib.reload(affichage_module)
    signature = inspect.signature(fresh_affichage.plot_suspended_response)
    if "show_geometry" in signature.parameters:
        return fresh_affichage.plot_suspended_response(data, show_geometry=show_geometry)
    return fresh_affichage.plot_suspended_response(data)


def run_with_progress(label: str, estimated_seconds: float, function, *args):
    progress = st.progress(0, text=f"{label} : démarrage")
    status = st.empty()
    start = time.perf_counter()
    estimated_seconds = max(0.1, float(estimated_seconds))

    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(function, *args)
        while not future.done():
            elapsed = time.perf_counter() - start
            fraction = min(0.95, 0.95 * elapsed / estimated_seconds)
            remaining = max(0.0, estimated_seconds - elapsed)
            progress.progress(
                int(100 * fraction),
                text=(
                    f"{label} : {format_seconds(elapsed)} écoulées, "
                    f"total estimé {format_seconds(estimated_seconds)}, "
                    f"~{format_seconds(remaining)} restantes"
                ),
            )
            time.sleep(0.2)

        try:
            result = future.result()
        except Exception as exc:
            progress.empty()
            status.error(f"{label} interrompu : {exc}")
            st.stop()

    elapsed = time.perf_counter() - start
    progress.progress(100, text=f"{label} : terminé en {format_seconds(elapsed)}")
    status.caption(f"Estimé : {format_seconds(estimated_seconds)} | réel : {format_seconds(elapsed)}")
    return result, elapsed


def run_model(settings: dict[str, SettingValue], measured_history=None):
    config = build_config(settings)
    if measured_history is None:
        pressure_time, pressure_mpa = make_pressure_history(config)
    else:
        pressure_time = np.asarray(measured_history["time"], dtype=float)
        pressure_mpa = np.asarray(measured_history["pressure_MPa"], dtype=float)
        measured_period = estimate_measured_period(pressure_time, pressure_mpa)
        config.measured_cycle_period_s = measured_period if measured_period is not None else float(pressure_time[-1])
        if measured_period is not None:
            config.n_cycles = max(1, int(round(float(pressure_time[-1]) / measured_period)))
        else:
            config.n_cycles = 1
    _, data = modele.run_blocked_actuation(config, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    return config, data, modele.summary(data)


def result_cycle_period(config) -> float:
    measured_period = getattr(config, "measured_cycle_period_s", None)
    return float(measured_period) if measured_period is not None else cycle_period_seconds(config)


def run_relaxation_model(settings: dict[str, SettingValue]):
    config = build_config(settings)
    ramp_time = float(settings["relaxation_ramp_time_s"])
    hold_time = float(settings["relaxation_hold_time_s"])
    pressure_time, pressure_mpa = modele.ramp_hold_pressure_history(
        P_hold=config.Pmax,
        ramp_time=ramp_time,
        hold_time=hold_time,
        dt=config.dt,
        nonlinear_ramp=bool(config.nonlinear_pressure),
    )
    _, data = modele.run_blocked_actuation(config, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    hold_start_index = int(np.searchsorted(data["time"], ramp_time, side="left"))
    hold_start_index = min(max(hold_start_index, 0), len(data["time"]) - 1)
    data["ramp_time"] = float(ramp_time)
    data["hold_time"] = float(hold_time)
    data["hold_start_index"] = hold_start_index
    data["force_hold_relax_mN"] = data["force_total_mN"] - data["force_total_mN"][hold_start_index]
    data["torque_hold_relax_microNm"] = data["torque_act_microNm"] - data["torque_act_microNm"][hold_start_index]
    return config, data, modele.summary(data)


def make_suspended_pressure_history(config, settings: dict[str, SettingValue]):
    duration = float(settings["suspended_duration_s"])
    rate = float(settings["suspended_pressure_rate_mpa_s"])
    if duration <= 0.0:
        raise ValueError("La durée de simulation masse suspendue doit être positive.")
    if rate <= 0.0:
        raise ValueError("La vitesse d'actionnement masse suspendue doit être positive.")

    ramp_time = config.Pmax / rate if config.Pmax > 0.0 else 0.0
    transition_times = [0.0, duration, ramp_time]
    if not bool(settings["suspended_hold_pressure"]):
        transition_times.append(2.0 * ramp_time)

    regular_time = np.arange(0.0, duration + config.dt, config.dt)
    time_values = np.unique(np.concatenate((regular_time, np.asarray(transition_times, dtype=float))))
    time_values = time_values[(time_values >= 0.0) & (time_values <= duration)]
    if time_values.size < 2:
        time_values = np.array([0.0, duration], dtype=float)

    if bool(settings["suspended_hold_pressure"]):
        pressure = np.minimum(rate * time_values, config.Pmax)
    else:
        pressure = np.where(
            time_values <= ramp_time,
            rate * time_values,
            np.maximum(config.Pmax - rate * (time_values - ramp_time), 0.0),
        )

    pressure = np.clip(pressure, 0.0, config.Pmax)
    hold_start_time = min(ramp_time, duration)
    return time_values, pressure, hold_start_time


def run_suspended_model(settings: dict[str, SettingValue]):
    config = build_config(settings)
    pressure_time, pressure_mpa, hold_start_time = make_suspended_pressure_history(config, settings)
    load_N = float(settings["suspended_mass_g"]) * 1.0e-3 * 9.80665
    _, data = modele.run_suspended_actuation(config, load_N=load_N, pressure_time=pressure_time, pressure_MPa=pressure_mpa)
    hold_start_index = int(np.searchsorted(data["time"], hold_start_time, side="left"))
    hold_start_index = min(max(hold_start_index, 0), len(data["time"]) - 1)
    data["suspended_duration_s"] = float(settings["suspended_duration_s"])
    data["suspended_pressure_rate_mpa_s"] = float(settings["suspended_pressure_rate_mpa_s"])
    data["suspended_hold_pressure"] = bool(settings["suspended_hold_pressure"])
    data["suspended_hold_start_time_s"] = float(hold_start_time)
    data["suspended_hold_start_index"] = hold_start_index
    data["free_contraction_hold_relax_mm"] = data["free_contraction_mm"] - data["free_contraction_mm"][hold_start_index]
    data["free_actuation_hold_relax_percent"] = (
        data["free_actuation_percent"] - data["free_actuation_percent"][hold_start_index]
    )
    data["axial_length_hold_relax_mm"] = data["axial_length_mm"] - data["axial_length_mm"][hold_start_index]
    return config, data, modele.summary(data)


def run_prestrain_study_with_progress(
    settings: dict[str, SettingValue],
    eps_values: np.ndarray,
    estimated_seconds: float,
    parallel_workers: int = 1,
):
    eps_values = np.asarray(eps_values, dtype=float)
    total = len(eps_values)
    progress = st.progress(0, text="Étude de précontrainte : démarrage")
    status = st.empty()
    start = time.perf_counter()
    worker_count = max(1, min(int(parallel_workers), total))

    def update_progress(done: int, label: str) -> None:
        elapsed = time.perf_counter() - start
        if done > 0:
            avg = elapsed / float(done)
            remaining = avg * float(total - done)
        else:
            remaining = max(0.0, estimated_seconds - elapsed)
        progress.progress(
            int(100 * done / max(1, total)),
            text=(
                f"Étude de précontrainte : {label} "
                f"({done}/{total}), ~{format_seconds(remaining)} restantes"
            ),
        )

    def run_sequential() -> list[dict[str, float]]:
        rows_seq = []
        for index, eps_value in enumerate(eps_values, start=1):
            update_progress(index - 1, f"précontrainte = {eps_value:.3f}")
            rows_seq.append(run_prestrain_case(dict(settings), float(eps_value)))
        return rows_seq

    if worker_count <= 1:
        rows = run_sequential()
    else:
        try:
            rows = []
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(run_prestrain_case, dict(settings), float(eps_value)): float(eps_value)
                    for eps_value in eps_values
                }
                for done, future in enumerate(as_completed(futures), start=1):
                    eps_value = futures[future]
                    rows.append(future.result())
                    update_progress(done, f"précontrainte = {eps_value:.3f} | {worker_count} cœurs")
            rows.sort(key=lambda row: row["eps"])
        except Exception as exc:
            st.warning(f"Calcul parallèle indisponible ({exc}). Repli en calcul séquentiel.")
            rows = run_sequential()

    elapsed = time.perf_counter() - start
    progress.progress(100, text=f"Étude de précontrainte : terminée en {format_seconds(elapsed)}")
    status.caption(
        f"Estimé : {format_seconds(estimated_seconds)} | réel : {format_seconds(elapsed)} | "
        f"cœurs utilisés : {worker_count}"
    )
    return {key: np.array([row[key] for row in rows], dtype=float) for key in rows[0]}, elapsed


def run_hysteresis_comparison_with_progress(
    settings: dict[str, SettingValue],
    mode: str,
    values: list[float],
    cycles: list[int],
    estimated_seconds: float,
    parallel_workers: int = 1,
):
    progress = st.progress(0, text="Comparaison d'hystérèse : démarrage")
    status = st.empty()
    start = time.perf_counter()
    total = max(1, len(values))
    max_cycle = max(cycles) if cycles else 1
    worker_count = max(1, min(int(parallel_workers), total))

    def update_progress(done: int, label: str) -> None:
        elapsed = time.perf_counter() - start
        if done > 0:
            avg = elapsed / float(done)
            remaining = avg * float(total - done)
        else:
            remaining = max(0.0, estimated_seconds - elapsed)
        progress.progress(
            int(100 * done / total),
            text=(
                f"Comparaison d'hystérèse : {label} "
                f"({done}/{total}), ~{format_seconds(remaining)} restantes"
            ),
        )

    def run_one(value: float) -> dict:
        if mode == "prestrain":
            return run_hysteresis_prestrain_case(dict(settings), float(value))
        if mode == "pressure_rate":
            return run_hysteresis_pressure_rate_case(dict(settings), float(value), max_cycle)
        raise ValueError("Mode de comparaison d'hystérèse inconnu.")

    if worker_count <= 1:
        cases = []
        for index, value in enumerate(values, start=1):
            label = f"précontrainte = {value:.3f}" if mode == "prestrain" else f"vitesse = {value:.3f} MPa/s"
            update_progress(index - 1, label)
            cases.append(run_one(float(value)))
    else:
        try:
            cases = []
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                if mode == "prestrain":
                    futures = {
                        executor.submit(run_hysteresis_prestrain_case, dict(settings), float(value)): float(value)
                        for value in values
                    }
                elif mode == "pressure_rate":
                    futures = {
                        executor.submit(
                            run_hysteresis_pressure_rate_case,
                            dict(settings),
                            float(value),
                            max_cycle,
                        ): float(value)
                        for value in values
                    }
                else:
                    raise ValueError("Mode de comparaison d'hystérèse inconnu.")
                for done, future in enumerate(as_completed(futures), start=1):
                    value = futures[future]
                    cases.append(future.result())
                    label = f"précontrainte = {value:.3f}" if mode == "prestrain" else f"vitesse = {value:.3f} MPa/s"
                    update_progress(done, f"{label} | {worker_count} cœurs")
            cases.sort(key=lambda case: float(case.get("value", 0.0)))
        except Exception as exc:
            st.warning(f"Calcul parallèle indisponible ({exc}). Repli en calcul séquentiel.")
            cases = []
            for index, value in enumerate(values, start=1):
                label = f"précontrainte = {value:.3f}" if mode == "prestrain" else f"vitesse = {value:.3f} MPa/s"
                update_progress(index - 1, label)
                cases.append(run_one(float(value)))

    elapsed = time.perf_counter() - start
    progress.progress(100, text=f"Comparaison d'hystérèse : terminée en {format_seconds(elapsed)}")
    status.caption(
        f"Estimé : {format_seconds(estimated_seconds)} | réel : {format_seconds(elapsed)} | "
        f"cœurs utilisés : {worker_count}"
    )
    return {"mode": mode, "values": values, "cycles": cycles, "cases": cases}, elapsed


def initialize_cached_results() -> None:
    cache_paths = {
        "calculator_result": BLOCKED_RESULT_PATH,
        "relaxation_result": RELAXATION_RESULT_PATH,
        "prestrain_study_result": PRESTRAIN_RESULT_PATH,
        "suspended_result": SUSPENDED_RESULT_PATH,
        "hysteresis_comparison_result": HYSTERESIS_RESULT_PATH,
    }
    for key, path in cache_paths.items():
        if key not in st.session_state:
            st.session_state[key] = load_result_cache(path)


def sidebar_help(lines: list[str]) -> None:
    with st.popover("Aide"):
        st.markdown("\n".join(f"- {line}" for line in lines))


settings = dict(DEFAULT_SETTINGS)
settings.update(load_settings())
settings.setdefault("parallel_workers", INTERFACE_DEFAULT_PARALLEL_WORKERS)
settings.setdefault("suspended_show_geometry_plot", INTERFACE_DEFAULT_SUSPENDED_SHOW_GEOMETRY_PLOT)
settings = apply_view_query_params(settings)
timing_profile = load_timing_profile(TIMING_PROFILE_PATH)

st.set_page_config(page_title="Calculateur Cavatappi - Beta", layout="wide")
st.title("Calculateur d'actionneur Cavatappi")
st.caption("Interface de simulation pour la géométrie, l'actionnement et la visualisation du modèle TCPA.")

st.sidebar.header("Paramètres d'entrée")
if st.sidebar.button("Paramètres par défaut", type="secondary", use_container_width=True):
    reset_settings_and_cache()
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    try:
        st.query_params.clear()
    except Exception:
        pass
    st.rerun()

with st.sidebar.expander("Géométrie", expanded=True):
    sidebar_help(
        [
            "Ces paramètres définissent la forme fabriquée de l'actionneur et la section tube/nylon.",
            "Rout et Rin règlent l'épaisseur du tube, donc sa raideur et la surface soumise à la pression.",
            "rho0, alpha0 et la longueur initiale définissent l'hélice de départ.",
            "theta_f est l'angle de biais des fibres/matière du tube utilisé pour l'anisotropie.",
            "La mise à jour hélicoïdale conserve la section initiale ; le mode évolutif actualise aussi les rayons et l'orientation du matériau.",
        ]
    )
    rout_mm = st.number_input("Rayon extérieur du tube Rout (mm)", 0.05, 5.0, float(settings["rout_mm"]), 0.05)
    rin_mm = st.number_input("Rayon intérieur du tube Rin (mm)", 0.01, 4.0, float(settings["rin_mm"]), 0.05)
    nylon_diameter_mm = st.number_input("Diamètre du nylon (mm)", 0.01, 4.0, float(settings["nylon_diameter_mm"]), 0.01)
    rho0_mm = st.number_input("Rayon de ligne centrale rho0 (mm)", 0.05, 10.0, float(settings["rho0_mm"]), 0.05)
    alpha0_deg = st.number_input("Angle hélicoïdal initial alpha0 (deg)", 0.1, 85.0, float(settings["alpha0_deg"]), 0.1)
    theta_f_deg = st.number_input("Angle de biais du tube theta_f (deg)", 0.0, 89.0, float(settings["theta_f_deg"]), 0.1)
    bias_angle_profile = st.selectbox(
        "Profil radial de l'angle de biais",
        BIAS_ANGLE_PROFILE_OPTIONS,
        index=option_index(BIAS_ANGLE_PROFILE_OPTIONS, settings.get("bias_angle_profile", "paper_linear")),
        format_func=lambda value: BIAS_ANGLE_PROFILE_LABELS.get(value, value),
        help=(
            "La variation linéaire fait évoluer l'angle de biais proportionnellement au rayon. "
            "La loi en tangente représente un taux de torsion uniforme dans le tube droit."
        ),
    )
    initial_length_mm = st.number_input(
        "Longueur initiale de l'actionneur (mm)", 1.0, 500.0, float(settings["initial_length_mm"]), 0.5
    )
    section_update_mode = st.selectbox(
        "Mise à jour de la section du tube",
        SECTION_UPDATE_OPTIONS,
        index=option_index(SECTION_UPDATE_OPTIONS, settings.get("section_update_mode", "fixed")),
        format_func=lambda value: SECTION_UPDATE_LABELS.get(value, value),
    )

uploaded_pressure_payload = None
measured_pressure_time_column = str(settings.get("measured_pressure_time_column", ""))
measured_pressure_column = str(settings.get("measured_pressure_column", ""))
measured_pressure_unit = str(settings.get("measured_pressure_unit", "MPa"))
measured_pressure_file_hash = str(settings.get("measured_pressure_file_hash", ""))
measured_pressure_subtract_initial = bool(settings.get("measured_pressure_subtract_initial", True))

with st.sidebar.expander("Pression et actionnement", expanded=True):
    sidebar_help(
        [
            "Ces paramètres définissent le chargement appliqué pendant l'actionnement bloqué.",
            "La précontrainte initiale étire l'actionneur avant l'injection de pression.",
            "La pression maximale fixe l'amplitude du cycle de pression.",
            "Le débit et le volume définissent uniquement la durée d'un demi-cycle : ils ne constituent pas un modèle hydraulique pression-volume.",
        ]
    )
    pressure_input_mode = st.selectbox(
        "Source de pression",
        PRESSURE_INPUT_OPTIONS,
        index=option_index(PRESSURE_INPUT_OPTIONS, settings.get("pressure_input_mode", "generated")),
        format_func=lambda value: PRESSURE_INPUT_LABELS.get(value, value),
    )
    eps = st.slider("Précontrainte initiale", 0.0, 1.5, float(settings["eps"]), 0.05)
    p_max_mpa = st.slider("Pression maximale (MPa)", 0.0, 1.5, min(float(settings["p_max_mpa"]), 1.5), 0.05)
    n_cycles = st.slider("Cycles", 1, 60, int(settings["n_cycles"]), 1)
    use_fixed_duration = st.checkbox("Utiliser une durée totale fixe", bool(settings["use_fixed_duration"]))
    duration_s = st.number_input("Durée totale (s)", 1.0, 5000.0, float(settings["duration_s"]), 10.0)
    flow_rate_mL_min = st.number_input("Débit (mL/min)", 0.01, 200.0, float(settings["flow_rate_mL_min"]), 0.5)
    volume_mL = st.number_input("Volume de commande par demi-cycle (mL)", 0.001, 100.0, float(settings["volume_mL"]), 0.05)
    nonlinear_pressure = st.checkbox(
        "Profil de pression phénoménologique non linéaire",
        bool(settings["nonlinear_pressure"]),
        help="Ce profil utilise des exposants empiriques. Pour un protocole contrôlé, conservez le profil linéaire.",
    )
    if use_fixed_duration:
        st.caption("La durée fixe remplace le débit demandé par un débit effectif calculé à partir du nombre de cycles.")
    if pressure_input_mode == "measured_csv":
        uploaded_pressure_file = st.file_uploader(
            "Historique de pression mesuré",
            type=["csv", "txt"],
            help="Le fichier doit contenir une colonne de temps et une colonne de pression.",
        )
        if uploaded_pressure_file is not None:
            raw_pressure_file = uploaded_pressure_file.getvalue()
            measured_pressure_file_hash = hashlib.sha256(raw_pressure_file).hexdigest()
            try:
                uploaded_columns = parse_uploaded_numeric_csv(raw_pressure_file)
                uploaded_headers = list(uploaded_columns)
                default_time = (
                    measured_pressure_time_column
                    if measured_pressure_time_column in uploaded_headers
                    else infer_column(uploaded_headers, ("temps", "time", "seconde", "second", " t"), 0)
                )
                default_pressure = (
                    measured_pressure_column
                    if measured_pressure_column in uploaded_headers
                    else infer_column(uploaded_headers, ("pression", "pressure", "press"), min(1, len(uploaded_headers) - 1))
                )
                measured_pressure_time_column = st.selectbox(
                    "Colonne de temps",
                    uploaded_headers,
                    index=option_index(uploaded_headers, default_time),
                )
                measured_pressure_column = st.selectbox(
                    "Colonne de pression",
                    uploaded_headers,
                    index=option_index(uploaded_headers, default_pressure),
                )
                unit_options = ["MPa", "bar", "kPa", "psi"]
                inferred_unit = infer_pressure_unit(measured_pressure_column)
                selected_unit = measured_pressure_unit if measured_pressure_unit in unit_options else inferred_unit
                if not str(settings.get("measured_pressure_column", "")):
                    selected_unit = inferred_unit
                measured_pressure_unit = st.selectbox(
                    "Unité de pression du fichier",
                    unit_options,
                    index=option_index(unit_options, selected_unit),
                )
                measured_pressure_subtract_initial = st.checkbox(
                    "Soustraire le zéro initial du capteur",
                    measured_pressure_subtract_initial,
                )
                uploaded_pressure_payload = measured_pressure_payload(
                    uploaded_columns,
                    measured_pressure_time_column,
                    measured_pressure_column,
                    measured_pressure_unit,
                    measured_pressure_subtract_initial,
                )
                st.caption(
                    f"{len(uploaded_pressure_payload['time'])} points | "
                    f"durée {uploaded_pressure_payload['time'][-1]:.3f} s | "
                    f"pression max {np.max(uploaded_pressure_payload['pressure_MPa']):.4f} MPa"
                )
            except (ValueError, KeyError) as exc:
                st.error(f"Historique de pression invalide : {exc}")
        else:
            st.info("Chargez un fichier CSV pour pouvoir relancer le calcul avec la pression mesurée.")

with st.sidebar.expander("Masse suspendue"):
    sidebar_help(
        [
            "Ce volet règle le mode d'actionnement libre avec une masse accrochée au bas de l'actionneur.",
            "La masse est convertie en charge F_load = m g.",
            "Le solveur cherche la géométrie libre qui équilibre la force, la flexion et la torsion.",
            "Ce mode calcule un déplacement et une déformation d'actionnement, pas une force bloquée.",
            "La pression peut monter à vitesse imposée puis redescendre, ou rester maintenue pour observer la relaxation libre.",
        ]
    )
    suspended_mass_g = st.number_input("Masse suspendue (g)", 0.1, 5000.0, max(float(settings["suspended_mass_g"]), 0.1), 10.0)
    suspended_duration_s = st.number_input(
        "Durée de simulation masse suspendue (s)",
        0.1,
        5000.0,
        float(settings["suspended_duration_s"]),
        5.0,
    )
    suspended_pressure_rate_mpa_s = st.number_input(
        "Vitesse d'actionnement masse suspendue (MPa/s)",
        0.001,
        10.0,
        float(settings["suspended_pressure_rate_mpa_s"]),
        0.01,
    )
    suspended_hold_pressure = st.checkbox(
        "Maintenir la pression après la rampe",
        bool(settings["suspended_hold_pressure"]),
    )
    suspended_show_geometry_plot = st.checkbox(
        "Ajouter le graphe Rh et beta_h",
        bool(settings.get("suspended_show_geometry_plot", INTERFACE_DEFAULT_SUSPENDED_SHOW_GEOMETRY_PLOT)),
    )

with st.sidebar.expander("Étude de précontrainte"):
    sidebar_help(
        [
            "Ce volet prépare une série de simulations avec plusieurs précontraintes initiales.",
            "La force maximale est ensuite tracée en fonction de la précontrainte.",
            "Plus le nombre de valeurs est grand, plus l'étude est précise mais longue à calculer.",
        ]
    )
    eps_study_min = st.number_input("Précontrainte minimale", 0.0, 3.0, float(settings["eps_study_min"]), 0.05)
    eps_study_max = st.number_input("Précontrainte maximale", 0.0, 3.0, float(settings["eps_study_max"]), 0.05)
    eps_study_points = st.slider("Nombre de valeurs de précontrainte", 2, 60, int(settings["eps_study_points"]), 1)

with st.sidebar.expander("Hystérèse"):
    sidebar_help(
        [
            "Ce volet choisit les courbes à superposer dans le graphe d'hystérèse.",
            "Cycles à afficher accepte plusieurs valeurs séparées par des virgules, par exemple 1, 2, 5.",
            "La comparaison peut garder le calcul courant ou relancer plusieurs cas avec différentes précontraintes ou vitesses de pression injectée.",
            "La vitesse de pression injectée est exprimée en MPa/s et impose une rampe de pression linéaire.",
        ]
    )
    hysteresis_cycles = st.text_input(
        "Cycles à afficher",
        str(settings.get("hysteresis_cycles", settings.get("hysteresis_cycle", 1))),
    )
    hysteresis_compare_mode = st.selectbox(
        "Comparaison",
        HYSTERESIS_COMPARE_OPTIONS,
        index=HYSTERESIS_COMPARE_OPTIONS.index(str(settings.get("hysteresis_compare_mode", "current")))
        if str(settings.get("hysteresis_compare_mode", "current")) in HYSTERESIS_COMPARE_OPTIONS
        else 0,
        format_func=lambda value: HYSTERESIS_COMPARE_LABELS.get(value, value),
    )
    if hysteresis_compare_mode == "prestrain":
        hysteresis_prestrain_values = st.text_input(
            "Précontraintes à comparer",
            str(settings.get("hysteresis_prestrain_values", "0.6, 0.8, 1.0")),
        )
        hysteresis_pressure_rates_mpa_s = str(settings.get("hysteresis_pressure_rates_mpa_s", "0.05, 0.10, 0.20"))
    elif hysteresis_compare_mode == "pressure_rate":
        hysteresis_pressure_rates_mpa_s = st.text_input(
            "Vitesses de pression injectée (MPa/s)",
            str(settings.get("hysteresis_pressure_rates_mpa_s", "0.05, 0.10, 0.20")),
        )
        hysteresis_prestrain_values = str(settings.get("hysteresis_prestrain_values", "0.6, 0.8, 1.0"))
    else:
        hysteresis_prestrain_values = str(settings.get("hysteresis_prestrain_values", "0.6, 0.8, 1.0"))
        hysteresis_pressure_rates_mpa_s = str(settings.get("hysteresis_pressure_rates_mpa_s", "0.05, 0.10, 0.20"))

with st.sidebar.expander("Relaxation"):
    sidebar_help(
        [
            "Ces paramètres servent au calcul de relaxation après actionnement bloqué.",
            "Le temps de montée en pression définit la rampe jusqu'à la pression maximale.",
            "Le temps de maintien définit la durée pendant laquelle la pression reste constante.",
            "La relaxation observée vient du modèle viscoélastique de Maxwell.",
        ]
    )
    relaxation_ramp_time_s = st.number_input(
        "Temps de montée en pression (s)", 0.01, 1000.0, float(settings["relaxation_ramp_time_s"]), 1.0
    )
    relaxation_hold_time_s = st.number_input(
        "Temps de maintien à pression constante (s)", 0.0, 5000.0, float(settings["relaxation_hold_time_s"]), 10.0
    )

with st.sidebar.expander("Matériau tube"):
    sidebar_help(
        [
            "Ces paramètres règlent la réponse élastique anisotrope du tube.",
            "E_axial agit principalement sur la raideur dans la direction de l'actionneur.",
            "E_radius et G12 influencent la réponse radiale et le couplage pression-torsion.",
            "Les coefficients de Poisson règlent le couplage entre les déformations transverses.",
        ]
    )
    E_axial_mpa = st.number_input("Module axial du tube E_axial (MPa)", 0.001, 10000.0, float(settings["E_axial_mpa"]), 0.1)
    axial_modulus_mode = st.selectbox(
        "Convention du module axial",
        AXIAL_MODULUS_OPTIONS,
        index=option_index(AXIAL_MODULUS_OPTIONS, settings.get("axial_modulus_mode", "maxwell_sum")),
        format_func=lambda value: AXIAL_MODULUS_LABELS.get(value, value),
    )
    maxwell_anisotropy_mode = st.selectbox(
        "Anisotropie de la relaxation",
        MAXWELL_ANISOTROPY_OPTIONS,
        index=option_index(
            MAXWELL_ANISOTROPY_OPTIONS,
            settings.get("maxwell_anisotropy_mode", "paper_equal"),
        ),
        format_func=lambda value: MAXWELL_ANISOTROPY_LABELS.get(value, value),
        help=(
            "Le mode axial applique la relaxation uniquement à la direction axiale. "
            "Le mode proportionnel applique les mêmes fractions de relaxation à toutes les directions du matériau."
        ),
    )
    E_radius_mpa = st.number_input("Module radial du tube E_radius (MPa)", 0.001, 10000.0, float(settings["E_radius_mpa"]), 0.1)
    G12_mpa = st.number_input("Module de cisaillement du tube G12 (MPa)", 0.001, 10000.0, float(settings["G12_mpa"]), 0.1)
    nu12 = st.number_input("Coefficient de Poisson nu12", -0.49, 0.49, float(np.clip(settings["nu12"], -0.49, 0.49)), 0.005)
    nu23 = st.number_input("Coefficient de Poisson nu23", -0.49, 0.49, float(np.clip(settings["nu23"], -0.49, 0.49)), 0.005)

with st.sidebar.expander("Maxwell généralisé"):
    sidebar_help(
        [
            "Ce volet règle la partie viscoélastique du tube.",
            "E0 est la raideur permanente qui ne relaxe pas.",
            "Chaque branche Ei, etai ajoute une relaxation avec un temps caractéristique proche de etai / Ei.",
            "Des viscosités plus grandes ralentissent la relaxation.",
        ]
    )
    constitutive_mode = st.selectbox(
        "Réponse constitutive",
        CONSTITUTIVE_OPTIONS,
        index=option_index(CONSTITUTIVE_OPTIONS, settings.get("constitutive_mode", "generalized_maxwell")),
        format_func=lambda value: CONSTITUTIVE_LABELS.get(value, value),
        help="Le mode élastique instantané désactive la relaxation et l'hystérésis sans changer la rigidité instantanée totale.",
    )
    prestrain_reference_mode = st.selectbox(
        "Traitement de la précontrainte",
        PRESTRAIN_REFERENCE_OPTIONS,
        index=option_index(
            PRESTRAIN_REFERENCE_OPTIONS,
            settings.get("prestrain_reference_mode", "elastic_tk_reference"),
        ),
        format_func=lambda value: PRESTRAIN_REFERENCE_LABELS.get(value, value),
        help=(
            "Précontrainte élastique conservée : la contrainte créée avant l'actionnement reste une base fixe. "
            "Précontrainte viscoélastique : la précontrainte est aussi traitée par les branches de Maxwell et peut relaxer."
        ),
    )
    maxwell_E0_mpa = st.number_input("Ressort permanent E0 (MPa)", 0.0, 10000.0, float(settings["maxwell_E0_mpa"]), 0.1)
    maxwell_E1_mpa = st.number_input("Branche E1 (MPa)", 0.0, 10000.0, float(settings["maxwell_E1_mpa"]), 0.1)
    maxwell_eta1_mpa_s = st.number_input("Branche eta1 (MPa s)", 1.0e-9, 1.0e9, float(settings["maxwell_eta1_mpa_s"]), 1.0)
    maxwell_E2_mpa = st.number_input("Branche E2 (MPa)", 0.0, 10000.0, float(settings["maxwell_E2_mpa"]), 0.1)
    maxwell_eta2_mpa_s = st.number_input("Branche eta2 (MPa s)", 1.0e-9, 1.0e9, float(settings["maxwell_eta2_mpa_s"]), 1.0)
    maxwell_E3_mpa = st.number_input("Branche E3 (MPa)", 0.0, 10000.0, float(settings["maxwell_E3_mpa"]), 0.1)
    maxwell_eta3_mpa_s = st.number_input("Branche eta3 (MPa s)", 1.0e-9, 1.0e9, float(settings["maxwell_eta3_mpa_s"]), 1.0)

with st.sidebar.expander("Nylon"):
    sidebar_help(
        [
            "Ces paramètres règlent la contribution du filament de nylon.",
            "E_nylon contrôle la force axiale liée à l'étirement du nylon.",
            "G_nylon contrôle sa contribution en torsion.",
            "Les facteurs de couplage indiquent dans quelle mesure le nylon participe à la précontrainte et à l'actionnement.",
        ]
    )
    E_nylon_mpa = st.number_input("Module axial du nylon E_nylon (MPa)", 0.001, 100000.0, float(settings["E_nylon_mpa"]), 10.0)
    G_nylon_mpa = st.number_input("Module de cisaillement du nylon G_nylon (MPa)", 0.001, 100000.0, float(settings["G_nylon_mpa"]), 10.0)
    nylon_condition_mode = st.selectbox(
        "Condition mécanique du nylon",
        NYLON_CONDITION_OPTIONS,
        index=option_index(NYLON_CONDITION_OPTIONS, settings.get("nylon_condition_mode", "bonded_linear")),
        format_func=lambda value: NYLON_CONDITION_LABELS.get(value, value),
        help=(
            "Linéaire bilatéral autorise le nylon à travailler en traction et en compression. Traction seulement "
            "représente un filament qui se détend au lieu de pousser. Glissant annule uniquement sa force axiale."
        ),
    )
    nylon_scale = 1.0
    nylon_axial_prestrain_coupling = 0.0 if nylon_condition_mode == "axially_sliding_confined" else 1.0
    nylon_axial_actuation_coupling = nylon_axial_prestrain_coupling
    st.caption("Les multiplicateurs continus du nylon sont verrouillés à leur valeur physique, sans calibration.")

with st.sidebar.expander("Force de fond pression"):
    sidebar_help(
        [
            "Ce volet ajoute éventuellement une force axiale directe due à la pression sur une surface de fond.",
            "Le mode choisi définit la surface utilisée pour cette force.",
            "La contribution vaut strictement pression multipliée par aire, sans facteur ajustable.",
            "Par défaut, cette contribution est nulle pour rester proche du modèle mécanique principal.",
            "Une valeur non nulle ajoute la poussée directe exercée par la pression sur un fond fermé.",
        ]
    )
    pressure_end_force_mode = st.selectbox(
        "Mode de force de fond",
        PRESSURE_END_FORCE_OPTIONS,
        index=option_index(PRESSURE_END_FORCE_OPTIONS, settings["pressure_end_force_mode"]),
        format_func=lambda value: PRESSURE_END_FORCE_LABELS.get(value, value),
    )
    pressure_end_force_scale = 0.0 if pressure_end_force_mode == "none" else 1.0
    st.caption(f"Facteur physique appliqué : {pressure_end_force_scale:.0f}")

with st.sidebar.expander("Solveur"):
    sidebar_help(
        [
            "Ces paramètres contrôlent la précision et le temps de calcul.",
            "dt est le pas de temps de la simulation.",
            "Les couches radiales et divisions angulaires définissent le maillage numérique de la section du tube.",
            "Les étapes de précontrainte divisent l'étirement initial en petits incréments.",
            "Les cœurs CPU parallèles sont utilisés pour les études avec plusieurs simulations indépendantes.",
            "Euler explicite : applique directement la loi incrémentale, mais impose dt < 2 fois le plus petit temps de relaxation.",
            "Exponentielle cohérente : intègre exactement la relaxation sur chaque pas à déformation affine et utilise le même opérateur pour le gonflement radial.",
        ]
    )
    dt = st.number_input("Pas de temps dt (s)", 0.01, 20.0, float(settings["dt"]), 0.05)
    n_layers = st.slider("Couches radiales du tube", 1, 30, int(settings["n_layers"]), 1)
    n_phi = st.slider("Divisions angulaires phi", 4, 120, int(settings["n_phi"]), 4)
    pre_steps = st.slider("Étapes de précontrainte", 1, 240, int(settings["pre_steps"]), 1)
    cpu_count = max(1, os.cpu_count() or 1)
    parallel_workers_default = int(settings.get("parallel_workers", INTERFACE_DEFAULT_PARALLEL_WORKERS))
    parallel_workers_default = max(1, min(parallel_workers_default, cpu_count))
    parallel_workers = st.slider("Cœurs CPU parallèles", 1, cpu_count, parallel_workers_default, 1)
    integration = st.selectbox(
        "Intégration temporelle",
        INTEGRATION_OPTIONS,
        index=option_index(INTEGRATION_OPTIONS, settings["integration"]),
        format_func=lambda value: INTEGRATION_LABELS.get(value, value),
    )

with st.sidebar.expander("Visualiseur"):
    sidebar_help(
        [
            "Ces paramètres ne modifient pas le calcul mécanique.",
            "Ils changent seulement l'angle de vue du visualiseur Cavatappi.",
            "Les flèches du clavier peuvent aussi être utilisées pour faire tourner la vue.",
        ]
    )
    view_elev_deg = st.slider("Élévation de vue (deg)", 0.0, 90.0, float(settings["view_elev_deg"]), 1.0)
    view_azim_deg = st.slider("Azimut de vue (deg)", -180.0, 180.0, float(settings["view_azim_deg"]), 1.0)

current_settings = {
    "_settings_schema_version": SETTINGS_SCHEMA_VERSION,
    "eps": float(eps),
    "eps_study_min": float(eps_study_min),
    "eps_study_max": float(eps_study_max),
    "eps_study_points": int(eps_study_points),
    "p_max_mpa": float(p_max_mpa),
    "rout_mm": float(rout_mm),
    "rin_mm": float(rin_mm),
    "nylon_diameter_mm": float(nylon_diameter_mm),
    "rho0_mm": float(rho0_mm),
    "alpha0_deg": float(alpha0_deg),
    "theta_f_deg": float(theta_f_deg),
    "bias_angle_profile": str(bias_angle_profile),
    "initial_length_mm": float(initial_length_mm),
    "section_update_mode": str(section_update_mode),
    "pressure_end_force_mode": str(pressure_end_force_mode),
    "pressure_end_force_scale": float(pressure_end_force_scale),
    "n_cycles": int(n_cycles),
    "hysteresis_cycle": parse_cycle_list(hysteresis_cycles, int(n_cycles))[0],
    "hysteresis_cycles": str(hysteresis_cycles),
    "hysteresis_compare_mode": str(hysteresis_compare_mode),
    "hysteresis_prestrain_values": str(hysteresis_prestrain_values),
    "hysteresis_pressure_rates_mpa_s": str(hysteresis_pressure_rates_mpa_s),
    "use_fixed_duration": bool(use_fixed_duration),
    "duration_s": float(duration_s),
    "flow_rate_mL_min": float(flow_rate_mL_min),
    "volume_mL": float(volume_mL),
    "nonlinear_pressure": bool(nonlinear_pressure),
    "pressure_input_mode": str(pressure_input_mode),
    "measured_pressure_time_column": str(measured_pressure_time_column),
    "measured_pressure_column": str(measured_pressure_column),
    "measured_pressure_unit": str(measured_pressure_unit),
    "measured_pressure_file_hash": str(measured_pressure_file_hash if pressure_input_mode == "measured_csv" else ""),
    "measured_pressure_subtract_initial": bool(measured_pressure_subtract_initial),
    "relaxation_ramp_time_s": float(relaxation_ramp_time_s),
    "relaxation_hold_time_s": float(relaxation_hold_time_s),
    "suspended_mass_g": float(suspended_mass_g),
    "suspended_duration_s": float(suspended_duration_s),
    "suspended_pressure_rate_mpa_s": float(suspended_pressure_rate_mpa_s),
    "suspended_hold_pressure": bool(suspended_hold_pressure),
    "suspended_show_geometry_plot": bool(suspended_show_geometry_plot),
    "dt": float(dt),
    "pre_steps": int(pre_steps),
    "integration": str(integration),
    "prestrain_reference_mode": str(prestrain_reference_mode),
    "constitutive_mode": str(constitutive_mode),
    "maxwell_anisotropy_mode": str(maxwell_anisotropy_mode),
    "axial_modulus_mode": str(axial_modulus_mode),
    "E_axial_mpa": float(E_axial_mpa),
    "E_radius_mpa": float(E_radius_mpa),
    "G12_mpa": float(G12_mpa),
    "nu12": float(nu12),
    "nu23": float(nu23),
    "maxwell_E0_mpa": float(maxwell_E0_mpa),
    "maxwell_E1_mpa": float(maxwell_E1_mpa),
    "maxwell_eta1_mpa_s": float(maxwell_eta1_mpa_s),
    "maxwell_E2_mpa": float(maxwell_E2_mpa),
    "maxwell_eta2_mpa_s": float(maxwell_eta2_mpa_s),
    "maxwell_E3_mpa": float(maxwell_E3_mpa),
    "maxwell_eta3_mpa_s": float(maxwell_eta3_mpa_s),
    "E_nylon_mpa": float(E_nylon_mpa),
    "G_nylon_mpa": float(G_nylon_mpa),
    "nylon_axial_prestrain_coupling": float(nylon_axial_prestrain_coupling),
    "nylon_axial_actuation_coupling": float(nylon_axial_actuation_coupling),
    "nylon_condition_mode": str(nylon_condition_mode),
    "nylon_scale": float(nylon_scale),
    "n_layers": int(n_layers),
    "n_phi": int(n_phi),
    "parallel_workers": int(parallel_workers),
    "view_elev_deg": float(view_elev_deg),
    "view_azim_deg": float(view_azim_deg),
}
save_settings(current_settings)
inject_keyboard_view_controls(current_settings)
initialize_cached_results()
if "hysteresis_comparison_result" not in st.session_state:
    st.session_state["hysteresis_comparison_result"] = None

error = settings_error(current_settings)
if error:
    st.error(error)
if nylon_condition_mode != "bonded_linear":
    st.info(
        "La condition sélectionnée représente un filament qui peut se détendre ou glisser axialement. "
        "Aucun coefficient d'ajustement intermédiaire n'est appliqué."
    )
maxwell_sum = float(maxwell_E0_mpa + maxwell_E1_mpa + maxwell_E2_mpa + maxwell_E3_mpa)
effective_axial_modulus = maxwell_sum if axial_modulus_mode == "maxwell_sum" else float(E_axial_mpa)
if maxwell_sum > 0.0 and abs(effective_axial_modulus - maxwell_sum) / maxwell_sum > 0.05:
    st.info(
        f"Le module axial ({float(E_axial_mpa):.2f} MPa) diffère de la somme des modules de Maxwell "
        f"({maxwell_sum:.2f} MPa). Le modèle conserve E_axial pour l'anisotropie et les modules de Maxwell pour leurs fractions relatives."
    )
if constitutive_mode == "instantaneous_elastic":
    st.info(
        "Le mode élastique instantané conserve la rigidité instantanée totale mais désactive la relaxation et l'hystérésis. "
        "Il sert de référence physique, pas de calibration."
    )
if maxwell_anisotropy_mode == "axial_test_only":
    st.caption(
        "Les paramètres de relaxation proviennent d'un essai de traction : seule la direction matérielle axiale "
        "reçoit ces branches. La raideur radiale et le cisaillement restent élastiques faute de mesures dédiées."
    )
if section_update_mode == "updated":
    st.caption(
        "La section radiale évolutive actualise les rayons et l'orientation du matériau à chaque incrément."
    )
if str(integration) == "paper_explicit":
    active_tau = [
        eta / modulus
        for modulus, eta in (
            (float(maxwell_E1_mpa), float(maxwell_eta1_mpa_s)),
            (float(maxwell_E2_mpa), float(maxwell_eta2_mpa_s)),
            (float(maxwell_E3_mpa), float(maxwell_eta3_mpa_s)),
        )
        if modulus > 0.0
    ]
    if active_tau and float(dt) > min(active_tau):
        st.warning(
            f"Le pas dépasse le plus petit temps de relaxation ({min(active_tau):.3g} s). "
            "Euler reste éventuellement stable sous 2 tau, mais peut osciller ; l'intégration exponentielle est recommandée."
        )

if pressure_input_mode == "measured_csv" and uploaded_pressure_payload is not None:
    estimated_duration_s = float(uploaded_pressure_payload["time"][-1])
    estimated_steps = max(1, len(uploaded_pressure_payload["time"]) - 1)
elif bool(use_fixed_duration):
    estimated_duration_s = float(duration_s)
    estimated_steps = int(np.ceil(estimated_duration_s / dt))
else:
    estimated_duration_s = int(n_cycles) * 2.0 * 60.0 * float(volume_mL) / float(flow_rate_mL_min)
    estimated_steps = int(np.ceil(estimated_duration_s / dt))
estimated_cost = model_cost_index(estimated_steps, n_layers, n_phi, pre_steps)
estimated_relaxation_steps = int(np.ceil((relaxation_ramp_time_s + relaxation_hold_time_s) / dt))
estimated_relaxation_cost = model_cost_index(estimated_relaxation_steps, n_layers, n_phi, pre_steps)
estimated_prestrain_cost = estimated_cost * int(eps_study_points)
estimated_suspended_steps = int(np.ceil(float(suspended_duration_s) / dt))
estimated_suspended_cost = model_cost_index(estimated_suspended_steps, n_layers, n_phi, pre_steps, step_multiplier=6.0)
parallel_worker_count = max(1, int(parallel_workers))
estimated_compute_s = estimate_compute_seconds(timing_profile, "blocked", estimated_cost)
estimated_relaxation_compute_s = estimate_compute_seconds(timing_profile, "relaxation", estimated_relaxation_cost)
estimated_prestrain_compute_s = estimate_compute_seconds(
    timing_profile,
    "prestrain_study",
    estimated_prestrain_cost,
    parallel_worker_count,
    int(eps_study_points),
)
estimated_suspended_compute_s = estimate_compute_seconds(timing_profile, "suspended", estimated_suspended_cost)
selected_hysteresis_cycles = parse_cycle_list(current_settings["hysteresis_cycles"], int(n_cycles))
hysteresis_prestrain_compare_values = parse_positive_float_list(current_settings["hysteresis_prestrain_values"])
hysteresis_pressure_rate_values = parse_positive_float_list(current_settings["hysteresis_pressure_rates_mpa_s"])
if hysteresis_compare_mode == "prestrain":
    estimated_hysteresis_compare_cost = estimated_cost * max(1, len(hysteresis_prestrain_compare_values))
elif hysteresis_compare_mode == "pressure_rate":
    estimated_hysteresis_compare_cost = sum(
        model_cost_index(
            int(np.ceil((max(selected_hysteresis_cycles) * 2.0 * max(float(p_max_mpa), 1e-12) / rate) / dt)),
            n_layers,
            n_phi,
            pre_steps,
        )
        for rate in hysteresis_pressure_rate_values
    )
else:
    estimated_hysteresis_compare_cost = 0
if hysteresis_compare_mode == "current":
    estimated_hysteresis_compare_s = estimate_compute_seconds(timing_profile, "hysteresis_compare", estimated_hysteresis_compare_cost)
else:
    estimated_hysteresis_parallel_factor = max(
        1,
        min(
            parallel_worker_count,
            len(hysteresis_prestrain_compare_values)
            if hysteresis_compare_mode == "prestrain"
            else len(hysteresis_pressure_rate_values),
        ),
    )
    estimated_hysteresis_compare_s = estimate_compute_seconds(
        timing_profile,
        "hysteresis_compare",
        estimated_hysteresis_compare_cost,
        parallel_worker_count,
        estimated_hysteresis_parallel_factor,
    )
hysteresis_compare_disabled = error is not None or estimated_hysteresis_compare_cost > 750_000
blocked_result_signature = settings_signature(current_settings, BLOCKED_RESULT_IGNORE_KEYS)
relaxation_result_signature = settings_signature(current_settings, RELAXATION_RESULT_IGNORE_KEYS)
prestrain_result_signature = settings_signature(current_settings, PRESTRAIN_RESULT_IGNORE_KEYS)
suspended_result_signature = settings_signature(current_settings, SUSPENDED_RESULT_IGNORE_KEYS)
hysteresis_settings_signature = settings_signature(current_settings, HYSTERESIS_COMPARE_IGNORE_KEYS)
hysteresis_comparison_signature = (
    hysteresis_compare_mode,
    tuple(selected_hysteresis_cycles),
    tuple(hysteresis_prestrain_compare_values)
    if hysteresis_compare_mode == "prestrain"
    else tuple(hysteresis_pressure_rate_values),
    hysteresis_settings_signature,
)

left, right = st.columns([1.3, 1.2], gap="large")

with left:
    visual_state = st.radio(
        "État du visualiseur",
        VISUAL_STATE_OPTIONS,
        horizontal=True,
        format_func=lambda value: VISUAL_STATE_LABELS.get(value, value),
    )
    fig_visual = make_cavatappi_figure(current_settings, visual_state)
    st.pyplot(fig_visual)
    plt.close(fig_visual)

with right:
    geom = derived_geometry(current_settings)
    st.subheader("Résultats géométriques")
    c1, c2 = st.columns(2)
    c1.metric("Nombre de spires", f"{geom['turns']:.2f}")
    c2.metric("Pas", f"{geom['pitch0_mm']:.2f} mm")
    c1.metric("Indice rho/Rout", f"{geom['spring_index']:.2f}")
    c2.metric("Mandrin estimé", f"{geom['equivalent_mandrel_diameter_mm']:.2f} mm")
    c1.metric("Aire de paroi du tube", f"{geom['wall_area_mm2']:.3f} mm²")
    c2.metric("Volume interne", f"{geom['tube_internal_volume_ml']:.4f} mL")
    c1.metric("Remplissage nylon", f"{100.0 * geom['nylon_fill_ratio']:.1f} %")
    c2.metric("Longueur précontrainte", f"{geom['prestrained_length_mm']:.2f} mm")

    fig_section = make_cross_section_figure(current_settings)
    st.pyplot(fig_section)
    plt.close(fig_section)

tabs = st.tabs(
    [
        "Sortie modèle",
        "Courbes temporelles",
        "Étude de l'hystérèse",
        "Étude de la relaxation",
        "Étude de précontrainte",
        "Masse suspendue",
    ]
)

missing_measured_pressure = pressure_input_mode == "measured_csv" and uploaded_pressure_payload is None
run_disabled = error is not None or estimated_cost > 250_000 or missing_measured_pressure
relaxation_run_disabled = error is not None or estimated_relaxation_cost > 250_000
suspended_run_disabled = error is not None or estimated_suspended_cost > 250_000
prestrain_range_error = float(eps_study_max) <= float(eps_study_min)
prestrain_study_run_disabled = error is not None or prestrain_range_error or estimated_prestrain_cost > 750_000

with tabs[0]:
    blocked_csv_payload = None
    st.caption(f"Temps de calcul estimé : ~{format_seconds(estimated_compute_s)}.")
    if estimated_cost > 120_000 and not run_disabled:
        st.warning("Cette simulation peut être lente. Augmentez le pas de temps ou réduisez le maillage pour l'interaction.")
    if run_disabled and error is None:
        if missing_measured_pressure:
            st.warning("Chargez l'historique CSV avant de relancer l'actionnement bloqué.")
        else:
            st.warning("Les réglages du solveur sont trop lourds pour l'interface interactive.")

    if st.button("Calculer l'actionnement bloqué", disabled=run_disabled):
        (config, data, summary), elapsed_s = run_with_progress(
            "Actionnement bloqué",
            estimated_compute_s,
            run_model,
            current_settings,
            uploaded_pressure_payload,
        )
        timing_profile = record_timing_sample(timing_profile, TIMING_PROFILE_PATH, "blocked", estimated_cost, elapsed_s)
        st.session_state["calculator_result"] = {
            "config": config,
            "data": data,
            "summary": summary,
            "elapsed_s": elapsed_s,
            "estimated_s": estimated_compute_s,
            "settings": dict(current_settings),
            "signature": blocked_result_signature,
        }
        save_result_cache(BLOCKED_RESULT_PATH, st.session_state["calculator_result"])

    stored_result = st.session_state["calculator_result"]
    result = stored_result if result_matches_settings(stored_result, blocked_result_signature, BLOCKED_RESULT_IGNORE_KEYS) else None
    if result is None:
        if stored_result is None:
            st.info("Cliquez sur 'Calculer l'actionnement bloqué' pour lancer le modèle.")
        else:
            st.warning("Les paramètres ont changé depuis le dernier calcul. Veuillez relancer l'actionnement bloqué pour mettre les résultats à jour.")
    else:
        summary = result["summary"]
        data = result["data"]
        config = result["config"]
        k1, k2, k3, k4 = st.columns(4)
        force_gain = float(np.nanmax(data["force_act_mN"]))
        k1.metric("Force min", f"{summary['force_min_mN']:.1f} mN")
        k2.metric("Force max", f"{summary['force_max_mN']:.1f} mN")
        k3.metric("Gain d'actionnement", f"{force_gain:.1f} mN")
        k4.metric("Couple max", f"{summary['torque_act_max_microNm']:.1f} microN m")
        st.caption(
            f"pression : {'CSV mesuré' if result['settings'].get('pressure_input_mode') == 'measured_csv' else 'profil généré'} | "
            f"durée : {data['time'][-1]:.1f} s | "
            f"période de cycle : {result_cycle_period(config):.2f} s | "
            f"échantillons : {len(data['time'])} | "
            f"résidu maximal : {summary['max_abs_residual_Nmm']:.2e} N·mm | "
            f"calcul : {format_seconds(result.get('elapsed_s', 0.0))} "
            f"(estimé {format_seconds(result.get('estimated_s', 0.0))})"
        )
        blocked_csv_payload = mapping_to_csv_bytes(data)

    render_csv_download(
        blocked_csv_payload,
        "resultats_actionnement_bloque.csv",
        "download_blocked_csv",
    )

with tabs[1]:
    temporal_csv_payload = None
    stored_result = st.session_state["calculator_result"]
    result = stored_result if result_matches_settings(stored_result, blocked_result_signature, BLOCKED_RESULT_IGNORE_KEYS) else None
    if result is None:
        if stored_result is None:
            st.info("Veuillez d'abord lancer le modèle dans l'onglet 'Sortie modèle' pour afficher les courbes temporelles.")
        else:
            st.warning("Les paramètres ont changé depuis le dernier calcul. Veuillez relancer l'actionnement bloqué pour mettre les courbes à jour.")
    else:
        fig_response = plot_time_response_fr(result["data"])
        st.pyplot(fig_response)
        plt.close(fig_response)
        temporal_csv_payload = mapping_to_csv_bytes(result["data"])

    render_csv_download(
        temporal_csv_payload,
        "courbes_temporelles_actionnement_bloque.csv",
        "download_temporal_csv",
    )

with tabs[2]:
    hysteresis_export_cases: list[dict] = []
    hysteresis_csv_payload = None
    st.caption(f"Cycles sélectionnés : {', '.join(str(cycle) for cycle in selected_hysteresis_cycles)}.")
    if hysteresis_compare_mode == "current":
        stored_result = st.session_state["calculator_result"]
        result = stored_result if result_matches_settings(stored_result, blocked_result_signature, BLOCKED_RESULT_IGNORE_KEYS) else None
        if result is None:
            if stored_result is None:
                st.info("Veuillez d'abord lancer le modèle dans l'onglet 'Sortie modèle' pour afficher l'hystérèse.")
            else:
                st.warning("Les paramètres ont changé depuis le dernier calcul. Veuillez relancer l'actionnement bloqué pour mettre l'hystérèse à jour.")
        else:
            try:
                config = result["config"]
                cycles_to_plot = parse_cycle_list(current_settings["hysteresis_cycles"], int(config.n_cycles))
                cases = [
                    {
                        "label": "calcul courant",
                        "data": result["data"],
                        "period": result_cycle_period(config),
                    }
                ]
                hysteresis_export_cases = cases
                fig_hyst = plot_hysteresis_overlay(cases, cycles_to_plot, show=False)
                st.pyplot(fig_hyst)
                plt.close(fig_hyst)
            except Exception as exc:
                st.error(f"Impossible de tracer l'hystérèse : {exc}")
    else:
        if hysteresis_compare_mode == "prestrain":
            comparison_values = hysteresis_prestrain_compare_values
            comparison_label = "précontraintes initiales"
        else:
            comparison_values = hysteresis_pressure_rate_values
            comparison_label = "vitesses de pression injectée"

        st.caption(f"Temps de calcul estimé : ~{format_seconds(estimated_hysteresis_compare_s)}.")
        if not comparison_values:
            st.warning(f"Veuillez indiquer au moins une valeur positive pour les {comparison_label}.")
        if estimated_hysteresis_compare_cost > 120_000 and not hysteresis_compare_disabled:
            st.warning("Cette comparaison peut être lente. Réduisez le nombre de valeurs, augmentez le pas de temps ou réduisez le maillage.")
        if hysteresis_compare_disabled and error is None:
            st.warning("Cette comparaison d'hystérèse est trop lourde pour l'interface interactive.")

        if st.button("Calculer la comparaison d'hystérèse", disabled=hysteresis_compare_disabled or not comparison_values):
            comparison_result, elapsed_s = run_hysteresis_comparison_with_progress(
                current_settings,
                hysteresis_compare_mode,
                comparison_values,
                selected_hysteresis_cycles,
                estimated_hysteresis_compare_s,
                parallel_worker_count,
            )
            timing_profile = record_timing_sample(
                timing_profile,
                TIMING_PROFILE_PATH,
                "hysteresis_compare",
                estimated_hysteresis_compare_cost,
                elapsed_s,
                parallel_worker_count,
                len(comparison_values),
            )
            comparison_result["elapsed_s"] = elapsed_s
            comparison_result["estimated_s"] = estimated_hysteresis_compare_s
            comparison_result["settings"] = dict(current_settings)
            comparison_result["signature"] = hysteresis_comparison_signature
            st.session_state["hysteresis_comparison_result"] = comparison_result
            save_result_cache(HYSTERESIS_RESULT_PATH, comparison_result)

        comparison_result = st.session_state["hysteresis_comparison_result"]
        if comparison_result is None or comparison_result.get("signature") != hysteresis_comparison_signature:
            st.info("Veuillez lancer la comparaison pour afficher les courbes d'hystérèse.")
        else:
            try:
                fig_hyst = plot_hysteresis_overlay(
                    comparison_result["cases"],
                    selected_hysteresis_cycles,
                    show=False,
                )
                st.pyplot(fig_hyst)
                plt.close(fig_hyst)
                hysteresis_export_cases = comparison_result["cases"]
                st.caption(
                    f"Calcul : {format_seconds(comparison_result.get('elapsed_s', 0.0))} "
                    f"(estimé {format_seconds(comparison_result.get('estimated_s', 0.0))})"
                )
            except Exception as exc:
                st.error(f"Impossible de tracer la comparaison d'hystérèse : {exc}")

    if hysteresis_export_cases:
        try:
            hysteresis_csv_payload = hysteresis_to_csv_bytes(
                hysteresis_export_cases,
                selected_hysteresis_cycles,
                str(hysteresis_compare_mode),
            )
        except ValueError:
            hysteresis_csv_payload = None
    render_csv_download(
        hysteresis_csv_payload,
        "resultats_hysterese.csv",
        "download_hysteresis_csv",
    )

with tabs[3]:
    relaxation_csv_payload = None
    st.caption(f"Temps de calcul estimé : ~{format_seconds(estimated_relaxation_compute_s)}.")
    if estimated_relaxation_cost > 120_000 and not relaxation_run_disabled:
        st.warning("La relaxation peut être lente. Augmentez le pas de temps, réduisez le maintien ou réduisez le maillage.")
    if relaxation_run_disabled and error is None:
        st.warning("Les réglages de relaxation sont trop lourds pour l'interface interactive.")

    if st.button("Calculer la relaxation à pression constante", disabled=relaxation_run_disabled):
        (config, data, summary), elapsed_s = run_with_progress(
            "Relaxation à pression constante",
            estimated_relaxation_compute_s,
            run_relaxation_model,
            current_settings,
        )
        timing_profile = record_timing_sample(
            timing_profile,
            TIMING_PROFILE_PATH,
            "relaxation",
            estimated_relaxation_cost,
            elapsed_s,
        )
        st.session_state["relaxation_result"] = {
            "config": config,
            "data": data,
            "summary": summary,
            "elapsed_s": elapsed_s,
            "estimated_s": estimated_relaxation_compute_s,
            "settings": dict(current_settings),
            "signature": relaxation_result_signature,
        }
        save_result_cache(RELAXATION_RESULT_PATH, st.session_state["relaxation_result"])

    stored_relaxation_result = st.session_state["relaxation_result"]
    relaxation_result = (
        stored_relaxation_result
        if result_matches_settings(stored_relaxation_result, relaxation_result_signature, RELAXATION_RESULT_IGNORE_KEYS)
        else None
    )
    if relaxation_result is None:
        if stored_relaxation_result is None:
            st.info("Veuillez lancer la relaxation pour afficher les courbes de maintien à pression constante.")
        else:
            st.warning("Les paramètres ont changé depuis le dernier calcul. Veuillez relancer la relaxation pour mettre les courbes à jour.")
    else:
        relaxation_data = relaxation_result["data"]
        hold_start_index = int(relaxation_data["hold_start_index"])
        r1, r2, r3, r4 = st.columns(4)
        r1.metric("Force au début du maintien", f"{relaxation_data['force_total_mN'][hold_start_index]:.1f} mN")
        r2.metric("Force en fin de maintien", f"{relaxation_data['force_total_mN'][-1]:.1f} mN")
        r3.metric("Relaxation de force", f"{relaxation_data['force_hold_relax_mN'][-1]:.1f} mN")
        r4.metric("Pression maintenue", f"{np.nanmax(relaxation_data['pressure_MPa']):.3f} MPa")
        st.caption(
            f"Calcul : {format_seconds(relaxation_result.get('elapsed_s', 0.0))} "
            f"(estimé {format_seconds(relaxation_result.get('estimated_s', 0.0))})"
        )
        fig_relax = plot_relaxation_response(relaxation_data)
        st.pyplot(fig_relax)
        plt.close(fig_relax)
        relaxation_csv_payload = mapping_to_csv_bytes(relaxation_data)

    render_csv_download(
        relaxation_csv_payload,
        "resultats_relaxation.csv",
        "download_relaxation_csv",
    )

with tabs[4]:
    prestrain_csv_payload = None
    st.caption(f"Temps de calcul estimé : ~{format_seconds(estimated_prestrain_compute_s)}.")
    if estimated_prestrain_cost > 120_000 and not prestrain_study_run_disabled:
        st.warning(
            "L'étude de précontrainte peut être lente. Réduisez le nombre de valeurs, augmentez le pas de temps ou réduisez le maillage."
        )
    if prestrain_range_error:
        st.warning("La précontrainte maximale doit être strictement supérieure à la précontrainte minimale.")
    elif prestrain_study_run_disabled and error is None:
        st.warning("Les réglages de l'étude de précontrainte sont trop lourds pour l'interface interactive.")

    if st.button("Calculer l'étude force max / précontrainte", disabled=prestrain_study_run_disabled):
        eps_values = np.linspace(float(eps_study_min), float(eps_study_max), int(eps_study_points))
        study_result, elapsed_s = run_prestrain_study_with_progress(
            current_settings,
            eps_values,
            estimated_prestrain_compute_s,
            parallel_worker_count,
        )
        timing_profile = record_timing_sample(
            timing_profile,
            TIMING_PROFILE_PATH,
            "prestrain_study",
            estimated_prestrain_cost,
            elapsed_s,
            parallel_worker_count,
            int(eps_study_points),
        )
        st.session_state["prestrain_study_result"] = {
            "data": study_result,
            "elapsed_s": elapsed_s,
            "estimated_s": estimated_prestrain_compute_s,
            "settings": dict(current_settings),
            "signature": prestrain_result_signature,
        }
        save_result_cache(PRESTRAIN_RESULT_PATH, st.session_state["prestrain_study_result"])

    stored_prestrain_result = st.session_state["prestrain_study_result"]
    prestrain_result = (
        stored_prestrain_result
        if result_matches_settings(stored_prestrain_result, prestrain_result_signature, PRESTRAIN_RESULT_IGNORE_KEYS)
        else None
    )
    if prestrain_result is None:
        if stored_prestrain_result is None:
            st.info("Veuillez lancer l'étude pour afficher la force maximale en fonction de la précontrainte initiale.")
        else:
            st.warning("Les paramètres ont changé depuis le dernier calcul. Veuillez relancer l'étude pour mettre la courbe à jour.")
    else:
        study_data = prestrain_result["data"]
        best_index = int(np.nanargmax(study_data["force_max_mN"]))
        s1, s2, s3 = st.columns(3)
        s1.metric("Force max", f"{study_data['force_max_mN'][best_index]:.1f} mN")
        s2.metric("Gain max", f"{study_data['force_gain_mN'][best_index]:.1f} mN")
        s3.metric("Points calculés", f"{len(study_data['eps'])}")
        st.caption(
            f"Calcul : {format_seconds(prestrain_result.get('elapsed_s', 0.0))} "
            f"(estimé {format_seconds(prestrain_result.get('estimated_s', 0.0))})"
        )
        fig_study = plot_prestrain_study(study_data)
        st.pyplot(fig_study)
        plt.close(fig_study)
        prestrain_csv_payload = mapping_to_csv_bytes(study_data)

    render_csv_download(
        prestrain_csv_payload,
        "resultats_etude_precontrainte.csv",
        "download_prestrain_csv",
    )

with tabs[5]:
    suspended_csv_payload = None
    load_N = float(suspended_mass_g) * 1.0e-3 * 9.80665
    suspended_ramp_time_s = float(p_max_mpa) / max(float(suspended_pressure_rate_mpa_s), 1e-12)
    suspended_profile_label = (
        "rampe puis maintien à pression constante"
        if bool(suspended_hold_pressure)
        else "rampe puis descente de pression"
    )
    st.caption(
        f"Charge appliquée : {load_N:.3f} N pour une masse de {float(suspended_mass_g):.1f} g | "
        f"temps de calcul estimé : ~{format_seconds(estimated_suspended_compute_s)}."
    )
    st.caption(
        f"Profil masse suspendue : {suspended_profile_label} | "
        f"durée : {float(suspended_duration_s):.1f} s | "
        f"vitesse : {float(suspended_pressure_rate_mpa_s):.3f} MPa/s | "
        f"temps de rampe jusqu'à Pmax : {format_seconds(suspended_ramp_time_s)}."
    )
    if float(suspended_duration_s) < suspended_ramp_time_s:
        st.warning(
            "La durée masse suspendue est plus courte que le temps nécessaire pour atteindre la pression maximale."
        )
    if bool(suspended_hold_pressure) and float(suspended_duration_s) <= suspended_ramp_time_s:
        st.warning("Le maintien à pression constante ne sera visible que si la durée dépasse le temps de rampe.")
    st.info(
        "Ce mode résout les équilibres mécaniques suivants : "
        "F_tube + F_nylon = F_load sin(beta_h), "
        "M_tube + M_nylon = -F_load Rh sin(beta_h), "
        "T_tube + T_nylon = F_load Rh cos(beta_h)."
    )
    if estimated_suspended_cost > 120_000 and not suspended_run_disabled:
        st.warning("Le mode masse suspendue peut être lent. Augmentez le pas de temps ou réduisez le maillage.")
    if suspended_run_disabled and error is None:
        st.warning("Les réglages du mode masse suspendue sont trop lourds pour l'interface interactive.")

    if st.button("Calculer l'actionnement avec masse suspendue", disabled=suspended_run_disabled):
        (config, data, summary), elapsed_s = run_with_progress(
            "Actionnement masse suspendue",
            estimated_suspended_compute_s,
            run_suspended_model,
            current_settings,
        )
        timing_profile = record_timing_sample(
            timing_profile,
            TIMING_PROFILE_PATH,
            "suspended",
            estimated_suspended_cost,
            elapsed_s,
        )
        st.session_state["suspended_result"] = {
            "config": config,
            "data": data,
            "summary": summary,
            "elapsed_s": elapsed_s,
            "estimated_s": estimated_suspended_compute_s,
            "settings": dict(current_settings),
            "signature": suspended_result_signature,
        }
        save_result_cache(SUSPENDED_RESULT_PATH, st.session_state["suspended_result"])

    stored_suspended_result = st.session_state["suspended_result"]
    suspended_result = (
        stored_suspended_result
        if result_matches_settings(stored_suspended_result, suspended_result_signature, SUSPENDED_RESULT_IGNORE_KEYS)
        else None
    )
    if suspended_result is None:
        if stored_suspended_result is None:
            st.info("Veuillez lancer le calcul pour afficher l'actionnement libre sous masse suspendue.")
        else:
            st.warning("Les paramètres ont changé depuis le dernier calcul. Veuillez relancer le calcul masse suspendue pour mettre les courbes à jour.")
    else:
        suspended_data = suspended_result["data"]
        a1, a2, a3, a4 = st.columns(4)
        a1.metric("Contraction max", f"{np.nanmax(suspended_data['free_contraction_mm']):.3f} mm")
        a2.metric("Actionnement max", f"{np.nanmax(suspended_data['free_actuation_percent']):.2f} %")
        a3.metric("Longueur finale", f"{suspended_data['axial_length_mm'][-1]:.2f} mm")
        a4.metric("Résidu max", f"{np.nanmax(suspended_data['residual']):.2e}")
        if bool(suspended_data.get("suspended_hold_pressure", False)):
            hold_index = int(suspended_data.get("suspended_hold_start_index", 0))
            hold_time = float(suspended_data["time"][hold_index])
            h1, h2, h3 = st.columns(3)
            h1.metric("Début du maintien", f"{hold_time:.1f} s")
            h2.metric("Relaxation contraction", f"{suspended_data['free_contraction_hold_relax_mm'][-1]:.4f} mm")
            h3.metric("Relaxation actionnement", f"{suspended_data['free_actuation_hold_relax_percent'][-1]:.3f} %")
        st.caption(
            f"Calcul : {format_seconds(suspended_result.get('elapsed_s', 0.0))} "
            f"(estimé {format_seconds(suspended_result.get('estimated_s', 0.0))})"
        )
        fig_suspended = make_suspended_response_figure(
            suspended_data,
            show_geometry=bool(current_settings["suspended_show_geometry_plot"]),
        )
        st.pyplot(fig_suspended)
        plt.close(fig_suspended)
        suspended_csv_payload = mapping_to_csv_bytes(suspended_data)

    render_csv_download(
        suspended_csv_payload,
        "resultats_masse_suspendue.csv",
        "download_suspended_csv",
    )
