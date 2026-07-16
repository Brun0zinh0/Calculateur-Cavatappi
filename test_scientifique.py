"""Tests numériques autonomes de l'interface Cavatappi Beta.

Exécution : python test_scientifique.py
"""

from __future__ import annotations

import csv
from io import StringIO
import sys
from pathlib import Path

import numpy as np


APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

import Base
import affichage
import parametres
import pression


def test_pressure_histories() -> None:
    time, pressure = Base.ramp_hold_pressure_history(P_hold=1.5, ramp_time=9.0, hold_time=100.0, dt=20.0)
    assert time[-1] == 109.0
    assert pressure[-1] == 1.5
    assert np.all(np.diff(time) > 0.0)

    time, pressure = Base.cyclic_pressure_history(n_cycles=2, Pmax=1.5, dt=7.0)
    assert time[-1] == 36.0
    assert pressure[-1] == 0.0
    assert np.isclose(np.max(pressure), 1.5)


def test_blocked_equilibrium() -> None:
    settings = dict(parametres.DEFAULT_SETTINGS)
    settings.update({"eps": 0.8, "n_cycles": 1, "dt": 1.0, "n_layers": 2, "n_phi": 8, "pre_steps": 4})
    config = parametres.build_config(settings)
    _, data = Base.run_blocked_actuation(config)
    assert data["time"][0] == 0.0
    assert data["pressure_MPa"][0] == 0.0
    assert np.max(np.abs(data["residual"])) < 1.0e-5
    force_error = np.max(
        np.abs(
            data["force_total_mN"]
            - data["force_tube_mN"]
            - data["force_nylon_mN"]
        )
    )
    torque_error = np.max(
        np.abs(
            data["torque_total_signed_microNm"]
            - data["torque_tube_microNm"]
            - data["torque_nylon_microNm"]
        )
    )
    assert force_error < 1.0e-8
    assert torque_error < 1.0e-8
    assert "force_pressure_end_mN" not in data
    assert "force_structural_mN" not in data
    assert not any("pressure_end" in key for key in parametres.DEFAULT_SETTINGS)
    assert np.all(data["Rin_mm"] > 0.0)
    assert np.all(data["Rout_mm"] > data["Rin_mm"])


def test_suspended_equation_24() -> None:
    settings = dict(parametres.DEFAULT_SETTINGS)
    settings.update({"eps": 0.8, "dt": 1.0, "n_layers": 2, "n_phi": 8, "pre_steps": 4})
    config = parametres.build_config(settings)
    time = np.array([0.0, 1.0, 2.0, 3.0])
    pressure = np.array([0.0, 0.05, 0.10, 0.15])
    model, data = Base.run_suspended_actuation(
        config,
        load_N=0.981,
        pressure_time=time,
        pressure_MPa=pressure,
    )
    expected_length = (
        config.geom.initial_length
        * data["axial_stretch"]
        * np.sin(data["alpha_rad"])
        / np.sin(model.alpha0)
    )
    assert np.max(np.abs(expected_length - data["axial_length_mm"])) < 1.0e-10
    assert data["free_contraction_mm"][0] == 0.0
    assert np.max(data["residual"]) < 1.0e-4


def test_input_guards() -> None:
    bad = dict(parametres.DEFAULT_SETTINGS)
    bad.update({"maxwell_E0_mpa": 0.0, "maxwell_E1_mpa": 0.0, "maxwell_E2_mpa": 0.0, "maxwell_E3_mpa": 0.0})
    assert parametres.material_error(bad) is not None

    model = Base.TCPAMaxwellBlockedModel(
        integration="paper_explicit",
        disc=Base.default_discretization(n_layers=1, n_phi=4, pre_steps=1),
    )
    try:
        model._validate_time_step(20.0)
    except ValueError:
        pass
    else:
        raise AssertionError("Un pas explicite instable aurait dû être refusé.")


def test_paper_physical_conventions() -> None:
    settings = dict(parametres.DEFAULT_SETTINGS)
    assert settings["bias_angle_profile"] == "paper_linear"
    assert settings["section_update_mode"] == "fixed"
    assert settings["axial_modulus_mode"] == "maxwell_sum"
    assert settings["prestrain_reference_mode"] == "elastic_tk_reference"
    assert settings["maxwell_anisotropy_mode"] == "paper_equal"
    assert settings["nylon_condition_mode"] == "bonded_linear"

    config = parametres.build_config(settings)
    assert np.isclose(config.mat.E_axial, config.mat.maxwell.E_total)
    model = Base.TCPAMaxwellBlockedModel(
        mat=config.mat,
        geom=config.geom,
        disc=Base.default_discretization(n_layers=3, n_phi=4, pre_steps=1),
    )
    expected = model.theta_f * model.R_centers / config.geom.Rout
    assert np.max(np.abs(model.theta_layers - expected)) < 1.0e-14


def test_tk_reference_state() -> None:
    settings = dict(parametres.DEFAULT_SETTINGS)
    settings.update({"eps": 0.8, "n_layers": 2, "n_phi": 8, "pre_steps": 4})
    config = parametres.build_config(settings)
    model = Base.TCPAMaxwellBlockedModel(
        mat=config.mat,
        geom=config.geom,
        disc=Base.default_discretization(n_layers=2, n_phi=8, pre_steps=4),
        integration=config.integration,
        prestrain_reference_mode=config.prestrain_reference_mode,
    )
    model.prestretch_to(config.eps)
    assert np.max(np.abs(model.sigma_reference)) > 0.0
    assert np.max(np.abs(model.sigma_i)) == 0.0
    reference_before = model.sigma_reference.copy()
    model.step(0.1, 1.0, h_target=model.h_blocked)
    assert np.max(np.abs(model.sigma_i)) > 0.0
    assert np.max(np.abs(model.sigma_reference - reference_before)) == 0.0


def test_nylon_physical_modes() -> None:
    settings = dict(parametres.DEFAULT_SETTINGS)
    settings.update({"nylon_condition_mode": "tension_only", "n_layers": 1, "n_phi": 4, "pre_steps": 1})
    config = parametres.build_config(settings)
    model = Base.TCPAMaxwellBlockedModel(
        mat=config.mat,
        geom=config.geom,
        disc=Base.default_discretization(n_layers=1, n_phi=4, pre_steps=1),
        prestrain_reference_mode=config.prestrain_reference_mode,
    )
    assert model._next_nylon_axial_force(-0.01, 1.0) == 0.0


def test_measured_pressure_csv() -> None:
    raw = "Temps (s);Pression (bar)\n0,0;-0,1\n1,0;0,9\n2,0;1,9\n".encode("utf-8")
    columns = pression.parse_uploaded_numeric_csv(raw)
    payload = pression.measured_pressure_payload(
        columns,
        "Temps (s)",
        "Pression (bar)",
        "bar",
        subtract_initial=True,
    )
    assert np.allclose(payload["time"], [0.0, 1.0, 2.0])
    assert np.allclose(payload["pressure_MPa"], [0.0, 0.1, 0.2])
    time = np.arange(0.0, 9.0)
    pressure_values = np.array([0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0])
    assert pression.estimate_measured_period(time, pressure_values) == 3.0


def test_result_csv_exports() -> None:
    data = {
        "time": np.array([0.0, 1.0, 2.0]),
        "pressure_MPa": np.array([0.0, 0.5, 1.0]),
        "force_total_mN": np.array([100.0, 120.0, 140.0]),
        "hold_start_index": 1,
    }
    csv_payload = affichage.mapping_to_csv_bytes(data)
    rows = list(csv.DictReader(StringIO(csv_payload.decode("utf-8-sig")), delimiter=";"))
    assert len(rows) == 3
    assert rows[1]["time"] == "1.0"
    assert rows[1]["force_total_mN"] == "120.0"
    assert rows[1]["hold_start_index"] == "1"

    cases = [{"label": "cas test", "data": data, "period": 2.0, "value": 0.8}]
    hysteresis_payload = affichage.hysteresis_to_csv_bytes(cases, [1], "prestrain")
    hysteresis_rows = list(
        csv.DictReader(StringIO(hysteresis_payload.decode("utf-8-sig")), delimiter=";")
    )
    assert len(hysteresis_rows) == 3
    assert hysteresis_rows[0]["case_label"] == "cas test"
    assert hysteresis_rows[0]["cycle"] == "1"
    assert np.isclose(float(hysteresis_rows[-1]["pressure_psi"]), 1.0 / parametres.PSI_TO_MPA)


def main() -> None:
    tests = (
        test_pressure_histories,
        test_blocked_equilibrium,
        test_suspended_equation_24,
        test_input_guards,
        test_paper_physical_conventions,
        test_tk_reference_state,
        test_nylon_physical_modes,
        test_measured_pressure_csv,
        test_result_csv_exports,
    )
    for test in tests:
        test()
        print(f"OK - {test.__name__}")
    print(f"Tous les tests scientifiques sont validés avec Base {Base.MODEL_VERSION}.")


if __name__ == "__main__":
    main()
