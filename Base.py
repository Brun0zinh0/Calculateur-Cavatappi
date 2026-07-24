"""
Base.py - module autonome final pour le modele TCPA corrige.

    - les valeurs par defaut materiau/geometrie/discretisation,
    - la rotation de raideur corrigee de la V6,
    - le solveur generalized-Maxwell bloque,
    - les conventions de sortie de la V6 corrected,
    - les graphes principaux.

Un code d'appel peut simplement faire :

    import Base
    model, data = Base.run_blocked_actuation(eps=0.8)
    Base.plot_all(data)
"""

from __future__ import annotations

from dataclasses import dataclass, is_dataclass, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar


MODEL_VERSION = "2026.07.24-fixed-maxwell-prestrain-nylon-8"


# ---------------------------------------------------------------------------
# Valeurs par defaut et resultats
# ---------------------------------------------------------------------------


def default_maxwell_tensile_params(**overrides):
    values = {
        "E0": 6.36,
        "E1": 20.67,
        "eta1": 154.57,
        "E2": 5.98,
        "eta2": 977.79,
        "E3": 4.75,
        "eta3": 11044.83,
    }
    values.update({k: v for k, v in overrides.items() if v is not None})
    return SimpleNamespace(**values)


def default_material_params(maxwell=None, **overrides):
    values = {
        "E_axial": 37.76,
        "E_radius": 8.82,
        "G12": 7.24,
        "nu12": 0.205,
        "nu23": 0.422,
        "maxwell": default_maxwell_tensile_params() if maxwell is None else maxwell,
        "E_nylon": 3.69e3,
        "G_nylon": 0.79e3,
        "maxwell_anisotropy_mode": "paper_equal",
        "nylon_condition_mode": "bonded_linear",
        "nylon_axial_prestrain_coupling": 1.0,
        "nylon_axial_actuation_coupling": 1.0,
    }
    values.update({k: v for k, v in overrides.items() if v is not None})
    return SimpleNamespace(**values)


def default_geometry_params(**overrides):
    values = {
        "Rout": 1.0,
        "Rin": 0.4,
        "r_nylon": 0.77 / 2.0,
        "rho0": 2.16,
        "alpha0_deg": 10.53,
        "theta_f_deg": 37.91,
        "initial_length": 32.45,
        "bias_angle_profile": "paper_linear",
        "section_update_mode": "fixed",
    }
    values.update({k: v for k, v in overrides.items() if v is not None})
    return SimpleNamespace(**values)


def default_discretization(**overrides):
    values = {
        "n_layers": 18,
        "n_phi": 72,
        "pre_steps": 120,
        "dw_bracket": (-0.08, 0.08),
    }
    values.update({k: v for k, v in overrides.items() if v is not None})
    return SimpleNamespace(**values)


def default_simulation_config(**overrides):
    values = {
        "eps": 0.8,
        "n_cycles": 3,
        "Pmax": 1.3,
        "dt": 0.25,
        "n_layers": 4,
        "n_phi": 24,
        "pre_steps": 24,
        "integration": "exponential",
        "prestrain_reference_mode": "elastic_tk_reference",
        "flow_rate_mL_min": 10.0,
        "volume_mL": 1.50,
        "nonlinear_pressure": False,
        "mat": default_material_params(),
        "geom": default_geometry_params(),
    }
    values.update({k: v for k, v in overrides.items() if v is not None})
    return SimpleNamespace(**values)


def _maxwell_E(maxwell) -> np.ndarray:
    if hasattr(maxwell, "E"):
        return np.asarray(maxwell.E, dtype=float)
    return np.array([maxwell.E1, maxwell.E2, maxwell.E3], dtype=float)


def _maxwell_eta(maxwell) -> np.ndarray:
    if hasattr(maxwell, "eta"):
        return np.asarray(maxwell.eta, dtype=float)
    return np.array([maxwell.eta1, maxwell.eta2, maxwell.eta3], dtype=float)


def _maxwell_E_total(maxwell) -> float:
    if hasattr(maxwell, "E_total"):
        return float(maxwell.E_total)
    return float(maxwell.E0 + maxwell.E1 + maxwell.E2 + maxwell.E3)


def _maxwell_rates(maxwell) -> np.ndarray:
    return _maxwell_E(maxwell) / _maxwell_eta(maxwell)


def _maxwell_branch_count(maxwell) -> int:
    return int(len(_maxwell_E(maxwell)))


@dataclass
class HelixState:
    rho: float
    alpha: float
    h: float
    pressure: float = 0.0
    time: float = 0.0


@dataclass
class StepResult:
    time: float
    pressure: float
    rho: float
    alpha: float
    dw: float
    dv: float
    dkappa: float
    Ft: float
    Tt: float
    residual: float
    Ftube: float
    Mtube: float
    Ttube: float
    Fnylon: float
    Mnylon: float
    Tnylon: float
    axial_stretch: float
    Rin: float
    Rout: float


# ---------------------------------------------------------------------------
# Raideur anisotrope et rotation corrigee
# ---------------------------------------------------------------------------


def ti_stiffness_from_paper(
    E_axial: float,
    E_radius: float,
    G12: float,
    nu12: float,
    nu23: float,
) -> np.ndarray:
    """Raideur transverse-isotrope locale, ordre Voigt engineering shear."""
    nu21 = nu12 * E_radius / E_axial
    A = np.array(
        [
            [1.0, -2.0 * nu12, 0.0, 0.0],
            [0.0, 1.0, -nu12, -nu12],
            [-nu21, 1.0 - nu23, 0.0, 0.0],
            [0.0, -nu21, 1.0, -nu23],
        ],
        dtype=float,
    )
    b = np.array([E_axial, 0.0, 0.0, E_radius], dtype=float)
    C11, C12, C22, C23 = np.linalg.solve(A, b)

    C = np.zeros((6, 6), dtype=float)
    C[0, 0] = C11
    C[0, 1] = C[1, 0] = C12
    C[0, 2] = C[2, 0] = C12
    C[1, 1] = C22
    C[1, 2] = C[2, 1] = C23
    C[2, 2] = C22
    C[3, 3] = 0.5 * (C22 - C23)
    C[4, 4] = G12
    C[5, 5] = G12
    return C


VOIGT_PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def voigt_to_tensor(Cv: np.ndarray) -> np.ndarray:
    """Conversion corrigee V6 : pas de facteur shear supplementaire."""
    C4 = np.zeros((3, 3, 3, 3), dtype=float)
    for I, (i, j) in enumerate(VOIGT_PAIRS):
        for J, (k, l) in enumerate(VOIGT_PAIRS):
            val = Cv[I, J]
            for a, b in ((i, j), (j, i)):
                for c, d in ((k, l), (l, k)):
                    C4[a, b, c, d] = val
    return C4


def tensor_to_voigt(C4: np.ndarray) -> np.ndarray:
    Cv = np.zeros((6, 6), dtype=float)
    for I, (i, j) in enumerate(VOIGT_PAIRS):
        for J, (k, l) in enumerate(VOIGT_PAIRS):
            Cv[I, J] = C4[i, j, k, l]
    return Cv


def rotate_stiffness_bias(C_local: np.ndarray, theta: float) -> np.ndarray:
    """Rotation tensorielle corrigee de la raideur locale vers [s, phi, r]."""
    c, s = np.cos(theta), np.sin(theta)
    q = np.array(
        [
            [c, -s, 0.0],
            [s, c, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    C4_local = voigt_to_tensor(C_local)
    C4_global = np.einsum(
        "iA,jB,kC,lD,ABCD->ijkl",
        q,
        q,
        q,
        q,
        C4_local,
        optimize=True,
    )
    Cbar = tensor_to_voigt(C4_global)
    return 0.5 * (Cbar + Cbar.T)


def effective_poissons(Cbar: np.ndarray) -> Tuple[float, float, float, float]:
    idx = [0, 1, 2, 5]
    S = np.linalg.inv(Cbar[np.ix_(idx, idx)])
    EL = 1.0 / S[0, 0]
    v = S @ np.array([EL, 0.0, 0.0, 0.0])
    return EL, -v[1], -v[2], -v[3]


def _normalize_integration_name(integration: str) -> str:
    aliases = {
        "paper_incremental": "paper_explicit",
        "paper": "paper_explicit",
        "explicit": "paper_explicit",
        "stable": "exponential",
    }
    normalized = aliases.get(str(integration), str(integration))
    if normalized not in {"paper_explicit", "exponential"}:
        raise ValueError("integration must be 'paper_explicit' or 'exponential'.")
    return normalized


def _time_grid_with_events(total_time: float, dt: float, events=()) -> np.ndarray:
    """Return a bounded grid containing the requested physical transitions."""
    if total_time < 0.0 or dt <= 0.0:
        raise ValueError("total_time must be non-negative and dt must be positive.")
    regular = np.arange(0.0, total_time, dt, dtype=float)
    values = [regular, np.array([0.0, total_time], dtype=float)]
    physical_events = np.asarray(list(events), dtype=float)
    if physical_events.size:
        physical_events = physical_events[np.isfinite(physical_events)]
        physical_events = physical_events[(physical_events >= 0.0) & (physical_events <= total_time)]
        values.append(physical_events)
    grid = np.unique(np.concatenate(values))
    grid[np.isclose(grid, total_time, rtol=0.0, atol=1e-12)] = total_time
    return grid[(grid >= 0.0) & (grid <= total_time)]


# ---------------------------------------------------------------------------
# Solveur principal
# ---------------------------------------------------------------------------


class TCPAMaxwellBlockedModel:
    """Modele d'actionnement bloque, autonome, avec corrections V6."""

    def __init__(
        self,
        mat=None,
        geom=None,
        disc=None,
        integration: str = "exponential",
        prestrain_reference_mode: str = "elastic_tk_reference",
    ):
        if mat is None:
            mat = default_material_params()
        if geom is None:
            geom = default_geometry_params()
        if disc is None:
            disc = default_discretization()
        self.mat = mat
        self.geom = geom
        self.disc = disc
        self.integration = _normalize_integration_name(integration)
        requested_prestrain_mode = str(prestrain_reference_mode)
        if requested_prestrain_mode != "elastic_tk_reference":
            raise ValueError("Beta supports only the conserved elastic prestrain reference.")
        self.prestrain_reference_mode = "elastic_tk_reference"
        self.maxwell_anisotropy_mode = str(getattr(mat, "maxwell_anisotropy_mode", "paper_equal"))
        requested_nylon_mode = str(getattr(mat, "nylon_condition_mode", "bonded_linear"))
        if requested_nylon_mode != "bonded_linear":
            raise ValueError("Beta supports only bonded bilateral linear nylon.")
        self.nylon_condition_mode = "bonded_linear"
        self._building_reference_state = False
        self.nylon_axial_prestrain_coupling = float(getattr(mat, "nylon_axial_prestrain_coupling", 1.0))
        self.nylon_axial_actuation_coupling = float(getattr(mat, "nylon_axial_actuation_coupling", 1.0))
        self.section_update_mode = str(getattr(geom, "section_update_mode", "fixed"))
        self.bias_angle_profile = str(getattr(geom, "bias_angle_profile", "paper_linear"))
        for name, value in (
            ("nylon_axial_prestrain_coupling", self.nylon_axial_prestrain_coupling),
            ("nylon_axial_actuation_coupling", self.nylon_axial_actuation_coupling),
        ):
            if not np.isclose(value, 1.0, rtol=0.0, atol=1.0e-12):
                raise ValueError(f"Beta fixes {name} to 1.0.")
        if self.section_update_mode not in {"fixed", "updated"}:
            raise ValueError("section_update_mode must be 'fixed' or 'updated'.")
        if self.bias_angle_profile not in {"paper_linear", "uniform_twist"}:
            raise ValueError("bias_angle_profile must be 'paper_linear' or 'uniform_twist'.")
        if self.maxwell_anisotropy_mode not in {"paper_equal", "axial_test_only"}:
            raise ValueError("maxwell_anisotropy_mode must be 'paper_equal' or 'axial_test_only'.")

        self._validate_inputs()

        self.alpha0 = np.deg2rad(geom.alpha0_deg)
        self.theta_f = np.deg2rad(geom.theta_f_deg)
        self.h0 = geom.rho0 * np.tan(self.alpha0)
        self.turns = geom.initial_length / (2.0 * np.pi * self.h0)
        self.helix = HelixState(rho=geom.rho0, alpha=self.alpha0, h=self.h0)
        self.h_blocked = self.h0

        self.R_edges = np.linspace(geom.Rin, geom.Rout, disc.n_layers + 1)
        self.R_centers = 0.5 * (self.R_edges[:-1] + self.R_edges[1:])
        self.dR = np.diff(self.R_edges)
        self.phi = np.linspace(0.0, 2.0 * np.pi, disc.n_phi, endpoint=False)
        self.dphi = 2.0 * np.pi / disc.n_phi

        self.C_local_total = ti_stiffness_from_paper(
            mat.E_axial,
            mat.E_radius,
            mat.G12,
            mat.nu12,
            mat.nu23,
        )
        self.theta_layers = self._bias_angles(self.R_centers, geom.Rout)
        self.C_total: List[np.ndarray] = []
        self.C0: List[np.ndarray] = []
        self.Ci: List[List[np.ndarray]] = []
        self.vbar: List[Tuple[float, float, float, float]] = []
        self.n_maxwell = _maxwell_branch_count(mat.maxwell)
        self._rebuild_section_properties()

        shape = (disc.n_layers, disc.n_phi, 6)
        self.sigma_reference = np.zeros(shape)
        self.sigma0 = np.zeros(shape)
        self.sigma_i = np.zeros((self.n_maxwell,) + shape)
        self.sigma_total = np.zeros(shape)
        self.Fnylon = 0.0
        self.Mnylon = 0.0
        self.Tnylon = 0.0
        self.axial_stretch = 1.0
        self.history: List[StepResult] = []

    def _validate_inputs(self) -> None:
        positive = {
            "E_axial": self.mat.E_axial,
            "E_radius": self.mat.E_radius,
            "G12": self.mat.G12,
            "E_nylon": self.mat.E_nylon,
            "G_nylon": self.mat.G_nylon,
            "Rin": self.geom.Rin,
            "Rout": self.geom.Rout,
            "rho0": self.geom.rho0,
            "initial_length": self.geom.initial_length,
        }
        for name, value in positive.items():
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and strictly positive.")
        if self.geom.Rin >= self.geom.Rout:
            raise ValueError("Rin must be smaller than Rout.")
        if self.geom.rho0 <= self.geom.Rout:
            raise ValueError("rho0 must be larger than Rout.")
        if self.geom.r_nylon < 0.0 or self.geom.r_nylon > self.geom.Rin:
            raise ValueError("The nylon radius must lie between zero and Rin.")
        if not 0.0 < self.geom.alpha0_deg < 90.0:
            raise ValueError("alpha0_deg must lie strictly between 0 and 90 degrees.")
        if not 0.0 <= self.geom.theta_f_deg < 90.0:
            raise ValueError("theta_f_deg must lie between 0 and 90 degrees.")
        if int(self.disc.n_layers) < 1 or int(self.disc.n_phi) < 4:
            raise ValueError("The mesh requires at least one radial layer and four angular divisions.")
        branch_E = _maxwell_E(self.mat.maxwell)
        branch_eta = _maxwell_eta(self.mat.maxwell)
        if self.mat.maxwell.E0 < 0.0 or np.any(branch_E < 0.0):
            raise ValueError("Maxwell moduli must be non-negative.")
        if _maxwell_E_total(self.mat.maxwell) <= 0.0:
            raise ValueError("At least one Maxwell modulus must be strictly positive.")
        if np.any((branch_E > 0.0) & (branch_eta <= 0.0)):
            raise ValueError("Each active Maxwell branch requires a strictly positive viscosity.")
        C_local = ti_stiffness_from_paper(
            self.mat.E_axial,
            self.mat.E_radius,
            self.mat.G12,
            self.mat.nu12,
            self.mat.nu23,
        )
        eig_min = float(np.linalg.eigvalsh(C_local).min())
        if not np.isfinite(eig_min) or eig_min <= 1e-10:
            raise ValueError("The elastic constants produce a non-physical stiffness matrix.")

    def _rebuild_section_properties(self) -> None:
        Etotal = _maxwell_E_total(self.mat.maxwell)
        self.C_total = []
        self.C0 = []
        self.Ci = []
        self.vbar = []
        if self.maxwell_anisotropy_mode == "axial_test_only":
            axial_projector = np.zeros((6, 6), dtype=float)
            axial_projector[0, 0] = 1.0
            Ci_local = [Ei * axial_projector for Ei in _maxwell_E(self.mat.maxwell)]
            C0_local = self.C_local_total - sum(Ci_local, np.zeros((6, 6), dtype=float))
            if float(np.linalg.eigvalsh(C0_local).min()) <= 1.0e-10:
                raise ValueError(
                    "The axial-only Maxwell decomposition leaves a non-physical permanent stiffness matrix."
                )
        else:
            C0_local = (self.mat.maxwell.E0 / Etotal) * self.C_local_total
            Ci_local = [(Ei / Etotal) * self.C_local_total for Ei in _maxwell_E(self.mat.maxwell)]
        for theta_j in self.theta_layers:
            Cbar = rotate_stiffness_bias(self.C_local_total, float(theta_j))
            self.C_total.append(Cbar)
            self.C0.append(rotate_stiffness_bias(C0_local, float(theta_j)))
            self.Ci.append([rotate_stiffness_bias(Ci, float(theta_j)) for Ci in Ci_local])
            self.vbar.append(effective_poissons(Cbar))

    def _bias_angles(self, radii: np.ndarray, outer_radius: float) -> np.ndarray:
        radial_fraction = np.asarray(radii, dtype=float) / float(outer_radius)
        if self.bias_angle_profile == "paper_linear":
            return radial_fraction * self.theta_f
        return np.arctan(radial_fraction * np.tan(self.theta_f))

    def _validate_time_step(self, dt: float) -> None:
        if dt < 0.0 or not np.isfinite(dt):
            raise ValueError("dt must be finite and non-negative.")
        if dt == 0.0 or self.integration != "paper_explicit":
            return
        active = _maxwell_E(self.mat.maxwell) > 0.0
        if not np.any(active):
            return
        tau_min = float(np.min(_maxwell_eta(self.mat.maxwell)[active] / _maxwell_E(self.mat.maxwell)[active]))
        if dt >= 2.0 * tau_min:
            raise ValueError(
                f"Explicit Maxwell integration is unstable for dt={dt:g} s; "
                f"use dt < {2.0 * tau_min:.6g} s or select exponential integration."
            )

    def _nylon_axial_coupling(self, h_target: float) -> float:
        if abs(h_target - self.h_blocked) > 1e-10:
            return self.nylon_axial_prestrain_coupling
        return self.nylon_axial_actuation_coupling

    def _next_nylon_axial_force(self, dw: float, axial_coupling: float) -> float:
        raw_force = self.Fnylon + axial_coupling * np.pi * self.mat.E_nylon * self.geom.r_nylon**2 * dw
        return float(raw_force)

    def _new_geometry_from_dw_and_h(self, dw: float, h_target: float) -> Tuple[float, float]:
        l_old = self.helix.rho / np.cos(self.helix.alpha)
        l_new = (1.0 + dw) * l_old
        if l_new <= abs(h_target):
            raise ValueError("Inadmissible geometry: centerline length <= pitch projection.")
        rho_new = np.sqrt(l_new * l_new - h_target * h_target)
        alpha_new = np.arctan2(h_target, rho_new)
        return rho_new, alpha_new

    def _kinematic_increments(self, dw: float, h_target: float) -> Tuple[float, float, float, float, float]:
        rho_old, alpha_old = self.helix.rho, self.helix.alpha
        rho_new, alpha_new = self._new_geometry_from_dw_and_h(dw, h_target)
        dv = np.sin(2.0 * alpha_new) / (2.0 * rho_new) - np.sin(2.0 * alpha_old) / (2.0 * rho_old)
        dkappa = (np.cos(alpha_new) ** 2) / rho_new - (np.cos(alpha_old) ** 2) / rho_old
        return rho_new, alpha_new, dw, dv, dkappa

    def _layer_AB_mu(self, j: int, C_layers: Optional[List[np.ndarray]] = None) -> Tuple[float, float, float]:
        C = self.C_total[j] if C_layers is None else C_layers[j]
        C12, C13 = C[0, 1], C[0, 2]
        C22, C26 = C[1, 1], C[1, 5]
        C33, C36 = C[2, 2], C[2, 5]
        mu = np.sqrt(max(C22 / C33, 1e-14))
        A = (C26 - 2.0 * C36) / (4.0 * C33 - C22)
        B = (C12 - C13) / (C33 - C22)
        return A, B, mu

    @staticmethod
    def _u_base(R: float, A: float, B: float, mu: float, dv: float, dw: float):
        coeff_u = np.array([R**mu, R ** (-mu)], dtype=float)
        known_u = A * dv * R * R + B * dw * R
        coeff_du = np.array([mu * R ** (mu - 1.0), -mu * R ** (-mu - 1.0)], dtype=float)
        known_du = 2.0 * A * dv * R + B * dw
        return coeff_u, known_u, coeff_du, known_du

    def _radial_stress_axisym_coeff(
        self,
        j: int,
        R: float,
        A: float,
        B: float,
        mu: float,
        dv: float,
        dw: float,
        Cvisc_r: float = 0.0,
        C_layers: Optional[List[np.ndarray]] = None,
    ) -> Tuple[np.ndarray, float]:
        C = self.C_total[j] if C_layers is None else C_layers[j]
        cu, ku, cdu, kdu = self._u_base(R, A, B, mu, dv, dw)
        coeff = np.zeros((2, 6), dtype=float)
        known = np.array([dw, ku / R, kdu, 0.0, 0.0, dv * R], dtype=float)
        for k in range(2):
            strain_coeff = np.array([0.0, cu[k] / R, cdu[k], 0.0, 0.0, 0.0], dtype=float)
            coeff[k] = C @ strain_coeff
        sig_known = C @ known
        return coeff[:, 2], sig_known[2] + Cvisc_r

    def _algorithmic_data(self, dt: float) -> Tuple[List[np.ndarray], np.ndarray]:
        """Return the consistent tangent and stress-history increment."""
        if self._building_reference_state:
            return [C.copy() for C in self.C_total], np.zeros_like(self.sigma_total)
        rates = _maxwell_rates(self.mat.maxwell)
        history = np.zeros_like(self.sigma_total)
        if self.integration == "paper_explicit":
            for ib in range(self.n_maxwell):
                history -= dt * rates[ib] * self.sigma_i[ib]
            return [C.copy() for C in self.C_total], history

        factors = np.ones(self.n_maxwell, dtype=float)
        decays = np.ones(self.n_maxwell, dtype=float)
        for ib, rate in enumerate(rates):
            if rate > 0.0 and dt > 0.0:
                decays[ib] = np.exp(-dt * rate)
                factors[ib] = (1.0 - decays[ib]) / (dt * rate)
            history += (decays[ib] - 1.0) * self.sigma_i[ib]
        tangents = [
            self.C0[j] + sum((factors[ib] * self.Ci[j][ib] for ib in range(self.n_maxwell)), np.zeros((6, 6)))
            for j in range(self.disc.n_layers)
        ]
        return tangents, history

    def _update_maxwell_branches(self, j: int, k: int, de: np.ndarray, dt: float) -> np.ndarray:
        rates = _maxwell_rates(self.mat.maxwell)
        sigma_i_point = np.empty_like(self.sigma_i[:, j, k])
        for ib in range(self.n_maxwell):
            if self.integration == "paper_explicit":
                ds = self.Ci[j][ib] @ de - dt * rates[ib] * self.sigma_i[ib, j, k]
                sigma_i_point[ib] = self.sigma_i[ib, j, k] + ds
            elif self.integration == "exponential":
                if rates[ib] <= 0.0:
                    sigma_i_point[ib] = self.sigma_i[ib, j, k] + self.Ci[j][ib] @ de
                else:
                    tau = 1.0 / rates[ib]
                    a = np.exp(-dt / tau)
                    factor = tau / dt * (1.0 - a) if dt > 0 else 1.0
                    sigma_i_point[ib] = a * self.sigma_i[ib, j, k] + factor * (self.Ci[j][ib] @ de)
            else:
                raise ValueError(
                    "integration must be 'paper_explicit' or 'exponential'."
                )
        return sigma_i_point

    def _solve_radial_constants(
        self,
        dw: float,
        dv: float,
        dP: float,
        dt: float,
        C_algorithmic: List[np.ndarray],
        history_increment: np.ndarray,
    ) -> np.ndarray:
        n = self.disc.n_layers
        A_mat = np.zeros((2 * n, 2 * n), dtype=float)
        b_vec = np.zeros(2 * n, dtype=float)
        CSmean = history_increment.mean(axis=1)

        def cols(j: int) -> Tuple[int, int]:
            return 2 * j, 2 * j + 1

        row = 0
        j = 0
        Aj, Bj, muj = self._layer_AB_mu(j, C_algorithmic)
        c_sig, k_sig = self._radial_stress_axisym_coeff(
            j, self.R_edges[0], Aj, Bj, muj, dv, dw, CSmean[j, 2], C_algorithmic
        )
        c0, c1 = cols(j)
        A_mat[row, c0 : c1 + 1] = c_sig
        b_vec[row] = -dP - k_sig
        row += 1

        for j in range(n - 1):
            Rb = self.R_edges[j + 1]
            Aj, Bj, muj = self._layer_AB_mu(j, C_algorithmic)
            Ak, Bk, muk = self._layer_AB_mu(j + 1, C_algorithmic)
            cu_j, ku_j, _, _ = self._u_base(Rb, Aj, Bj, muj, dv, dw)
            cu_k, ku_k, _, _ = self._u_base(Rb, Ak, Bk, muk, dv, dw)
            j0, j1 = cols(j)
            k0, k1 = cols(j + 1)
            A_mat[row, j0 : j1 + 1] = cu_j
            A_mat[row, k0 : k1 + 1] = -cu_k
            b_vec[row] = ku_k - ku_j
            row += 1

            cs_j, ks_j = self._radial_stress_axisym_coeff(
                j, Rb, Aj, Bj, muj, dv, dw, CSmean[j, 2], C_algorithmic
            )
            cs_k, ks_k = self._radial_stress_axisym_coeff(
                j + 1, Rb, Ak, Bk, muk, dv, dw, CSmean[j + 1, 2], C_algorithmic
            )
            A_mat[row, j0 : j1 + 1] = cs_j
            A_mat[row, k0 : k1 + 1] = -cs_k
            b_vec[row] = ks_k - ks_j
            row += 1

        j = n - 1
        Aj, Bj, muj = self._layer_AB_mu(j, C_algorithmic)
        c_sig, k_sig = self._radial_stress_axisym_coeff(
            j, self.R_edges[-1], Aj, Bj, muj, dv, dw, CSmean[j, 2], C_algorithmic
        )
        c0, c1 = cols(j)
        A_mat[row, c0 : c1 + 1] = c_sig
        b_vec[row] = -k_sig
        return np.linalg.solve(A_mat, b_vec).reshape(n, 2)

    def _u_du_layer(
        self,
        j: int,
        R: float,
        coeffs: np.ndarray,
        dv: float,
        dw: float,
        C_layers: Optional[List[np.ndarray]] = None,
    ) -> Tuple[float, float]:
        A, B, mu = self._layer_AB_mu(j, C_layers)
        C1, C2 = coeffs[j]
        u = C1 * R**mu + C2 * R ** (-mu) + A * dv * R * R + B * dw * R
        du = C1 * mu * R ** (mu - 1.0) - C2 * mu * R ** (-mu - 1.0) + 2.0 * A * dv * R + B * dw
        return u, du

    def _trial_section_geometry(
        self,
        coeffs: np.ndarray,
        dv: float,
        dw: float,
        C_algorithmic: List[np.ndarray],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        if self.section_update_mode == "fixed":
            return (
                self.R_edges.copy(),
                self.R_centers.copy(),
                self.dR.copy(),
                self.theta_layers.copy(),
            )
        u_edges = np.empty_like(self.R_edges)
        for edge in range(len(self.R_edges)):
            layers = [0] if edge == 0 else [self.disc.n_layers - 1] if edge == self.disc.n_layers else [edge - 1, edge]
            values = [
                self._u_du_layer(j, self.R_edges[edge], coeffs, dv, dw, C_algorithmic)[0]
                for j in layers
            ]
            u_edges[edge] = float(np.mean(values))
        edges_new = self.R_edges + u_edges
        if not np.all(np.isfinite(edges_new)) or edges_new[0] <= 0.0 or np.any(np.diff(edges_new) <= 1e-9):
            raise ValueError("The radial update produced an inadmissible tube cross-section.")
        centers_new = 0.5 * (edges_new[:-1] + edges_new[1:])
        dR_new = np.diff(edges_new)

        if self.bias_angle_profile == "paper_linear":
            theta_new = self._bias_angles(centers_new, edges_new[-1])
            return edges_new, centers_new, dR_new, theta_new

        theta_new = np.empty_like(self.theta_layers)
        for j, R in enumerate(self.R_centers):
            u, _ = self._u_du_layer(j, R, coeffs, dv, dw, C_algorithmic)
            theta = self.theta_layers[j]
            axial_component = (1.0 + dw) * np.cos(theta)
            hoop_component = dv * R * np.cos(theta) + (1.0 + u / R) * np.sin(theta)
            theta_new[j] = np.arctan2(hoop_component, axial_component)
        if np.any(np.abs(theta_new) >= np.deg2rad(89.9)):
            raise ValueError("The material bias angle left the admissible range.")
        return edges_new, centers_new, dR_new, theta_new

    def _trial_state(self, dw: float, dP: float, dt: float, h_target: float) -> Dict[str, object]:
        rho_new, alpha_new, dw, dv, dkappa = self._kinematic_increments(dw, h_target)
        C_algorithmic, history_increment = self._algorithmic_data(dt)
        coeffs = self._solve_radial_constants(dw, dv, dP, dt, C_algorithmic, history_increment)
        R_edges_new, R_centers_new, dR_new, theta_layers_new = self._trial_section_geometry(
            coeffs, dv, dw, C_algorithmic
        )
        if R_edges_new[-1] >= rho_new:
            raise ValueError("The deformed tube cross-section intersects the helix axis.")
        sigma_reference_new = np.empty_like(self.sigma_reference)
        sigma0_new = np.empty_like(self.sigma0)
        sigma_i_new = np.empty_like(self.sigma_i)
        sigma_total_new = np.empty_like(self.sigma_total)
        K_old = (np.cos(self.helix.alpha) ** 2) / self.helix.rho

        for j, R in enumerate(self.R_centers):
            C0 = self.C0[j]
            _, v12b, v13b, v14b = self.vbar[j]
            u, du = self._u_du_layer(j, R, coeffs, dv, dw, C_algorithmic)
            for k, Phi in enumerate(self.phi):
                denom = 1.0 + K_old * R * np.cos(Phi)
                curv = (dkappa * R * np.cos(Phi) + u * K_old * np.cos(Phi)) / denom
                eps_r = du - v12b * curv
                eps_phi = u / R - v13b * curv
                eps_s = dw + curv
                gamma_sphi = dv * R / denom - v14b * curv
                de = np.array([eps_s, eps_phi, eps_r, 0.0, 0.0, gamma_sphi], dtype=float)
                if self._building_reference_state:
                    sigma_reference_new[j, k] = self.sigma_reference[j, k] + self.C_total[j] @ de
                    sigma0_new[j, k] = self.sigma0[j, k]
                    sigma_i_new[:, j, k] = self.sigma_i[:, j, k]
                else:
                    sigma_reference_new[j, k] = self.sigma_reference[j, k]
                    sigma0_new[j, k] = self.sigma0[j, k] + C0 @ de
                    sigma_i_new[:, j, k] = self._update_maxwell_branches(j, k, de, dt)
                sigma_total_new[j, k] = (
                    sigma_reference_new[j, k] + sigma0_new[j, k] + sigma_i_new[:, j, k].sum(axis=0)
                )

        Ftube = 0.0
        Mtube = 0.0
        Ttube = 0.0
        for j, R in enumerate(R_centers_new):
            weight_R = R * dR_new[j] * self.dphi
            for k, Phi in enumerate(self.phi):
                sig = sigma_total_new[j, k]
                Ftube += sig[0] * weight_R
                Mtube += sig[0] * (R * np.cos(Phi)) * weight_R
                Ttube += sig[5] * R * weight_R

        axial_coupling = self._nylon_axial_coupling(h_target)
        rn = self.geom.r_nylon
        Fny = self._next_nylon_axial_force(dw, axial_coupling)
        Mny = self.Mnylon + 0.25 * np.pi * self.mat.E_nylon * rn**4 * dkappa
        Tny = self.Tnylon + 0.5 * np.pi * self.mat.G_nylon * rn**4 * dv

        Ares = Ftube + Fny
        Bres = Mtube + Mny
        Cres = Ttube + Tny
        sin_a, cos_a = np.sin(alpha_new), np.cos(alpha_new)
        if abs(sin_a) < 1e-8 or abs(cos_a) < 1e-8:
            Ft = np.nan
            Tt = np.nan
            residual = np.inf
        else:
            Ft = Ares / sin_a
            Tt = (Bres + Ft * rho_new * sin_a) / cos_a
            residual = (Tt * sin_a + Ft * rho_new * cos_a) - Cres

        return {
            "rho_new": rho_new,
            "alpha_new": alpha_new,
            "dw": dw,
            "dv": dv,
            "dkappa": dkappa,
            "sigma_reference_new": sigma_reference_new,
            "sigma0_new": sigma0_new,
            "sigma_i_new": sigma_i_new,
            "sigma_total_new": sigma_total_new,
            "Fny": Fny,
            "Mny": Mny,
            "Tny": Tny,
            "Ftube": Ftube,
            "Mtube": Mtube,
            "Ttube": Ttube,
            "Ft": Ft,
            "Tt": Tt,
            "residual": residual,
            "R_edges_new": R_edges_new,
            "R_centers_new": R_centers_new,
            "dR_new": dR_new,
            "theta_layers_new": theta_layers_new,
            "axial_stretch_new": self.axial_stretch * (1.0 + dw),
        }

    def _trial_state_from_geometry(
        self,
        dw: float,
        dP: float,
        dt: float,
        rho_new: float,
        alpha_new: float,
        axial_coupling: float = 1.0,
    ) -> Dict[str, object]:
        if rho_new <= self.R_edges[-1]:
            raise ValueError("Inadmissible geometry: helix radius <= tube outer radius.")
        if not 0.0 < alpha_new < 0.5 * np.pi:
            raise ValueError("Inadmissible geometry: helix angle outside (0, 90 deg).")

        rho_old, alpha_old = self.helix.rho, self.helix.alpha
        dv = np.sin(2.0 * alpha_new) / (2.0 * rho_new) - np.sin(2.0 * alpha_old) / (2.0 * rho_old)
        dkappa = (np.cos(alpha_new) ** 2) / rho_new - (np.cos(alpha_old) ** 2) / rho_old

        C_algorithmic, history_increment = self._algorithmic_data(dt)
        coeffs = self._solve_radial_constants(dw, dv, dP, dt, C_algorithmic, history_increment)
        R_edges_new, R_centers_new, dR_new, theta_layers_new = self._trial_section_geometry(
            coeffs, dv, dw, C_algorithmic
        )
        if R_edges_new[-1] >= rho_new:
            raise ValueError("The deformed tube cross-section intersects the helix axis.")
        sigma_reference_new = np.empty_like(self.sigma_reference)
        sigma0_new = np.empty_like(self.sigma0)
        sigma_i_new = np.empty_like(self.sigma_i)
        sigma_total_new = np.empty_like(self.sigma_total)
        K_old = (np.cos(self.helix.alpha) ** 2) / self.helix.rho

        for j, R in enumerate(self.R_centers):
            C0 = self.C0[j]
            _, v12b, v13b, v14b = self.vbar[j]
            u, du = self._u_du_layer(j, R, coeffs, dv, dw, C_algorithmic)
            for k, Phi in enumerate(self.phi):
                denom = 1.0 + K_old * R * np.cos(Phi)
                curv = (dkappa * R * np.cos(Phi) + u * K_old * np.cos(Phi)) / denom
                eps_r = du - v12b * curv
                eps_phi = u / R - v13b * curv
                eps_s = dw + curv
                gamma_sphi = dv * R / denom - v14b * curv
                de = np.array([eps_s, eps_phi, eps_r, 0.0, 0.0, gamma_sphi], dtype=float)
                if self._building_reference_state:
                    sigma_reference_new[j, k] = self.sigma_reference[j, k] + self.C_total[j] @ de
                    sigma0_new[j, k] = self.sigma0[j, k]
                    sigma_i_new[:, j, k] = self.sigma_i[:, j, k]
                else:
                    sigma_reference_new[j, k] = self.sigma_reference[j, k]
                    sigma0_new[j, k] = self.sigma0[j, k] + C0 @ de
                    sigma_i_new[:, j, k] = self._update_maxwell_branches(j, k, de, dt)
                sigma_total_new[j, k] = (
                    sigma_reference_new[j, k] + sigma0_new[j, k] + sigma_i_new[:, j, k].sum(axis=0)
                )

        Ftube = 0.0
        Mtube = 0.0
        Ttube = 0.0
        for j, R in enumerate(R_centers_new):
            weight_R = R * dR_new[j] * self.dphi
            for k, Phi in enumerate(self.phi):
                sig = sigma_total_new[j, k]
                Ftube += sig[0] * weight_R
                Mtube += sig[0] * (R * np.cos(Phi)) * weight_R
                Ttube += sig[5] * R * weight_R

        rn = self.geom.r_nylon
        Fny = self._next_nylon_axial_force(dw, axial_coupling)
        Mny = self.Mnylon + 0.25 * np.pi * self.mat.E_nylon * rn**4 * dkappa
        Tny = self.Tnylon + 0.5 * np.pi * self.mat.G_nylon * rn**4 * dv

        return {
            "rho_new": rho_new,
            "alpha_new": alpha_new,
            "dw": dw,
            "dv": dv,
            "dkappa": dkappa,
            "sigma_reference_new": sigma_reference_new,
            "sigma0_new": sigma0_new,
            "sigma_i_new": sigma_i_new,
            "sigma_total_new": sigma_total_new,
            "Fny": Fny,
            "Mny": Mny,
            "Tny": Tny,
            "Ftube": Ftube,
            "Mtube": Mtube,
            "Ttube": Ttube,
            "Ft": np.nan,
            "Tt": np.nan,
            "residual": np.nan,
            "R_edges_new": R_edges_new,
            "R_centers_new": R_centers_new,
            "dR_new": dR_new,
            "theta_layers_new": theta_layers_new,
            "axial_stretch_new": self.axial_stretch * (1.0 + dw),
        }

    def _find_dw(self, dP: float, dt: float, h_target: float) -> float:
        lo, hi = self.disc.dw_bracket

        def f(x: float) -> float:
            try:
                return float(self._trial_state(x, dP, dt, h_target)["residual"])
            except Exception:
                return np.nan

        grid = np.linspace(lo, hi, 65)
        vals = np.array([f(x) for x in grid])
        roots = []
        for x, value in zip(grid, vals):
            if np.isfinite(value) and abs(value) <= 1e-10:
                roots.append(float(x))
        for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if np.isfinite(fa) and np.isfinite(fb) and fa * fb <= 0:
                roots.append(float(brentq(lambda z: f(z), a, b, xtol=1e-10, rtol=1e-9, maxiter=100)))
        if roots:
            return min(roots, key=abs)

        def obj(x: float) -> float:
            y = f(x)
            return 1e100 if not np.isfinite(y) else y * y

        res = minimize_scalar(obj, bounds=(lo, hi), method="bounded", options={"xatol": 1e-8})
        residual = abs(f(float(res.x))) if res.success else np.inf
        span = hi - lo
        at_bound = min(float(res.x) - lo, hi - float(res.x)) <= 1e-5 * span
        if not res.success or not np.isfinite(residual) or residual > 1e-6 or at_bound:
            raise RuntimeError(
                f"Could not solve blocked equilibrium: residual={residual:.3e} N mm, dw={float(res.x):.6g}."
            )
        return float(res.x)

    def _commit_trial(self, trial: Dict[str, object]) -> None:
        self.sigma_reference = trial["sigma_reference_new"]
        self.sigma0 = trial["sigma0_new"]
        self.sigma_i = trial["sigma_i_new"]
        self.sigma_total = trial["sigma_total_new"]
        self.Fnylon = float(trial["Fny"])
        self.Mnylon = float(trial["Mny"])
        self.Tnylon = float(trial["Tny"])
        self.axial_stretch = float(trial["axial_stretch_new"])
        self.R_edges = np.asarray(trial["R_edges_new"], dtype=float)
        self.R_centers = np.asarray(trial["R_centers_new"], dtype=float)
        self.dR = np.asarray(trial["dR_new"], dtype=float)
        self.theta_layers = np.asarray(trial["theta_layers_new"], dtype=float)
        self._rebuild_section_properties()

    def step(self, pressure_new: float, dt: float, h_target: Optional[float] = None) -> StepResult:
        self._validate_time_step(dt)
        if not np.isfinite(pressure_new) or pressure_new < 0.0:
            raise ValueError("pressure must be finite and non-negative.")
        if h_target is None:
            h_target = self.h_blocked
        dP = pressure_new - self.helix.pressure
        dw = self._find_dw(dP, dt, h_target)
        trial = self._trial_state(dw, dP, dt, h_target)

        self._commit_trial(trial)
        self.helix = HelixState(
            rho=float(trial["rho_new"]),
            alpha=float(trial["alpha_new"]),
            h=h_target,
            pressure=pressure_new,
            time=self.helix.time + dt,
        )

        out = StepResult(
            time=self.helix.time,
            pressure=pressure_new,
            rho=float(trial["rho_new"]),
            alpha=float(trial["alpha_new"]),
            dw=float(trial["dw"]),
            dv=float(trial["dv"]),
            dkappa=float(trial["dkappa"]),
            Ft=float(trial["Ft"]),
            Tt=float(trial["Tt"]),
            residual=float(trial["residual"]),
            Ftube=float(trial["Ftube"]),
            Mtube=float(trial["Mtube"]),
            Ttube=float(trial["Ttube"]),
            Fnylon=float(trial["Fny"]),
            Mnylon=float(trial["Mny"]),
            Tnylon=float(trial["Tny"]),
            axial_stretch=self.axial_stretch,
            Rin=float(self.R_edges[0]),
            Rout=float(self.R_edges[-1]),
        )
        self.history.append(out)
        return out

    def _suspended_residuals(self, trial: Dict[str, object], load_N: float) -> np.ndarray:
        rho = float(trial["rho_new"])
        alpha = float(trial["alpha_new"])
        sin_a = np.sin(alpha)
        cos_a = np.cos(alpha)
        force_scale = max(abs(load_N), 0.05)
        moment_scale = max(abs(load_N * rho), 0.05)
        return np.array(
            [
                (float(trial["Ftube"]) + float(trial["Fny"]) - load_N * sin_a) / force_scale,
                (float(trial["Mtube"]) + float(trial["Mny"]) + load_N * rho * sin_a) / moment_scale,
                (float(trial["Ttube"]) + float(trial["Tny"]) - load_N * rho * cos_a) / moment_scale,
            ],
            dtype=float,
        )

    def _find_suspended_state(self, dP: float, dt: float, load_N: float) -> Dict[str, object]:
        rho0 = float(self.helix.rho)
        alpha0 = float(self.helix.alpha)
        lower = np.array([-0.25, max(1.001 * self.R_edges[-1], 1e-6), np.deg2rad(0.5)], dtype=float)
        upper = np.array([0.25, max(3.0 * rho0, 2.0 * self.R_edges[-1]), np.deg2rad(85.0)], dtype=float)

        guesses = [
            np.array([0.0, rho0, alpha0], dtype=float),
            np.array([0.0, rho0, np.clip(0.95 * alpha0, lower[2], upper[2])], dtype=float),
            np.array([0.0, rho0, np.clip(1.05 * alpha0, lower[2], upper[2])], dtype=float),
            np.array([0.02, rho0, alpha0], dtype=float),
            np.array([-0.02, rho0, alpha0], dtype=float),
        ]

        best = None
        best_cost = np.inf
        best_trial = None

        def residual_from_x(x: np.ndarray) -> np.ndarray:
            try:
                trial = self._trial_state_from_geometry(
                    float(x[0]),
                    dP,
                    dt,
                    float(x[1]),
                    float(x[2]),
                    axial_coupling=1.0,
                )
                r = self._suspended_residuals(trial, load_N)
                if np.all(np.isfinite(r)):
                    return r
            except Exception:
                pass
            return np.array([1e6, 1e6, 1e6], dtype=float)

        for guess in guesses:
            x0 = np.clip(guess, lower, upper)
            res = least_squares(
                residual_from_x,
                x0,
                bounds=(lower, upper),
                x_scale=np.array([0.05, max(rho0, 1.0), 0.1], dtype=float),
                xtol=1e-8,
                ftol=1e-8,
                gtol=1e-8,
                max_nfev=120,
            )
            cost = float(2.0 * res.cost)
            if cost < best_cost:
                try:
                    trial = self._trial_state_from_geometry(
                        float(res.x[0]),
                        dP,
                        dt,
                        float(res.x[1]),
                        float(res.x[2]),
                        axial_coupling=1.0,
                    )
                except Exception:
                    trial = None
                best = res
                best_cost = cost
                best_trial = trial

        normalized_residual = float(np.sqrt(best_cost))
        bound_margin = np.minimum(best.x - lower, upper - best.x) if best is not None else np.zeros(3)
        scaled_span = np.maximum(upper - lower, 1e-12)
        at_bound = bool(np.any(bound_margin <= 1e-5 * scaled_span))
        if (
            best is None
            or best_trial is None
            or not best.success
            or not np.isfinite(normalized_residual)
            or normalized_residual > 1e-4
            or at_bound
        ):
            message = "no finite solution" if best is None else f"residual={normalized_residual:.3e}, x={best.x}"
            raise RuntimeError(f"Could not solve suspended-mass equilibrium: {message}.")
        best_trial["residual"] = normalized_residual
        best_trial["Ft"] = float(load_N)
        best_trial["Tt"] = float(load_N * float(best_trial["rho_new"]) * np.cos(float(best_trial["alpha_new"])))
        return best_trial

    def step_suspended(self, pressure_new: float, dt: float, load_N: float) -> StepResult:
        self._validate_time_step(dt)
        if load_N <= 0.0:
            raise ValueError("load_N must be strictly positive for suspended-mass equilibrium.")
        if not np.isfinite(pressure_new) or pressure_new < 0.0:
            raise ValueError("pressure must be finite and non-negative.")
        dP = pressure_new - self.helix.pressure
        trial = self._find_suspended_state(dP, dt, load_N)

        self._commit_trial(trial)
        rho_new = float(trial["rho_new"])
        alpha_new = float(trial["alpha_new"])
        h_new = rho_new * np.tan(alpha_new)
        self.helix = HelixState(
            rho=rho_new,
            alpha=alpha_new,
            h=h_new,
            pressure=pressure_new,
            time=self.helix.time + dt,
        )

        out = StepResult(
            time=self.helix.time,
            pressure=pressure_new,
            rho=rho_new,
            alpha=alpha_new,
            dw=float(trial["dw"]),
            dv=float(trial["dv"]),
            dkappa=float(trial["dkappa"]),
            Ft=float(trial["Ft"]),
            Tt=float(trial["Tt"]),
            residual=float(trial["residual"]),
            Ftube=float(trial["Ftube"]),
            Mtube=float(trial["Mtube"]),
            Ttube=float(trial["Ttube"]),
            Fnylon=float(trial["Fny"]),
            Mnylon=float(trial["Mny"]),
            Tnylon=float(trial["Tny"]),
            axial_stretch=self.axial_stretch,
            Rin=float(self.R_edges[0]),
            Rout=float(self.R_edges[-1]),
        )
        self.history.append(out)
        return out

    def prestretch_to(self, eps_tk: float, strain_rate_mm_min: float = 20.0) -> None:
        if eps_tk < 0.0 or not np.isfinite(eps_tk):
            raise ValueError("eps_tk must be finite and non-negative.")
        if strain_rate_mm_min <= 0.0:
            raise ValueError("strain_rate_mm_min must be strictly positive.")
        h_end = (1.0 + eps_tk) * self.h0
        if eps_tk == 0.0:
            self.h_blocked = h_end
            return
        L0 = 2.0 * np.pi * self.turns * self.h0
        total_time = 60.0 * eps_tk * L0 / strain_rate_mm_min
        dt = total_time / self.disc.pre_steps if self.disc.pre_steps > 0 else 1.0
        self._building_reference_state = self.prestrain_reference_mode == "elastic_tk_reference"
        try:
            for h in np.linspace(self.h0, h_end, self.disc.pre_steps + 1)[1:]:
                self.step(0.0, dt, h_target=h)
        finally:
            self._building_reference_state = False
        self.h_blocked = h_end

    def run_pressure_history(self, time: np.ndarray, pressure: np.ndarray) -> List[StepResult]:
        if len(time) != len(pressure):
            raise ValueError("time and pressure must have the same length.")
        for k in range(1, len(time)):
            dt = float(time[k] - time[k - 1])
            if dt <= 0.0:
                raise ValueError("time must be strictly increasing.")
            self.step(float(pressure[k]), dt, h_target=self.h_blocked)
        return self.history

    def run_pressure_history_suspended(self, time: np.ndarray, pressure: np.ndarray, load_N: float) -> List[StepResult]:
        if len(time) != len(pressure):
            raise ValueError("time and pressure must have the same length.")
        for k in range(1, len(time)):
            dt = float(time[k] - time[k - 1])
            if dt <= 0.0:
                raise ValueError("time must be strictly increasing.")
            self.step_suspended(float(pressure[k]), dt, load_N)
        return self.history

    def history_arrays(self) -> Dict[str, np.ndarray]:
        h = self.history
        arr = {
            "time": np.array([x.time for x in h], dtype=float),
            "pressure_MPa": np.array([x.pressure for x in h], dtype=float),
            "rho_mm": np.array([x.rho for x in h], dtype=float),
            "alpha_rad": np.array([x.alpha for x in h], dtype=float),
            "alpha_deg": np.rad2deg(np.array([x.alpha for x in h], dtype=float)),
            "dw": np.array([x.dw for x in h], dtype=float),
            "dv_invmm": np.array([x.dv for x in h], dtype=float),
            "dkappa_invmm": np.array([x.dkappa for x in h], dtype=float),
            "force_N": np.array([x.Ft for x in h], dtype=float),
            "force_mN": 1000.0 * np.array([x.Ft for x in h], dtype=float),
            "torque_Nmm": np.array([x.Tt for x in h], dtype=float),
            "torque_microNm": 1000.0 * np.array([x.Tt for x in h], dtype=float),
            "residual": np.array([x.residual for x in h], dtype=float),
            "axial_stretch": np.array([x.axial_stretch for x in h], dtype=float),
            "Rin_mm": np.array([x.Rin for x in h], dtype=float),
            "Rout_mm": np.array([x.Rout for x in h], dtype=float),
        }

        alpha = np.array([x.alpha for x in h], dtype=float)
        rho = np.array([x.rho for x in h], dtype=float)
        h_per_rad = rho * np.tan(alpha)
        axial_length_geometry = 2.0 * np.pi * self.turns * h_per_rad
        axial_stretch = np.asarray(arr["axial_stretch"], dtype=float)
        axial_length = self.geom.initial_length * axial_stretch * np.sin(alpha) / np.sin(self.alpha0)
        centerline_length = self.geom.initial_length * axial_stretch / np.sin(self.alpha0)
        arr.update(
            {
                "h_mm_per_rad": h_per_rad,
                "axial_length_mm": axial_length,
                "axial_length_geometry_mm": axial_length_geometry,
                "centerline_length_mm": centerline_length,
            }
        )
        sin_a = np.sin(alpha)
        cos_a = np.cos(alpha)
        Ftube = np.array([x.Ftube for x in h], dtype=float)
        Fny = np.array([x.Fnylon for x in h], dtype=float)
        Mtube = np.array([x.Mtube for x in h], dtype=float)
        Mny = np.array([x.Mnylon for x in h], dtype=float)
        Ttube = np.array([x.Ttube for x in h], dtype=float)
        Tny = np.array([x.Tnylon for x in h], dtype=float)
        force_tube = np.divide(Ftube, sin_a, out=np.full_like(Ftube, np.nan), where=np.abs(sin_a) > 1e-12)
        force_nylon = np.divide(Fny, sin_a, out=np.full_like(Fny, np.nan), where=np.abs(sin_a) > 1e-12)
        torque_tube = np.divide(
            Mtube + force_tube * rho * sin_a,
            cos_a,
            out=np.full_like(Mtube, np.nan),
            where=np.abs(cos_a) > 1e-12,
        )
        torque_nylon = np.divide(
            Mny + force_nylon * rho * sin_a,
            cos_a,
            out=np.full_like(Mny, np.nan),
            where=np.abs(cos_a) > 1e-12,
        )
        arr.update(
            {
                "Ftube_axis_N": Ftube,
                "Fnylon_axis_N": Fny,
                "Mtube_Nmm": Mtube,
                "Mnylon_Nmm": Mny,
                "Ttube_axis_Nmm": Ttube,
                "Tnylon_axis_Nmm": Tny,
                "force_tube_mN": 1000.0 * force_tube,
                "force_nylon_mN": 1000.0 * force_nylon,
                "torque_tube_microNm": 1000.0 * torque_tube,
                "torque_nylon_microNm": 1000.0 * torque_nylon,
                "torque_axis_tube_microNm": 1000.0 * Ttube,
                "torque_axis_nylon_microNm": 1000.0 * Tny,
            }
        )
        return arr


# ---------------------------------------------------------------------------
# Historiques, runners et conventions de sortie
# ---------------------------------------------------------------------------


def cyclic_pressure_history(
    n_cycles: int = 3,
    Pmax: float = 1.3,
    flow_rate_mL_min: float = 10.0,
    volume_mL: float = 1.50,
    dt: float = 0.15,
    nonlinear: bool = True,
    gamma_load: float = 3.5,
    gamma_unload: float = 2.8,
) -> Tuple[np.ndarray, np.ndarray]:
    if int(n_cycles) < 1 or flow_rate_mL_min <= 0.0 or volume_mL <= 0.0 or dt <= 0.0:
        raise ValueError("n_cycles, flow_rate_mL_min, volume_mL, and dt must be positive.")
    if Pmax < 0.0 or not np.isfinite(Pmax):
        raise ValueError("Pmax must be finite and non-negative.")
    half_period = 60.0 * volume_mL / flow_rate_mL_min
    period = 2.0 * half_period
    total_time = int(n_cycles) * period
    transitions = np.arange(0.0, total_time + 0.5 * half_period, half_period)
    t = _time_grid_with_events(total_time, dt, transitions)
    phase = (t % period) / period
    loading = phase < 0.5
    P = np.zeros_like(t)
    if nonlinear:
        x = phase[loading] / 0.5
        y = (phase[~loading] - 0.5) / 0.5
        P[loading] = Pmax * x**gamma_load
        P[~loading] = Pmax * (1.0 - y) ** gamma_unload
    else:
        x = phase[loading] / 0.5
        y = (phase[~loading] - 0.5) / 0.5
        P[loading] = Pmax * x
        P[~loading] = Pmax * (1.0 - y)
    P[np.isclose(t, total_time, rtol=0.0, atol=1e-12)] = 0.0
    return t, np.clip(P, 0.0, Pmax)


def ramp_hold_pressure_history(
    P_hold: float = 1.3,
    ramp_time: float = 9.0,
    hold_time: float = 300.0,
    dt: float = 0.50,
    nonlinear_ramp: bool = True,
    gamma_ramp: float = 3.5,
    unload: bool = False,
    unload_time: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    if ramp_time <= 0 or dt <= 0:
        raise ValueError("ramp_time and dt must be positive.")
    if hold_time < 0:
        raise ValueError("hold_time must be non-negative.")
    if unload_time is None:
        unload_time = ramp_time
    if P_hold < 0.0 or not np.isfinite(P_hold):
        raise ValueError("P_hold must be finite and non-negative.")
    if unload and unload_time <= 0.0:
        raise ValueError("unload_time must be positive when unloading is enabled.")
    total_time = ramp_time + hold_time + (unload_time if unload else 0.0)
    events = [ramp_time, ramp_time + hold_time, total_time]
    t = _time_grid_with_events(total_time, dt, events)
    P = np.zeros_like(t)
    ramp_mask = t <= ramp_time
    x = np.clip(t[ramp_mask] / ramp_time, 0.0, 1.0)
    P[ramp_mask] = P_hold * (x**gamma_ramp if nonlinear_ramp else x)
    hold_mask = (t > ramp_time) & (t <= ramp_time + hold_time)
    P[hold_mask] = P_hold
    if unload:
        unload_mask = t > ramp_time + hold_time
        y = np.clip((t[unload_mask] - ramp_time - hold_time) / unload_time, 0.0, 1.0)
        P[unload_mask] = P_hold * (1.0 - y)
    else:
        P[t > ramp_time + hold_time] = P_hold
    return t, P


def add_corrected_output_conventions(arr: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    if len(arr.get("torque_microNm", [])) > 0:
        raw = np.asarray(arr["torque_microNm"], dtype=float)
        arr["torque_total_signed_microNm"] = raw.copy()
        arr["torque_total_magnitude_microNm"] = np.abs(raw)
        arr["torque_act_signed_microNm"] = raw - raw[0]
        arr["torque_act_microNm"] = arr["torque_act_signed_microNm"].copy()
        arr["torque_act_magnitude_microNm"] = np.abs(arr["torque_act_signed_microNm"])
    if len(arr.get("force_mN", [])) > 0:
        force = np.asarray(arr["force_mN"], dtype=float)
        arr["force_total_mN"] = force.copy()
        arr["force_act_mN"] = force - force[0]
    for key in (
        "force_tube_mN",
        "force_nylon_mN",
        "torque_tube_microNm",
        "torque_nylon_microNm",
    ):
        if key in arr and len(arr[key]) > 0:
            act_key = key.replace("_mN", "_act_mN").replace("_microNm", "_act_microNm")
            arr[act_key] = np.asarray(arr[key], dtype=float) - float(arr[key][0])
    return arr


def _prepare_actuation_history(
    n_cycles: int,
    Pmax: float,
    dt: float,
    pressure_time: Optional[np.ndarray],
    pressure_MPa: Optional[np.ndarray],
    flow_rate_mL_min: float,
    volume_mL: float,
    nonlinear_pressure: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    if pressure_time is None and pressure_MPa is None:
        return cyclic_pressure_history(n_cycles, Pmax, flow_rate_mL_min, volume_mL, dt, nonlinear_pressure)
    if pressure_time is None or pressure_MPa is None:
        raise ValueError("pressure_time and pressure_MPa must be provided together.")
    t = np.asarray(pressure_time, dtype=float)
    p = np.asarray(pressure_MPa, dtype=float)
    if t.ndim != 1 or p.ndim != 1 or len(t) != len(p):
        raise ValueError("pressure_time and pressure_MPa must be 1D arrays with equal length.")
    if len(t) < 2 or not np.all(np.diff(t) > 0.0):
        raise ValueError("pressure history must contain at least two strictly increasing samples.")
    if not np.all(np.isfinite(p)) or np.any(p < 0.0):
        raise ValueError("pressure history must contain finite, non-negative pressures.")
    return t, p


def _config_with_overrides(config: Optional[object], **overrides):
    cfg = config if config is not None else default_simulation_config()
    clean = {k: v for k, v in overrides.items() if v is not None}
    if not clean:
        return cfg
    if is_dataclass(cfg):
        usable = {k: v for k, v in clean.items() if hasattr(cfg, k)}
        return replace(cfg, **usable) if usable else cfg

    defaults = vars(default_simulation_config())
    current = dict(defaults)
    current.update(getattr(cfg, "__dict__", {}))
    current.update(clean)
    return SimpleNamespace(**current)


def run_blocked_actuation(
    config: Optional[object] = None,
    pressure_time: Optional[np.ndarray] = None,
    pressure_MPa: Optional[np.ndarray] = None,
    **overrides,
) -> Tuple[TCPAMaxwellBlockedModel, Dict[str, np.ndarray]]:
    """Lancer la simulation bloquee corrigee.

    Exemples :
        run_blocked_actuation(eps=0.8)
        run_blocked_actuation(eps=0.8, n_cycles=3, Pmax=1.3)
    """
    cfg = _config_with_overrides(config, **overrides)
    disc = default_discretization(
        n_layers=cfg.n_layers,
        n_phi=cfg.n_phi,
        pre_steps=cfg.pre_steps,
        dw_bracket=(-0.05, 0.05),
    )
    model = TCPAMaxwellBlockedModel(
        mat=cfg.mat,
        geom=cfg.geom,
        disc=disc,
        integration=cfg.integration,
        prestrain_reference_mode=getattr(cfg, "prestrain_reference_mode", "elastic_tk_reference"),
    )
    model.prestretch_to(cfg.eps)
    t_start = model.helix.time
    t_local, pressure = _prepare_actuation_history(
        cfg.n_cycles,
        cfg.Pmax,
        cfg.dt,
        pressure_time,
        pressure_MPa,
        cfg.flow_rate_mL_min,
        cfg.volume_mL,
        cfg.nonlinear_pressure,
    )
    model.step(float(pressure[0]), 0.0, h_target=model.h_blocked)
    i_act0 = len(model.history) - 1
    model.run_pressure_history(t_start + t_local, pressure)
    full = model.history_arrays()
    arr = {key: value[i_act0:].copy() for key, value in full.items()}
    arr["time"] = arr["time"] - t_start
    add_corrected_output_conventions(arr)
    return model, arr


def _equilibrate_suspended_load(model: TCPAMaxwellBlockedModel, load_N: float) -> float:
    """Place the suspended load and relax its Maxwell transients before actuation."""
    initial_time = float(model.helix.time)
    model.step_suspended(0.0, 0.0, float(load_N))

    branch_E = _maxwell_E(model.mat.maxwell)
    branch_eta = _maxwell_eta(model.mat.maxwell)
    active = (branch_E > 0.0) & (branch_eta > 0.0)
    equivalent_settling_time = 0.0
    if np.any(active):
        relaxation_times = np.sort(branch_eta[active] / branch_E[active])
        settling_steps = np.unique(
            np.concatenate((relaxation_times, [5.0 * relaxation_times[-1], 20.0 * relaxation_times[-1]]))
        )
        original_integration = model.integration
        model.integration = "exponential"
        try:
            for settling_dt in settling_steps:
                model.step_suspended(0.0, float(settling_dt), float(load_N))
                equivalent_settling_time += float(settling_dt)
        finally:
            model.integration = original_integration

    # Stabilization defines the initial state; its clock and intermediate
    # records are not part of the pressure experiment.
    model.helix.time = initial_time
    model.history.clear()
    return equivalent_settling_time


def run_suspended_actuation(
    config: Optional[object] = None,
    load_N: float = 1.0,
    pressure_time: Optional[np.ndarray] = None,
    pressure_MPa: Optional[np.ndarray] = None,
    equilibrate_load_before_pressure: bool = True,
    **overrides,
) -> Tuple[TCPAMaxwellBlockedModel, Dict[str, np.ndarray]]:
    """Lancer la simulation en actionnement libre avec une masse suspendue.

    La geometrie libre est obtenue en resolvant, a chaque pas de temps, les
    equilibres de l'article :

        F_tube + F_nylon = F_load sin(beta_h)
        M_tube + M_nylon = -F_load Rh sin(beta_h)
        T_tube + T_nylon = F_load Rh cos(beta_h)
    """
    cfg = _config_with_overrides(config, **overrides)
    disc = default_discretization(
        n_layers=cfg.n_layers,
        n_phi=cfg.n_phi,
        pre_steps=cfg.pre_steps,
        dw_bracket=(-0.05, 0.05),
    )
    model = TCPAMaxwellBlockedModel(
        mat=cfg.mat,
        geom=cfg.geom,
        disc=disc,
        integration=cfg.integration,
        prestrain_reference_mode=getattr(cfg, "prestrain_reference_mode", "elastic_tk_reference"),
    )
    model.prestretch_to(cfg.eps)
    t_start = model.helix.time
    equivalent_settling_time = 0.0
    if equilibrate_load_before_pressure:
        equivalent_settling_time = _equilibrate_suspended_load(model, float(load_N))
    t_local, pressure = _prepare_actuation_history(
        cfg.n_cycles,
        cfg.Pmax,
        cfg.dt,
        pressure_time,
        pressure_MPa,
        cfg.flow_rate_mL_min,
        cfg.volume_mL,
        cfg.nonlinear_pressure,
    )
    model.step_suspended(float(pressure[0]), 0.0, float(load_N))
    i_act0 = len(model.history) - 1
    reference_length = (
        cfg.geom.initial_length
        * model.axial_stretch
        * np.sin(model.helix.alpha)
        / np.sin(model.alpha0)
    )
    model.run_pressure_history_suspended(t_start + t_local, pressure, float(load_N))
    full = model.history_arrays()
    arr = {key: value[i_act0:].copy() for key, value in full.items()}
    arr["time"] = arr["time"] - t_start
    add_corrected_output_conventions(arr)
    arr["load_N"] = np.full_like(arr["time"], float(load_N), dtype=float)
    arr["load_mN"] = 1000.0 * arr["load_N"]
    if len(arr.get("axial_length_mm", [])) > 0:
        axial = np.asarray(arr["axial_length_mm"], dtype=float)
        contraction_from_initial = reference_length - axial
        arr["free_displacement_mm"] = axial - reference_length
        arr["free_contraction_mm"] = contraction_from_initial
        arr["free_actuation_strain"] = contraction_from_initial / max(abs(float(reference_length)), 1e-12)
        arr["free_actuation_percent"] = 100.0 * arr["free_actuation_strain"]
        arr["reference_axial_length_mm"] = np.full_like(axial, float(reference_length))
    arr["suspended_load_equilibrated"] = bool(equilibrate_load_before_pressure)
    arr["suspended_equivalent_settling_time_s"] = float(equivalent_settling_time)
    arr["mode"] = np.array(["suspended_mass"] * len(arr["time"]), dtype=object)
    return model, arr


def run_hold_relaxation(
    eps: float = 0.8,
    P_hold: float = 1.3,
    hold_time: float = 300.0,
    ramp_time: float = 9.0,
    dt: float = 0.5,
    n_layers: int = 4,
    n_phi: int = 24,
    pre_steps: int = 24,
    integration: str = "exponential",
):
    t, p = ramp_hold_pressure_history(P_hold=P_hold, ramp_time=ramp_time, hold_time=hold_time, dt=dt)
    config = default_simulation_config(
        eps=eps,
        Pmax=P_hold,
        dt=dt,
        n_layers=n_layers,
        n_phi=n_phi,
        pre_steps=pre_steps,
        integration=integration,
    )
    model, arr = run_blocked_actuation(config, pressure_time=t, pressure_MPa=p)
    i0 = int(np.searchsorted(arr["time"], ramp_time, side="left"))
    arr["ramp_time"] = float(ramp_time)
    arr["hold_time"] = float(hold_time)
    arr["hold_start_index"] = min(max(i0, 0), len(arr["time"]) - 1)
    arr["force_hold_relax_mN"] = arr["force_total_mN"] - arr["force_total_mN"][arr["hold_start_index"]]
    arr["torque_hold_relax_microNm"] = arr["torque_act_microNm"] - arr["torque_act_microNm"][
        arr["hold_start_index"]
    ]
    return model, arr


# ---------------------------------------------------------------------------
# Graphes
# ---------------------------------------------------------------------------


def _pyplot(show: bool):
    import matplotlib

    if not show:
        matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def plot_response(
    arr: Dict[str, np.ndarray],
    show: bool = True,
    save_path: Optional[str | Path] = None,
    title: str = "Base - modele TCPA corrige",
):
    plt = _pyplot(show)
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(arr["time"], arr["force_total_mN"], lw=1.5, color="#1f77b4")
    axes[0].set_ylabel("Blocked force (mN)")
    axes[1].plot(arr["time"], arr["torque_act_microNm"], lw=1.5, color="#d62728")
    axes[1].set_ylabel("Actuation torque (microN m)")
    axes[2].plot(arr["time"], arr["pressure_MPa"], lw=1.5, color="#2ca02c")
    axes[2].set_ylabel("Pressure (MPa)")
    axes[2].set_xlabel("Time since pressure start (s)")
    axes[0].set_title(title)
    for ax in axes:
        ax.grid(True, alpha=0.3)
    if save_path is not None:
        fig.savefig(save_path, dpi=180)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_decomposition(
    arr: Dict[str, np.ndarray],
    show: bool = True,
    save_path: Optional[str | Path] = None,
    title: str = "Base - decomposition tube / nylon",
    force_mode: str = "total",
    torque_mode: str = "actuation",
):
    if force_mode not in {"total", "actuation"}:
        raise ValueError("force_mode must be 'total' or 'actuation'.")
    if torque_mode not in {"total", "actuation"}:
        raise ValueError("torque_mode must be 'total' or 'actuation'.")
    plt = _pyplot(show)
    fig, axes = plt.subplots(3, 1, figsize=(9.4, 7.4), sharex=True, constrained_layout=True)
    t = arr["time"]
    f_total = "force_total_mN" if force_mode == "total" else "force_act_mN"
    f_tube = "force_tube_mN" if force_mode == "total" else "force_tube_act_mN"
    f_nylon = "force_nylon_mN" if force_mode == "total" else "force_nylon_act_mN"
    tq_total = "torque_total_signed_microNm" if torque_mode == "total" else "torque_act_microNm"
    tq_tube = "torque_tube_microNm" if torque_mode == "total" else "torque_tube_act_microNm"
    tq_nylon = "torque_nylon_microNm" if torque_mode == "total" else "torque_nylon_act_microNm"

    axes[0].plot(t, arr[f_total], lw=1.6, label="total")
    axes[0].plot(t, arr[f_tube], lw=1.2, label="tube")
    axes[0].plot(t, arr[f_nylon], lw=1.2, label="nylon")
    axes[0].set_ylabel("Force (mN)")
    axes[0].set_title(title)
    axes[0].legend(loc="best")
    axes[1].plot(t, arr[tq_total], lw=1.6, label="total")
    axes[1].plot(t, arr[tq_tube], lw=1.2, label="tube")
    axes[1].plot(t, arr[tq_nylon], lw=1.2, label="nylon")
    axes[1].set_ylabel("Torque (microN m)")
    axes[1].legend(loc="best")
    axes[2].plot(t, arr["pressure_MPa"], lw=1.5, color="#2ca02c")
    axes[2].set_ylabel("Pressure (MPa)")
    axes[2].set_xlabel("Time since pressure start (s)")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    if save_path is not None:
        fig.savefig(save_path, dpi=180)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_hysteresis(
    arr: Dict[str, np.ndarray],
    cycle: int = 1,
    show: bool = True,
    save_path: Optional[str | Path] = None,
    period: Optional[float] = None,
):
    if cycle < 1:
        raise ValueError("cycle must be 1-based.")
    if period is None:
        period = 2.0 * 60.0 * 1.50 / 10.0
    time = np.asarray(arr["time"], dtype=float)
    mask = (time >= (cycle - 1) * period) & (time <= cycle * period)
    if mask.sum() < 3:
        raise ValueError("Selected cycle contains too few points.")
    P = arr["pressure_MPa"][mask]
    F = arr["force_total_mN"][mask]
    T = arr["torque_act_microNm"][mask]
    plt = _pyplot(show)
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6), constrained_layout=True)
    axes[0].plot(P, F, lw=1.8)
    axes[0].set_xlabel("Pressure (MPa)")
    axes[0].set_ylabel("Blocked force (mN)")
    axes[0].set_title(f"Pressure-force hysteresis - cycle {cycle}")
    axes[1].plot(P, T, lw=1.8, color="#d62728")
    axes[1].set_xlabel("Pressure (MPa)")
    axes[1].set_ylabel("Actuation torque (microN m)")
    axes[1].set_title(f"Pressure-torque hysteresis - cycle {cycle}")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    if save_path is not None:
        fig.savefig(save_path, dpi=180)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_all(arr: Dict[str, np.ndarray], cycle: int = 1, show: bool = True):
    return {
        "response": plot_response(arr, show=show),
        "decomposition": plot_decomposition(arr, show=show),
        "hysteresis": plot_hysteresis(arr, cycle=cycle, show=show),
    }


# ---------------------------------------------------------------------------
# Validation et resume
# ---------------------------------------------------------------------------


@dataclass
class SanityReport:
    isotropic_rotation_error: float
    rotated_stiffness_symmetry_error: float
    local_stiffness_min_eigenvalue: float


def stiffness_sanity_report(theta: float = 0.7) -> SanityReport:
    E = 10.0
    nu = 0.3
    G = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    Ciso = np.zeros((6, 6), dtype=float)
    Ciso[:3, :3] = lam
    np.fill_diagonal(Ciso[:3, :3], lam + 2.0 * G)
    Ciso[3, 3] = Ciso[4, 4] = Ciso[5, 5] = G
    Crot_iso = rotate_stiffness_bias(Ciso, theta)
    C_local = ti_stiffness_from_paper(31.24, 8.82, 7.24, 0.205, 0.422)
    Crot_local = rotate_stiffness_bias(C_local, theta)
    return SanityReport(
        isotropic_rotation_error=float(np.max(np.abs(Crot_iso - Ciso))),
        rotated_stiffness_symmetry_error=float(np.max(np.abs(Crot_local - Crot_local.T))),
        local_stiffness_min_eigenvalue=float(np.linalg.eigvalsh(C_local).min()),
    )


def quick_validation() -> Dict[str, float]:
    report = stiffness_sanity_report()
    _, arr = run_blocked_actuation(
        eps=0.5,
        n_cycles=1,
        Pmax=1.2,
        dt=1.0,
        n_layers=2,
        n_phi=8,
        pre_steps=3,
        integration="exponential",
    )
    return {
        "isotropic_rotation_error": report.isotropic_rotation_error,
        "rotated_stiffness_symmetry_error": report.rotated_stiffness_symmetry_error,
        "local_stiffness_min_eigenvalue": report.local_stiffness_min_eigenvalue,
        "n_steps": float(len(arr["time"])),
        "max_abs_residual_Nmm": float(np.max(np.abs(arr["residual"]))),
        "force_total_max_mN": float(np.max(arr["force_total_mN"])),
        "torque_act_max_microNm": float(np.max(arr["torque_act_microNm"])),
        "force_decomposition_error_mN": float(
            np.nanmax(
                np.abs(
                    arr["force_total_mN"]
                    - arr["force_tube_mN"]
                    - arr["force_nylon_mN"]
                )
            )
        ),
        "torque_decomposition_error_microNm": float(
            np.nanmax(
                np.abs(
                    arr["torque_total_signed_microNm"]
                    - arr["torque_tube_microNm"]
                    - arr["torque_nylon_microNm"]
                )
            )
        ),
    }


def summary(arr: Dict[str, np.ndarray]) -> Dict[str, float]:
    return {
        "force_min_mN": float(np.nanmin(arr["force_total_mN"])),
        "force_max_mN": float(np.nanmax(arr["force_total_mN"])),
        "torque_act_min_microNm": float(np.nanmin(arr["torque_act_microNm"])),
        "torque_act_max_microNm": float(np.nanmax(arr["torque_act_microNm"])),
        "pressure_max_MPa": float(np.nanmax(arr["pressure_MPa"])),
        "max_abs_residual_Nmm": float(np.nanmax(np.abs(arr["residual"]))),
    }


if __name__ == "__main__":
    _, data = run_blocked_actuation(eps=0.8, n_cycles=3, Pmax=1.3)
    print("Resume Base:")
    for key, value in summary(data).items():
        print(f"  {key}: {value:.6g}")
    plot_all(data, cycle=1, show=True)
