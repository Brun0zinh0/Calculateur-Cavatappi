"""
Base.py - module autonome final pour le modele TCPA corrige.

Ce fichier ne depend ni de v5 ni de v6. Il contient directement :

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
        "E_axial": 31.24,
        "E_radius": 8.82,
        "G12": 7.24,
        "nu12": 0.205,
        "nu23": 0.422,
        "maxwell": default_maxwell_tensile_params() if maxwell is None else maxwell,
        "E_nylon": 3.69e3,
        "G_nylon": 0.79e3,
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
        "pressure_end_force_mode": "none",
        "pressure_end_force_scale": 0.0,
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
        "integration": "paper_incremental",
        "flow_rate_mL_min": 10.0,
        "volume_mL": 1.50,
        "nonlinear_pressure": True,
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
    Fpressure: float


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
        integration: str = "paper_incremental",
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
        self.integration = integration
        self.nylon_axial_prestrain_coupling = float(getattr(mat, "nylon_axial_prestrain_coupling", 1.0))
        self.nylon_axial_actuation_coupling = float(getattr(mat, "nylon_axial_actuation_coupling", 1.0))
        self.pressure_end_force_mode = str(getattr(geom, "pressure_end_force_mode", "none"))
        self.pressure_end_force_scale = float(getattr(geom, "pressure_end_force_scale", 0.0))
        for name, value in (
            ("nylon_axial_prestrain_coupling", self.nylon_axial_prestrain_coupling),
            ("nylon_axial_actuation_coupling", self.nylon_axial_actuation_coupling),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1.")
        if self.pressure_end_force_scale < 0.0:
            raise ValueError("pressure_end_force_scale must be non-negative.")
        if self.pressure_end_force_mode not in {"none", "projected_inner_area", "axial_inner_area"}:
            raise ValueError("pressure_end_force_mode must be 'none', 'projected_inner_area', or 'axial_inner_area'.")

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

        C_local_total = ti_stiffness_from_paper(
            mat.E_axial,
            mat.E_radius,
            mat.G12,
            mat.nu12,
            mat.nu23,
        )
        self.C_total: List[np.ndarray] = []
        self.C0: List[np.ndarray] = []
        self.Ci: List[List[np.ndarray]] = []
        self.vbar: List[Tuple[float, float, float, float]] = []
        Etotal = _maxwell_E_total(mat.maxwell)
        self.n_maxwell = _maxwell_branch_count(mat.maxwell)
        for R in self.R_centers:
            theta_j = np.arctan((R / geom.Rout) * np.tan(self.theta_f))
            Cbar = rotate_stiffness_bias(C_local_total, theta_j)
            self.C_total.append(Cbar)
            self.C0.append((mat.maxwell.E0 / Etotal) * Cbar)
            self.Ci.append([(Ei / Etotal) * Cbar for Ei in _maxwell_E(mat.maxwell)])
            self.vbar.append(effective_poissons(Cbar))

        shape = (disc.n_layers, disc.n_phi, 6)
        self.sigma0 = np.zeros(shape)
        self.sigma_i = np.zeros((self.n_maxwell,) + shape)
        self.sigma_total = np.zeros(shape)
        self.Fnylon = 0.0
        self.Mnylon = 0.0
        self.Tnylon = 0.0
        self.history: List[StepResult] = []

    def _nylon_axial_coupling(self, h_target: float) -> float:
        if abs(h_target - self.h_blocked) > 1e-10:
            return self.nylon_axial_prestrain_coupling
        return self.nylon_axial_actuation_coupling

    def _pressure_end_force(self, pressure: float, alpha: float) -> float:
        if self.pressure_end_force_mode == "none" or self.pressure_end_force_scale == 0.0:
            return 0.0
        thrust = self.pressure_end_force_scale * pressure * np.pi * self.geom.Rin**2
        if self.pressure_end_force_mode == "projected_inner_area":
            return thrust * np.sin(alpha)
        return thrust

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

    def _layer_AB_mu(self, j: int) -> Tuple[float, float, float]:
        C = self.C_total[j]
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
    ) -> Tuple[np.ndarray, float]:
        C = self.C_total[j]
        cu, ku, cdu, kdu = self._u_base(R, A, B, mu, dv, dw)
        coeff = np.zeros((2, 6), dtype=float)
        known = np.array([dw, ku / R, kdu, 0.0, 0.0, dv * R], dtype=float)
        for k in range(2):
            strain_coeff = np.array([0.0, cu[k] / R, cdu[k], 0.0, 0.0, 0.0], dtype=float)
            coeff[k] = C @ strain_coeff
        sig_known = C @ known
        return coeff[:, 2], sig_known[2] + Cvisc_r

    def _viscous_CS_layer_mean(self, dt: float) -> np.ndarray:
        rates = _maxwell_rates(self.mat.maxwell)
        CS = np.zeros((self.disc.n_layers, 6), dtype=float)
        for i in range(self.n_maxwell):
            CS -= dt * rates[i] * self.sigma_i[i].mean(axis=1)
        return CS

    def _viscous_CS_point(self, j: int, k: int, dt: float) -> np.ndarray:
        rates = _maxwell_rates(self.mat.maxwell)
        CS = np.zeros(6, dtype=float)
        for i in range(self.n_maxwell):
            CS -= dt * rates[i] * self.sigma_i[i, j, k]
        return CS

    def _update_maxwell_branches(self, j: int, k: int, de: np.ndarray, dt: float) -> np.ndarray:
        rates = _maxwell_rates(self.mat.maxwell)
        sigma_i_point = np.empty_like(self.sigma_i[:, j, k])
        for ib in range(self.n_maxwell):
            if self.integration in {"paper_explicit", "explicit", "paper_incremental", "paper"}:
                ds = self.Ci[j][ib] @ de - dt * rates[ib] * self.sigma_i[ib, j, k]
                sigma_i_point[ib] = self.sigma_i[ib, j, k] + ds
            elif self.integration in {"exponential", "stable"}:
                if rates[ib] <= 0.0:
                    sigma_i_point[ib] = self.sigma_i[ib, j, k] + self.Ci[j][ib] @ de
                else:
                    tau = 1.0 / rates[ib]
                    a = np.exp(-dt / tau)
                    factor = tau / dt * (1.0 - a) if dt > 0 else 1.0
                    sigma_i_point[ib] = a * self.sigma_i[ib, j, k] + factor * (self.Ci[j][ib] @ de)
            else:
                raise ValueError(
                    "integration must be 'paper_explicit', 'paper_incremental', or 'exponential'."
                )
        return sigma_i_point

    def _solve_radial_constants(self, dw: float, dv: float, dP: float, dt: float) -> np.ndarray:
        n = self.disc.n_layers
        A_mat = np.zeros((2 * n, 2 * n), dtype=float)
        b_vec = np.zeros(2 * n, dtype=float)
        CSmean = self._viscous_CS_layer_mean(dt)

        def cols(j: int) -> Tuple[int, int]:
            return 2 * j, 2 * j + 1

        row = 0
        j = 0
        Aj, Bj, muj = self._layer_AB_mu(j)
        c_sig, k_sig = self._radial_stress_axisym_coeff(j, self.R_edges[0], Aj, Bj, muj, dv, dw, CSmean[j, 2])
        c0, c1 = cols(j)
        A_mat[row, c0 : c1 + 1] = c_sig
        b_vec[row] = -dP - k_sig
        row += 1

        for j in range(n - 1):
            Rb = self.R_edges[j + 1]
            Aj, Bj, muj = self._layer_AB_mu(j)
            Ak, Bk, muk = self._layer_AB_mu(j + 1)
            cu_j, ku_j, _, _ = self._u_base(Rb, Aj, Bj, muj, dv, dw)
            cu_k, ku_k, _, _ = self._u_base(Rb, Ak, Bk, muk, dv, dw)
            j0, j1 = cols(j)
            k0, k1 = cols(j + 1)
            A_mat[row, j0 : j1 + 1] = cu_j
            A_mat[row, k0 : k1 + 1] = -cu_k
            b_vec[row] = ku_k - ku_j
            row += 1

            cs_j, ks_j = self._radial_stress_axisym_coeff(j, Rb, Aj, Bj, muj, dv, dw, CSmean[j, 2])
            cs_k, ks_k = self._radial_stress_axisym_coeff(j + 1, Rb, Ak, Bk, muk, dv, dw, CSmean[j + 1, 2])
            A_mat[row, j0 : j1 + 1] = cs_j
            A_mat[row, k0 : k1 + 1] = -cs_k
            b_vec[row] = ks_k - ks_j
            row += 1

        j = n - 1
        Aj, Bj, muj = self._layer_AB_mu(j)
        c_sig, k_sig = self._radial_stress_axisym_coeff(j, self.R_edges[-1], Aj, Bj, muj, dv, dw, CSmean[j, 2])
        c0, c1 = cols(j)
        A_mat[row, c0 : c1 + 1] = c_sig
        b_vec[row] = -k_sig
        return np.linalg.solve(A_mat, b_vec).reshape(n, 2)

    def _u_du_layer(self, j: int, R: float, coeffs: np.ndarray, dv: float, dw: float) -> Tuple[float, float]:
        A, B, mu = self._layer_AB_mu(j)
        C1, C2 = coeffs[j]
        u = C1 * R**mu + C2 * R ** (-mu) + A * dv * R * R + B * dw * R
        du = C1 * mu * R ** (mu - 1.0) - C2 * mu * R ** (-mu - 1.0) + 2.0 * A * dv * R + B * dw
        return u, du

    def _trial_state(self, dw: float, dP: float, dt: float, h_target: float) -> Dict[str, object]:
        rho_new, alpha_new, dw, dv, dkappa = self._kinematic_increments(dw, h_target)
        coeffs = self._solve_radial_constants(dw, dv, dP, dt)
        sigma0_new = np.empty_like(self.sigma0)
        sigma_i_new = np.empty_like(self.sigma_i)
        sigma_total_new = np.empty_like(self.sigma_total)
        K_old = (np.cos(self.helix.alpha) ** 2) / self.helix.rho

        for j, R in enumerate(self.R_centers):
            C0 = self.C0[j]
            _, v12b, v13b, v14b = self.vbar[j]
            u, du = self._u_du_layer(j, R, coeffs, dv, dw)
            for k, Phi in enumerate(self.phi):
                denom = 1.0 + K_old * R * np.cos(Phi)
                curv = (dkappa * R * np.cos(Phi) + u * K_old * np.cos(Phi)) / denom
                eps_r = du - v12b * curv
                eps_phi = u / R - v13b * curv
                eps_s = dw + curv
                gamma_sphi = dv * R / denom - v14b * curv
                de = np.array([eps_s, eps_phi, eps_r, 0.0, 0.0, gamma_sphi], dtype=float)
                sigma0_new[j, k] = self.sigma0[j, k] + C0 @ de
                sigma_i_new[:, j, k] = self._update_maxwell_branches(j, k, de, dt)
                if self.integration in {"paper_incremental", "paper"}:
                    # Appendix A form: Delta sigma = C^E Delta epsilon + C^S,
                    # with C^S = -dt * sum_i (E_i/eta_i) sigma_i(t).
                    sigma_total_new[j, k] = self.sigma_total[j, k] + self.C_total[j] @ de + self._viscous_CS_point(
                        j, k, dt
                    )
                else:
                    sigma_total_new[j, k] = sigma0_new[j, k] + sigma_i_new[:, j, k].sum(axis=0)

        Ftube = 0.0
        Mtube = 0.0
        Ttube = 0.0
        for j, R in enumerate(self.R_centers):
            weight_R = R * self.dR[j] * self.dphi
            for k, Phi in enumerate(self.phi):
                sig = sigma_total_new[j, k]
                Ftube += sig[0] * weight_R
                Mtube += sig[0] * (R * np.cos(Phi)) * weight_R
                Ttube += sig[5] * R * weight_R

        rn = self.geom.r_nylon
        axial_coupling = self._nylon_axial_coupling(h_target)
        Fny = self.Fnylon + axial_coupling * np.pi * self.mat.E_nylon * rn**2 * dw
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
            Fpressure = np.nan
        else:
            Ft_structural = Ares / sin_a
            Tt = (Bres + Ft_structural * rho_new * sin_a) / cos_a
            residual = (Tt * sin_a + Ft_structural * rho_new * cos_a) - Cres
            Fpressure = self._pressure_end_force(self.helix.pressure + dP, alpha_new)
            Ft = Ft_structural + Fpressure

        return {
            "rho_new": rho_new,
            "alpha_new": alpha_new,
            "dw": dw,
            "dv": dv,
            "dkappa": dkappa,
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
            "Fpressure": Fpressure,
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
        if rho_new <= self.geom.Rout:
            raise ValueError("Inadmissible geometry: helix radius <= tube outer radius.")
        if not 0.0 < alpha_new < 0.5 * np.pi:
            raise ValueError("Inadmissible geometry: helix angle outside (0, 90 deg).")

        rho_old, alpha_old = self.helix.rho, self.helix.alpha
        dv = np.sin(2.0 * alpha_new) / (2.0 * rho_new) - np.sin(2.0 * alpha_old) / (2.0 * rho_old)
        dkappa = (np.cos(alpha_new) ** 2) / rho_new - (np.cos(alpha_old) ** 2) / rho_old

        coeffs = self._solve_radial_constants(dw, dv, dP, dt)
        sigma0_new = np.empty_like(self.sigma0)
        sigma_i_new = np.empty_like(self.sigma_i)
        sigma_total_new = np.empty_like(self.sigma_total)
        K_old = (np.cos(self.helix.alpha) ** 2) / self.helix.rho

        for j, R in enumerate(self.R_centers):
            C0 = self.C0[j]
            _, v12b, v13b, v14b = self.vbar[j]
            u, du = self._u_du_layer(j, R, coeffs, dv, dw)
            for k, Phi in enumerate(self.phi):
                denom = 1.0 + K_old * R * np.cos(Phi)
                curv = (dkappa * R * np.cos(Phi) + u * K_old * np.cos(Phi)) / denom
                eps_r = du - v12b * curv
                eps_phi = u / R - v13b * curv
                eps_s = dw + curv
                gamma_sphi = dv * R / denom - v14b * curv
                de = np.array([eps_s, eps_phi, eps_r, 0.0, 0.0, gamma_sphi], dtype=float)
                sigma0_new[j, k] = self.sigma0[j, k] + C0 @ de
                sigma_i_new[:, j, k] = self._update_maxwell_branches(j, k, de, dt)
                if self.integration in {"paper_incremental", "paper"}:
                    sigma_total_new[j, k] = self.sigma_total[j, k] + self.C_total[j] @ de + self._viscous_CS_point(
                        j, k, dt
                    )
                else:
                    sigma_total_new[j, k] = sigma0_new[j, k] + sigma_i_new[:, j, k].sum(axis=0)

        Ftube = 0.0
        Mtube = 0.0
        Ttube = 0.0
        for j, R in enumerate(self.R_centers):
            weight_R = R * self.dR[j] * self.dphi
            for k, Phi in enumerate(self.phi):
                sig = sigma_total_new[j, k]
                Ftube += sig[0] * weight_R
                Mtube += sig[0] * (R * np.cos(Phi)) * weight_R
                Ttube += sig[5] * R * weight_R

        rn = self.geom.r_nylon
        Fny = self.Fnylon + axial_coupling * np.pi * self.mat.E_nylon * rn**2 * dw
        Mny = self.Mnylon + 0.25 * np.pi * self.mat.E_nylon * rn**4 * dkappa
        Tny = self.Tnylon + 0.5 * np.pi * self.mat.G_nylon * rn**4 * dv
        Fpressure = self._pressure_end_force(self.helix.pressure + dP, alpha_new)

        return {
            "rho_new": rho_new,
            "alpha_new": alpha_new,
            "dw": dw,
            "dv": dv,
            "dkappa": dkappa,
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
            "Fpressure": Fpressure,
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
        for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if np.isfinite(fa) and np.isfinite(fb) and fa * fb <= 0:
                return float(brentq(lambda z: f(z), a, b, xtol=1e-9, rtol=1e-8, maxiter=80))

        def obj(x: float) -> float:
            y = f(x)
            return 1e100 if not np.isfinite(y) else y * y

        res = minimize_scalar(obj, bounds=(lo, hi), method="bounded", options={"xatol": 1e-8})
        if not res.success:
            raise RuntimeError("Could not solve for dw.")
        return float(res.x)

    def step(self, pressure_new: float, dt: float, h_target: Optional[float] = None) -> StepResult:
        if h_target is None:
            h_target = self.h_blocked
        dP = pressure_new - self.helix.pressure
        dw = self._find_dw(dP, dt, h_target)
        trial = self._trial_state(dw, dP, dt, h_target)

        self.sigma0 = trial["sigma0_new"]
        self.sigma_i = trial["sigma_i_new"]
        self.sigma_total = trial["sigma_total_new"]
        self.Fnylon = float(trial["Fny"])
        self.Mnylon = float(trial["Mny"])
        self.Tnylon = float(trial["Tny"])
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
            Fpressure=float(trial["Fpressure"]),
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
        lower = np.array([-0.25, max(1.001 * self.geom.Rout, 1e-6), np.deg2rad(0.5)], dtype=float)
        upper = np.array([0.25, max(3.0 * rho0, 2.0 * self.geom.Rout), np.deg2rad(85.0)], dtype=float)

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

        if best is None or best_trial is None or not np.isfinite(best_cost):
            raise RuntimeError("Could not solve suspended-mass equilibrium.")
        best_trial["residual"] = float(np.sqrt(best_cost))
        best_trial["Ft"] = float(load_N)
        best_trial["Tt"] = float(load_N * float(best_trial["rho_new"]) * np.cos(float(best_trial["alpha_new"])))
        return best_trial

    def step_suspended(self, pressure_new: float, dt: float, load_N: float) -> StepResult:
        if load_N < 0.0:
            raise ValueError("load_N must be non-negative.")
        dP = pressure_new - self.helix.pressure
        trial = self._find_suspended_state(dP, dt, load_N)

        self.sigma0 = trial["sigma0_new"]
        self.sigma_i = trial["sigma_i_new"]
        self.sigma_total = trial["sigma_total_new"]
        self.Fnylon = float(trial["Fny"])
        self.Mnylon = float(trial["Mny"])
        self.Tnylon = float(trial["Tny"])
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
            Fpressure=float(trial["Fpressure"]),
        )
        self.history.append(out)
        return out

    def prestretch_to(self, eps_tk: float, strain_rate_mm_min: float = 20.0) -> None:
        h_end = (1.0 + eps_tk) * self.h0
        L0 = 2.0 * np.pi * self.turns * self.h0
        total_time = 60.0 * eps_tk * L0 / strain_rate_mm_min
        dt = total_time / self.disc.pre_steps if self.disc.pre_steps > 0 else 1.0
        for h in np.linspace(self.h0, h_end, self.disc.pre_steps + 1)[1:]:
            self.step(0.0, dt, h_target=h)
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
        }

        alpha = np.array([x.alpha for x in h], dtype=float)
        rho = np.array([x.rho for x in h], dtype=float)
        h_per_rad = rho * np.tan(alpha)
        axial_length = 2.0 * np.pi * self.turns * h_per_rad
        centerline_length = 2.0 * np.pi * self.turns * rho / np.cos(alpha)
        arr.update(
            {
                "h_mm_per_rad": h_per_rad,
                "axial_length_mm": axial_length,
                "centerline_length_mm": centerline_length,
            }
        )
        sin_a = np.sin(alpha)
        cos_a = np.cos(alpha)
        Ftube = np.array([x.Ftube for x in h], dtype=float)
        Fny = np.array([x.Fnylon for x in h], dtype=float)
        Fpressure = np.array([x.Fpressure for x in h], dtype=float)
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
                "force_pressure_end_mN": 1000.0 * Fpressure,
                "force_structural_mN": 1000.0 * (np.asarray([x.Ft for x in h], dtype=float) - Fpressure),
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
    half_period = 60.0 * volume_mL / flow_rate_mL_min
    period = 2.0 * half_period
    t = np.arange(0.0, n_cycles * period + dt, dt)
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
    return t, P


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
    total_time = ramp_time + hold_time + (unload_time if unload else 0.0)
    t = np.arange(0.0, total_time + dt, dt)
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
        "force_pressure_end_mN",
        "force_structural_mN",
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
    model = TCPAMaxwellBlockedModel(mat=cfg.mat, geom=cfg.geom, disc=disc, integration=cfg.integration)
    model.prestretch_to(cfg.eps)
    i_act0 = len(model.history)
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
    model.run_pressure_history(t_start + t_local, pressure)
    full = model.history_arrays()
    arr = {key: value[i_act0:].copy() for key, value in full.items()}
    arr["time"] = arr["time"] - t_start
    add_corrected_output_conventions(arr)
    return model, arr


def run_suspended_actuation(
    config: Optional[object] = None,
    load_N: float = 1.0,
    pressure_time: Optional[np.ndarray] = None,
    pressure_MPa: Optional[np.ndarray] = None,
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
    model = TCPAMaxwellBlockedModel(mat=cfg.mat, geom=cfg.geom, disc=disc, integration=cfg.integration)
    model.prestretch_to(cfg.eps)
    i_act0 = len(model.history)
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
    model.run_pressure_history_suspended(t_start + t_local, pressure, float(load_N))
    full = model.history_arrays()
    arr = {key: value[i_act0:].copy() for key, value in full.items()}
    arr["time"] = arr["time"] - t_start
    add_corrected_output_conventions(arr)
    arr["load_N"] = np.full_like(arr["time"], float(load_N), dtype=float)
    arr["load_mN"] = 1000.0 * arr["load_N"]
    if len(arr.get("axial_length_mm", [])) > 0:
        axial = np.asarray(arr["axial_length_mm"], dtype=float)
        contraction = axial[0] - axial
        arr["free_displacement_mm"] = axial - axial[0]
        arr["free_contraction_mm"] = contraction
        arr["free_actuation_strain"] = contraction / max(float(cfg.geom.initial_length), 1e-12)
        arr["free_actuation_percent"] = 100.0 * arr["free_actuation_strain"]
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
    integration: str = "paper_incremental",
):
    disc = default_discretization(n_layers=n_layers, n_phi=n_phi, pre_steps=pre_steps, dw_bracket=(-0.05, 0.05))
    model = TCPAMaxwellBlockedModel(disc=disc, integration=integration)
    model.prestretch_to(eps)
    i_act0 = len(model.history)
    t_start = model.helix.time
    t, p = ramp_hold_pressure_history(P_hold=P_hold, ramp_time=ramp_time, hold_time=hold_time, dt=dt)
    model.run_pressure_history(t_start + t, p)
    full = model.history_arrays()
    arr = {key: value[i_act0:].copy() for key, value in full.items()}
    arr["time"] = arr["time"] - t_start
    add_corrected_output_conventions(arr)
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
    f_pressure = "force_pressure_end_mN" if force_mode == "total" else "force_pressure_end_act_mN"
    tq_total = "torque_total_signed_microNm" if torque_mode == "total" else "torque_act_microNm"
    tq_tube = "torque_tube_microNm" if torque_mode == "total" else "torque_tube_act_microNm"
    tq_nylon = "torque_nylon_microNm" if torque_mode == "total" else "torque_nylon_act_microNm"

    axes[0].plot(t, arr[f_total], lw=1.6, label="total")
    axes[0].plot(t, arr[f_tube], lw=1.2, label="tube")
    axes[0].plot(t, arr[f_nylon], lw=1.2, label="nylon")
    if f_pressure in arr and np.nanmax(np.abs(arr[f_pressure])) > 1e-9:
        axes[0].plot(t, arr[f_pressure], lw=1.2, label="pressure end")
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
        integration="paper_incremental",
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
                    - arr.get("force_pressure_end_mN", 0.0)
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
