from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from livrable_parametres import PSI_TO_MPA, SettingValue, VISUAL_STATE_LABELS, derived_geometry


def format_seconds(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    if seconds < 60.0:
        return f"{seconds:.1f} s"
    minutes, rem = divmod(seconds, 60.0)
    if minutes < 60.0:
        return f"{int(minutes)} min {rem:04.1f} s"
    hours, minutes = divmod(minutes, 60.0)
    return f"{int(hours)} h {int(minutes):02d} min"


def estimate_compute_seconds(cost_index: int) -> float:
    return min(3600.0, max(0.3, 0.35 + 1.55e-3 * max(1, int(cost_index))))


def draw_dimension_3d(ax, start, end, label: str, text_offset=(0.0, 0.0, 0.0), color="#222222") -> None:
    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)
    vec = end - start
    length = float(np.linalg.norm(vec))
    if length <= 1e-12:
        return

    unit = vec / length
    head = min(0.18 * length, max(0.25, 0.05 * length))
    ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=color, lw=1.2)
    ax.quiver(*(start + unit * head), *(-unit * head), color=color, linewidth=1.1, arrow_length_ratio=0.65)
    ax.quiver(*(end - unit * head), *(unit * head), color=color, linewidth=1.1, arrow_length_ratio=0.65)
    midpoint = 0.5 * (start + end) + np.asarray(text_offset, dtype=float)
    ax.text(midpoint[0], midpoint[1], midpoint[2], label, color=color, fontsize=9, ha="center", va="center")


def make_cavatappi_figure(settings: dict[str, SettingValue], state: str):
    geom = derived_geometry(settings)
    rho = float(settings["rho0_mm"])
    rout = float(settings["rout_mm"])
    length = float(settings["initial_length_mm"])
    pitch = geom["pitch0_mm"]
    view_elev = float(settings.get("view_elev_deg", 22.0))
    view_azim = float(settings.get("view_azim_deg", -58.0))
    if state == "prestrained":
        length = geom["prestrained_length_mm"]
        pitch = geom["prestrained_pitch_mm"]

    turns = max(length / pitch, 0.1)
    n = max(250, int(80 * turns))
    theta = np.linspace(0.0, 2.0 * np.pi * turns, n)
    z = pitch * theta / (2.0 * np.pi)
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)

    fig = plt.figure(figsize=(7.6, 5.6))
    ax = fig.add_subplot(111, projection="3d")
    color = "#1f77b4" if state == "fabricated" else "#d62728"
    ax.plot(x, y, z, color=color, lw=4.0, solid_capstyle="round")
    ax.plot([0, 0], [0, 0], [0, length], color="0.45", lw=1.2, ls="--")

    for zi in (0.0, length):
        circle_theta = np.linspace(0.0, 2.0 * np.pi, 160)
        ax.plot(
            rho * np.cos(circle_theta),
            rho * np.sin(circle_theta),
            np.full_like(circle_theta, zi),
            color="0.78",
            lw=0.8,
        )

    outer_radius = rho + rout
    dim_pad = max(0.7, 1.8 * rout)
    radius_limit = outer_radius + 2.2 * dim_pad
    z_pad = max(1.0, 0.08 * length)

    x_len = outer_radius + 0.95 * dim_pad
    y_len = -outer_radius - 0.85 * dim_pad
    draw_dimension_3d(
        ax,
        (x_len, y_len, 0.0),
        (x_len, y_len, length),
        f"L = {length:.1f} mm",
        text_offset=(0.45 * dim_pad, 0.0, 0.0),
    )
    ax.plot([0.0, x_len], [0.0, y_len], [0.0, 0.0], color="0.55", lw=0.8)
    ax.plot([0.0, x_len], [0.0, y_len], [length, length], color="0.55", lw=0.8)

    y_diam = -outer_radius - 1.55 * dim_pad
    z_diam = -0.25 * z_pad
    draw_dimension_3d(
        ax,
        (-outer_radius, y_diam, z_diam),
        (outer_radius, y_diam, z_diam),
        f"diamètre ext. spire = {2.0 * outer_radius:.2f} mm",
        text_offset=(0.0, -0.28 * dim_pad, 0.0),
    )
    ax.plot([-outer_radius, -outer_radius], [0.0, y_diam], [z_diam, z_diam], color="0.55", lw=0.8)
    ax.plot([outer_radius, outer_radius], [0.0, y_diam], [z_diam, z_diam], color="0.55", lw=0.8)

    visible_pitch = min(pitch, length)
    x_pitch = -outer_radius - 0.95 * dim_pad
    y_pitch = outer_radius + 0.75 * dim_pad
    draw_dimension_3d(
        ax,
        (x_pitch, y_pitch, 0.0),
        (x_pitch, y_pitch, visible_pitch),
        f"pas = {pitch:.2f} mm",
        text_offset=(-0.45 * dim_pad, 0.0, 0.0),
    )
    ax.plot([0.0, x_pitch], [0.0, y_pitch], [0.0, 0.0], color="0.55", lw=0.8)
    ax.plot([0.0, x_pitch], [0.0, y_pitch], [visible_pitch, visible_pitch], color="0.55", lw=0.8)

    ax.set_xlim(-radius_limit, radius_limit)
    ax.set_ylim(-radius_limit, radius_limit)
    ax.set_zlim(-z_pad, max(length + z_pad, 1e-6))
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("axe (mm)")
    ax.set_title(f"Géométrie du Cavatappi - {VISUAL_STATE_LABELS.get(state, state)}")
    ax.view_init(elev=view_elev, azim=view_azim)
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_alpha(0.0)
        axis._axinfo["grid"]["linewidth"] = 0.0
    try:
        ax.set_box_aspect((1, 1, max((length + 2.0 * z_pad) / (2.0 * radius_limit), 0.5)))
    except AttributeError:
        pass
    return fig


def make_cross_section_figure(settings: dict[str, SettingValue]):
    rout = float(settings["rout_mm"])
    rin = float(settings["rin_mm"])
    nylon_radius = 0.5 * float(settings["nylon_diameter_mm"])

    fig, ax = plt.subplots(figsize=(4.9, 4.9), constrained_layout=True)
    ax.add_patch(Circle((0.0, 0.0), rout, facecolor="#d7e8ff", edgecolor="#1f77b4", lw=2.0, label="tube PVC"))
    ax.add_patch(Circle((0.0, 0.0), rin, facecolor="white", edgecolor="#1f77b4", lw=1.5, label="alesage interne"))
    ax.add_patch(Circle((0.0, 0.0), nylon_radius, facecolor="#ffd8a8", edgecolor="#ff7f0e", lw=1.8, label="nylon"))
    limit = 1.55 * rout
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("mm")
    ax.set_ylabel("mm")
    ax.set_title("Section du tube")
    ax.grid(False)
    arrow = {"arrowstyle": "<->", "color": "#222222", "lw": 1.1, "shrinkA": 0.0, "shrinkB": 0.0}
    ax.annotate("", xy=(-rout, -1.22 * rout), xytext=(rout, -1.22 * rout), arrowprops=arrow)
    ax.text(0.0, -1.34 * rout, f"2Rout = {2.0 * rout:.2f} mm", ha="center", va="top", fontsize=9)
    ax.annotate("", xy=(-rin, 1.17 * rout), xytext=(rin, 1.17 * rout), arrowprops=arrow)
    ax.text(0.0, 1.28 * rout, f"2Rin = {2.0 * rin:.2f} mm", ha="center", va="bottom", fontsize=9)
    if nylon_radius > 0.0:
        ax.annotate(
            "",
            xy=(-nylon_radius, 0.0),
            xytext=(nylon_radius, 0.0),
            arrowprops={"arrowstyle": "<->", "color": "#8a4b08", "lw": 1.0, "shrinkA": 0.0, "shrinkB": 0.0},
        )
        ax.text(0.0, 0.12 * rout, f"nylon = {2.0 * nylon_radius:.2f} mm", ha="center", fontsize=8.5)
    ax.legend(loc="upper right")
    return fig


def plot_time_response_fr(data: dict[str, np.ndarray]):
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(data["time"], data["force_total_mN"], lw=1.5, color="#1f77b4")
    axes[0].set_ylabel("Force bloquée (mN)")
    axes[1].plot(data["time"], data["torque_act_microNm"], lw=1.5, color="#d62728")
    axes[1].set_ylabel("Couple d'actionnement (microN m)")
    axes[2].plot(data["time"], data["pressure_MPa"], lw=1.5, color="#2ca02c")
    axes[2].set_ylabel("Pression (MPa)")
    axes[2].set_xlabel("Temps depuis le debut de pression (s)")
    axes[2].set_xlabel("Temps depuis le début de pression (s)")
    axes[0].set_title("Réponse du modèle Cavatappi")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    return fig


def _add_direction_arrows(ax, x, y, color, n_arrows: int = 6) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if x.size < 3:
        return

    segment_step = max(1, x.size // 80)
    arrow_positions = np.linspace(0, x.size - 2, n_arrows + 2, dtype=int)[1:-1]
    for i in np.unique(arrow_positions):
        j = min(x.size - 1, i + segment_step)
        while j < x.size - 1 and np.hypot(x[j] - x[i], y[j] - y[i]) < 1e-12:
            j = min(x.size - 1, j + segment_step)
            if j == x.size - 1:
                break
        if j <= i or np.hypot(x[j] - x[i], y[j] - y[i]) < 1e-12:
            continue
        ax.annotate(
            "",
            xy=(x[j], y[j]),
            xytext=(x[i], y[i]),
            arrowprops={
                "arrowstyle": "->",
                "color": color,
                "lw": 1.4,
                "shrinkA": 0.0,
                "shrinkB": 0.0,
                "mutation_scale": 13,
            },
        )


def plot_hysteresis_with_arrows(data, cycle: int = 1, period: float | None = None, show: bool = False):
    if cycle < 1:
        raise ValueError("cycle must be 1-based.")
    if period is None:
        period = 2.0 * 60.0 * 1.50 / 10.0

    time = np.asarray(data["time"], dtype=float)
    mask = (time >= (cycle - 1) * period) & (time <= cycle * period)
    if mask.sum() < 3:
        raise ValueError("Le cycle selectionne contient trop peu de points.")

    pressure_psi = np.asarray(data["pressure_MPa"], dtype=float)[mask] / PSI_TO_MPA
    force_mN = np.asarray(data["force_total_mN"], dtype=float)[mask]
    torque_microNm = np.asarray(data["torque_act_microNm"], dtype=float)[mask]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6), constrained_layout=True)
    force_color = "#1f77b4"
    torque_color = "#d62728"

    axes[0].plot(pressure_psi, force_mN, lw=1.8, color=force_color)
    _add_direction_arrows(axes[0], pressure_psi, force_mN, force_color)
    axes[0].scatter(pressure_psi[0], force_mN[0], s=28, color=force_color, zorder=3, label="début")
    axes[0].scatter(pressure_psi[-1], force_mN[-1], s=28, facecolor="white", edgecolor=force_color, zorder=3, label="fin")
    axes[0].set_xlabel("Pression (psi)")
    axes[0].set_ylabel("Force bloquée (mN)")
    axes[0].set_title(f"Hystérèse pression-force - cycle {cycle}")
    axes[0].legend()

    axes[1].plot(pressure_psi, torque_microNm, lw=1.8, color=torque_color)
    _add_direction_arrows(axes[1], pressure_psi, torque_microNm, torque_color)
    axes[1].scatter(pressure_psi[0], torque_microNm[0], s=28, color=torque_color, zorder=3, label="début")
    axes[1].scatter(
        pressure_psi[-1], torque_microNm[-1], s=28, facecolor="white", edgecolor=torque_color, zorder=3, label="fin"
    )
    axes[1].set_xlabel("Pression (psi)")
    axes[1].set_ylabel("Couple d'actionnement (microN m)")
    axes[1].set_title(f"Hystérèse pression-couple - cycle {cycle}")
    axes[1].legend()

    for ax in axes:
        ax.grid(True, alpha=0.3)
    if show:
        plt.show()
    return fig


def plot_hysteresis_overlay(cases: list[dict], cycles: list[int], show: bool = False):
    """Trace plusieurs cycles ou plusieurs variantes de paramètres sur le même graphe."""
    if not cases:
        raise ValueError("Aucune donnée d'hystérèse à afficher.")
    if not cycles:
        raise ValueError("Aucun cycle d'hystérèse sélectionné.")

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), constrained_layout=True)
    color_count = max(1, len(cases) * len(cycles))
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, min(color_count, 10)))
    plotted = 0
    show_arrows = color_count <= 6

    for case in cases:
        data = case["data"]
        label = str(case.get("label", "cas"))
        period = float(case["period"])
        time = np.asarray(data["time"], dtype=float)

        for cycle in cycles:
            if cycle < 1:
                continue
            mask = (time >= (cycle - 1) * period) & (time <= cycle * period)
            if mask.sum() < 3:
                continue

            pressure_psi = np.asarray(data["pressure_MPa"], dtype=float)[mask] / PSI_TO_MPA
            force_mN = np.asarray(data["force_total_mN"], dtype=float)[mask]
            torque_microNm = np.asarray(data["torque_act_microNm"], dtype=float)[mask]
            color = colors[plotted % len(colors)]
            if len(cases) == 1:
                curve_label = f"cycle {cycle}"
            elif len(cycles) == 1:
                curve_label = label
            else:
                curve_label = f"{label} - cycle {cycle}"

            axes[0].plot(pressure_psi, force_mN, lw=1.7, color=color, label=curve_label)
            axes[1].plot(pressure_psi, torque_microNm, lw=1.7, color=color, label=curve_label)
            if show_arrows:
                _add_direction_arrows(axes[0], pressure_psi, force_mN, color, n_arrows=3)
                _add_direction_arrows(axes[1], pressure_psi, torque_microNm, color, n_arrows=3)
            plotted += 1

    if plotted == 0:
        raise ValueError("Les cycles sélectionnés ne contiennent pas assez de points.")

    axes[0].set_xlabel("Pression (psi)")
    axes[0].set_ylabel("Force bloquée (mN)")
    axes[0].set_title("Hystérèse pression-force")
    axes[1].set_xlabel("Pression (psi)")
    axes[1].set_ylabel("Couple d'actionnement (microN m)")
    axes[1].set_title("Hystérèse pression-couple")

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)
    if show:
        plt.show()
    return fig


def plot_relaxation_response(data: dict[str, np.ndarray]):
    time = np.asarray(data["time"], dtype=float)
    ramp_time = float(data["ramp_time"])
    hold_start_index = int(data["hold_start_index"])
    t_hold = time[hold_start_index:] - ramp_time

    fig, axes = plt.subplots(4, 1, figsize=(9.2, 8.8), sharex=False, constrained_layout=True)
    fig.suptitle("Relaxation en actionnement bloqué à pression constante")

    axes[0].plot(time, data["force_total_mN"], color="#1f77b4", lw=1.6)
    axes[0].axvline(ramp_time, color="0.35", lw=1.0, ls="--")
    axes[0].set_ylabel("Force (mN)")

    axes[1].plot(t_hold, data["force_hold_relax_mN"][hold_start_index:], color="#1f77b4", lw=1.6)
    axes[1].axhline(0.0, color="0.35", lw=0.9)
    axes[1].set_ylabel("Variation de force (mN)")

    axes[2].plot(t_hold, data["torque_hold_relax_microNm"][hold_start_index:], color="#d62728", lw=1.6)
    axes[2].axhline(0.0, color="0.35", lw=0.9)
    axes[2].set_ylabel("Variation de couple (microN m)")

    axes[3].plot(time, data["pressure_MPa"], color="#2ca02c", lw=1.6)
    axes[3].axvline(ramp_time, color="0.35", lw=1.0, ls="--")
    axes[3].set_ylabel("Pression (MPa)")
    axes[3].set_xlabel("Temps depuis le début de pression (s)")

    for ax in axes:
        ax.grid(True, alpha=0.28)
    axes[1].set_xlabel("Temps de maintien à pression constante (s)")
    axes[2].set_xlabel("Temps de maintien à pression constante (s)")
    return fig


def plot_suspended_response(data: dict[str, np.ndarray], show_geometry: bool = False):
    time = np.asarray(data["time"], dtype=float)
    n_axes = 4 if show_geometry else 3
    fig_height = 9.2 if show_geometry else 7.0
    fig, axes = plt.subplots(n_axes, 1, figsize=(9.2, fig_height), sharex=True, constrained_layout=True)
    fig.suptitle("Actionnement libre avec masse suspendue")
    axes = np.asarray(axes).ravel()
    hold_start_time = None
    if bool(data.get("suspended_hold_pressure", False)) and "suspended_hold_start_index" in data:
        hold_index = int(data["suspended_hold_start_index"])
        if 0 <= hold_index < len(time):
            hold_start_time = float(time[hold_index])

    axes[0].plot(time, data["free_actuation_percent"], color="#1f77b4", lw=1.6)
    axes[0].set_ylabel("Actionnement (%)")

    contraction_mm = np.asarray(data["free_contraction_mm"], dtype=float)
    axes[1].plot(time, contraction_mm, color="#ff7f0e", lw=1.7)
    axes[1].axhline(0.0, color="0.35", lw=0.8)
    axes[1].set_ylabel("Contraction (mm)")

    axes[2].plot(time, data["pressure_MPa"], color="#2ca02c", lw=1.6)
    axes[2].set_ylabel("Pression (MPa)")

    axes[2].set_xlabel("Temps depuis le debut de pression (s)")
    if not show_geometry:
        for ax in axes:
            if hold_start_time is not None:
                ax.axvline(hold_start_time, color="0.35", lw=1.0, ls="--")
            ax.grid(True, alpha=0.28)
        return fig

    axes[3].set_ylabel("Rh (mm)")
    ax_angle = axes[3].twinx()
    ax_angle.set_ylabel("beta_h (deg)")
    axes[3].set_xlabel("Temps depuis le début de pression (s)")

    if "rho_mm" in data and "alpha_deg" in data:
        axes[3].plot(time, data["rho_mm"], color="#d62728", lw=1.5, label="Rh")
        ax_angle.plot(time, data["alpha_deg"], color="#17becf", lw=1.2, label="beta_h")
        lines, labels = axes[3].get_legend_handles_labels()
        lines2, labels2 = ax_angle.get_legend_handles_labels()
        axes[3].legend(lines + lines2, labels + labels2, loc="best")
    else:
        axes[3].text(
            0.5,
            0.5,
            "Données Rh / beta_h absentes : relancez le calcul masse suspendue.",
            transform=axes[3].transAxes,
            ha="center",
            va="center",
        )

    for ax in axes:
        if hold_start_time is not None:
            ax.axvline(hold_start_time, color="0.35", lw=1.0, ls="--")
        ax.grid(True, alpha=0.28)
    return fig


def plot_prestrain_study(result: dict[str, np.ndarray]):
    eps_values = result["eps"]
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 7.2), sharex=True, constrained_layout=True)
    fig.suptitle("Étude de la force en fonction de la précontrainte initiale")

    axes[0].plot(eps_values, result["force_initial_mN"], marker="o", lw=1.6, label="force initiale")
    axes[0].plot(eps_values, result["force_max_mN"], marker="o", lw=1.8, label="force maximale")
    axes[0].set_ylabel("Force bloquée (mN)")
    axes[0].legend(loc="best")

    axes[1].plot(eps_values, result["force_gain_mN"], marker="o", color="#d62728", lw=1.8)
    axes[1].set_ylabel("Gain de force max (mN)")
    axes[1].set_xlabel("Précontrainte initiale")

    for ax in axes:
        ax.grid(True, alpha=0.3)
    return fig
