# -*- coding: utf-8 -*-
"""Confrontation modèle / essais en fonction de la vitesse de sollicitation (MPa/s).

Le banc n'est pas asservi : chaque rampe de pression a sa propre pente. On
segmente donc chaque essai (1 O, 3 C, 3 P) en rampes, on mesure la vitesse
de chacune (ajustement linéaire sur 10-90 % de la course) et l'on compare,
rampe par rampe et cycle par cycle, la réponse mesurée et la réponse simulée
sous la MÊME pression mesurée :
  - gain de rampe  ΔF/ΔP (mN/MPa) en montée et en descente ;
  - par cycle (montée puis descente) : aire de boucle brute et corrigée de
    la dérive (définition 7.5.1 du rapport, mN·bar), largeur à mi-course,
    retard force/pression, vitesse moyenne du cycle.
Puis un balayage synthétique du moteur à vitesse imposée (nouveau paramètre
pressure_rate_mpa_s) montre la dépendance intrinsèque du modèle.

Sorties : vitesse/confrontation_vitesse.json, vitesse/tables.md, PNG.
"""
import csv
import glob
import io
import json
import os
import re
import sys
import warnings
from multiprocessing import Pool

import numpy as np

# Chemins surchargeables par variables d'environnement : TCPA_APP (dossier du
# moteur alpha V4), TCPA_ESSAIS (dossier des CSV de la campagne), TCPA_OUT
# (dossier de sortie). Par défaut : l'arbre de travail de l'auteur.
W = r"C:\Users\b.pereiraazevedo\AppData\Local\Temp\tcpa_v3"
APP = os.environ.get("TCPA_APP", os.path.join(W, "alpha V4"))
ESSAIS = os.environ.get(
    "TCPA_ESSAIS",
    r"C:\Users\b.pereiraazevedo\OneDrive - House Of HR NV\Documents\Stage muscle artificièle"
    r"\modèle\modèle chinois\Espace de travail\essais",
)
OUT = os.environ.get("TCPA_OUT", os.path.join(W, "vitesse"))
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, APP)

import Base  # noqa: E402
import parametres  # noqa: E402
from scipy.ndimage import uniform_filter1d  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402
from scipy.stats import spearmanr, theilslopes, wilcoxon  # noqa: E402

PSI_TO_BAR = 0.0689476
E_TOTAL = 37.76
COMMON = dict(rout_mm=0.96, rin_mm=0.56, nylon_diameter_mm=0.85, rho0_mm=2.29, theta_f_deg=42.0)
# (longueur helice mm, angle initial deg, longueur desenroulee mm) — Dataset.xlsx
FICHES = {
    ("J", 0.8): (49.2, 7.788, 20.0), ("J", 1.0): (49.2, 8.109, 30.0), ("J", 1.2): (49.2, 8.456, 38.0),
    ("L", 0.8): (49.66, 7.860, 12.0), ("L", 1.0): (49.66, 7.860, 18.0), ("L", 1.2): (49.66, 7.860, 22.0),
    ("Q", 0.8): (56.82, 8.321, 8.0), ("Q", 1.0): (56.82, 8.321, 8.0), ("Q", 1.2): (56.82, 8.321, 8.0),
    ("R", 0.8): (53.0, 7.769, 20.0), ("R", 1.0): (53.0, 7.769, 20.0), ("R", 1.2): (53.0, 7.769, 20.0),
    ("G", 0.8): (43.6, 8.211, 18.0), ("G", 1.0): (43.6, 8.211, 20.0), ("G", 1.2): (43.6, 8.615, 22.0),
    ("A", 0.8): (55.61, 7.591, 20.0), ("A", 1.0): (55.61, 7.859, 22.0), ("A", 1.2): (55.61, 7.859, 24.0),
    ("E", 0.8): (47.37, 8.146, 16.0), ("E", 1.0): (47.37, 8.146, 18.0), ("E", 1.2): (47.37, 8.511, 20.0),
    ("F", 0.8): (48.0, 8.253, 20.0), ("F", 1.0): (48.0, 8.622, 22.0), ("F", 1.2): (48.0, 9.026, 24.0),
}
SPECTRES = {
    # médiane des 18 essais 10 N (spectre_identifie.json)
    "10N": dict(E0=33.4809, E1=1.6273, eta1=31.86, E2=2.6518, eta2=1014.66, E3=0.0, eta3=1.0),
    # maintiens longs 30 L (ANALYSE_PRELIMINAIRE § F)
    "30L": dict(E0=31.3, E1=2.4, eta1=115.0, E2=4.3, eta2=5270.0, E3=0.0, eta3=1.0),
}
V4_DEMO = dict(engagement_reform_pressure_mpa=0.257, engagement_unload_ratio=0.7, friction_pressure_coulomb_mpa=0.02)
PRESTRETCH_RATE_MM_MIN = 300.0
DT_SIM = 0.25
WORKERS = 8
NAME_RE = re.compile(r"Mus[a-z]+ (\w) (\d+) ([A-Z]) ([\d.]+)\.csv$")


# ---------------------------------------------------------------------------
# Lecture des essais, pression reconstruite (7.3.1 / recalcul_analyse.py)
# ---------------------------------------------------------------------------
def psi(adc):
    return (adc * 5.0 / 1023.0 - 0.5) * 500.0 / 4.0


def load_test(path):
    rows = list(csv.reader(io.StringIO(io.open(path, encoding="utf-8-sig").read())))
    idx = {k: i for i, k in enumerate(rows[0])}
    data = [r for r in rows[1:] if r and r[0]]
    g = lambda k, f=float: np.array([f(r[idx[k]]) for r in data])
    t, F, adc = g("time_s"), g("force_unfiltered_mN"), g("pressure_adc", int)
    adc_rest = int(np.bincount(adc[:15]).argmax())
    P = (psi(adc.astype(float)) - psi(adc_rest)) * PSI_TO_BAR * 0.1  # MPa, référencée au repos
    keep = np.concatenate(([True], np.diff(t) > 0.0))
    return t[keep], np.clip(P[keep], 0.0, None), F[keep]


def list_tests():
    out = []
    for folder in sorted(os.listdir(ESSAIS)):
        d = os.path.join(ESSAIS, folder)
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            m = NAME_RE.match(fn)
            if not m:
                continue
            muscle, minutes, kind, eps = m.group(1), int(m.group(2)), m.group(3), float(m.group(4))
            if kind not in ("O", "C", "P") or (muscle, eps) not in FICHES:
                continue
            out.append(dict(muscle=muscle, eps=eps, kind=kind, minutes=minutes, path=os.path.join(d, fn), name=fn))
    return out


# ---------------------------------------------------------------------------
# Segmentation en rampes et métriques
# ---------------------------------------------------------------------------
def smooth(x, t, window_s=1.0):
    dt = float(np.median(np.diff(t)))
    k = max(3, int(round(window_s / dt)))
    return uniform_filter1d(x, size=k, mode="nearest")


def ramps(t, P, min_dp=0.08):
    """Liste de (i_debut, i_fin, montee: bool) sur les extrema de la pression lissée."""
    Ps = smooth(P, t)
    peaks, _ = find_peaks(Ps, prominence=0.05)
    troughs, _ = find_peaks(-Ps, prominence=0.05)
    ext = sorted(set([0, len(t) - 1] + peaks.tolist() + troughs.tolist()))
    # fusion : on ne garde que les extrema qui produisent une excursion suffisante
    kept = [ext[0]]
    for i in ext[1:]:
        if abs(Ps[i] - Ps[kept[-1]]) >= min_dp:
            kept.append(i)
        elif (len(kept) >= 2) and ((Ps[i] - Ps[kept[-2]]) * (Ps[kept[-1]] - Ps[kept[-2]]) > 0) and (
            abs(Ps[i] - Ps[kept[-2]]) > abs(Ps[kept[-1]] - Ps[kept[-2]])
        ):
            kept[-1] = i  # même sens, plus loin : remplace l'extremum
    out = []
    for a, b in zip(kept[:-1], kept[1:]):
        dp = Ps[b] - Ps[a]
        if abs(dp) >= min_dp:
            out.append((a, b, dp > 0))
    return out


def level(F, t, i, half_window_s=0.5):
    m = (t >= t[i] - half_window_s) & (t <= t[i] + half_window_s)
    return float(np.median(F[m])) if m.sum() >= 3 else float(F[i])


def ramp_metrics(t, P, F, a, b):
    """Vitesse (ajustement linéaire de P(t) sur 10-90 % de la course), ΔP, ΔF, gain."""
    Pa, Pb = level(P, t, a), level(P, t, b)
    dP = Pb - Pa
    lo, hi = Pa + 0.1 * dP, Pa + 0.9 * dP
    seg = slice(a, b + 1)
    Pseg, tseg = smooth(P[seg], t[seg]), t[seg]
    if dP > 0:
        inside = (Pseg >= min(lo, hi)) & (Pseg <= max(lo, hi))
    else:
        inside = (Pseg <= max(lo, hi)) & (Pseg >= min(lo, hi))
    if inside.sum() < 4:
        return None
    i0, i1 = int(np.argmax(inside)), int(len(inside) - 1 - np.argmax(inside[::-1]))
    tt, pp = tseg[i0:i1 + 1], P[seg][i0:i1 + 1]
    slope = float(np.polyfit(tt, pp, 1)[0])
    dF = level(F, t, b) - level(F, t, a)
    return dict(
        i_start=int(a), i_end=int(b), up=bool(dP > 0), t_start=float(t[a]), t_end=float(t[b]),
        duration_s=float(t[b] - t[a]), duration_1090_s=float(tt[-1] - tt[0]),
        dP_MPa=float(dP), rate_MPa_s=slope, abs_rate_MPa_s=abs(slope), P_start=float(Pa), P_end=float(Pb),
        dF_mN=float(dF), gain_mN_per_MPa=float(dF / dP) if abs(dP) > 1e-9 else float("nan"),
    )


def loop_metrics(t, P, F, a, m, b, n_bins=40):
    """Cycle montée (a→m) puis descente (m→b) : aire sur plage commune de
    pression (7.5.1), brute et corrigée d'une dérive linéaire en temps entre
    les deux extrémités à basse pression ; largeur à mi-course ; retard."""
    up, down = slice(a, m + 1), slice(m, b + 1)
    P_lo = max(level(P, t, a), level(P, t, b))
    P_hi = min(P[up].max(), P[down].max())
    if P_hi - P_lo < 0.08:
        return None
    grid = np.linspace(P_lo, P_hi, n_bins + 1)

    def branch(sl, Fv):
        Pb, Fb = P[sl], Fv[sl]
        centers, vals = [], []
        for g0, g1 in zip(grid[:-1], grid[1:]):
            sel = (Pb >= g0) & (Pb < g1)
            if sel.sum() >= 1:
                centers.append(0.5 * (g0 + g1))
                vals.append(float(np.mean(Fb[sel])))
        return np.array(centers), np.array(vals)

    def area(Fv):
        cu, fu = branch(up, Fv)
        cd, fd = branch(down, Fv)
        if len(cu) < 5 or len(cd) < 5:
            return float("nan"), float("nan")
        c = np.linspace(P_lo, P_hi, n_bins + 1)
        fu_i, fd_i = np.interp(c, cu, fu), np.interp(c, cd, fd)
        A = float(np.trapezoid(fd_i - fu_i, c)) * 10.0  # mN·MPa -> mN·bar
        mid = 0.5 * (P_lo + P_hi)
        width = float(np.interp(mid, cd, fd) - np.interp(mid, cu, fu))
        return A, width

    A_raw, width_raw = area(F)
    F_a, F_b = level(F, t, a), level(F, t, b)
    drift = (F_b - F_a) / max(t[b] - t[a], 1e-9)
    F_cor = F - drift * (t - t[a])
    A_cor, width_cor = area(F_cor)
    i_pmax = a + int(np.argmax(P[a:b + 1]))
    i_fmax = a + int(np.argmax(F[a:b + 1]))
    return dict(
        i_start=int(a), i_peak=int(m), i_end=int(b), t_start=float(t[a]), t_end=float(t[b]),
        duration_s=float(t[b] - t[a]), P_lo=float(P_lo), P_hi=float(P_hi),
        area_raw_mNbar=A_raw, area_cor_mNbar=A_cor, width_mid_raw_mN=width_raw, width_mid_cor_mN=width_cor,
        drift_mN_s=float(drift), peak_force_mN=float(F[a:b + 1].max()), lag_s=float(t[i_fmax] - t[i_pmax]),
        gain_cycle_mN_per_MPa=float((F[a:b + 1].max() - F_a) / max(P_hi - P_lo, 1e-9)),
    )


def analyse_signal(t, P, F, rps, kind="C"):
    """Métriques par rampe et par cycle (montée suivie d'une descente).
    Pour un essai à palier (3 P) seule la rampe initiale de montée est une
    sollicitation : la « descente » est la fuite du circuit, pas une rampe."""
    if kind == "P":
        rps = [r for r in rps[:1] if r[2]]
    rm = [ramp_metrics(t, P, F, a, b) for a, b, _ in rps]
    rm = [r for r in rm if r is not None]
    cycles = []
    for k in range(len(rps) - 1):
        a, m, up1 = rps[k]
        m2, b, up2 = rps[k + 1]
        if up1 and not up2 and m2 == m:
            lm = loop_metrics(t, P, F, a, m, b)
            if lm is not None:
                r_up = ramp_metrics(t, P, F, a, m)
                r_dn = ramp_metrics(t, P, F, m, b)
                if r_up and r_dn:
                    lm["rate_up_MPa_s"] = r_up["rate_MPa_s"]
                    lm["rate_down_MPa_s"] = r_dn["rate_MPa_s"]
                    lm["rate_mean_MPa_s"] = 0.5 * (abs(r_up["rate_MPa_s"]) + abs(r_dn["rate_MPa_s"]))
                    lm["rate_cycle_MPa_s"] = 2.0 * (lm["P_hi"] - lm["P_lo"]) / lm["duration_s"]
                    lm["gain_up_mN_per_MPa"] = r_up["gain_mN_per_MPa"]
                    lm["gain_down_mN_per_MPa"] = r_dn["gain_mN_per_MPa"]
                    cycles.append(lm)
    return rm, cycles


# ---------------------------------------------------------------------------
# Simulation sous pression mesurée
# ---------------------------------------------------------------------------
def build_settings(muscle, eps, spectre, v4=False):
    L0, alpha0, Ld = FICHES[(muscle, eps)]
    s = dict(parametres.DEFAULT_SETTINGS)
    s.update(COMMON)
    s.update(dict(initial_length_mm=L0, alpha0_deg=alpha0, uncoiled_length_mm=Ld, eps=eps,
                  n_layers=3, n_phi=8, pre_steps=24, dt=DT_SIM,
                  prestrain_reference_mode="viscoelastic_history",
                  maxwell_anisotropy_mode="axial_test_only"))
    sp = SPECTRES[spectre]
    s.update(dict(maxwell_E0_mpa=sp["E0"], maxwell_E1_mpa=sp["E1"], maxwell_eta1_mpa_s=sp["eta1"],
                  maxwell_E2_mpa=sp["E2"], maxwell_eta2_mpa_s=sp["eta2"],
                  maxwell_E3_mpa=sp["E3"], maxwell_eta3_mpa_s=sp["eta3"]))
    if v4:
        s.update(V4_DEMO)
    return s


def simulate(settings, pressure_time, pressure_mpa):
    cfg = parametres.build_config(settings)
    disc = Base.default_discretization(n_layers=cfg.n_layers, n_phi=cfg.n_phi, pre_steps=cfg.pre_steps,
                                       dw_bracket=(-0.05, 0.05))
    model = Base.TCPAMaxwellBlockedModel(mat=cfg.mat, geom=cfg.geom, disc=disc, integration=cfg.integration,
                                         prestrain_reference_mode=cfg.prestrain_reference_mode)
    model.prestretch_to(cfg.eps, strain_rate_mm_min=PRESTRETCH_RATE_MM_MIN)
    t_start = model.helix.time
    model.step(float(pressure_mpa[0]), 0.0, h_target=model.h_blocked)
    model.lock_blocked_series_reference()
    i0 = len(model.history) - 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.run_pressure_history(t_start + np.asarray(pressure_time, float), np.asarray(pressure_mpa, float))
    arr = model.history_arrays()
    return np.asarray(arr["time"], float)[i0:] - t_start, np.asarray(arr["force_mN"], float)[i0:]


def job_measured(task):
    """Un essai × une variante de modèle : renvoie la force simulée sur la grille de l'essai."""
    test, variant = task
    spectre, v4 = variant
    t, P, F = load_test(test["path"])
    t_ds = np.arange(0.0, float(t[-1]), DT_SIM)
    p_ds = np.clip(np.interp(t_ds, t, P), 0.0, None)
    try:
        settings = build_settings(test["muscle"], test["eps"], spectre, v4)
        t_sim, f_sim = simulate(settings, t_ds, p_ds)
        f_on = np.interp(t, t_sim, f_sim)
        return dict(name=test["name"], variant=f"{spectre}{'+V4' if v4 else ''}", ok=True, f_sim=f_on.tolist())
    except Exception as exc:  # noqa: BLE001
        return dict(name=test["name"], variant=f"{spectre}{'+V4' if v4 else ''}", ok=False, error=f"{type(exc).__name__}: {exc}")


def job_sweep(task):
    """Balayage synthétique : 3 cycles triangulaires à vitesse imposée."""
    muscle, eps, spectre, v4, rate, pmax = task
    settings = build_settings(muscle, eps, spectre, v4)
    half = pmax / rate
    dt = min(DT_SIM, half / 60.0)
    settings.update(dict(p_max_mpa=pmax, pressure_rate_mpa_s=rate, n_cycles=3, dt=dt, use_fixed_duration=False))
    try:
        t, p = Base.cyclic_pressure_history(3, pmax, dt=dt, pressure_rate_mpa_s=rate)
        t_sim, f_sim = simulate(settings, t, p)
        f_on = np.interp(t, t_sim, f_sim)
        rps = [(int(k * len(t) // 6), int((k + 1) * len(t) // 6), k % 2 == 0) for k in range(6)]
        # bornes exactes des rampes : indices des transitions
        edges = [int(np.argmin(np.abs(t - k * half))) for k in range(7)]
        rps = [(edges[k], edges[k + 1], k % 2 == 0) for k in range(6)]
        _, cycles = analyse_signal(t, p, f_on, rps)
        out = dict(muscle=muscle, eps=eps, variant=f"{spectre}{'+V4' if v4 else ''}", rate=rate, pmax=pmax,
                   half_period_s=half, dt=dt, ok=True,
                   cycles=[{k: v for k, v in c.items() if not k.startswith("i_")} for c in cycles],
                   F0=float(f_on[0]), gain_first_mN=float(f_on[:len(t) // 6 + 1].max() - f_on[0]))
        return out
    except Exception as exc:  # noqa: BLE001
        return dict(muscle=muscle, eps=eps, variant=f"{spectre}{'+V4' if v4 else ''}", rate=rate, pmax=pmax, ok=False,
                    error=f"{type(exc).__name__}: {exc}")


# ---------------------------------------------------------------------------
# Statistiques
# ---------------------------------------------------------------------------
def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return dict(n=int(m.sum()), rho=float("nan"), p=float("nan"))
    r, p = spearmanr(x[m], y[m])
    return dict(n=int(m.sum()), rho=float(r), p=float(p))


def theil(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return dict(n=int(m.sum()), slope=float("nan"), lo=float("nan"), hi=float("nan"), intercept=float("nan"))
    s, b, lo, hi = theilslopes(y[m], x[m])
    return dict(n=int(m.sum()), slope=float(s), lo=float(lo), hi=float(hi), intercept=float(b))


def within_test_slopes(records, xkey, ykey, min_points=3):
    """Pente de Theil-Sen de y(x) à l'intérieur de chaque essai (élimine l'effet spécimen)."""
    out = []
    by = {}
    for r in records:
        by.setdefault(r["test"], []).append(r)
    for name, recs in by.items():
        x = [r[xkey] for r in recs]
        y = [r[ykey] for r in recs]
        if len(recs) >= min_points:
            th = theil(x, y)
            out.append(dict(test=name, n=len(recs), slope=th["slope"], x_range=[float(min(x)), float(max(x))]))
    return out


def sign_summary(slopes):
    s = np.array([v for v in slopes if np.isfinite(v)])
    if len(s) == 0:
        return dict(n=0)
    res = dict(n=int(len(s)), median=float(np.median(s)), n_pos=int((s > 0).sum()), n_neg=int((s < 0).sum()))
    if len(s) >= 5 and np.any(s != 0):
        try:
            res["wilcoxon_p"] = float(wilcoxon(s).pvalue)
        except ValueError:
            pass
    return res


# ---------------------------------------------------------------------------
def main():
    tests = list_tests()
    print(f"{len(tests)} essais 1 O / 3 C / 3 P avec fiche géométrique", flush=True)
    variants = [("10N", False), ("30L", False), ("10N", True)]
    cache_path = os.path.join(OUT, "sims_cache.json")
    cached = []
    if os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8") as fh:
            cached = [c for c in json.load(fh) if c.get("ok")]
    done = {(c["name"], c["variant"]) for c in cached}
    tasks = [(tst, v) for tst in tests for v in variants
             if (tst["name"], f"{v[0]}{'+V4' if v[1] else ''}") not in done]
    print(f"{len(cached)} simulations en cache, {len(tasks)} à calculer", flush=True)
    sims = list(cached)
    if tasks:
        with Pool(WORKERS) as pool:
            for k, res in enumerate(pool.imap_unordered(job_measured, tasks, chunksize=1), 1):
                sims.append(res)
                if k % 20 == 0 or k == len(tasks):
                    print(f"  {k}/{len(tasks)} simulations", flush=True)
                    with open(cache_path, "w", encoding="utf-8") as fh:
                        json.dump(sims, fh)
    sim_by = {(s["name"], s["variant"]): s for s in sims}
    failures = [s for s in sims if not s["ok"]]
    print(f"{len(sims)} simulations, {len(failures)} échecs", flush=True)
    for f in failures[:10]:
        print("  ", f["name"], f["variant"], f["error"])

    ramp_records, cycle_records, test_summaries = [], [], []
    overlay = {}
    for tst in tests:
        t, P, F = load_test(tst["path"])
        rps = ramps(t, P)
        rm_exp, cy_exp = analyse_signal(t, P, F, rps, tst["kind"])
        summary = dict(name=tst["name"], muscle=tst["muscle"], eps=tst["eps"], kind=tst["kind"],
                       n_ramps=len(rm_exp), n_cycles=len(cy_exp), P_max=float(P.max()),
                       rates=[r["abs_rate_MPa_s"] for r in rm_exp])
        for r in rm_exp:
            rec = dict(test=tst["name"], muscle=tst["muscle"], eps=tst["eps"], kind=tst["kind"], source="exp", **r)
            ramp_records.append(rec)
        for c in cy_exp:
            rec = dict(test=tst["name"], muscle=tst["muscle"], eps=tst["eps"], kind=tst["kind"], source="exp", **c)
            cycle_records.append(rec)
        for spectre, v4 in variants:
            key = (tst["name"], f"{spectre}{'+V4' if v4 else ''}")
            s = sim_by.get(key)
            if not s or not s["ok"]:
                continue
            f_sim = np.asarray(s["f_sim"])
            rm_sim, cy_sim = analyse_signal(t, P, f_sim, rps, tst["kind"])
            for r in rm_sim:
                ramp_records.append(dict(test=tst["name"], muscle=tst["muscle"], eps=tst["eps"], kind=tst["kind"],
                                         source=key[1], **r))
            for c in cy_sim:
                cycle_records.append(dict(test=tst["name"], muscle=tst["muscle"], eps=tst["eps"], kind=tst["kind"],
                                          source=key[1], **c))
            if tst["kind"] == "C" and tst["muscle"] == "J" and tst["eps"] == 1.0:
                overlay[key[1]] = f_sim[::3].tolist()
        if tst["kind"] == "C" and tst["muscle"] == "J" and tst["eps"] == 1.0:
            overlay["t"], overlay["P"], overlay["F"] = t[::3].tolist(), P[::3].tolist(), F[::3].tolist()
        test_summaries.append(summary)

    # --- balayage synthétique --------------------------------------------------
    rates = [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.1667, 0.5]
    sweep_tasks = [("J", 1.0, sp, v4, r, 0.45) for sp, v4 in variants for r in rates]
    sweep_tasks += [("J", 1.0, "10N", False, r, 1.5) for r in rates]  # amplitude de l'article
    sweep_cache = os.path.join(OUT, "sweep_cache.json")
    if os.path.exists(sweep_cache):
        with open(sweep_cache, encoding="utf-8") as fh:
            sweep = json.load(fh)
    else:
        with Pool(WORKERS) as pool:
            sweep = pool.map(job_sweep, sweep_tasks, chunksize=1)
        with open(sweep_cache, "w", encoding="utf-8") as fh:
            json.dump(sweep, fh)
    print(f"balayage : {sum(1 for s in sweep if s['ok'])}/{len(sweep)} réussis", flush=True)

    # --- statistiques ------------------------------------------------------------
    stats = {}
    sources = ["exp"] + [f"{sp}{'+V4' if v4 else ''}" for sp, v4 in variants]
    for src in sources:
        rr = [r for r in ramp_records if r["source"] == src]
        cc = [c for c in cycle_records if c["source"] == src]
        st = {}
        for kind in ("O", "C", "P", "all"):
            sel = [r for r in rr if kind == "all" or r["kind"] == kind]
            up = [r for r in sel if r["up"]]
            dn = [r for r in sel if not r["up"]]
            st[f"gain_up_vs_rate_{kind}"] = dict(spearman=spearman([r["abs_rate_MPa_s"] for r in up], [r["gain_mN_per_MPa"] for r in up]),
                                                theil=theil([r["abs_rate_MPa_s"] for r in up], [r["gain_mN_per_MPa"] for r in up]),
                                                median_gain=float(np.nanmedian([r["gain_mN_per_MPa"] for r in up])) if up else float("nan"))
            st[f"gain_down_vs_rate_{kind}"] = dict(spearman=spearman([r["abs_rate_MPa_s"] for r in dn], [r["gain_mN_per_MPa"] for r in dn]),
                                                  median_gain=float(np.nanmedian([r["gain_mN_per_MPa"] for r in dn])) if dn else float("nan"))
        for kind in ("O", "C", "all"):
            sel = [c for c in cc if kind == "all" or c["kind"] == kind]
            st[f"area_cor_vs_rate_{kind}"] = dict(spearman=spearman([c["rate_mean_MPa_s"] for c in sel], [c["area_cor_mNbar"] for c in sel]),
                                                 theil=theil([c["rate_mean_MPa_s"] for c in sel], [c["area_cor_mNbar"] for c in sel]),
                                                 median_area=float(np.nanmedian([c["area_cor_mNbar"] for c in sel])) if sel else float("nan"),
                                                 median_area_raw=float(np.nanmedian([c["area_raw_mNbar"] for c in sel])) if sel else float("nan"),
                                                 median_width=float(np.nanmedian([c["width_mid_cor_mN"] for c in sel])) if sel else float("nan"),
                                                 median_lag=float(np.nanmedian([c["lag_s"] for c in sel])) if sel else float("nan"),
                                                 n=len(sel))
            st[f"lag_vs_rate_{kind}"] = spearman([c["rate_mean_MPa_s"] for c in sel], [c["lag_s"] for c in sel])
        # pentes intra-essai (3 C : plusieurs cycles par essai, vitesse variable)
        cC = [c for c in cc if c["kind"] == "C"]
        st["within_test_area_slopes"] = within_test_slopes(cC, "rate_mean_MPa_s", "area_cor_mNbar")
        st["within_test_area_sign"] = sign_summary([w["slope"] for w in st["within_test_area_slopes"]])
        st["within_test_gain_slopes"] = within_test_slopes([r for r in rr if r["kind"] == "C" and r["up"]], "abs_rate_MPa_s", "gain_mN_per_MPa")
        st["within_test_gain_sign"] = sign_summary([w["slope"] for w in st["within_test_gain_slopes"]])
        stats[src] = st

    # rapport gain sim/exp par rampe (même rampe, même pression)
    exp_by = {(r["test"], r["i_start"]): r for r in ramp_records if r["source"] == "exp"}
    ratio_records = []
    for r in ramp_records:
        if r["source"] == "exp":
            continue
        e = exp_by.get((r["test"], r["i_start"]))
        if e and np.isfinite(e["gain_mN_per_MPa"]) and abs(e["gain_mN_per_MPa"]) > 1e-9:
            ratio_records.append(dict(test=r["test"], source=r["source"], kind=r["kind"], up=r["up"],
                                      rate=r["abs_rate_MPa_s"], ratio=r["gain_mN_per_MPa"] / e["gain_mN_per_MPa"]))
    ratio_stats = {}
    for src in sources[1:]:
        for kind in ("O", "C", "P", "all"):
            sel = [x for x in ratio_records if x["source"] == src and x["up"] and (kind == "all" or x["kind"] == kind)]
            ratio_stats[f"{src}_{kind}"] = dict(n=len(sel), median_ratio=float(np.median([x["ratio"] for x in sel])) if sel else float("nan"),
                                                 spearman_ratio_vs_rate=spearman([x["rate"] for x in sel], [x["ratio"] for x in sel]))

    rate_summary = {}
    for kind in ("O", "C", "P"):
        v = np.array([r["abs_rate_MPa_s"] for r in ramp_records if r["source"] == "exp" and r["kind"] == kind])
        vu = np.array([r["abs_rate_MPa_s"] for r in ramp_records if r["source"] == "exp" and r["kind"] == kind and r["up"]])
        if len(v):
            rate_summary[kind] = dict(n=int(len(v)), n_up=int(len(vu)), min=float(v.min()), q1=float(np.percentile(v, 25)),
                                      median=float(np.median(v)), q3=float(np.percentile(v, 75)), max=float(v.max()),
                                      median_up=float(np.median(vu)) if len(vu) else float("nan"))

    result = dict(engine=Base.MODEL_VERSION, n_tests=len(tests), tests=test_summaries, rate_summary=rate_summary,
                  ramps=ramp_records, cycles=cycle_records, stats=stats, ratio_stats=ratio_stats,
                  sweep=sweep, overlay_J_3C_1_0=overlay, failures=failures, spectres=SPECTRES, v4_demo=V4_DEMO)
    with open(os.path.join(OUT, "confrontation_vitesse.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=1)
    print("JSON écrit", flush=True)
    write_tables(result)
    make_figures(result)
    print("TERMINE", flush=True)


def fmt(x, d=2):
    return "—" if x is None or not np.isfinite(x) else f"{x:.{d}f}"


def write_tables(res):
    L = ["# Tables générées par confrontation_vitesse.py", "", f"Moteur {res['engine']} ; {res['n_tests']} essais.", ""]
    L += ["## Vitesses mesurées (|dP/dt| sur 10-90 % de la course, MPa/s)", "", "| Essai | n rampes | min | Q1 | médiane | Q3 | max |", "|---|---|---|---|---|---|---|"]
    for k, v in res["rate_summary"].items():
        L.append(f"| {k} | {v['n']} | {fmt(v['min'],3)} | {fmt(v['q1'],3)} | {fmt(v['median'],3)} | {fmt(v['q3'],3)} | {fmt(v['max'],3)} |")
    L += ["", "## Gain de rampe en montée ΔF/ΔP (mN/MPa) et dépendance à la vitesse (Spearman)", "",
          "| Source | Essai | n | gain médian | ρ(gain, vitesse) | p |", "|---|---|---|---|---|---|"]
    for src, st in res["stats"].items():
        for kind in ("O", "C", "P", "all"):
            g = st[f"gain_up_vs_rate_{kind}"]
            L.append(f"| {src} | {kind} | {g['spearman']['n']} | {fmt(g['median_gain'],0)} | {fmt(g['spearman']['rho'])} | {fmt(g['spearman']['p'],3)} |")
    L += ["", "## Aire de boucle corrigée (mN·bar) et dépendance à la vitesse", "",
          "| Source | Essai | n cycles | aire médiane (brute) | largeur mi-course | retard F/P (s) | ρ(aire, vitesse) | p | pente Theil (mN·bar par MPa/s) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for src, st in res["stats"].items():
        for kind in ("O", "C", "all"):
            a = st[f"area_cor_vs_rate_{kind}"]
            L.append(f"| {src} | {kind} | {a['n']} | {fmt(a['median_area'],1)} ({fmt(a['median_area_raw'],1)}) | {fmt(a['median_width'],1)} | {fmt(a['median_lag'],2)} | {fmt(a['spearman']['rho'])} | {fmt(a['spearman']['p'],3)} | {fmt(a['theil']['slope'],0)} [{fmt(a['theil']['lo'],0)} ; {fmt(a['theil']['hi'],0)}] |")
    L += ["", "## Pentes intra-essai (3 C) de l'aire corrigée vs vitesse moyenne du cycle", "",
          "| Source | n essais | pente médiane (mN·bar par MPa/s) | n > 0 | n < 0 | Wilcoxon p |", "|---|---|---|---|---|---|"]
    for src, st in res["stats"].items():
        s = st["within_test_area_sign"]
        L.append(f"| {src} | {s.get('n',0)} | {fmt(s.get('median',float('nan')),0)} | {s.get('n_pos','—')} | {s.get('n_neg','—')} | {fmt(s.get('wilcoxon_p',float('nan')),3)} |")
    L += ["", "## Rapport gain simulé / gain mesuré, rampes de montée", "", "| Variante | Essai | n | rapport médian | ρ(rapport, vitesse) | p |", "|---|---|---|---|---|---|"]
    for k, v in res["ratio_stats"].items():
        src, kind = k.rsplit("_", 1)
        L.append(f"| {src} | {kind} | {v['n']} | {fmt(v['median_ratio'])} | {fmt(v['spearman_ratio_vs_rate']['rho'])} | {fmt(v['spearman_ratio_vs_rate']['p'],3)} |")
    L += ["", "## Balayage synthétique (muscle J, ε = 1,0, 3 cycles triangulaires, cycle n° 2)", "",
          "| Variante | Pmax | vitesse (MPa/s) | demi-cycle (s) | gain 1er cycle (mN) | aire brute | aire corrigée | largeur mi-course | retard (s) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for s in res["sweep"]:
        if not s["ok"] or len(s["cycles"]) < 2:
            L.append(f"| {s['variant']} | {s['pmax']} | {s['rate']} | — | — | échec ou cycle absent | | | |")
            continue
        c = s["cycles"][1]
        L.append(f"| {s['variant']} | {s['pmax']} | {s['rate']} | {fmt(s['half_period_s'],1)} | {fmt(s['gain_first_mN'],0)} | {fmt(c['area_raw_mNbar'],1)} | {fmt(c['area_cor_mNbar'],1)} | {fmt(c['width_mid_cor_mN'],1)} | {fmt(c['lag_s'],2)} |")
    with open(os.path.join(OUT, "tables.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("tables.md écrit", flush=True)


def make_figures(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    C = {"exp": "#f3a43b", "10N": "#67a4ff", "30L": "#3fbf7f", "10N+V4": "#c678dd"}
    # 1. distribution des vitesses mesurées
    fig, ax = plt.subplots(figsize=(8, 3.6), constrained_layout=True)
    for kind, lab in (("O", "1 O"), ("C", "3 C"), ("P", "3 P")):
        v = [r["abs_rate_MPa_s"] for r in res["ramps"] if r["source"] == "exp" and r["kind"] == kind]
        if v:
            ax.hist(v, bins=np.linspace(0, 0.12, 25), alpha=0.6, label=f"{lab} ({len(v)} rampes)")
    ax.axvline(1.5 / 9.0, color="k", ls="--", lw=1, label="article : 1,5 MPa / 9 s")
    ax.set_xlabel("|dP/dt| mesuré (MPa/s)")
    ax.set_ylabel("rampes")
    ax.legend(fontsize=8)
    ax.set_title("Vitesses de sollicitation de la campagne (banc non asservi)")
    fig.savefig(os.path.join(OUT, "fig1_vitesses_mesurees.png"), dpi=130)

    # 2. aire corrigée vs vitesse : exp et modèles (cycles 3 C + 1 O)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    for ax, kind, title in zip(axes, ("C", "O"), ("Cyclage 3 C (un point par cycle)", "Rampe lente 1 O (un point par essai)")):
        for src in ("exp", "10N", "30L", "10N+V4"):
            cc = [c for c in res["cycles"] if c["source"] == src and c["kind"] == kind]
            if cc:
                ax.scatter([c["rate_mean_MPa_s"] for c in cc], [c["area_cor_mNbar"] for c in cc], s=16, alpha=0.75,
                           color=C[src], label={"exp": "mesure", "10N": "modèle, spectre 10 N", "30L": "modèle, spectre 30 L", "10N+V4": "modèle + V4 (P_r0, P_c, r du muscle J)"}[src])
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel("vitesse moyenne du cycle (MPa/s)")
        ax.set_ylabel("aire de boucle corrigée (mN·bar)")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.savefig(os.path.join(OUT, "fig2_aire_vs_vitesse.png"), dpi=130)

    # 3. gain de montée vs vitesse
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for ax, kind, title in zip(axes, ("O", "C", "P"), ("1 O", "3 C", "3 P (rampe initiale)")):
        for src in ("exp", "10N", "30L", "10N+V4"):
            rr = [r for r in res["ramps"] if r["source"] == src and r["kind"] == kind and r["up"]]
            if rr:
                ax.scatter([r["abs_rate_MPa_s"] for r in rr], [r["gain_mN_per_MPa"] for r in rr], s=16, alpha=0.75, color=C[src],
                           label={"exp": "mesure", "10N": "modèle 10 N", "30L": "modèle 30 L", "10N+V4": "modèle + V4"}[src])
        ax.set_xlabel("|dP/dt| de la montée (MPa/s)")
        ax.set_ylabel("gain de montée ΔF/ΔP (mN/MPa)")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.savefig(os.path.join(OUT, "fig3_gain_vs_vitesse.png"), dpi=130)

    # 4. balayage synthétique
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    for variant, col in (("10N", C["10N"]), ("30L", C["30L"]), ("10N+V4", C["10N+V4"])):
        ss = sorted([s for s in res["sweep"] if s["ok"] and s["variant"] == variant and s["pmax"] == 0.45 and len(s["cycles"]) >= 2], key=lambda s: s["rate"])
        if ss:
            axes[0].plot([s["rate"] for s in ss], [s["cycles"][1]["area_cor_mNbar"] for s in ss], "o-", color=col, label=variant + " (Pmax 0,45 MPa)")
            axes[1].plot([s["rate"] for s in ss], [s["gain_first_mN"] for s in ss], "o-", color=col, label=variant)
    ss = sorted([s for s in res["sweep"] if s["ok"] and s["variant"] == "10N" and s["pmax"] == 1.5 and len(s["cycles"]) >= 2], key=lambda s: s["rate"])
    if ss:
        axes[0].plot([s["rate"] for s in ss], [s["cycles"][1]["area_cor_mNbar"] for s in ss], "s--", color="#888", label="10N (Pmax 1,5 MPa)")
        axes[1].plot([s["rate"] for s in ss], [s["gain_first_mN"] for s in ss], "s--", color="#888", label="10N (Pmax 1,5 MPa)")
    exp_rates = [c["rate_mean_MPa_s"] for c in res["cycles"] if c["source"] == "exp"]
    for ax in axes:
        ax.set_xscale("log")
        if exp_rates:
            ax.axvspan(min(exp_rates), max(exp_rates), color="#f3a43b", alpha=0.15, label="plage des essais")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7)
        ax.set_xlabel("vitesse de pression imposée (MPa/s)")
    axes[0].axhline(0, color="k", lw=0.6)
    axes[0].set_ylabel("aire de boucle corrigée, cycle 2 (mN·bar)")
    axes[1].set_ylabel("gain du 1er cycle (mN)")
    fig.suptitle("Balayage synthétique du moteur à vitesse imposée (muscle J, ε = 1,0)")
    fig.savefig(os.path.join(OUT, "fig4_balayage_synthetique.png"), dpi=130)

    # 5. superposition J 3 C 1.0
    ov = res["overlay_J_3C_1_0"]
    if "t" in ov:
        fig, ax = plt.subplots(figsize=(11, 4), constrained_layout=True)
        ax.plot(ov["t"], ov["F"], color=C["exp"], lw=1.2, label="force mesurée")
        for src in ("10N", "30L", "10N+V4"):
            if src in ov:
                ax.plot(ov["t"], ov[src], color=C[src], lw=1.1, label=f"simulée ({src})")
        ax2 = ax.twinx()
        ax2.plot(ov["t"], ov["P"], color="#ff5c68", lw=0.8, alpha=0.7)
        ax2.set_ylabel("P (MPa)", color="#ff5c68")
        ax.set_xlabel("t (s)")
        ax.set_ylabel("F (mN)")
        ax.set_title("Muscle J, 3 C, ε = 1,0 — force mesurée et simulée sous la pression mesurée")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.savefig(os.path.join(OUT, "fig5_superposition_J_3C.png"), dpi=130)
    print("figures écrites", flush=True)


if __name__ == "__main__":
    main()
