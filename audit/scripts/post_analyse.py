# -*- coding: utf-8 -*-
"""Post-analyse : hysteresis (charge vs decharge) et zone morte basse pression."""
import numpy as np

OUT = r"C:/Users/B3DCB~1.PER/AppData/Local/Temp/claude/C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail/c5852e10-6d38-4703-b218-a1649a5b12c9/scratchpad"


def hysteresis_at(t, p_bar, df, P0, t_load, t_unload, tol=0.12):
    """ΔF moyen a P0 pendant la charge et pendant la decharge."""
    out = []
    for (ta, tb) in (t_load, t_unload):
        m = (t >= ta) & (t <= tb) & (np.abs(p_bar - P0) < tol)
        out.append(float(np.mean(df[m])) if m.sum() >= 2 else np.nan)
    return out  # [charge, decharge]


for name, spans in {
    "essai_153353": dict(load=(5.0, 57.0), unload=(64.0, 98.0), Ps=[3.3, 3.5, 3.7]),
    "essai_155135": dict(load=(5.0, 31.0), unload=(76.0, 100.0), Ps=[3.0, 3.5, 4.0, 4.5]),
}.items():
    d = np.load(f"{OUT}/sim_{name}.npz")
    t, p, dfm, dfs = d["t"], d["p_mpa"] * 10.0, d["df_meas"], d["f_sim_act"]
    print("=" * 66)
    print(name)
    for P0 in spans["Ps"]:
        cm, dm = hysteresis_at(t, p, dfm, P0, spans["load"], spans["unload"])
        cs, ds = hysteresis_at(t, p, dfs, P0, spans["load"], spans["unload"])
        print(f"  P={P0:4.1f} bar | mesure: charge {cm:7.1f} / decharge {dm:7.1f} "
              f"(ecart {dm-cm:+6.1f}) | modele: charge {cs:6.2f} / decharge {ds:6.2f} "
              f"(ecart {ds-cs:+6.2f}) mN")
    # zone morte : reponse a basse pression (<= 2 bar) pendant la premiere charge
    m = (t <= spans["load"][1]) & (p > 0.2) & (p <= 2.0)
    if m.sum() > 3:
        sm = np.polyfit(p[m] / 10.0, dfm[m], 1)[0]
        ss = np.polyfit(p[m] / 10.0, dfs[m], 1)[0]
        print(f"  pente basse pression (0.02-0.2 MPa, 1re charge) : mesure {sm:.0f} "
              f"/ modele {ss:.1f} mN/MPa")
    m2 = (t <= spans["load"][1]) & (p > 2.0)
    if m2.sum() > 3:
        sm2 = np.polyfit(p[m2] / 10.0, dfm[m2], 1)[0]
        ss2 = np.polyfit(p[m2] / 10.0, dfs[m2], 1)[0]
        print(f"  pente haute pression (>0.2 MPa, 1re charge)     : mesure {sm2:.0f} "
              f"/ modele {ss2:.1f} mN/MPa")
