"""Corrélations simulation / mesure des essais "1 O", sous les deux échelles
de pression, toutes choses égales par ailleurs.

Seule l'échelle de pression injectée dans le moteur change d'une colonne à
l'autre : même géométrie, même spectre, mêmes réglages. Ce sont donc les
écarts entre colonnes qui sont exploitables, non les valeurs absolues — le
spectre employé ici est celui de la Table A1 et non le spectre ajusté sur les
essais « 10 N ».

À lancer depuis le dossier "Espace de travail" :
    python essais/recalcul_correlations.py
"""
import csv, importlib.util, io, json, os, sys
import numpy as np

V3 = os.path.abspath('alpha V3')
sys.path.insert(0, V3)
spec = importlib.util.spec_from_file_location('_base', os.path.join(V3, 'Base.py'))
Base = importlib.util.module_from_spec(spec)
sys.modules['_base'] = Base
spec.loader.exec_module(Base)

PSI_TO_BAR = 0.0689476
psi = lambda a: (a * 5.0 / 1023.0 - 0.5) * 500.0 / 4.0


def load_csv(path):
    rows = list(csv.reader(io.StringIO(io.open(path, encoding='utf-8-sig').read())))
    idx = {k: i for i, k in enumerate(rows[0])}
    data = [r for r in rows[1:] if r and r[0]]
    g = lambda k, f=float: np.array([f(r[idx[k]]) for r in data])
    return g('time_s'), g('force_unfiltered_mN'), g('pressure_adc', int), g('pressure_bar')


SPECTRE = os.path.join('essais', 'spectre_identifie.json')
if os.path.exists(SPECTRE):
    s = json.load(open(SPECTRE, encoding='utf-8'))
    SPEC = dict(E0=s['maxwell_E0_mpa'],
                E1=s['maxwell_E1_mpa'], eta1=s['maxwell_eta1_mpa_s'],
                E2=s['maxwell_E2_mpa'], eta2=s['maxwell_eta2_mpa_s'],
                E3=s['maxwell_E3_mpa'], eta3=s['maxwell_eta3_mpa_s'])
    ORIGINE = 'spectre identifie sur les essais 10 N'
else:
    SPEC = None
    ORIGINE = 'spectre Table A1 de l article'


def simule(fiche, t, P_mpa):
    """Actionnement bloqué sous la pression mesurée. Retourne (t, F_act, masque)."""
    mx = Base.default_maxwell_tensile_params(**SPEC) if SPEC else Base.default_maxwell_tensile_params()
    mat = Base.default_material_params(
        maxwell=mx, E_axial=mx.E0 + mx.E1 + mx.E2 + mx.E3,
        maxwell_anisotropy_mode='axial_test_only')
    geom = Base.default_geometry_params(
        Rout=fiche['rout_mm'], Rin=fiche['rin_mm'],
        r_nylon=0.5 * fiche['nylon_diameter_mm'], rho0=fiche['rho0_mm'],
        alpha0_deg=fiche['alpha0_deg'], theta_f_deg=fiche['theta_f_deg'],
        initial_length=fiche['initial_length_mm'],
        uncoiled_length=fiche['uncoiled_length_mm'])
    cfg = Base.default_simulation_config(
        mat=mat, geom=geom, eps=fiche['eps'],
        n_layers=3, n_phi=8, pre_steps=24, dt=0.5,
        prestrain_reference_mode='viscoelastic_history')
    keep = np.concatenate(([True], np.diff(t) > 1e-9))
    tt, pp = t[keep], np.maximum(P_mpa[keep], 0.0)
    _, arr = Base.run_blocked_actuation(cfg, pressure_time=tt, pressure_MPa=pp)
    return tt, np.interp(tt, arr['time'], arr['force_act_mN']), keep


FICHES = os.path.join(V3, 'identification', 'fiches')
print('Materiau : %s' % ORIGINE)
if SPEC:
    print('  E0=%.3f  E1=%.3f (eta=%.1f)  E2=%.3f (eta=%.1f)' %
          (SPEC['E0'], SPEC['E1'], SPEC['eta1'], SPEC['E2'], SPEC['eta2']))
print()
print('%-8s | %8s %8s %7s | %7s %7s %7s %7s' %
      ('cas', 'corr_anc', 'corr_cor', 'ecart', 'gain_me', 'gain_an', 'gain_co', 'ratio'))
print('-' * 74)
res = []
for m in ('J', 'L', 'Q'):
    for eps, tag in (('0.8', 'eps08'), ('1.0', 'eps10'), ('1.2', 'eps12')):
        fpath = os.path.join(FICHES, 'muscle_%s_%s.json' % (m, tag))
        cpath = os.path.join('essais', 'Muscle %s' % m, 'Muscle %s 1 O %s.csv' % (m, eps))
        if not (os.path.exists(fpath) and os.path.exists(cpath)):
            continue
        fiche = json.load(open(fpath))
        t, F, adc, P_old = load_csv(cpath)
        rest = int(np.bincount(adc[:15]).argmax())
        P_cor = (psi(adc.astype(float)) - psi(rest)) * PSI_TO_BAR * 0.1
        P_anc = P_old * 0.1
        dF = F - F[:10].mean()

        gain_mesure = float(np.nanmax(dF))
        corr, gain = [], []
        for P in (P_anc, P_cor):
            try:
                tt, Fsim, keep = simule(fiche, t, P)
                corr.append(float(np.corrcoef(dF[keep], Fsim)[0, 1]))
                gain.append(float(np.nanmax(Fsim)))
            except Exception as exc:
                corr.append(np.nan)
                gain.append(np.nan)
                print('   (%s %s : %s)' % (m, eps, str(exc)[:60]))
        print('%-8s | %8.4f %8.4f %+7.4f | %7.0f %7.0f %7.0f %6.2f' %
              ('%s %s' % (m, eps), corr[0], corr[1], corr[1] - corr[0],
               gain_mesure, gain[0], gain[1], gain[1] / gain_mesure))
        res.append([corr[0], corr[1], gain_mesure, gain[0], gain[1]])

r = np.array(res)
print()
print('corr ancienne echelle     : %.3f a %.3f  (mediane %.3f)' % (np.nanmin(r[:, 0]), np.nanmax(r[:, 0]), np.nanmedian(r[:, 0])))
print('corr echelle reconstruite : %.3f a %.3f  (mediane %.3f)' % (np.nanmin(r[:, 1]), np.nanmax(r[:, 1]), np.nanmedian(r[:, 1])))
print('ecart median : %+.4f' % np.nanmedian(r[:, 1] - r[:, 0]))
print()
ra, rc = r[:, 3] / r[:, 2], r[:, 4] / r[:, 2]
print('gain simule / gain mesure, ancienne echelle     : %.0f a %.0f %%  (mediane %.0f %%)'
      % (100 * np.nanmin(ra), 100 * np.nanmax(ra), 100 * np.nanmedian(ra)))
print('gain simule / gain mesure, echelle reconstruite : %.0f a %.0f %%  (mediane %.0f %%)'
      % (100 * np.nanmin(rc), 100 * np.nanmax(rc), 100 * np.nanmedian(rc)))
print()
n = 3
print('a eps = 0,8 seulement (%d cas) : %.0f a %.0f %% (reconstruite)'
      % (n, 100 * np.nanmin(rc[::3]), 100 * np.nanmax(rc[::3])))
