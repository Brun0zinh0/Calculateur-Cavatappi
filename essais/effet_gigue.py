"""Effet de la gigue d'horodatage sur les temps caractéristiques identifiés.

Les fichiers de la campagne sont datés par le PC à réception de la trame. La
date enregistrée vaut donc t_i = t_vrai_i + latence_i, où la latence varie
d'un échantillon à l'autre. Le PC n'ayant pas de dérive d'horloge, ces erreurs
ne s'accumulent pas : ce sont des perturbations indépendantes de chaque date.

Si l'échantillonnage réel est périodique de période T, l'intervalle observé
vaut Dt_i = T + eps_i - eps_{i-1}, d'où Var(Dt) = 2 Var(eps) et

    sigma_eps = sigma_Dt / sqrt(2)

On propage ensuite cette incertitude par Monte-Carlo : on reperturbe la grille
temporelle, on réajuste, et on observe la dispersion de tau1 et tau2.

À lancer depuis le dossier "Espace de travail" :
    python essais/effet_gigue.py
"""
import csv, io, os
import numpy as np
from scipy.optimize import curve_fit

N_TIRAGES = 40
RNG = np.random.default_rng(12345)


def load(path):
    rows = list(csv.reader(io.StringIO(io.open(path, encoding='utf-8-sig').read())))
    idx = {k: i for i, k in enumerate(rows[0])}
    data = [r for r in rows[1:] if r and r[0]]
    g = lambda k: np.array([float(r[idx[k]]) for r in data])
    return g('time_s'), g('force_unfiltered_mN')


def biexp(t, Finf, A1, tau1, A2, tau2):
    return Finf + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)


def ajuste(t, F):
    F0 = F[:10].mean()
    p0 = [0.90 * F0, 0.04 * F0, 20.0, 0.05 * F0, 400.0]
    bornes = ([0.0, 0.0, 2.0, 0.0, 80.0], [2.0 * F0, 0.6 * F0, 120.0, 0.6 * F0, 4000.0])
    popt, _ = curve_fit(biexp, t, F, p0=p0, bounds=bornes, maxfev=40000)
    Finf, A1, tau1, A2, tau2 = popt
    if tau1 > tau2:
        A1, tau1, A2, tau2 = A2, tau2, A1, tau1
    return tau1, tau2


# ---- 1. gigue reellement presente dans les fichiers -------------------------
fichiers = []
base = 'essais'
for muscle in sorted(os.listdir(base)):
    d = os.path.join(base, muscle)
    if not os.path.isdir(d):
        continue
    for fn in sorted(os.listdir(d)):
        if ' 10 N ' in fn and fn.endswith('.csv'):
            fichiers.append(os.path.join(d, fn))

sig_dt, medianes = [], []
for p in fichiers:
    t, _ = load(p)
    dt = np.diff(t)
    dt = dt[(dt > 0) & (dt < 1.0)]           # ecarte les trous eventuels
    sig_dt.append(dt.std())
    medianes.append(np.median(dt))

sd = float(np.median(sig_dt))
se = sd / np.sqrt(2.0)
print('Gigue mesuree sur les %d essais « 10 N »' % len(fichiers))
print('  periode mediane      : %.1f ms' % (1000 * np.median(medianes)))
print('  ecart-type des Dt    : %.1f ms  (mediane sur les essais)' % (1000 * sd))
print('  -> sigma sur la date : %.1f ms' % (1000 * se))
print('  a comparer a tau1 ~ 20 s : rapport %.0e' % (se / 20.0))
print()

# ---- 2. propagation Monte-Carlo -------------------------------------------
print('Propagation Monte-Carlo (%d tirages par essai)' % N_TIRAGES)
print()
print('%-9s | %8s %9s | %8s %9s' % ('cas', 'tau1', 'ecart-type', 'tau2', 'ecart-type'))
print('-' * 54)

rel1, rel2 = [], []
for p in fichiers:
    t, F = load(p)
    t = t - t[0]
    try:
        tau1_ref, tau2_ref = ajuste(t, F)
    except Exception:
        continue
    T1, T2 = [], []
    for _ in range(N_TIRAGES):
        tp = t + RNG.normal(0.0, se, size=t.size)
        o = np.argsort(tp)                    # la gigue peut inverser deux dates
        try:
            a, b = ajuste(tp[o] - tp[o][0], F[o])
            T1.append(a)
            T2.append(b)
        except Exception:
            pass
    if len(T1) < 5:
        continue
    s1, s2 = np.std(T1), np.std(T2)
    nom = os.path.basename(p).replace('Muscle ', '').replace('.csv', '').replace(' 10 N ', ' ')
    print('%-9s | %8.2f %8.4f s | %8.1f %8.3f s' % (nom, tau1_ref, s1, tau2_ref, s2))
    rel1.append(s1 / tau1_ref)
    rel2.append(s2 / tau2_ref)

r1, r2 = np.array(rel1), np.array(rel2)
print()
print('=' * 62)
print('EFFET DE LA GIGUE SUR LES TEMPS CARACTERISTIQUES')
print('=' * 62)
print('  tau1 : ecart-type relatif median %.2e  (max %.2e)' % (np.median(r1), r1.max()))
print('  tau2 : ecart-type relatif median %.2e  (max %.2e)' % (np.median(r2), r2.max()))
print()
print('  soit %.4f %% sur tau1 et %.4f %% sur tau2 en median.'
      % (100 * np.median(r1), 100 * np.median(r2)))
print()
print('A comparer a la dispersion entre specimens, qui est de l ordre de')
print('50 %% sur tau1 (12 a 32 s) et d un facteur 15 sur tau2 (108 a 4000 s).')
print('La gigue d horodatage est donc sans effet mesurable sur le spectre.')
