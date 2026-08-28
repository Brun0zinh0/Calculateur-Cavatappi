"""Identification du spectre de Maxwell à partir des essais « 10 N ».

Un essai « 10 N » maintient la précontrainte à pression nulle pendant 600 s.
À déformation bloquée et en régime linéaire, la force est proportionnelle au
module de relaxation :

    F(t)/F(0) = E(t)/E(0) = [E0 + somme(Ei * exp(-t/tau_i))] / somme(E)

Un ajustement bi-exponentiel F(t) = Finf + A1 exp(-t/tau1) + A2 exp(-t/tau2)
donne donc directement les fractions du spectre :

    E0 / somme(E) = Finf / F(0)        Ei / somme(E) = Ai / F(0)
    eta_i = Ei * tau_i

La somme des modules est fixée à celle de l'article (37,76 MPa) : c'est la
seule grandeur qui ne provient pas des essais.

À lancer depuis le dossier "Espace de travail" :
    python essais/identifier_spectre.py
"""
import csv, io, os
import numpy as np
from scipy.optimize import curve_fit

E_SIGMA = 37.76  # MPa, somme des modules (article, annexe A)


def load(path):
    rows = list(csv.reader(io.StringIO(io.open(path, encoding='utf-8-sig').read())))
    idx = {k: i for i, k in enumerate(rows[0])}
    data = [r for r in rows[1:] if r and r[0]]
    g = lambda k, f=float: np.array([f(r[idx[k]]) for r in data])
    return g('time_s'), g('force_unfiltered_mN')


def biexp(t, Finf, A1, tau1, A2, tau2):
    return Finf + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)


def ajuste(t, F):
    F0 = F[:10].mean()
    p0 = [0.90 * F0, 0.04 * F0, 20.0, 0.05 * F0, 400.0]
    bornes = ([0.0, 0.0, 2.0, 0.0, 80.0],
              [2.0 * F0, 0.6 * F0, 120.0, 0.6 * F0, 4000.0])
    popt, _ = curve_fit(biexp, t, F, p0=p0, bounds=bornes, maxfev=40000)
    Finf, A1, tau1, A2, tau2 = popt
    if tau1 > tau2:  # on ordonne : branche rapide d'abord
        A1, tau1, A2, tau2 = A2, tau2, A1, tau1
    resid = F - biexp(t, Finf, A1, tau1, A2, tau2)
    r2 = 1.0 - resid.var() / F.var() if F.var() > 0 else np.nan
    return F0, Finf, A1, tau1, A2, tau2, r2


base = 'essais'
print('Ajustement bi-exponentiel des essais « 10 N » (600 s, P = 0)')
print()
print('%-8s %8s %7s | %7s %7s | %7s %7s | %6s %6s' %
      ('cas', 'F0 (mN)', 'relax%', 'f1 (%)', 'tau1', 'f2 (%)', 'tau2', 'E0/SE', 'R2'))
print('-' * 82)

lignes = []
for muscle in sorted(os.listdir(base)):
    d = os.path.join(base, muscle)
    if not os.path.isdir(d):
        continue
    for fn in sorted(os.listdir(d)):
        if ' 10 N ' not in fn or not fn.endswith('.csv'):
            continue
        t, F = load(os.path.join(d, fn))
        t = t - t[0]
        try:
            F0, Finf, A1, tau1, A2, tau2, r2 = ajuste(t, F)
        except Exception as exc:
            print('%-8s  echec : %s' % (fn, str(exc)[:40]))
            continue
        f1, f2, f0 = A1 / F0, A2 / F0, Finf / F0
        relax = 100.0 * (1.0 - Finf / F0)
        nom = fn.replace('Muscle ', '').replace('.csv', '').replace(' 10 N ', ' ')
        print('%-8s %8.0f %7.1f | %7.2f %7.1f | %7.2f %7.1f | %6.3f %6.3f' %
              (nom, F0, relax, 100 * f1, tau1, 100 * f2, tau2, f0, r2))
        lignes.append((f0, f1, tau1, f2, tau2, relax, r2, tau2 > 3900))

a = np.array([l[:7] for l in lignes], dtype=float)
sature = sum(1 for l in lignes if l[7])

print()
print('%d essais ajustes, %d avec tau2 en butee de borne' % (len(a), sature))
print()
print('                    mediane      etendue')
noms = ['E0/SE      ', 'f1 (rapide)', 'tau1 (s)   ', 'f2 (lente) ', 'tau2 (s)   ', 'relax (%)  ', 'R2         ']
for i, n in enumerate(noms):
    print('%s %10.3f   %.3f - %.3f' % (n, np.median(a[:, i]), a[:, i].min(), a[:, i].max()))

f0, f1, tau1, f2, tau2 = (np.median(a[:, i]) for i in range(5))
# renormalisation : les fractions doivent sommer a 1
s = f0 + f1 + f2
f0, f1, f2 = f0 / s, f1 / s, f2 / s

E0, E1, E2 = E_SIGMA * f0, E_SIGMA * f1, E_SIGMA * f2
eta1, eta2 = E1 * tau1, E2 * tau2

print()
print('=' * 60)
print('SPECTRE IDENTIFIE  (somme des modules fixee a %.2f MPa)' % E_SIGMA)
print('=' * 60)
print('  E0    = %8.3f MPa      (ressort permanent)' % E0)
print('  E1    = %8.3f MPa      eta1 = %9.1f MPa.s   tau1 = %6.1f s' % (E1, eta1, tau1))
print('  E2    = %8.3f MPa      eta2 = %9.1f MPa.s   tau2 = %6.1f s' % (E2, eta2, tau2))
print('  E3    = %8.3f MPa      (non contrainte par 600 s)' % 0.0)
print()
print('  E0/somme(E) = %.3f      (article : 0,168)' % f0)
print('  somme       = %.3f MPa' % (E0 + E1 + E2))
print()
print('A reporter dans les champs maxwell_* de l interface :')
print('  maxwell_E0_mpa      = %.3f' % E0)
print('  maxwell_E1_mpa      = %.3f     maxwell_eta1_mpa_s = %.1f' % (E1, eta1))
print('  maxwell_E2_mpa      = %.3f     maxwell_eta2_mpa_s = %.1f' % (E2, eta2))
print('  maxwell_E3_mpa      = 0.0       maxwell_eta3_mpa_s = 1.0')
import json
# chemin derive de __file__ : le repertoire courant peut etre mal decode selon
# le shell, alors que le chemin du script lui-meme est toujours exploitable.
sortie = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'spectre_identifie.json')
payload = {
        'source': 'essais 10 N, ajustement bi-exponentiel, mediane sur %d essais' % len(a),
        'E_sigma_mpa': E_SIGMA,
        'maxwell_E0_mpa': round(E0, 4),
        'maxwell_E1_mpa': round(E1, 4), 'maxwell_eta1_mpa_s': round(eta1, 2),
        'maxwell_E2_mpa': round(E2, 4), 'maxwell_eta2_mpa_s': round(eta2, 2),
        'maxwell_E3_mpa': 0.0, 'maxwell_eta3_mpa_s': 1.0,
        'tau1_s': round(tau1, 2), 'tau2_s': round(tau2, 2),
        'E0_sur_somme': round(f0, 4),
        'relaxation_600s_mediane_pct': round(float(np.median(a[:, 5])), 2),
}
try:
    with io.open(sortie, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    print('spectre ecrit dans %s' % sortie)
except OSError as exc:
    print('spectre non ecrit (%s) ; contenu :' % exc.__class__.__name__)
    print(json.dumps(payload, indent=2, ensure_ascii=False))

print()
print('Reserves :')
print(' - la force mesuree est celle du composite tube + nylon + extremites ;')
print('   le nylon etant elastique et portant une large part de la raideur,')
print('   ces fractions minorent la relaxation propre du PVC ;')
print(' - la derive de la cellule sur 600 s n est pas soustraite (essais a')
print('   blanc requis) : elle est comptee ici comme de la relaxation ;')
print(' - tau2 est mal contrainte par une fenetre de 600 s, et E3 ne l est pas')
print('   du tout : un essai de 30 a 60 min est necessaire pour la queue lente.')
