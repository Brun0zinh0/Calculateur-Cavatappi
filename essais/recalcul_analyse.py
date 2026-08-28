"""Recalcul de la section A de ANALYSE_PRELIMINAIRE.md avec la pression
reconstruite depuis les comptes ADC.

La pression est référencée au repos propre à chaque essai, sans écrêtage :

    P(bar) = [psi(ADC) - psi(ADC_repos)] * 0.0689476
    psi(ADC) = (ADC * 5 / 1023 - 0.5) * 500 / 4

Section A : seuils de dead-band des essais "1 O" (montée et descente).
Section B : vérifie que les essais "10 N" sont insensibles au référencement.

À lancer depuis le dossier "Espace de travail" :
    python essais/recalcul_analyse.py
"""
import csv, io, os
import numpy as np

PSI_TO_BAR = 0.0689476


def psi(adc):
    return (adc * 5.0 / 1023.0 - 0.5) * 500.0 / 4.0


def load(path):
    rows = list(csv.reader(io.StringIO(io.open(path, encoding='utf-8-sig').read())))
    idx = {k: i for i, k in enumerate(rows[0])}
    data = [r for r in rows[1:] if r and r[0]]
    g = lambda k, f=float: np.array([f(r[idx[k]]) for r in data])
    return g('time_s'), g('force_unfiltered_mN'), g('pressure_adc', int), g('pressure_bar')


def pressions(adc, P_old_bar):
    """Retourne (ancienne, reconstruite) en MPa, la zone morte en bar, le repos."""
    adc_rest = int(np.bincount(adc[:15]).argmax())
    P_cor = (psi(adc.astype(float)) - psi(adc_rest)) * PSI_TO_BAR * 0.1
    return P_old_bar * 0.1, P_cor, -psi(adc_rest) * PSI_TO_BAR, adc_rest


def seuil_montee(t, P, dF, lim):
    """Premiere pression ou dF depasse le seuil et s'y maintient 2 s."""
    for i in range(len(dF)):
        if dF[i] <= lim:
            continue
        j = np.searchsorted(t, t[i] + 2.0)
        if j > i and np.all(dF[i:j] > lim):
            return P[i]
    return np.nan


def seuil_descente(t, P, dF, lim):
    """Pression a laquelle la force retombe sous le seuil, branche descendante."""
    for i in range(len(dF)):
        if dF[i] > lim:
            continue
        j = np.searchsorted(t, t[i] + 2.0)
        if j > i and np.all(dF[i:j] <= lim):
            return P[i]
    return np.nan


base = 'essais'
print('=' * 78)
print('SECTION A - seuils de dead-band, essais "1 O"')
print('=' * 78)
print('%-8s %5s %6s | %7s %7s | %7s %7s | %s' %
      ('cas', 'repos', 'zone', 'Ps_anc', 'Ps_cor', 'Pd_cor', 'ecart', 'hysterese'))
print('-' * 78)

table, ecarts = {}, []
for muscle in sorted(os.listdir(base)):
    d = os.path.join(base, muscle)
    if not os.path.isdir(d):
        continue
    for fn in sorted(os.listdir(d)):
        if ' 1 O ' not in fn or not fn.endswith('.csv'):
            continue
        t, F, adc, P_old = load(os.path.join(d, fn))
        P_anc, P_cor, zone, rest = pressions(adc, P_old)

        imax = int(np.argmax(P_cor))
        au_repos = P_cor <= P_cor[0] + 1e-12
        F0 = F[au_repos].mean() if au_repos.sum() >= 3 else F[:3].mean()
        sig = F[au_repos].std() if au_repos.sum() >= 5 else F[:5].std()
        lim = max(10.0, 3.0 * sig)
        dF = F - F0

        Ps_a = seuil_montee(t[:imax + 1], P_anc[:imax + 1], dF[:imax + 1], lim)
        Ps_c = seuil_montee(t[:imax + 1], P_cor[:imax + 1], dF[:imax + 1], lim)
        Pd_c = seuil_descente(t[imax:], P_cor[imax:], dF[imax:], lim)

        hyst = '-' if np.isnan(Pd_c) else ('oui (%+.3f)' % (Pd_c - Ps_c) if Pd_c < Ps_c - 0.005 else 'quasi nulle')
        nom = fn.replace('Muscle ', '').replace('.csv', '').replace(' 1 O ', ' ')
        print('%-8s %5d %6.3f | %7.3f %7.3f | %7.3f %+6.3f | %s' %
              (nom, rest, zone, Ps_a, Ps_c, Pd_c, Ps_c - Ps_a, hyst))
        m, eps = nom.split()
        table.setdefault(m, {})[eps] = (Ps_c, Pd_c)
        ecarts.append(Ps_c - Ps_a)

print()
print('| Muscle | P_s a 0,8 | a 1,0 | a 1,2 | P_d a 0,8 | Hysterese |')
print('|---|---|---|---|---|---|')
for m in sorted(table):
    v = table[m]
    g = lambda e, i: ('%.3f' % v[e][i]) if e in v and not np.isnan(v[e][i]) else '-'
    d08 = v.get('0.8')
    h = 'oui (%+.3f)' % (d08[1] - d08[0]) if d08 and not np.isnan(d08[1]) and d08[1] < d08[0] - 0.005 else 'quasi nulle'
    print('| %s | %s | %s | %s | %s | %s |' % (m, g('0.8', 0), g('1.0', 0), g('1.2', 0), g('0.8', 1), h))

e = np.array(ecarts)
print()
print('effet du referencement : +%.3f a +%.3f MPa, mediane +%.3f' % (e.min(), e.max(), np.median(e)))

print()
print('=' * 78)
print('SECTION B - essais "10 N" : sensibilite au referencement')
print('=' * 78)
n, pa, pc = 0, [], []
for muscle in sorted(os.listdir(base)):
    d = os.path.join(base, muscle)
    if not os.path.isdir(d):
        continue
    for fn in sorted(os.listdir(d)):
        if ' 10 N ' not in fn or not fn.endswith('.csv'):
            continue
        t, F, adc, P_old = load(os.path.join(d, fn))
        P_anc, P_cor, zone, rest = pressions(adc, P_old)
        n += 1
        pa.append(np.abs(P_anc).max())
        pc.append(np.abs(P_cor).max())
print('%d essais examines' % n)
print('|P| max, echelle ancienne     : %.4f MPa' % max(pa))
print('|P| max, echelle reconstruite : %.4f MPa  (un compte ADC)' % max(pc))
print()
print('Les ajustements bi-exponentiels de la section B portent sur F(t) a P = 0')
print('et ne font intervenir la pression a aucun moment : ils sont inchanges.')
