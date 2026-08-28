"""Spectre du PVC seul, déconvolué de la mesure composite des essais « 10 N ».

Un essai « 10 N » mesure la force du muscle monté, c'est-à-dire la somme de la
contribution du tube de PVC et de celle du filament de nylon :

    F(t) = F_tube(t) + F_nylon

Le nylon est élastique et ne relaxe pas : F_nylon est constant. Toute la
relaxation observée provient donc du tube, mais elle est rapportée à une force
totale plus grande que celle du tube. La fraction relaxante mesurée est ainsi
diluée d'un facteur

    d = F_total / F_tube        (>= 1)

et les fractions du spectre propre au PVC valent f_i^PVC = d * f_i^mesure.
Les temps caractéristiques, eux, ne sont pas affectés par la dilution.

Le partage F_tube / F_nylon est fourni par le moteur, qui projette séparément
les deux contributions sur la géométrie hélicoïdale (colonnes force_tube_mN et
force_nylon_mN).

À lancer depuis le dossier "Espace de travail" :
    python essais/spectre_pvc_seul.py
"""
import importlib.util, io, json, os, sys
import numpy as np

V3 = os.path.abspath('alpha V3')
sys.path.insert(0, V3)
spec = importlib.util.spec_from_file_location('_base', os.path.join(V3, 'Base.py'))
Base = importlib.util.module_from_spec(spec)
sys.modules['_base'] = Base
spec.loader.exec_module(Base)

ICI = os.path.dirname(os.path.abspath(__file__))
SPECTRE = json.load(io.open(os.path.join(ICI, 'spectre_identifie.json'), encoding='utf-8'))
FICHES = os.path.join(V3, 'identification', 'fiches')


def partage(fiche):
    """Retourne (F_tube, F_nylon) en mN a l'etat precontraint, pression nulle."""
    mx = Base.default_maxwell_tensile_params(
        E0=SPECTRE['maxwell_E0_mpa'],
        E1=SPECTRE['maxwell_E1_mpa'], eta1=SPECTRE['maxwell_eta1_mpa_s'],
        E2=SPECTRE['maxwell_E2_mpa'], eta2=SPECTRE['maxwell_eta2_mpa_s'],
        E3=SPECTRE['maxwell_E3_mpa'], eta3=SPECTRE['maxwell_eta3_mpa_s'])
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
        n_layers=3, n_phi=8, pre_steps=24, dt=0.5)
    t = np.array([0.0, 1.0])
    p = np.array([0.0, 0.0])
    _, arr = Base.run_blocked_actuation(cfg, pressure_time=t, pressure_MPa=p)
    return float(arr['force_tube_mN'][0]), float(arr['force_nylon_mN'][0])


print('Partage de la force de precontrainte entre tube et nylon')
print()
print('%-10s %10s %10s %8s %8s' % ('fiche', 'F_tube', 'F_nylon', 'part_PVC', 'dilution'))
print('-' * 50)
dilutions = []
for fn in sorted(os.listdir(FICHES)):
    if not fn.endswith('.json'):
        continue
    fiche = json.load(io.open(os.path.join(FICHES, fn), encoding='utf-8'))
    try:
        Ft, Fn = partage(fiche)
    except Exception as exc:
        print('%-10s  echec : %s' % (fn, str(exc)[:40]))
        continue
    tot = Ft + Fn
    if abs(tot) < 1e-9 or Ft <= 0:
        print('%-10s  partage inexploitable (F_tube = %.1f)' % (fn, Ft))
        continue
    d = tot / Ft
    print('%-10s %10.1f %10.1f %7.1f %% %8.2f' %
          (fn.replace('muscle_', '').replace('.json', ''), Ft, Fn, 100 * Ft / tot, d))
    dilutions.append(d)

d = float(np.median(dilutions))
print()
print('facteur de dilution median : %.2f  (etendue %.2f - %.2f)' % (d, min(dilutions), max(dilutions)))

# --- application au spectre -------------------------------------------------
ES = SPECTRE['E_sigma_mpa']
f1 = SPECTRE['maxwell_E1_mpa'] / ES
f2 = SPECTRE['maxwell_E2_mpa'] / ES
f1p, f2p = d * f1, d * f2
f0p = 1.0 - f1p - f2p

print()
print('=' * 62)
print('SPECTRE DU PVC SEUL  (deconvolue, somme des modules %.2f MPa)' % ES)
print('=' * 62)
print('%-22s %12s %12s' % ('', 'composite', 'PVC seul'))
print('%-22s %11.3f %12.3f' % ('E0 (MPa)', ES * (1 - f1 - f2), ES * f0p))
print('%-22s %11.3f %12.3f' % ('E1 (MPa)', ES * f1, ES * f1p))
print('%-22s %11.3f %12.3f' % ('E2 (MPa)', ES * f2, ES * f2p))
print('%-22s %11.1f %12.1f' % ('eta1 (MPa.s)', ES * f1 * SPECTRE['tau1_s'], ES * f1p * SPECTRE['tau1_s']))
print('%-22s %11.1f %12.1f' % ('eta2 (MPa.s)', ES * f2 * SPECTRE['tau2_s'], ES * f2p * SPECTRE['tau2_s']))
print('%-22s %11.3f %12.3f' % ('E0 / somme(E)', 1 - f1 - f2, f0p))
print()
print('tau1 = %.1f s et tau2 = %.1f s : inchanges, la dilution ne porte que'
      % (SPECTRE['tau1_s'], SPECTRE['tau2_s']))
print('sur les amplitudes.')
print()
print('Rappel : article, E0/somme(E) = 0,168.')

payload = {
    'source': 'deconvolution du spectre composite par le partage tube/nylon du moteur',
    'facteur_dilution_median': round(d, 4),
    'E_sigma_mpa': ES,
    'maxwell_E0_mpa': round(ES * f0p, 4),
    'maxwell_E1_mpa': round(ES * f1p, 4), 'maxwell_eta1_mpa_s': round(ES * f1p * SPECTRE['tau1_s'], 2),
    'maxwell_E2_mpa': round(ES * f2p, 4), 'maxwell_eta2_mpa_s': round(ES * f2p * SPECTRE['tau2_s'], 2),
    'maxwell_E3_mpa': 0.0, 'maxwell_eta3_mpa_s': 1.0,
    'tau1_s': SPECTRE['tau1_s'], 'tau2_s': SPECTRE['tau2_s'],
    'E0_sur_somme': round(f0p, 4),
}
sortie = os.path.join(ICI, 'spectre_pvc_seul.json')
try:
    with io.open(sortie, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    print('\nspectre ecrit dans %s' % sortie)
except OSError:
    print('\nspectre non ecrit ; contenu :')
    print(json.dumps(payload, indent=2, ensure_ascii=False))

print()
print('Limites de cette deconvolution :')
print(' - le partage tube/nylon vient du modele, non de la mesure : il herite')
print('   de ses hypotheses (nylon lie bilateralement, geometrie ideale) ;')
print(' - la derive de la cellule reste comptee comme de la relaxation, et se')
print('   trouve amplifiee par la dilution ;')
print(' - seul un essai de relaxation sur tube nu (E2) leve ces deux reserves.')
