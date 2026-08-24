# Validation figure 7 — alpha V3

Valide le moteur `alpha V3/Base.py` contre la figure 7 de l'article « Blocked
actuation » (numérisation de l'image embarquée du PDF, pages 180-380 s).

## Usage

```powershell
python validation_figure7.py             # run complet (6 couches, ~5-10 min)
python validation_figure7.py --quick     # discrétisation réduite (~2 min)
python validation_figure7.py --show      # affiche les figures à l'écran
python validation_figure7.py --baseline  # (ré)génère la baseline de non-régression
```

Sorties : figures PNG et métriques JSON dans `out/`.

## Corrections issues de l'audit 2026-08 (plan, item 0.1)

- moteur alpha V3 importé directement (l'ancien script de la racine validait le
  moteur du dossier `livrable`, version 2026.07.17) ;
- chemin du PDF résolu vers `../../Article support/` (l'ancien chemin était
  cassé) ; la variable d'environnement `TCPA_FIGURE7_PDF` permet un chemin
  explicite ;
- calibration de la numérisation corrigée : la conversion pixel→pression est
  ajustée sur les positions réelles des étiquettes d'axe (l'ancienne conversion
  sous-lisait les pics : 1,32-1,37 MPa au lieu de 1,41-1,43 — ~40 % du déficit
  apparent de couple) ;
- mode `rescale` par défaut : les pics de pression numérisés sont recalés sur
  la valeur du texte de l'article (1,42 MPa).

## Baseline de non-régression

`--baseline` écrit `../audit/baseline/baseline_figure7.json` (pics/vallées de
force et de couple à deux discrétisations, réduite et production). Le test
`test_validation_figure7_non_regression` de `../test_scientifique.py` compare
le moteur courant à cette baseline (tolérance ±2 %) à chaque exécution de la
suite ; `TCPA_SLOW_VALIDATION=1` ajoute la discrétisation de production,
`TCPA_SKIP_VALIDATION=1` saute le test.

Ne régénérer la baseline (`--baseline`) qu'après une modification **assumée**
de la physique, en notant la raison dans le plan de correction (item 0.3).

## Lecture des résultats

Les écarts résiduels connus au moment du gel de la baseline (audit 2026-08) :
force −4 à −8 % (ε=0,8) et −11 à −13 % (ε=1,0) vs théorie de l'article, couple
−10 % (ε=0,8) et −25 % (ε=1,0). Pistes documentées dans
`../PLAN_DE_CORRECTION.md` (items 3.1/3.2). La cible « théorie » numérisée est
tronquée aux pics de couple (occlusion par les marqueurs expérimentaux) : les
RMSE sont des ordres de grandeur, pas des mesures de précision.
