# Notes de version

## 2026.09.02 — Alpha V4 (moteur `2026.09.02-v4-15`)

- Six mécanismes physiques optionnels, off par défaut, moteur bit-identique
  à l'alpha V3 sinon : pression d'engagement (ovalité, variable d'état,
  hystérésis du seuil), frottement de Coulomb sur la pression motrice,
  convention de pré-étirement entre mors (solveur à deux inconnues par
  incrément), viscosité d'Eyring par couche, fluage d'ancrage logarithmique,
  outil `identification/identifier_spectre_moteur.py`.
- Réglages : schéma 17 (sept clés `prestretch_convention`,
  `engagement_*`, `friction_*`, `eyring_sigma_star_mpa`, `anchor_creep_*`) ;
  expander « Mécanismes Alpha V4 » dans l'interface ; rappel des mécanismes
  actifs sous les réglages.
- Sorties nouvelles : `pressure_effective_MPa`, `pressure_friction_MPa`,
  `ovality`, `anchor_creep_mm`,
  `prestretch_end_extension_mm`.
- Tests : 9 tests V4 ajoutés (35 au total), baseline figure 7 inchangée.

### Corrections de la contre-expertise adversariale du code (04/09/2026)

15 constats confirmés, tous traités : migration de schéma 16→17 non
destructive (les réglages V3 survivent) ; V4-3 en formulation totale
δ = C(α)·F (fonction d'état, plus de dépendance au chemin) et longueurs
entre mors exportées avec cet allongement ; V4-1 paramétrée par P_r0 seule
(e0 et k n'étaient pas identifiables séparément) ; V4-4 sur l'invariant de
von Mises, x borné ; garde-fous (fluage d'ancrage vs longueur active,
2·P_c vs P_max, σ* ≥ 1e-6, fluage inactif en suspendu signalé) ; règle du
nylon par phase ; outil d'identification : vitesse de pré-étirement 20 mm/min
par défaut et tracée, CSV robuste (délimiteur, virgule décimale), fenêtre t0.
