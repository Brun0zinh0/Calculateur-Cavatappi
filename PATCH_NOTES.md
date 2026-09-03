# Notes de version

## 2026.09.02 — Alpha V4 (moteur `2026.09.02-v4-15`)

- Six mécanismes physiques optionnels, off par défaut, moteur bit-identique
  à l'alpha V3 sinon : pression d'engagement (ovalité, variable d'état,
  hystérésis du seuil), frottement de Coulomb sur la pression motrice,
  convention de pré-étirement entre mors (solveur à deux inconnues par
  incrément), viscosité d'Eyring par couche, fluage d'ancrage logarithmique,
  outil `identification/identifier_spectre_moteur.py`.
- Réglages : schéma 17 (neuf clés `prestretch_convention`,
  `engagement_*`, `friction_*`, `eyring_sigma_star_mpa`, `anchor_creep_*`) ;
  expander « Mécanismes Alpha V4 » dans l'interface ; rappel des mécanismes
  actifs sous les réglages.
- Sorties nouvelles : `pressure_effective_MPa`, `pressure_friction_MPa`,
  `ovality`, `anchor_creep_mm`,
  `prestretch_end_extension_mm`.
- Tests : 7 tests V4 ajoutés (33 au total), baseline figure 7 inchangée.
