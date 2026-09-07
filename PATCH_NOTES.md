# Notes de version

## 2026.09.07 — Vitesse de pression en MPa/s (moteur `2026.09.07-v4-16`)

- Le profil de pression généré est défini par la vitesse de pression des
  rampes (`pressure_rate_mpa_s`, MPa/s ; demi-cycle = Pmax / vitesse) et non
  plus par le couple débit (mL/min) / volume (mL) de seringue, qui n'est pas
  une consigne du banc (seringue manuelle). Valeur par défaut Pmax / 9 s :
  profil de l'article et figure 7 inchangés.
- Réglages : schéma 18 ; migration non destructive débit/volume → vitesse
  équivalente p_max·Q/(60·V), fichiers et exports ; clés historiques encore
  honorées (prioritaires) par `build_config` et `Base.cyclic_pressure_history`.
- Interface : champ « Vitesse de pression (MPa/s) », demi-cycle et cycle
  affichés, durée estimée recalculée, vitesse effective affichée en durée
  fixe. API `parametres` : `effective_pressure_rate_mpa_s` remplace
  `effective_flow_rate_mL_min`.
- Revue (07/09) : priorité unique des clés historiques débit/volume sur tous
  les chemins, demi-période exacte 60·V/Q transportée par
  `SimulationParams.half_period_s` (bit-identité même à p_max = 0 ou pour une
  demi-période non représentable), contradiction vitesse explicite / clés
  historiques refusée, écrêtage aux bornes signalé, vitesse contrôlée par
  `settings_error` (message, pas d'exception) ; scripts d'audit
  `contre_pretension`, `diag_31_instrument`, `audit_stiffness_rotation`,
  `ce_minors_check2` remis en cohérence (demi-période historique conservée) ;
  `audit/code_map.md` mis à jour.
- Sémantique : la vitesse étant constante, la durée d'un demi-cycle suit
  désormais Pmax (Pmax / vitesse) au lieu d'être fixée à 9 s. À la pression
  maximale par défaut le profil est bit-identique à l'alpha V3 ; à une autre
  Pmax, retrouver exactement l'ancien profil demande vitesse = Pmax / 9 s.
- Tests : 2 tests ajoutés (37 au total), baseline figure 7 inchangée.

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
