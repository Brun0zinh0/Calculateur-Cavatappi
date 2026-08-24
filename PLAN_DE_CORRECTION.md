# Plan de correction — alpha V3

> **🏁 PLAN EXÉCUTÉ EN TOTALITÉ le 21/08/2026** (moteur final
> `2026.08.21-audit-phase4-14`, suite de 25 tests, baseline de non-régression
> intacte de bout en bout). Seule reste ouverte l'identification matériau
> finale (item 3.3), conditionnée aux essais E1-E4 de l'équipe — l'outil est
> prêt. Voir la section « Suites de l'audit » du rapport
> (`audit/rapport_audit.html`) pour le bilan complet.

Issu de l'audit complet du 17-18/08/2026 (rapport : `audit/rapport_audit.html` ; constats bruts : `audit/*.json` ; matrice figure 7 : `audit/matrix_summary.md`). Version de départ : copie conforme d'`alpha V2` (moteur `Base.py` v2026.07.30-series-compliant-ends-10), **aucune modification appliquée** — ce dossier est la ligne de base, le présent document est la feuille de route. *(Si ce dossier ne contient que le présent fichier, la copie de dépôt — commande robocopy fournie en séance — n'a pas encore été exécutée : la faire avant la Phase 0. Ce plan a lui-même été contre-vérifié par un agent indépendant : 16/16 références de lignes et l'ensemble des chiffres confirmés contre les fichiers d'audit.)*

Références : **M1-M8** = constats majeurs, **[m]** = constats mineurs (détail dans le rapport). Effort : **S** < 1 h, **M** = ½ journée, **L** = plusieurs jours. Tous les numéros de ligne se réfèrent à `Base.py` v2026.07.30 et ont été vérifiés en contre-expertise.

**Principe directeur : ne toucher à aucune physique avant d'avoir un filet de sécurité exécutable (Phase 0).** Chaque correction est petite, isolée, et porte un critère d'acceptation chiffré.

---

## Décisions à prendre avant la Phase 2

**✅ Décisions prises le 21/08/2026 (recommandations de l'audit retenues)** :
D1 = garder `paper_linear` par défaut et documenter l'écart à l'éq. (11) ;
D2 = dossier `livrable` déprécié (notice `livrable/AVERTISSEMENT_DEPRECIE.md`
déposée), régénération après phases 0-2 ; D3 = exposer les deux normalisations
avec L_T0 par défaut (à implémenter avec 1.4/4.1).
**✅ D3 IMPLÉMENTÉE le 21/08/2026 (Phase 4)** : `free_actuation_percent` est
désormais normalisé par la longueur non chargée L_T0 (éq. 25 d'EXP), capturée
après précontrainte et avant application de la masse ; l'ancienne
normalisation reste exportée en `_loaded_ref` ; sortie
`reference_unloaded_length_mm` ajoutée ; test `test_suspended_normalization_LT0`.

| # | Décision | Options | Recommandation de l'audit |
|---|----------|---------|---------------------------|
| D1 | Défaut du profil d'angle de biais (M1) | garder `paper_linear` / basculer sur `uniform_twist` (éq. 11) | Garder `paper_linear` pour la reproduction figure 7 (la matrice montre un effet marginal et de signe mixte : force +0,5 pt, couple −4 pts) **mais** documenter honnêtement l'écart à l'équation (item 1.1). Rebasculer éventuellement après 3.1. |
| D2 | Sort du dossier `livrable` | déprécier avec notice / régénérer depuis alpha V3 / supprimer | Déprécier avec notice tant que 0.1-0.2 ne sont pas faits, puis régénérer. Il contient un moteur plus ancien (2026.07.17) — source de confusion documentée par l'audit. |
| D3 | Normalisation de l'actionnement en mode suspendu | garder L(0) chargée / passer à L_T0 (éq. 25 EXP) / exposer les deux | Exposer les deux avec L_T0 par défaut pour toute comparaison aux figures d'EXP (l'écart est un facteur ≈ L(0)/L_T0 : 2,11 % vs 3,71 % sur le cas test de l'audit). |

---

## Phase 0 — Filet de sécurité (à faire en premier, avant toute autre modification)

> **✅ PHASE 0 EXÉCUTÉE le 21/08/2026.**
> 0.1 : `validation/validation_figure7.py` porté sur le moteur alpha V3 (chemin
> PDF corrigé + variable `TCPA_FIGURE7_PDF`, calibration des ticks, recalage
> des pics sur 1,42 MPa, sorties dans `validation/out/`) — critère
> d'acceptation atteint : pics production 1705,44 / 2145,01 mN et
> 595,88 / 590,91 µN·m vs cibles matrice 1705/2145 et 596/591 (< 0,05 %).
> 0.2 : `test_validation_figure7_non_regression` ajouté à la suite (baseline
> ±2 % sur pics/vallées force, couple, pression ; `TCPA_SLOW_VALIDATION=1`
> pour la discrétisation de production ; sauts propres si PDF/baseline
> absents). Suite complète verte : 15/15.
> 0.3 : baseline gelée dans `audit/baseline/baseline_figure7.json`
> (discrétisations réduite et production, empreintes de pression, ancrages
> texte). Écarts résiduels au gel (connus de l'audit, items 3.1/3.2) :
> force −8,1 %/−12,9 %, couple −32,9 %/−33,2 % vs ancrages texte (réduite).

### 0.1 Porter la validation figure 7 vers alpha V3 — [M8] — effort M
`validation/validation_figure7.py` (copié tel quel depuis la racine) ne peut pas s'exécuter contre ce moteur. Trois corrections :
1. **Chemin du PDF** (l.37-38) : pointer vers `..\..\Article support\Blocked actuation....pdf` (le chemin actuel suppose le PDF à la racine — FileNotFoundError reproduit par l'audit).
2. **Calibration de la numérisation** (`PANEL_SPECS`, l.58-90) : les valeurs de ticks sont supposées aux bords des crops alors que les ticks réels sont en retrait — ce biais sous-lisait la pression aux pics (1,32-1,37 MPa au lieu de 1,41-1,43) et expliquait **~40 % du déficit apparent de couple** à ε=0,8. La correction est prototypée dans `audit/scripts/ce_fig7_measure2.py` et `ce_fig7_rerun.py`.
3. **Imports** : remplacer `livrable_parametres` par le `parametres` d'alpha V3 et retirer les clés `pressure_end_force_*` (ces clés n'existent pas dans ce moteur).
4. **Recalage des pics de pression** : après la calibration des ticks (pics ≈ 1,397 MPa), recaler la pression numérisée sur les pics du texte de l'article (1,41-1,43 MPa) — sans ce recalage, les cibles d'acceptation ci-dessous ne sont pas atteignables (variante « corr » seule : 1692/2129 mN et 583/577 µNm, ~2 % sous les cibles de couple).

**Acceptation** : le script s'exécute depuis `alpha V3/validation/` et reproduit à ±1 % les chiffres de la matrice d'audit (config lin/fixed, pression **recalée**, 6 couches : pics force ≈ 1705 / 2145 mN, pics couple ≈ 596 / 591 µNm à ε = 0,8 / 1,0).

### 0.2 Test de non-régression quantitatif — [M8] — effort M
Ajouter à `test_scientifique.py` une classe « validation » qui :
- exécute le protocole figure 7 en discrétisation réduite et compare les pics force/couple à une **baseline enregistrée** (tolérance ±2 %) ;
- affiche (à titre informatif, sans assertion) l'écart aux ancrages du texte de l'article : 1824,60 / 2424,67 mN et 877,01 / 873,41 µNm ;
- comporte un test « lent » optionnel aux réglages de production (dt = 0,5, n_layers = 6, 11 cycles), jamais exercés par la suite actuelle.

**Motivation** : les 14 tests actuels ne vérifient que la cohérence interne (3 assertions tautologiques identifiées) ; aucune chaîne automatisée ne protège le moteur contre une régression physique.

### 0.3 Geler la baseline — effort S
Avant toute modification de code : exécuter 0.1 + 0.2 et archiver les sorties de référence dans `audit/baseline/` (npz + chiffres). Toute correction ultérieure se mesure contre cette baseline.

---

## Phase 1 — Documentation (aucun risque de régression)

> **✅ PHASE 1 EXÉCUTÉE le 21/08/2026.** Les cinq items sont appliqués, chaque
> passage corrigé étant marqué « (audit 2026-08) » dans le README :
> 1.1 profil de biais (texte vs équation 11, option `uniform_twist`, D1) ;
> 1.2 section figée par défaut + correction de « Transmission de la pression »
> (le défaut ne réoriente pas les directions matérielles) + avertissement sur
> la dérive du mode `updated` ; 1.3 précontrainte à branches vierges avec les
> conséquences chiffrées (−18/−22 %, training ~17 % vs ~2/3, dead-band et
> fig. 11 hors de portée) ; 1.4 encart « reproduire le protocole EXP »
> (README § masse suspendue + expander de l'interface) ; 1.5 commentaire
> 31,24/37,76 dans `stiffness_sanity_report`. Le README porte aussi la
> nouvelle identité alpha V3 (contenu du dossier, validation/, audit/, plan).
> Compilation vérifiée, suite de tests complète relancée après édition.

| Item | Correction | Constat | Effort |
|------|-----------|---------|--------|
| 1.1 | `README.md` l.158 : écrire que `paper_linear` suit la **phrase descriptive** de la section 3.3 de l'article mais **pas son équation (11)** (forme arctan, disponible via l'option `uniform_twist`). L'article se contredit lui-même texte/équation — le dire aussi. | M1 | S |
| 1.2 | `README.md` l.180-182 : le défaut (`section_update_mode="fixed"`) **ne réoriente pas** les directions matérielles ni la section — corriger la phrase, documenter l'option `updated` et son état réel (non conforme sous cyclage long tant que 3.1 n'est pas résolu). | M2 | S |
| 1.3 | `README.md` § précontrainte : expliciter l'écart au schéma de l'article (branches de Maxwell vierges à t_k) **avec ses conséquences chiffrées** : −18/−22 % sur les niveaux absolus de force ; sur l'atténuation de « training » de −10,3 % mesurée par l'article, le schéma actuel n'en reproduit que ~17 % (−1,8 %) là où le schéma viscoélastique (3.2) en récupère ~2/3 (−6,3 %) ; dead-band par fluage et Fig. 11 d'EXP hors de portée. Renvoyer vers l'option 3.2. | M3 | S |
| 1.4 | README + interface : encart « reproduire le protocole EXP en mode suspendu » : `eps = 0`, équilibrage préalable **désactivé**, normalisation par L_T0 (selon D3). Aujourd'hui cette combinaison n'est indiquée nulle part. | [m] | S |
| 1.5 | Commentaire dans `stiffness_sanity_report` (l.1994-2010) : 31,24 = valeur mesurée Table B1 ; 37,76 = ΣE_k annexe A ; les métriques du rapport sont identiques dans les deux conventions (vérifié). | [m] | S |

---

## Phase 2 — Corrections localisées à faible risque

> **✅ PHASE 2 EXÉCUTÉE le 21/08/2026** (moteur `2026.08.21-audit-phase2-11`).
> 2.1 inférence d'unité par mots entiers ; 2.2 θ_f = 0 rejeté (message
> explicite) ; 2.3 référence série verrouillée sur P = 0 post-précontrainte
> (test : référence identique que l'historique démarre à 0 ou P₀ > 0, et
> extension du premier saut conservée) ; 2.4 profil linéaire par défaut sur
> tous les points d'entrée, `run_hold_relaxation` transmet `nonlinear_ramp`,
> γ centralisés dans `Base.NONLINEAR_GAMMA_*` ; 2.5 compliance des extrémités
> réévaluée à α_tk au verrou ; 2.6/2.11 conventions documentées (README
> § Extrémités désenroulées) ; 2.7 `axial_length_geometry_mm` = NaN en mode
> suspendu ; 2.8 avertissement multi-racines + garde NaN dans `_find_dw`,
> résidu revérifié dans `step()` (RuntimeError > 1e-4 N·mm) ; 2.9 mise en
> charge de la masse en rampe (8 incréments par défaut, convergence testée
> ±1 % à N×4) ; 2.10 = D1 appliquée (défaut conservé, rien à changer).
> **Vérification : suite 20/20 (5 nouveaux tests) et baseline figure 7
> régénérée avec le moteur modifié strictement identique à l'ancienne
> (écart relatif max 0,00) — zéro régression.**

Chaque item : modification + test. Ordre libre, sauf 2.10 (dépend de D1).

| Item | Correction (fichier : lignes) | Critère d'acceptation | Constat | Effort |
|------|-------------------------------|------------------------|---------|--------|
| 2.1 | `pression.py` l.58-66 : l'inférence d'unité de force par sous-chaîne « mn » piège les en-têtes génériques (`column1` → mN). Matcher des mots entiers sur en-tête normalisé. | `column1` → aucune inférence ; `Force (mN)` → mN ; `force (N)` → N | [m] | S |
| 2.2 | `Base.py` l.434-435 + 576-584 : θ_f = 0 accepté par la validation puis B = 0/0 = NaN dans `_layer_AB_mu` (cas résonant μ=1 de Lekhnitskii non implémenté). Minimal : rejeter θ_f = 0 avec un message clair. | Erreur explicite au lieu de `RuntimeError: residual=nan` | [m] | S |
| 2.3 | `Base.py` l.1696-1698 : le verrou de la référence série est posé **après** le pas initial à P(0) — si un CSV utilisateur démarre à P(0) > 0, la référence est contaminée et l'extension du premier saut est perdue. Verrouiller sur l'état P = 0 post-précontrainte. | Extension identique que l'historique commence à 0 ou à P₀ > 0 | [m] | S |
| 2.4 | Harmoniser les défauts des historiques de pression : `Base.py` l.1525 et 1560 (`nonlinear=True`) vs `default_simulation_config` l.111 et `parametres.py` l.150 (`False`) ; `run_hold_relaxation` l.1824 ne transmet pas `nonlinear_ramp` (rampe γ=3,5 imposée, ~24 % d'écart d'amplitude de relaxation vs le chemin interface) ; γ = 3,5/2,8 dupliqués en dur (`parametres.py` l.876-887). Défaut **linéaire partout** (conforme au README), γ centralisés. | API `Base.run_hold_relaxation` et chemin interface produisent la même sortie | [m] | M |
| 2.5 | `Base.py` l.400 + 468-482 : compliance des extrémités calculée une fois à α₀ et figée. La recalculer au verrou de référence (α_tk), là où elle sert. | Écart attendu ≈ 7,7 % vs version figée (chiffre de l'audit) ; nul si eps = 0 | [m] | S |
| 2.6 | `Base.py` l.1432-1441 : longueur des extrémités comptée pleinement axiale alors que la raideur les modélise quasi transversales (`tangent_beam`). Choisir une convention (« longueur entre mors » assumée, ou projection tangente) et la documenter. | Convention écrite au README ; sorties de longueur cohérentes avec elle | [m] | S |
| 2.7 | `Base.py` l.1425-1440 : `axial_length_geometry_mm` (2πN·h, N figé) est invalide en mode suspendu (écart 1,32 mm ≈ 2,3 % mesuré) mais exportée sans distinction. NaN/omission en suspendu, ou recalcul avec N libre. | Export cohérent dans les deux modes | [m] | S |
| 2.8 | `Base.py` l.973-1006 (`_find_dw`) : avertissement en cas de racines multiples (actuellement `min(roots, key=abs)` silencieux), garde NaN explicite (le comportement propre actuel dépend de scipy ≥ 1.12), revérification du résidu après `step()` (aujourd'hui seulement journalisé). | Log/exception explicites sur les 3 chemins ; aucun changement de résultat en régime nominal | [m] | M |
| 2.9 | `Base.py` l.1710 + 1783 : mise en charge de la masse en un seul incrément dt = 0 (+55 % de longueur d'un coup, hors hypothèse des incréments mineurs). Rampe en N pas (paramètre, défaut ≈ 10). | Résultat stable à ±1 % quand N double | [m] | S |
| 2.10 | Appliquer D1 (défaut `bias_angle_profile`). Si bascule vers `uniform_twist` : mettre à jour la baseline 0.3 (effet attendu de la bascule : **~−5 à −6 % sur le couple crête, ~+7 % sur l'amplitude de force d'actionnement**). | Baseline mise à jour et tracée | M1 | S |
| 2.11 | Documenter (ou corriger) la convention de précontrainte vis-à-vis des extrémités : `prestretch_to` (Base.py l.1361-1379) applique eps à la spire seule, sans compatibilité série — la correspondance eps ↔ ε_tk expérimental n'est exacte que pour des extrémités rigides. Convention absente du README. Minimal : l'écrire au README ; optionnel : offrir eps « série complète ». | Convention documentée ; nul par défaut (uncoiled_length = 0) | [m] | S |

---

## Phase 3 — Chantiers physiques majeurs

### 3.1 Comprendre et corriger la dérive du mode `updated` sous cyclage long — [M2, M4, M5] — effort L

> **✅ ITEM 3.1 EXÉCUTÉ le 21/08/2026** (moteur `2026.08.21-audit-phase31-12`).
> **Diagnostic** (instrumentation en 3 volets) : la dérive n'était pas du
> fluage viscoélastique — un run quasi élastique (η×10⁹) dérivait autant que
> le Maxwell complet (Rin −6,5 %/11 cycles, taux constant) — mais le **biais
> du premier ordre de l'intégrateur géométrique explicite** : chaque pas de
> décharge est évalué sur une configuration plus dilatée (plus complaisante)
> que son homologue de charge ; erreur de fermeture de cycle vérifiée
> strictement **O(h)** (÷2 quand dt ÷2), ratchet net −1,33·10⁻³ mm/cycle sur
> Rin porté par le terme homogène du BVP radial.
> **Correctif** : schéma de **point milieu complet** — le BVP radial est
> re-résolu sur la configuration à mi-incrément ET le champ u est évalué aux
> positions à mi-incrément (géométrie, angles matériels, déformations), dans
> les deux chemins (bloqué et suspendu). Actif uniquement en mode `updated`
> (`fixed` bit-identique, baseline vérifiée strictement inchangée, 0,00).
> **Résultats** : fermeture de cycle quasi élastique ×28 (aire +2,8·10⁻⁴/cycle,
> soit ~0,01 %) ; 11 cycles Maxwell stables ; figure 7 en `updated` corrigé :
> force −0,7/+0,7 % (ε=0,8) et −8/−9 % (ε=1,0) vs théorie, couple **−6 à
> −12 %** des ancrages expérimentaux (contre −33 % en `fixed` et −60/−85 %
> en `updated` non corrigé) — mieux que la courbe théorique publiée, laquelle
> sous-estime ses propres ancrages (~−20 % à ε=0,8), vraisemblablement par le
> même biais d'intégration. **Décision de défaut** : `fixed` conservé
> (référence de baseline) ; `updated` corrigé recommandé pour les études de
> fidélité. Test permanent ajouté : `test_updated_section_cycle_closure`.
> **Résidu attribué** : la force à ε=1,0 (−8/−9 %) reste le candidat de
> l'item 3.2 (précontrainte). Dérive résiduelle O(h) (~0,01 %/cycle)
> documentée au README. Scripts de diagnostic : `diag_31*.py` (archivés dans
> audit/scripts après dépôt).
La matrice d'audit a montré que `section_update_mode="updated"` (conforme à la lettre de l'article) **s'effondre sous le protocole figure 7 complet** (couple −60/−85 % vs théorie article) alors qu'il améliore l'amplitude sur protocole court (+55 %). C'est très probablement un défaut d'implémentation de la réactualisation (composition des incréments r = R + Δu au fil des cycles), pas une fatalité.
1. Instrumenter : tracer `R_edges`, `theta_layers`, C̄ par couche cycle après cycle (dérive monotone ? oscillation ?).
2. Comparer à la composition attendue du lagrangien réactualisé de l'article (éq. 16-17).
3. Corriger, puis réévaluer : la configuration conforme referme-t-elle une partie du résidu figure 7 (force −4/−12 %, couple −10/−25 %) ?

**Acceptation** : `updated` stable sur 11 cycles (pas d'effondrement du couple) ; décision de défaut documentée ; baseline mise à jour.

### 3.2 Option de précontrainte viscoélastique (schéma de l'article) — [M3, M4] — effort L

> **✅ ITEM 3.2 EXÉCUTÉ le 21/08/2026** (moteur `2026.08.21-audit-phase32-13`).
> **Implémentation** : mode `prestrain_reference_mode = "viscoelastic_history"`
> (branches de Maxwell actives dès l'élongation à 20 mm/min ; aucune référence
> élastique conservée), exposé dans l'interface (sélecteur « Précontrainte »),
> avec coercition des anciens modes retirés et test permanent
> `test_viscoelastic_prestress_mode` (σ_réf = 0 et branches chargées à t_k,
> force à t_k ~0,8× l'élastique, relaxation de prétension > 1 % à P = 0).
> Suite 22/22, baseline du défaut strictement identique (0,00).
> **Acceptation atteinte** : training reproduit en ordre de grandeur
> (atténuation des pics c1→c11 = −5,0 % contre −10,3 % article et −0,9 %
> élastique) ; fig. 11 d'EXP simulable (fluage sous poids du bon signe :
> +0,19 mm/60 s, là où le mode élastique donne un signe opposé).
> **Paradoxe arbitré** (figure 7, pics ET vallées, avec et sans section
> réactualisée) : les niveaux publiés par l'article coïncident avec le mode
> élastique conservé (vallée théorie ~1450 mN vs 1448 simulé à ε=0,8) — la
> simulation des auteurs n'a vraisemblablement pas relaxé la prétension. Le
> mode viscoélastique amène le couple à −2/−1 % des ancrages expérimentaux
> mais abaisse la force de ~22 % (−25/−30 %), conséquence directe du spectre
> Table A1 qui sur-relaxe (M6). **Décision : défaut inchangé
> (`elastic_tk_reference`)** ; `viscoelastic_history` recommandé pour
> training/dead-band/relaxation de prétension et pour tout spectre identifié
> sur le matériau réel (3.3). Le résidu de force à ε=1,0 n'est donc pas
> refermable par le seul schéma de précontrainte avec ce spectre : il est
> couplé à l'identification matériau (3.3). Étude archivée :
> `audit/scripts/diag_32_studies.py`.
Ajouter un mode `prestress_mode="viscoelastic_history"` à côté d'`elastic_tk_reference` (actuellement verrouillé seul, l.327-330) : rampe d'élongation à vitesse finie (20 mm/min du protocole) avec branches de Maxwell **actives**. Le prototype de l'audit (`audit/scripts/contre_pretension.py`) montre que le moteur le permet déjà en pilotant `step()` sans le drapeau `_building_reference_state`.
**Acceptation** : atténuation des pics sur 11 cycles reproduite en ordre de grandeur (l'article mesure −10,3 % ; le schéma viscoélastique prototype donne −6,4 % ; l'actuel : −1,8 %) ; Fig. 11 d'EXP (relaxation 60 s post-étirement) simulable. Garder `elastic_tk_reference` par défaut tant que 3.1/3.2 ne sont pas validés ensemble — noter que l'audit a relevé un paradoxe (les vallées théoriques publiées collent mieux au schéma élastique), donc trancher par comparaison aux figures, pas par principe.

### 3.3 Identification matériau pour les muscles de l'équipe — [M6, M7] — effort L, dépend des essais E1-E4

> **✅ ITEM 3.3 EXÉCUTÉ le 21/08/2026** — outil livré et dégrossissage fait ;
> l'identification finale reste conditionnée aux essais E1-E4.
> **Outil** : `identification/identifier_materiau.py` (+ README et fiche
> muscle type) — trois volets : A) détection de paliers et fit
> multi-exponentiel (τ_k et fractions de relaxation, sans géométrie ; un
> palier à force croissante est signalé, jamais forcé) ; B) simulation de
> l'essai avec la pression mesurée et facteur d'échelle en forme fermée
> (absorbe la géométrie inconnue tant que E1 manque) ; C) verdict
> « spectre vs structure » sur la courbure F(P) de la première montée.
> Formats : CSV Bruno, CSV français, xlsx. Test permanent :
> `test_identification_tools`.
> **Dégrossissage (Bruno 153353/155135, Sacha muscle B)** : relaxation aux
> paliers de **4 à 9 %** (souvent nulle, paliers à force croissante) contre
> **~83 %** à long terme pour la Table A1 → spectre publié **non
> transférable, chiffré** (M6 clos) ; les données pointent vers E0/ΣE ≈ 0,9
> (article : 0,17). Courbure F(P) de Bruno : b/a mesuré **+100/+118 /MPa**
> contre −3/−1 simulé, quel que soit le facteur d'échelle → **limite
> STRUCTURELLE** du modèle de l'article (M7 clos : documenter, ne pas
> compenser par le spectre). Échelles 5,6-6,2 = géométrie inconnue (E1).
> Sacha muscle B : corrélation 0,87, courbure compatible — cadence 1 s
> supposée (E3 pour la confirmer).
> **Reste à faire quand les essais arrivent** : remplir une fiche par muscle
> (E1), identifier le vrai spectre sur relaxation uniaxiale propre (E2) ou
> maintien sans fuite (E3), rejouer l'outil, et reporter le spectre identifié
> dans les réglages (les champs maxwell_* de l'interface).
Le spectre de la Table A1 n'est **pas transférable** (relaxation simulée de signe opposé aux essais de Bruno ; sensibilité dF/dP ×2 chez Sacha). Construire un petit outil d'identification :
1. Ingestion d'une fiche géométrique par muscle (E1) ;
2. Fit moindres carrés du spectre {E_k, τ_k} sur des paliers de relaxation propres (E2/E3) — commencer avec les données Bruno/Sacha existantes pour dégrossir ;
3. Trancher « spectre mal paramétré » vs « Maxwell structurellement insuffisant » pour la courbure F(P) quasi quadratique (E4) : si aucun spectre ne reproduit la courbure en première montée, elle est structurelle (dead-band/raidissement en pression) et sort du cadre de l'article — le documenter comme limite du modèle plutôt que de forcer un fit. Même traitement pour le **fluage positif de la force aux paliers** observé chez Bruno (~1 % du pic, de signe interdit pour un Maxwell généralisé bloqué) : l'essai E3 tranchera entre phénomène du muscle et dérive du capteur HX711 (même ordre et même signe mesurés à P = 0) avant toute conclusion sur le modèle.

### 3.4 Appariement des coefficients de Poisson effectifs — [m] — effort S

> **✅ ITEM 3.4 EXÉCUTÉ le 21/08/2026** (moteur `2026.08.21-audit-phase4-14`).
> Option `poisson_pairing` ajoutée (réglage + moteur, validée dans
> normalize_settings/settings_error) : `paper_crossed` (défaut, lettre des
> articles — physique par défaut inchangée au bit près) / `physical`
> (échange ν̄₁₂↔ν̄₁₃ dans les deux chemins d'essai). Test permanent
> `test_poisson_pairing_option` (impact non nul et < 5 % sur le couple,
> < 2 % force ; rejet des valeurs inconnues). Documenté au README
> (« Domaine d'utilisation ») avec la recommandation de poser la question
> aux auteurs.
Le code suit la lettre (probablement coquillée) des deux articles : ν̄(s→φ) porté par ε_r et ν̄(s→r) par ε_φ (facteur ~9 entre les deux). Ajouter une option d'échange derrière un drapeau (impact ≤ 1,5 % sur le couple) et, idéalement, poser la question aux auteurs.

---

## Phase 4 — Validation étendue et rangement

| Item | Contenu | Effort |
|------|---------|--------|
| 4.1 | Valider le mode suspendu contre EXP Fig. 11 (numériser la courbe de contraction, protocole eps=0 + équilibrage off — items 1.4/D3). C'est le seul mode du moteur sans validation numérique. **✅ Fait le 21/08/2026** — `validation/validation_fig11_exp.py` (numérisation calibrée sur les graduations, calée à +0,8 % du texte ; référence de la figure = état chargé à t=0+, démontré par le pas d'hélice 2R et les photos). **Fluage (fig. 11a) : VALIDÉ** — forme, échelles de temps et amplitude à ~10 %, RMSE 0,023 vs théorie article, surestimation du même signe que celle que l'article reconnaît pour sa propre théorie ; les deux modes de précontrainte coïncident exactement à eps=0 (vérifié). **Cycles (11b) : partiels** — cycle 1 reproduit (0,306 vs 0,308) puis ~50 % d'amplitude : **limite structurelle nouvelle documentée, absence de butée d'auto-contact des spires** (l'article s'appuie sur les spires jointives, pas = 2R, comme état de contraction maximale) + état « entraîné » jamais spécifié. Candidat de travaux futurs hors plan : butée d'auto-contact. | M |
| 4.2 | Étude de convergence documentée (n_layers 3→10, n_phi, dt) et réglages de production recommandés au README. L'audit a mesuré ≲1 % entre L3 et L6 sur la figure 7 — à confirmer après 3.1. **✅ Fait le 21/08/2026** — `validation/etude_convergence.py` + `validation/CONVERGENCE.md` (30 runs, 7 balayages incl. updated et suspendu, validation croisée des réglages combinés). Genoux : 3-4 couches ; n_phi neutre à partir de 8 (mais nφ4 : +1,9 % couple) ; **dt = paramètre limitant via le couple** (1er ordre, +0,5 % de biais résiduel à dt 0,1, extrapolation documentée) ; surcoût updated ×1,2 (pas ×2). Recommandés : interactif L3/nφ8/dt0,5 (< 2 %, ~4,5 s) ; production L6/nφ8/dt0,1 (< 0,15 %, ~40 s). | M |
| 4.3 | Appliquer D2 au dossier `livrable` (notice de dépréciation puis régénération depuis alpha V3). **✅ Fait le 21/08/2026** — décision finale : `livrable` archivé définitivement (notice mise à jour), **alpha V3 déclaré livrable courant** ; pas de régénération de copie, pour ne pas recréer la divergence de versions relevée par l'audit. | S |
| 4.4 | Clôture : cocher ici chaque constat corrigé, mettre à jour `audit/rapport_audit.html` ou l'artifact en ligne. **✅ Fait le 21/08/2026** — section « Suites de l'audit » ajoutée au rapport (copie `audit/rapport_audit.html` et artifact en ligne republié, même lien). Bilan des constats majeurs : M1 documenté (D1), M2/M4/M5 refermés pour l'essentiel (item 3.1), M3 clos (item 3.2, paradoxe arbitré), M6/M7 clos et chiffrés (item 3.3), M8 clos (Phase 0). Suite de tests finale : 25 tests. | S |

---

## Plan d'essais (équipe — préalable au quantitatif)

| # | Essai | Débloque |
|---|-------|----------|
| E1 | **Fiche géométrique par muscle** : rayons interne/externe du tube, rayon et pas d'hélice, angle de biais ou tours insérés, longueur active, extrémités non enroulées, prédéformation réellement appliquée. | Toute comparaison quantitative (3.3) |
| E2 | **Relaxation uniaxiale (ou DMA) sur le tube seul** → spectre de Maxwell du matériau réel. | 3.3 — le verrou principal |
| E3 | **Maintien bloqué sans fuite** (vanne fermée, capteur au plus près du muscle, ou pression asservie). L'essai « 20 min » actuel mesure la fuite du circuit (P chute de 37 %, corr(F,P) = 0,997). | Validation des constantes de temps ; arbitrage muscle vs capteur du fluage positif aux paliers (Bruno) |
| E4 | **Rampe quasi-statique lente 0→0,6→0 MPa**, pression mesurée au muscle. | Caractériser la courbure F(P) (M7) séparément de la viscoélasticité |
| E5 | **Essai suspendu avec mesure de déplacement** sous masse connue (aucune donnée de déplacement n'existe). | Validation du mode libre (4.1) |

---

## Correspondance constats → items

| Constat | Items |
|---------|-------|
| M1 profil de biais | 1.1, D1, 2.10 |
| M2 section figée | 1.2, 3.1 |
| M3 précontrainte élastique t_k | 1.3, 3.2 |
| M4 déficit force figure 7 | 0.1, 3.1, 3.2 |
| M5 déficit couple figure 7 | 0.1 (calibration ≈ 40 % du déficit), 3.1 |
| M6 relaxation non transférable (Bruno) | 3.3, E2 |
| M7 non-linéarité F(P) (Bruno) | 3.3, E4 |
| M8 validation absente de la suite | 0.1, 0.2, 0.3 |
| Mineurs code | 2.1 → 2.9, 2.11, 3.4 |
| Fluage positif aux paliers (Bruno) | 3.3 (pt 3), E3 |
| Mineurs documentation/protocole | 1.4, 1.5, D3, 4.1 |
| Divergence livrable | D2, 4.3 |

**Ordre recommandé** : Phase 0 → Phase 1 → décisions D1-D3 → Phase 2 → 3.1 → 3.2 → (E1-E4 en parallèle côté équipe) → 3.3 → Phase 4.
