# Protocole expérimental — discrimination du dead-band (préparation alpha V4)

Objectif : identifier le mécanisme physique responsable du retard initial de
la force par rapport à la pression (dead-band), afin de choisir la voie
d'implémentation V4 (voir `PLAN_V4.md`, § mécanismes M-A / M-B1 / M-B2).
Durée estimée : **une demi-journée de banc** (essais P0-P5) + **1 h** de
photos/vidéo (P6). Configuration : actionnement **bloqué** (montage de Bruno),
sauf mention contraire.

---

## 1. Préparation (avant tout essai)

**Checklist matériel :**

- [ ] Capteur de pression taré à l'atmosphère (noter l'offset) ; re-vérifier
      le zéro à la fin de chaque essai (dérive à consigner).
- [ ] Cellule de force : enregistrer **30 s de repos** avant chaque essai pour
      mesurer le bruit (σ) et la dérive du HX711 — ces 30 s font partie du
      fichier.
- [ ] Acquisition à **10 Hz minimum**, colonnes CSV du format de Bruno
      (`time_s`, `pressure_bar`, `force_mN`, `force_unfiltered_mN`, …) — la
      colonne **non filtrée est obligatoire** : le filtre ajoute ~1-2 s de
      retard pur qui fausserait le seuil.
- [ ] **Test d'étanchéité (P0, prérequis éliminatoire)** : monter à 0,3 MPa,
      fermer la vanne, attendre 5 min. Chute de pression < 2 % → OK.
      Sinon, corriger la fuite avant de continuer (l'essai « maintien
      20 min » de Sacha perdait 37 % : inutilisable pour un seuil).
- [ ] **Fiche par muscle testé** (c'est aussi l'essai E1 du plan de
      correction) : rayons du tube, rayon et pas d'hélice, longueur active,
      **longueur désenroulée à chaque extrémité (au pied à coulisse)**, angle
      de sortie approximatif des extrémités, prédéformation appliquée.
- [ ] Températures ambiantes notées ; **10 min de repos à P = 0** entre deux
      essais du même muscle (recouvrance viscoélastique).

**Définition opérationnelle du seuil P_s (la même pour tous les essais) :**
sur la force **non filtrée**, F₀ = médiane des 30 s de repos, σ = écart-type
de ces 30 s. P_s = pression au premier instant où F − F₀ > max(10 mN, 3σ)
**maintenu au moins 2 s**. Reporter aussi P à F − F₀ = 30 mN (seuil « franc »)
pour contrôle.

---

## 2. Les essais, dans l'ordre

### P1 — Rampe lente montée/descente (essai E4⁺) — le test d'hystérèse

- Muscle : celui de Bruno (référence), précharge habituelle (~600 mN).
- Rampe **0 → 0,3 → 0 MPa à ≤ 0,005 MPa/s** (≥ 60 s par sens ; réduire le
  débit du pousse-seringue en conséquence). **2 répétitions.**
- Mesures : P_s à la montée ; à la descente, pression P_d à laquelle la force
  retombe à F₀ + 10 mN.
- **Ce qui est tranché** : P_s ≈ P_d (pas d'hystérèse) → mécanisme de type
  jeu/engagement élastique (M-B). P_s > P_d nettement (> 0,02 MPa) →
  décollement/recollement interne (M-A).

### P2 — Seuil vs précharge (essai E7) — le test discriminant principal

- Même muscle, **3 précharges : ~300 / ~600 / ~1200 mN** (≈ 30/60/120 g de
  tension de montage), mesurées au repos avant chaque rampe.
- À chaque précharge : rampe lente 0 → 0,3 MPa (comme P1), **2 répétitions**.
- **Ce qui est tranché** : P_s **décroît** quand la précharge augmente →
  engagement en série (extrémités M-B1 ou jeu de montage M-B2 : la précharge
  pré-tend la chaîne). P_s **invariant** (± 0,01 MPa) → mécanisme interne à
  la section (M-A, regonflage).

### P3 — Effet de l'entraînement (variante « affaissement » de M-A)

- Même muscle, précharge ~600 mN : appliquer **10 cycles 0 → 0,5 → 0 MPa**
  au débit habituel, puis refaire immédiatement une rampe lente P1.
- **Ce qui est tranché** : P_s chute nettement après entraînement (comme le
  0,38 → 0,23 MPa de l'article) → l'état d'affaissement initial de la
  section participe (M-A confirmé dans sa variante « état entraînable »).

### P4 — Volume-pression à basse pression (essai E6) — la signature du regonflage

- Pousse-seringue à **débit constant connu** (le plus faible disponible),
  montée 0 → 0,15 MPa. Le volume injecté = débit × temps (noter le volume
  mort du circuit, estimé en répétant l'essai vanne fermée côté muscle si
  possible).
- Tracer V(P). **Ce qui est tranché** : compliance initiale anormalement
  élevée (beaucoup de volume avalé avant 0,05 MPa) puis coude → la basse
  pression regonfle la section (M-A), et la valeur du volume excédentaire
  donne directement le paramètre δ_c du modèle V4.1.

### P5 — Seuil vs longueur d'extrémités (essai E9) — le test des extrémités

- **Au moins 3 muscles** à longueurs d'extrémités désenroulées différentes
  (mesurées à la checklist : panel B/I/J/D de Sacha + muscle(s) de Bruno).
  Idéalement inclure un muscle à extrémités quasi nulles.
- Pour chacun : précharge ~600 mN, rampe lente 0 → 0,3 MPa, 2 répétitions.
- Tracer P_s en fonction de la longueur d'extrémités totale.
- **Ce qui est tranché** : P_s croît avec la longueur d'extrémités → les
  extrémités désenroulées sont le mécanisme (M-B1). P_s identique, y compris
  sur le muscle sans extrémités → M-A (ou M-B2 si P2 a montré la dépendance
  à la précharge).

### P6 — Observation directe (essais E8/E8⁺) — 1 h, hors banc de force si besoin

- Photos calibrées (règle dans le champ) du **diamètre externe** du tube à
  0 / 0,05 / 0,10 MPa (M-A : le regonflage se mesure).
- **Vidéo des extrémités désenroulées** pendant une rampe lente 0 → 0,1 MPa :
  le redressement des tangentes de sortie (effet Bourdon, M-B1) se voit à
  l'œil nu s'il est en jeu.

---

## 3. Feuille de résultats (à remplir)

| Essai | Muscle | Précharge (mN) | Sens | P_s (MPa) | P_d (MPa) | Notes (T°, dérive zéro…) |
|---|---|---|---|---|---|---|
| P1-a | | ~600 | ↑↓ | | | |
| P1-b | | ~600 | ↑↓ | | | |
| P2 à 300 ×2 | | ~300 | ↑ | | — | |
| P2 à 600 ×2 | | ~600 | ↑ | | — | |
| P2 à 1200 ×2 | | ~1200 | ↑ | | — | |
| P3 (après 10 cycles) | | ~600 | ↑ | | — | |
| P5 — muscle 1 (L_ext = … mm) | | ~600 | ↑ | | — | |
| P5 — muscle 2 (L_ext = … mm) | | ~600 | ↑ | | — | |
| P5 — muscle 3 (L_ext = … mm) | | ~600 | ↑ | | — | |
| P4 | | — | ↑ | volume excédentaire : … mL | | |

Critère de validité : les 2 répétitions d'un même point doivent donner P_s à
± 0,01 MPa ; sinon, chercher la cause (fuite, dérive, repos insuffisant) et
refaire.

---

## 4. Arbre de décision → démarche alpha V4

Lire les résultats dans cet ordre :

1. **P2 : P_s décroît avec la précharge ?**
   - **OUI** → mécanisme d'engagement en série. Puis **P5** :
     - P_s croît avec la longueur d'extrémités → **M-B1** :
       V4 = élément d'extrémité non linéaire (transition
       `tangent_beam` → `axial_rod` + redressement de Bourdon),
       paramètres = géométrie mesurée des extrémités. → **voie V4.0-enrichie**.
     - P_s indépendant des extrémités → **M-B2** (jeu de montage) :
       V4 = simple jeu série δ₀ (mesuré) ; le dead-band est un artefact de
       banc, pas une propriété du muscle — le documenter comme tel.
       → **voie V4.0 minimale**.
   - **NON (P_s invariant)** → mécanisme interne. Puis **P4 + P6 + P1** :
     - compliance initiale élevée en P4, regonflage visible en P6,
       hystérèse en P1 → **M-A** : V4 = contact paroi-nylon dans le BVP
       radial, paramètre δ_c tiré de P4. → **voie V4.1**.
     - P3 positif en plus → variante « état d'affaissement entraînable »
       (reproduire aussi le 0,38 → 0,23 MPa de l'article).
2. **Cas mixte** (précharge ET extrémités jouent) : implémenter M-B1 d'abord
   (le moins coûteux), re-mesurer le seuil résiduel, puis décider M-A sur ce
   résidu.
3. **Aucun signal net** (P_s erratique, non répétable, insensible à tout) :
   ne rien implémenter — documenter le dead-band comme non modélisé
   (falsifiabilité assumée, PLAN_V4 § garde-fous) et vérifier l'instrumentation
   (filtre, zéro, fuites).

---

## 5. Après la séance

- Envoyer les CSV bruts (avec les 30 s de repos et la colonne non filtrée) —
  ils passeront dans `identification/identifier_materiau.py` pour un
  dépouillement automatique (paliers, seuils, courbure), et la fiche
  géométrique de chaque muscle dans `identification/fiche_muscle_exemple.json`.
- Les rampes lentes de P1/P2 servent AUSSI à l'essai E4 du plan de
  correction (courbure F(P) séparée de la viscoélasticité) : deux plans
  alimentés par la même séance.

---

## 6. E2 — relaxation du tube nu (ajouté le 28/08, après E10)

**But.** Mesurer la relaxation du PVC étiré seul, sans aucun nœud ni scotch
dans le chemin d'effort, à plusieurs amplitudes. C'est le juge de paix de
l'incompatibilité d'amplitude (résultat n° 3 de
`essais/MODELISATION_QUANTITATIVE.md`) : E10 a montré que le stick-slip
d'ancrage existe mais ne porte que ~7 % de la relaxation — le reste est
lisse et E2 dit si c'est le matériau. Un tube tiré selon son axe mesure
exactement E_axial(t), la grandeur que le mode `axial_test_only` fait
relaxer : le résultat se reporte tel quel dans les champs maxwell_*.

**Éprouvettes.** Tube de PVC étiré du même lot que les muscles, nylon
retiré (le faire coulisser hors d'un segment droit). 4 éprouvettes de
~210 mm (150 de longueur libre + 2×30 pour les mors) : une par amplitude
(1 %, 3 %, 10 %), une pour l'essai long. Tracer deux repères fins à
100,0 mm d'écart sur la longueur libre (déformation vraie au pied à
coulisse — ne pas se fier au déplacement des mors).

**Mors sans nœud** (la leçon d'E10) — deux options :
- *Broche interne + collier* (rapide, matériel existant) : enfiler 10 mm de
  corde à piano Ø ≈ 1,1 mm (= diamètre interne) dans chaque extrémité, puis
  serrer le tube sur la broche dans l'accouplement / un collier à vis. La
  broche empêche l'écrasement et le glissement.
- *Collage* (le plus sûr) : percer deux embouts alu Ø 2,0 mm, dégraisser,
  cyanoacrylate, insérer 8-10 mm, 24 h de prise.
Dans les deux cas : trait blanc à cheval tube/mors de chaque côté + photos
caméra FIXE avant/après (E2 s'auto-valide).

**Mise en charge.** Échelon rapide (< 2 s) : préparer la position d'accroche
à la longueur cible et accrocher d'un geste, comme la précontrainte des
muscles — ne PAS visser l'écrou pendant l'essai (30 tours = rampe lente).
Enregistrer 30 s de repos avant l'échelon. Aire de section
A = π(0,96² − 0,56²) = 1,91 mm² → forces attendues (ΣE ≈ 37,8 MPa) :

| ε | allongement (L₀ = 150 mm) | F attendue | équivalent |
|---|---|---|---|
| 1 % (régime actionnement) | 1,5 mm | ~0,7 N | 74 g |
| 3 % (≈ contrainte du tube en précontrainte, 1,12 MPa) | 4,5 mm | ~2,2 N | 220 g |
| 10 % (grande amplitude) | 15 mm | ~7,2 N | 735 g — vérifier la capacité de la cellule ; sinon plafonner à 5 % |

**Séquence.** Banc sous tension 10-15 min avant (StabCell) ; tare ; par
éprouvette vierge : repos 30 s → échelon → maintien 600 s sans toucher au
banc → photos. Ordre 1 % → 3 % → 10 %. Puis l'essai long : 3 %, maintien
45-60 min (contraint enfin τ₂ et E₃). Noter la température de la pièce.

**Dépouillement.** E(t) = F(t)/(A·ε) ; passer les CSV dans
`essais/identifier_spectre.py` — cette fois les modules sont absolus
(E₀ = F∞/(A·ε), ΣE = F(0⁺)/(A·ε)), à comparer à ΣE = 37,76 (article) et
E₀/ΣE = 0,930 (déconvolution). Lecture :
- fractions relaxantes **faibles à toutes les amplitudes** → le matériau est
  quasi élastique : la composante lisse des « 10 N » est du fluage d'ancrage
  continu ; la configuration actionnement vaut partout, le spectre des
  « 10 N » n'est pas matériel ;
- fractions **croissantes avec ε** (petites à 1 %, ~9 % à 3-10 %) →
  dépendance d'amplitude (type Payne) confirmée : limite documentée du cadre
  Maxwell linéaire, spectres par régime (recommandation actuelle) ;
- fractions **fortes et égales partout** → matériau viscoélastique linéaire :
  la quasi-élasticité de l'actionnement vient du trajet de chargement
  multiaxial — renforce `axial_test_only`, spectre mesuré utilisable partout.

*Optionnel : E2b, même essai sur le filament de nylon seul (rôle secondaire,
il travaille en compression dans le muscle).*
