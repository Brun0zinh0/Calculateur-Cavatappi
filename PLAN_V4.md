# Plan V4 — prédire physiquement le retard initial de la force (dead-band)

Objectif : que le modèle **prédise** le seuil de pression avant montée de la
force (dead-band) et la super-linéarité F(P) qui l'accompagne, à partir de
mécanismes physiques à paramètres **mesurables indépendamment** — jamais par
un coefficient de calage destiné à reproduire artificiellement la courbe
(principe déjà posé par les PATCH_NOTES de la V2). Chantier complémentaire :
la butée d'auto-contact des spires, l'autre limite structurelle identifiée
(validation fig. 11).

Base de départ : alpha V3, moteur `2026.08.21-audit-phase4-14`, baseline de
non-régression figure 7 en place. **Règle absolue : tous les mécanismes V4
sont des options désactivées par défaut** — les modes actuels restent
bit-identiques, la baseline le vérifie à chaque étape.

---

## 1. Ce que l'on sait déjà (acquis de l'audit 2026-08)

| Fait | Source |
|---|---|
| Dead-band mesuré ~0,05-0,1 MPa sur les essais de Bruno (précharge ~580 mN) | essais 153353/155135 |
| Dead-band de l'article : 0,38 MPa (cycle 1) → 0,23 MPa (cycle 11) — attribué par les auteurs à « l'expansion du diamètre interne » | BLOCKED p. 5, fig. 3(b) |
| Courbure F(P) quasi quadratique : b/a ≈ +100 à +118 /MPa, hors de portée du modèle quel que soit le facteur d'échelle (simulé : −3 à −1) | outil identification, M7 |
| Le modèle actuel répond linéairement dès P = 0⁺ (aucun mécanisme de contact ou d'engagement dans le cadre de l'article) | audit, M7 |
| La forme hors seuil est quasi élastique (relaxation 4-9 % vs 83 % Table A1) — le seuil n'est PAS un effet viscoélastique | outil identification, M6 |
| En suspendu, l'article s'appuie sur les spires jointives (pas = 2R_TO0) comme état de contraction maximale ; le moteur n'a pas cette butée (cycles fig. 11 à ~50 %) | validation fig. 11 |
| Le tube est fabriqué autour d'un monofilament nylon anti-collapse : à P = 0 la paroi étirée/torsadée/recuite peut être partiellement affaissée sur le nylon | articles, protocole de fabrication |

## 2. Mécanismes physiques candidats

### M-A. Regonflage de section / contact paroi-nylon (candidat principal)

À P = 0, la section n'est pas parfaitement circulaire et porteuse : la paroi
repose en partie sur le monofilament. La basse pression sert d'abord à
**asseoir et regonfler la section** (travail radial) avant que le mécanisme
biais → détorsion → force ne s'engage. Prédit naturellement : un seuil P_s
qui dépend de l'état d'affaissement initial (donc du spécimen et de
l'entraînement — cohérent avec le 0,38 → 0,23 MPa de l'article), et une
montée progressive (engagement de la circonférence portante croissant avec
P) → super-linéarité.

Formulation : condition de contact unilatéral interne dans le BVP radial
existant — une fraction de la face interne s'appuie sur le nylon (pression
reprise par le contact) tant que le décollement n'est pas atteint ; la
fraction engagée χ(P) ∈ [0,1] multiplie la contribution de la couche à la
génération de force. Paramètre physique : le défaut de circularité /
l'interférence initiale δ_c (mm), mesurable (E6/E8), pas un coefficient
libre.

### M-B. Engagement des extrémités désenroulées / jeu série (candidat banc + spécimen)

Deux variantes du même mécanisme d'engagement en série :

- **M-B1 — redressement des extrémités désenroulées** : les portions non
  hélicoïdales quittent la spire courbées/coudées. À basse pression, deux
  effets absorbent les premiers incréments : le **redressement de type
  Bourdon** du tube courbé pressurisé (piloté par la pression → produit un
  seuil en pression), et la **transition flexion → traction** de l'extrémité
  (très souple tant qu'elle fléchit — asymptote `tangent_beam` —, quasi
  rigide une fois alignée — asymptote `axial_rod`). L'élément V4 interpole
  entre ces deux asymptotes déjà présentes dans alpha V3, avec des paramètres
  purement géométriques et mesurables : longueur désenroulée, angle de
  sortie, rayon du coude. Prédit un seuil croissant avec la longueur
  d'extrémités et décroissant avec la précharge, et la super-linéarité
  (pente faible puis raidissement à l'alignement).
- **M-B2 — jeu de montage pur** : un jeu δ₀ dans la chaîne
  actionneur-capteur qui doit se fermer avant transmission :
  F = k_ext·max(0, Δ − δ₀). Seuil dépendant de la précharge, identique
  montée/descente. Paramètre : δ₀ (mm) au comparateur, ou identifié sur E7.

Note : le ressort série **linéaire** actuel d'alpha V3 ne peut PAS produire
de seuil (il adoucit la pente sans la retarder) — c'est bien la non-linéarité
d'engagement qui fait le dead-band.

### M-C. Raidissement géométrique grande déformation (candidat courbure)

À ~20 % de déformation radiale sous 0,6-1,4 MPa, la raideur tangente réelle
du PVC étiré croît (hyperélasticité) et la géométrie porteuse évolue — le
modèle actuel est à modules constants. Contribue à la courbure au-dessus du
seuil plus qu'au seuil lui-même. À n'ouvrir que si le résidu de courbure
persiste après M-A/M-B (dépend de la reprise du mode `updated`, item 3.1).

### M-D. Butée d'auto-contact des spires (mode suspendu)

Contrainte unilatérale pas ≥ 2·R_TO0 dans le solveur suspendu (les spires ne
s'interpénètrent pas) : referme le déficit des cycles de la fig. 11
(amplitudes à ~50 %) en donnant l'état de contraction maximale sur lequel
l'article s'appuie. Paramètre : aucun nouveau (géométrie existante).

## 3. D'abord DISCRIMINER, ensuite implémenter

Principe : ne coder que le mécanisme que les essais désignent. Programme de
discrimination (une séance de banc) :

| Essai | Question tranchée |
|---|---|
| **E7 — seuil vs précharge** : rampe lente 0→0,3 MPa à 3 précharges (300/600/1200 mN) | Seuil qui diminue avec la précharge → M-B (jeu) ; seuil invariant → M-A (interne) |
| **E4⁺ — montée/descente lente** (0→0,3→0 MPa, < 0,005 MPa/s) | Seuil identique dans les deux sens → M-B ; hystérèse du seuil → M-A (décollement/recollement) |
| **E6 — volume-pression à basse pression** (pousse-seringue, V(P) de 0 à 0,15 MPa) | Compliance initiale anormalement élevée (regonflage) → M-A, et mesure directe de δ_c |
| **E8 — observation directe** (photo/mesure du diamètre externe à 0 / 0,05 / 0,1 MPa) | Confirmation visuelle du regonflage ; entrée quantitative pour M-A |
| Rejouer E7 après 10 cycles d'entraînement | Reproduire le 0,38 → 0,23 MPa de l'article → valide la variante « état d'affaissement » de M-A |
| **E9 — seuil vs longueur d'extrémités** : comparer le seuil sur des muscles à extrémités désenroulées différentes (B/I/J/D de Sacha, muscles de Bruno) | Seuil croissant avec la longueur d'extrémités → M-B1 ; identique sur un muscle à extrémités quasi nulles → M-A |
| **E8⁺ — vidéo des extrémités** entre 0 et 0,1 MPa | Redressement visible des tangentes de sortie → M-B1 confirmé visuellement |

## 4. Implémentation par étapes

### V4.0 — cadre « contacts et engagements » + les deux contraintes unilatérales simples — effort M

- Infrastructure : gestion propre de contraintes unilatérales dans les deux
  solveurs (résidus avec complémentarité régularisée, garde de convergence,
  journalisation de l'état engagé/décollé).
- **M-B (jeu série δ₀)** : extension de `lock_blocked_series_reference` /
  `step_blocked_series` — nouveau réglage `series_gap_mm` (défaut 0, inactif).
- **M-D (butée de spires)** : contrainte pas ≥ 2·R_TO0 dans
  `_find_suspended_state` — option `coil_contact_stop` (défaut off).
- Tests : seuil exact reproduit sur cas synthétique (M-B) ; re-run de
  `validation_fig11_exp.py` avec butée → critère : amplitudes des cycles 2-10
  recouvrées à mieux que −20 % (contre −50 % aujourd'hui) ; baseline figure 7
  strictement identique (options off).

### V4.1 — regonflage de section / contact paroi-nylon (M-A) — effort L

- État initial : section porteuse partielle paramétrée par δ_c (interférence
  paroi-nylon) issue de E6/E8.
- BVP radial : condition interne mixte (portion en contact : continuité avec
  le nylon rigide ; portion libre : −P), fraction engagée résolue par
  complémentarité ; la force axiale du tube est intégrée sur la fraction
  engagée.
- Sorties : P_seuil prédit, χ(P), courbe F(P) complète.
- Validation : **un seul jeu (δ_c, δ₀) mesuré doit reproduire À LA FOIS** le
  seuil de Bruno (~0,05-0,1 MPa à sa précharge) et, avec l'état d'affaissement
  « non entraîné », l'ordre de grandeur de l'article (0,2-0,4 MPa). Critère :
  seuil prédit à ±30 %, b/a simulé ≥ 50 % du mesuré sur 0-0,3 MPa.

### V4.2 — raidissement grande déformation (M-C) — effort L, conditionnel

- Uniquement si le résidu de courbure au-dessus du seuil reste > 30 % après
  V4.1. Pré-requis : élucider la composition long-terme du mode `updated`
  (reprise de l'item 3.1 sur protocole cyclique) ; puis tangente
  hyperélastique (1 paramètre mesuré par E2, pas ajusté).

### V4.3 — clôture — effort S

- Documentation README (mécanismes, domaines de validité, paramètres et leur
  essai de mesure), entrée PATCH_NOTES, mise à jour de la baseline UNIQUEMENT
  si un défaut change (décision explicite), rapport de validation V4.

## 5. Garde-fous

- **Aucun coefficient sans essai** : chaque paramètre nouveau (δ_c, δ₀) a son
  protocole de mesure dans le tableau §3 ; s'il n'est pas mesuré, le
  mécanisme reste off.
- **Non-régression** : options off par défaut ; baseline figure 7 et suite de
  tests (26) inchangées ; chaque mécanisme ajoute ses propres tests.
- **Falsifiabilité** : si E7/E4⁺ ne désignent aucun mécanisme (seuil ni
  dépendant de la précharge ni hystérétique ni lié au volume), documenter le
  dead-band comme non modélisé et s'arrêter — ne pas forcer.
- Attention à l'artefact connu : la colonne de force filtrée ajoute ~1-2 s de
  retard temporel pur — toute identification de seuil se fait sur
  `force_unfiltered` (ou avec le retard du filtre caractérisé).

## 6. Ordre recommandé et dépendances

1. Séance de discrimination (§3 : E7, E4⁺, E6, E8 — une demi-journée de banc) ;
2. V4.0 (indépendant des essais pour M-D ; M-B activable dès δ₀ mesuré) ;
3. V4.1 si M-A désigné ; V4.2 si résidu ; V4.3.
4. Synergie avec le plan de correction : E1 (fiche géométrique) et E2 (spectre
   réel) restent les préalables de tout accord quantitatif — la V4 traite la
   *forme* à basse pression, pas l'échelle ni la relaxation.
