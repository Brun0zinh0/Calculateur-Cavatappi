# Confrontation modèle / essais selon la vitesse de sollicitation (07/09/2026)

Objet : mesurer la vitesse de pression réellement appliquée sur le banc (MPa/s), comparer
la réponse mesurée et la réponse simulée **sous la même pression mesurée** en fonction de
cette vitesse, et exprimer le profil de pression du modèle en MPa/s plutôt qu'en débit de
seringue (mL/min). Moteur `2026.09.07-v4-16`, réglages schéma 18.

Fichiers : `confrontation_vitesse.py` (script, réexécutable), `vitesse/tables.md` (tables
brutes générées), `vitesse/fig1…fig5.png`. Les données brutes ne sont pas dans le dépôt.

## 1. Pourquoi une vitesse en MPa/s

Le profil généré par l'outil était défini par un débit (10 mL/min) et un volume de
seringue (1,5 mL), hérités du protocole de l'article : demi-période 60·V/Q = 9 s quelle
que soit la pression maximale. Sur le banc de la campagne, la seringue est comprimée à la
main (§ 6.2 du rapport), le débit n'est pas une consigne, et la pression n'est pas
asservie : chaque rampe a sa propre pente. La seule grandeur de sollicitation mesurable
est donc la vitesse de pression dP/dt.

Depuis le moteur v4-16 le profil généré est paramétré par `pressure_rate_mpa_s`
(demi-cycle = Pmax / vitesse). Valeur par défaut Pmax / 9 s (0,167 MPa/s à 1,5 MPa,
0,144 MPa/s à 1,3 MPa) : figure 7 de l'article et suite de tests inchangées. Les fichiers
de réglages et exports anciens sont migrés (vitesse équivalente p_max·Q/(60·V)) sans
toucher au reste ; les scripts qui passent encore débit et volume obtiennent exactement
l'ancien profil. Détails dans `README.md` (section « Vitesse de pression ») et
`PATCH_NOTES.md`.

## 2. Méthode

- 60 essais des huit muscles disposant d'une fiche géométrique (A, E, F, G, J, L, Q, R) :
  21 rampes lentes « 1 O », 15 cyclages « 3 C », 24 paliers « 3 P » (rampe initiale
  seulement, la « descente » d'un palier est la fuite du circuit, pas une rampe).
- Pression reconstruite depuis les comptes ADC et référencée au repos de chaque essai
  (§ 7.3.1), force non filtrée.
- Segmentation en rampes sur les extrema de la pression lissée (1 s, proéminence
  0,05 MPa, excursion ≥ 0,08 MPa) : 137 montées et 112 descentes exploitables. Vitesse
  d'une rampe = pente de l'ajustement linéaire de P(t) entre 10 et 90 % de la course.
- Par rampe : gain ΔF/ΔP (mN/MPa, niveaux médians sur ±0,5 s). Par cycle (montée puis
  descente) : aire de boucle sur plage commune de pression, brute et corrigée d'une
  dérive linéaire entre les deux extrémités (définition § 7.5.1), largeur à mi-course,
  retard entre le pic de force et le pic de pression, vitesse moyenne du cycle.
- Modèle : mêmes fiches géométriques et même configuration que la modélisation
  quantitative (précontrainte viscoélastique à 300 mm/min, relaxation axiale, ΣE 37,76 MPa,
  L3/nφ8/dt 0,25 s), pression mesurée injectée. Trois variantes : spectre médian des
  essais 10 N (τ₁ 20 s, τ₂ 380 s), spectre des maintiens longs 30 L (τ₁ 48 s, τ₂ 1 200 s),
  et spectre 10 N + mécanismes V4 avec les valeurs de la démonstration du muscle J
  (P_r0 0,257 MPa, r 0,7, P_c 0,02 MPa), appliquées telles quelles à tous les muscles.
  180 simulations, aucun échec. Les métriques sont calculées sur la force simulée avec
  exactement la même segmentation que sur la force mesurée.
- Balayage synthétique : muscle J à ε = 1,0, trois cycles triangulaires à vitesse imposée
  de 0,005 à 0,5 MPa/s (Pmax 0,45 MPa, et 1,5 MPa pour l'amplitude de l'article).

## 3. Les vitesses de la campagne

| Essai | rampes | min | Q1 | médiane | Q3 | max (MPa/s) |
|---|---|---|---|---|---|---|
| 1 O | 42 | 0,023 | 0,036 | 0,045 | 0,049 | 0,061 |
| 3 C | 183 | 0,010 | 0,034 | 0,040 | 0,047 | 0,060 |
| 3 P (montée) | 24 | 0,010 | 0,035 | 0,049 | 0,055 | 0,063 |

Toute la campagne se situe entre 0,01 et 0,06 MPa/s, soit trois à quatre fois plus lent
que le protocole de l'article (0,167 MPa/s). À l'intérieur d'un cyclage, la vitesse
varie d'un cycle à l'autre de ±20 % (banc manuel), ce qui fournit une plage étroite mais
réelle pour une comparaison intra-essai. Un 3 C dure 7 à 10 s par rampe, un 1 O 8 à
10 s sur la partie utile.

## 4. Résultats

### 4.1 Le modèle est presque insensible à la vitesse dans la plage des essais

Balayage synthétique (muscle J, Pmax 0,45 MPa) : entre 0,02 et 0,05 MPa/s, le gain du
premier cycle passe de 79 à 82 mN (+4 %) et l'aire corrigée de −2,8 à −1,4 mN·bar ; sur
deux décennies (0,005 à 0,5 MPa/s) le gain ne varie que de 74 à 84 mN. Les deux spectres
donnent la même chose. Explication : les constantes de temps identifiées (20 s et
380 s, ou 48 s et 1 200 s) sont longues devant les rampes de 7 à 10 s ; pendant une
rampe le tube répond quasi élastiquement, ce que la modélisation quantitative avait déjà
constaté sur les 1 O. À l'amplitude de l'article (1,5 MPa) l'aire simulée est franchement
négative (−55 à −14 mN·bar de 0,005 à 0,17 MPa/s) : c'est la relaxation de la
précontrainte, plus forte à haute pression, qui abaisse la branche de descente.

### 4.2 Gain d'actionnement : déficit ×0,7, très faiblement dépendant de la vitesse

| Rampes de montée | n | gain médian mesuré (mN/MPa) | gain simulé 10 N | rapport sim/exp |
|---|---|---|---|---|
| 1 O | 21 | 334 | 246 | 0,69 |
| 3 C | 92 | 369 | 261 | 0,68 |
| 3 P | 24 | 324 | 248 | 0,74 |
| ensemble | 137 | 359 | 255 | 0,69 |

Le rapport simulé/mesuré vaut 0,69 en médiane, identique pour les deux spectres (0,70
pour 30 L) et légèrement plus bas avec les mécanismes V4 (0,65, la pression d'engagement
retranchant une fraction de la pression à basse pression). Il croît faiblement avec la
vitesse : 0,66 sur le tiers des rampes les plus lentes, 0,73 sur le tiers le plus rapide
(Spearman ρ = 0,21, p = 0,012, n = 137). La cause est du côté du modèle, dont le gain
augmente avec la vitesse (ρ = +0,25, p = 0,016 sur les 3 C), alors que le gain mesuré
n'en dépend pas sur les 3 C (ρ = 0,06) et décroît sur les 1 O et 3 P (ρ = −0,39 et −0,49,
n = 21 et 24). Cette décroissance mesurée est fragile : peu de points, plage de vitesse
étroite, et une part vient du retard force/pression du § 4.4 (la force n'a pas fini de
monter quand la pression culmine, d'autant plus que la rampe est rapide). Elle ne doit pas
être lue comme un résultat établi.

### 4.3 Hystérésis : un plancher indépendant de la vitesse, plus une part qui croît avec elle

| Aire de boucle corrigée (mN·bar) | 3 C (91 cycles) | 1 O (21 boucles) | ρ(aire, vitesse) sur 3 C | pente (mN·bar par MPa/s) |
|---|---|---|---|---|
| mesure | +57,5 (brute 53,8) | +79,1 (brute 55,1) | +0,48, p < 0,001 | +970 [630 ; 1 320] |
| modèle, spectre 10 N | −3,4 | −10,0 | +0,23, p = 0,03 | +67 |
| modèle, spectre 30 L | −2,4 | −7,9 | +0,22, p = 0,04 | +51 |
| modèle + V4 (valeurs du muscle J) | +42,9 | +41,3 | +0,21, p = 0,05 | +257 [−12 ; 499] |

Trois constats.

- **Le cadre de Maxwell seul a le mauvais signe**, à toute vitesse de la plage : l'aire
  simulée reste négative (descente sous la montée) même après correction de dérive. Ce
  résultat étend à 112 cycles ce que le § 7.8 du rapport avait établi sur les rampes lentes.
- **Les mécanismes V4 restituent un plancher.** Avec les trois valeurs du muscle J
  appliquées à tous les muscles, l'aire simulée vaut +43 mN·bar sur les cyclages (75 % de
  la médiane mesurée) et +41 sur les rampes lentes (52 %), et elle est indépendante de la
  vitesse sur deux décennies (balayage : 31 à 37 mN·bar de 0,005 à 0,5 MPa/s), comme il se
  doit pour un frottement de Coulomb et un seuil d'engagement.
- **La mesure contient en plus une part qui croît avec la vitesse.** Sur l'ensemble des
  cycles, l'aire corrigée augmente d'environ 10 mN·bar par 0,01 MPa/s (Theil–Sen
  +970 mN·bar par MPa/s), avec une ordonnée à vitesse nulle de 19 mN·bar ; la corrélation
  tient à chaque pré-étirement (ρ = 0,58, 0,34 et 0,47 à ε = 0,8, 1,0 et 1,2). Elle ne
  tient pas à la correction de dérive : l'aire brute (ρ = 0,24, p = 0,02) et la largeur
  à mi-course (ρ = 0,31, p = 0,003) croissent aussi avec la vitesse, et la tendance
  subsiste en écartant le premier cycle de chaque essai, qui porte l'essentiel de la
  dérive (ρ = 0,50 corrigée, 0,33 brute, n = 76). Le rapport avait relevé la même
  tendance sur la durée des cycles (§ 7.5.5.1). Elle est
  cependant **pour l'essentiel inter-essais** : à l'intérieur d'un même cyclage, où la
  vitesse ne varie que de ±20 %, la pente est positive sur 9 essais sur 15 (médiane
  +540 mN·bar par MPa/s) mais la dispersion est telle que le signe n'est pas établi
  (Wilcoxon p = 0,85). Une part de la corrélation groupée peut donc être un effet
  spécimen (les muscles cyclés plus vite ont-ils aussi de plus grandes boucles ?), et
  seul un cyclage à vitesse imposée peut trancher.

Ce qu'une dissipation croissante avec la vitesse impliquerait pour le modèle : une
branche de Maxwell ne produit une boucle qui grandit avec la vitesse que tant que
ωτ < 1, soit, pour des rampes de 7 à 10 s, une constante de temps inférieure à 2 ou 3 s.
Aucun des deux spectres n'en contient : les essais de relaxation à pression nulle ne
résolvent pas les premières secondes après l'étirement manuel, et l'outil
d'identification en boucle fermée hérite de cette limite. Le retard du § 4.4 va dans le
même sens.

### 4.4 Un retard force/pression que le modèle ne produit pas

Sur les cyclages, le pic de force survient **0,56 s après** le pic de pression en médiane
(intervalle interquartile 0,32 à 0,83 s), 0,75 s sur les rampes lentes. Le modèle donne
0,07 s (10 N) et 0,13 s (V4), le balayage synthétique 0,00 s à toutes les vitesses. Le
retard mesuré diminue quand la vitesse augmente (ρ = −0,28, p = 0,007) et son produit par
la vitesse vaut 22 kPa en médiane : la force continue de monter pendant que la pression
lue redescend d'environ 0,02 MPa.

Un frottement de Coulomb ne l'explique pas : il donne un plateau de force au sommet
(la pression effective reste constante tant que l'élément de Jenkins se retourne), non
une montée qui se poursuit. Deux candidats restent : un retard hydraulique entre le
capteur et le muscle (résistance du raccord et compliance du muscle, la pression au
muscle culminant après celle lue au capteur), ou une retardation rapide du matériau ou de
la structure, de constante de temps voisine de 0,5 s, invisible dans les essais de
relaxation. Le premier est un artefact de mesure qui biaiserait toutes les boucles dans le
sens dissipatif ; le second est physique. Un échelon de pression avec un capteur au
raccord du muscle les sépare en une heure de banc.

## 5. Ce que cela change pour le modèle et pour le rapport

- L'outil simule désormais le banc à sa vitesse réelle (0,01 à 0,06 MPa/s) au lieu de
  la vitesse implicite de l'article (0,167 MPa/s). Dans cette plage, le moteur actuel est
  quasi insensible à la vitesse : les écarts établis au chapitre 7 (déficit de gain,
  absence de seuil, signe de la boucle) ne dépendent pas de la vitesse à laquelle les
  essais ont été conduits, et l'écart entre protocole de l'article et protocole du banc ne
  les explique pas.
- La confrontation en vitesse ajoute une signature nouvelle aux écarts du chapitre 7 :
  une dissipation qui croît avec la vitesse (plancher ≈ 20 mN·bar, +10 mN·bar par
  0,01 MPa/s, résultat groupé) et un retard force/pression de 0,5 s. Ni le spectre de
  Maxwell identifié, ni les mécanismes V4 (indépendants de la vitesse) ne les produisent.
  Le candidat le plus économe est une branche rapide (τ ≲ 2 s), à condition d'exclure
  d'abord le retard hydraulique.
- Pour le § 7.5.5.1 du rapport (« les cycles lents dissipent moins ») : la tendance est
  confirmée sur 91 cycles avec une vitesse mesurée rampe par rampe, mais elle n'est pas
  établie intra-essai. La conclusion reste la même : cyclage à vitesse imposée.

## 6. Essais qui trancheraient

| Essai | Ce qu'il tranche | Coût |
|---|---|---|
| Cyclage à vitesse imposée, 0,01 / 0,03 / 0,1 / 0,2 MPa/s, 10 cycles, sur un même muscle | dissipation intra-spécimen vs vitesse ; existence et τ de la branche rapide ; accommodation vs vitesse (§ 7.5.5) | 1 h |
| Échelon de pression avec capteur au raccord du muscle | retard hydraulique vs retardation matériau | 1 h |
| Relaxation à pression nulle échantillonnée dès la fin de l'étirement (premières 5 s) | branche τ < 3 s dans le spectre | 30 min |

Le premier essai se simule directement avec le nouveau paramètre (onglet hystérésis,
« plusieurs vitesses de pression injectée », ou profil généré à la vitesse voulue).

## 7. Réserves

- Les fiches géométriques de A, E, F et G viennent de Dataset.xlsx sans confrontation
  antérieure ; les niveaux absolus de force de ces muscles n'ont pas été vérifiés ici (la
  note porte sur les dépendances à la vitesse, calculées rampe par rampe).
- Les valeurs V4 (P_r0, r, P_c) sont celles ajustées sur un seul essai du muscle J ; leur
  application aux autres muscles n'est qu'une illustration du caractère indépendant de la
  vitesse de ces mécanismes.
- Le retard force/pression est mesuré sur l'instant du maximum de deux signaux bruités
  (quantification de la pression à un compte ADC) : la médiane est robuste, les valeurs
  individuelles ne le sont pas.
- La vitesse d'une rampe est une pente moyenne ; la commande manuelle produit des rampes
  dont la pente varie de 10 à 20 % en cours de course.
