# Modélisation quantitative des muscles — première confrontation sans facteur d'échelle (27/08/2026)

Pipeline : géométrie réelle du `Dataset.xlsx` (fiches déposées dans
`alpha V3/identification/fiches/`) + extrémités désenroulées actives en
compliance série + spectre de Maxwell ajusté sur les essais « 10 N »
(relaxation à P = 0) + précontrainte viscoélastique (rampe rapide
300 mm/min ≈ étirement manuel). **Seule grandeur empruntée à l'article :
ΣE = 37,76 MPa** (même famille de PVC étiré). Moteur alpha V3
`2026.08.21-audit-phase4-14`, réglages L3/nφ8 (CONVERGENCE.md).

> **Échelle de pression.** La pression injectée est reconstruite depuis les
> comptes ADC et référencée au repos propre à chaque essai (§ 6.4.2 du
> rapport). Le compte de repos variant de 91 à 103 selon la séance, l'échelle
> nominale sous-estimait la pression de façon variable. Scripts :
> `identifier_spectre.py`, `recalcul_correlations.py`, `spectre_identifie.json`.

**Spectre consolidé** (médianes sur les 18 essais « 10 N », `identifier_spectre.py`) :

| Paramètre | Valeur | Article (Table A1) |
|---|---|---|
| E₀ | **33,481 MPa** | 6,36 |
| E₁ / η₁ | **1,627** / **31,9** | 20,67 / 154,57 |
| E₂ / η₂ | **2,652** / **1 014,7** | 5,98 / 977,79 |
| E₃ | **0** (non contraint à 600 s) | 4,75 |
| τ₁ / τ₂ | **19,6 s** / **382,6 s** | 7,5 s / 163,5 s |
| E₀/ΣE | **0,887** | 0,168 |

> **Mise à jour 01/09 — essais longs « 30 L »** (1 800 s, ε = 1,0, muscles
> C/I/K/O ; `ANALYSE_PRELIMINAIRE.md` § F) : τ₂ vaut en réalité **≈ 1 200 s**
> (l'estimation 382,6 s était bornée par la fenêtre de 600 s) et E₀/ΣE
> descend à **≈ 0,83** — borne haute, la queue n'étant pas close à 30 min.
> L'écart à l'article demeure (facteur ~5). Branche rapide inchangée : rien
> ne change pour l'actionnement ; pour un maintien ≥ 10 min, utiliser
> E₀ ≈ 31,3, E₁ ≈ 2,4 / η₁ ≈ 115, E₂ ≈ 4,3 / η₂ ≈ 5 270 MPa.

Ce spectre relève de la famille « passe 4 » du résultat n° 3 : fractions de
force brutes, relaxation axiale seule. Les chiffres ci-dessous l'emploient
tel quel pour les neuf cas J/L/Q, là où les valeurs par muscle du § suivant
restent préférables pour un cas isolé.

C'est le spectre du **muscle monté**, non du matériau. Déconvolué du partage
tube/nylon (`spectre_pvc_seul.json`), le PVC seul donne E₀/ΣE = **0,930** au
lieu de 0,887 : la cellule mesure une force totale plus petite que celle du
tube, puisque le nylon y travaille en compression, si bien que les fractions
relaxantes du muscle monté **surestiment** celles du matériau. Voir
`ANALYSE_PRELIMINAIRE.md` § B.

## Résultat n° 1 — les niveaux absolus tombent juste (à ε = 0,8)

| Cas | F₀ simulé | F₀ mesuré | Ratio |
|---|---|---|---|
| J 0,8 | 1157 | 1067 | **1,08** |
| L 0,8 | 1170 | 1088 | **1,08** |
| Q 0,8 | 1250 | 1430 | **0,87** |
| R 0,8 | 1155 | 906 | 1,27 |
| D 0,8 | 1280 | 984 | 1,30 |

Avec la géométrie mesurée, la force de précontrainte absolue est prédite à
**±13 % pour J/L/Q** (contre un facteur ~6 d'écart avec la géométrie de
l'article) — le module ΣE du PVC de l'article est donc transférable en
première approximation. Le ratio dérive avec la précontrainte
(J : 1,08 → 1,36 → 1,53) : cohérent avec la convention d'application du
pre-strain (le moteur applique ε à la spire seule, la mesure porte sur la
longueur totale — item 2.11 de l'audit, amplifié chez J dont les extrémités
passent de 20 à 38 mm) et avec un adoucissement du matériau à grande
déformation.

## Résultat n° 2 — la forme d'actionnement est reproduite (corr 0,94-0,98)

Sur les « 1 O » (montée-descente lente, pression mesurée injectée), la
configuration finale donne des corrélations de **0,94 à 0,98, médiane 0,964**
sur les neuf cas J/L/Q × {0,8 ; 1,0 ; 1,2} — contre 0,44-0,69 avec le spectre
Table A1 de l'article lors de l'audit.

| Cas | corr | gain mesuré | gain simulé | ratio |
|---|---|---|---|---|
| J 0,8 / 1,0 / 1,2 | 0,973 / 0,975 / 0,953 | 180 / 185 / 153 | 124 / 88 / 65 | 0,69 / 0,47 / 0,43 |
| L 0,8 / 1,0 / 1,2 | 0,974 / 0,964 / 0,959 | 161 / 175 / 154 | 147 / 124 / 107 | 0,92 / 0,71 / 0,69 |
| Q 0,8 / 1,0 / 1,2 | 0,972 / 0,963 / 0,941 | 225 / 208 / 142 | 158 / 141 / 144 | 0,70 / 0,68 / 1,02 |

Le gain d'actionnement simulé vaut **69 % du mesuré en médiane** (43 à 102 %
selon les cas) ; le déficit porte la signature du **dead-band/super-linéarité**
absent du modèle (chantier V4) et des effets du résultat n° 1.

Le ratio se dégrade nettement avec la précontrainte chez J (0,69 → 0,43),
moins chez L et Q. C'est le même sens que la dérive du ratio F₀ du résultat
n° 1, et la même explication candidate : la convention ε spire seule contre
longueur totale, amplifiée chez J dont les extrémités passent de 20 à 38 mm.

**Effet du référencement de la pression.** Comparé à l'échelle nominale,
toutes choses égales par ailleurs : la corrélation perd 0,006 en médiane
(0,973 → 0,964) et le ratio de gain gagne 2 points (67 % → 69 %). Les deux
plus fortes dégradations de corrélation, J 1,2 (−0,013) et Q 1,2 (−0,017),
sont les deux cas de plus grand dead-band mesuré (0,253 et 0,312 MPa). Le
référencement ne change donc pas les conclusions : il déplace légèrement les
chiffres dans le sens attendu, et rend l'écart au modèle un peu plus visible
là où le dead-band est le plus marqué.

## Résultat n° 3 — une incompatibilité d'amplitude, découverte structurante

Quatre passes ont été nécessaires, et leur échec partiel est le résultat :

| Passe | Spectre | corr « 1 O » | relax 600 s simulée (mesurée : −3 à −16 %) |
|---|---|---|---|
| 1 | fractions de force brutes, anisotropie article | 0,96-0,98 | −0,5 à −2,3 % (trop faible) |
| 2 | fractions ×5 (amortissement composite corrigé) | 0,77-0,94 | −3 à −10 % ✓ |
| 3 | fractions ×5, relaxation **axiale seule** | 0,83-0,95 | −2 à −7 % |
| **4 (retenue)** | fractions brutes, relaxation axiale seule | **0,96-0,98** | −0,5 à −2,4 % |

*(Les corrélations de ce tableau sont établies sur l'échelle de pression
nominale. Elles ne sont pas directement comparables à celles du résultat n° 2,
qui emploient l'échelle référencée ; seule la comparaison entre passes est ici
en jeu, et le référencement les décale toutes du même ordre — environ −0,006.)*

**Aucun spectre de Maxwell linéaire ne reproduit à la fois** la relaxation de
précontrainte (~−9 % en 10 min, à ~100 % de déformation) **et** la réponse
d'actionnement quasi élastique (à ~1 % d'amplitude). Deux explications
candidates, à départager :

1. viscoélasticité **dépendante de l'amplitude** (type effet Payne) — hors du
   cadre Maxwell linéaire de l'article, à documenter comme limite ;
2. la « relaxation » des 10 N n'est pas matérielle mais du **fluage
   d'ancrage** (glissement des nœuds/fixations sous tension) — auquel cas la
   configuration actionnement est la bonne partout. **Essai discriminant
   E10** : marquer nœuds et fixations au feutre, photographier avant/après un
   « 10 N » — un glissement se voit.

Une troisième explication candidate — une **dérive de la cellule de charge**
comptée comme relaxation — est **exclue depuis le 28/08** : cinq essais à
blanc de 600 s (`StabCell/`, `ANALYSE_PRELIMINAIRE.md` § D) bornent la dérive
à +11 mN, de signe opposé aux −30 à −210 mN mesurés, avec un bruit de ±2 mN.
L'incompatibilité d'amplitude est donc bien physique (matériau ou montage),
et le départage se joue entre les deux candidates ci-dessus.

**Premier essai E10 réalisé (28/08, `test E10/`)** : relaxation −6,7 %
reproduite sous marquage ; photos non concluantes (caméra déplacée entre T0
et T10, à refaire caméra fixe) ; mais l'analyse des décrochements persistants
de la courbe — validée sur les blancs StabCell — montre que le **stick-slip
d'ancrage est réel et ne porte que ~7 % de la relaxation en médiane** (jusqu'à
33 % sur E10). L'essentiel de la chute est lisse : matériau, ou fluage
d'ancrage continu. Le juge de paix reste **E2** (relaxation sur tube nu, sans
nœud dans le chemin d'effort). Détail : `ANALYSE_PRELIMINAIRE.md` § E.

## Configuration recommandée par usage (champs de l'interface)

- **Prédire l'actionnement** (gains, formes, cycles) : fiche du muscle
  (`identification/fiches/`) + `maxwell_anisotropy_mode = axial_test_only` +
  fractions de force brutes du « 10 N » du muscle. Exemple J 0,8 :
  E0 = 32,29, E1 = 2,04 (η1 = 24,5), E2 = 1,55 (η2 = 420), E3 = 1,89
  (η3 = 9440) MPa. À défaut d'un « 10 N » propre sur le muscle considéré,
  le spectre consolidé de l'en-tête (`spectre_identifie.json`) donne des
  corrélations de 0,94 à 0,98 sur J/L/Q — l'écart au spectre par muscle
  reste inférieur à l'écart au modèle lui-même.
- **Prédire un maintien long** (dérive de précontrainte) : mêmes fiches,
  fractions ×5 (plafonnées à 35 %) — en acceptant la dégradation de la forme
  d'actionnement.
- Dans les deux cas : `prestrain_reference_mode = viscoelastic_history`.

## Limites et suites

- Pas de dead-band dans le modèle (V4) : la mesure démarre après **0,15 MPa
  en médiane** (0,097 à 0,312 selon muscle et précontrainte, cf.
  `ANALYSE_PRELIMINAIRE.md` § A), la simulation immédiatement.
- Convention ε spire seule vs longueur totale (2.11) : à trancher en
  renseignant la longueur réellement étirée ; corrigerait la dérive du ratio
  F₀ avec ε.
- ~~Queue lente du spectre non résolue (> 600 s)~~ **fait le 01/09** : six
  « 30 L » de 1 800 s → τ₂ ≈ 1 200 s, E₀/ΣE ≈ 0,83 (borne haute, queue non
  close à 30 min ; deux essais écartés pour dérive suspecte — logger la
  température ambiante). Voir l'encadré en tête et § F de l'analyse.
- E10 : premier essai fait (28/08) — à-coups d'ancrage confirmés mais
  minoritaires (~7 %) ; refaire les photos caméra fixe, puis E2 (tube nu)
  pour attribuer la composante lisse.
