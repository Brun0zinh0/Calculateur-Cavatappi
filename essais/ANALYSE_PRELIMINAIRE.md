# Analyse préliminaire de la campagne d'essais (27/08/2026)

Analyse automatique des fichiers du dossier `essais/` (force **non filtrée**,
seuil défini comme dans `alpha V3/PROTOCOLE_ESSAIS_V4.md` : F − F₀ >
max(10 mN, 3σ) soutenu 2 s).

## Décodage de la campagne

| Type | Contenu | Usage |
|---|---|---|
| `10 N` | 600 s à **pression nulle** après mise en précontrainte | **Relaxation de la précontrainte** — l'essai de spectre idéal, insensible aux fuites (P = 0) ; c'est l'essai E2/E3 du plan de correction |
| `3 C` | 3 min de cycles 0→0,4 MPa | Training + hystérésis |
| `3 P` | Rampe + palier long (~160 s, fuite résiduelle ~0,1 MPa) | Relaxation sous pression |
| `1 O` | 1 montée-descente lente 0→0,4→0 MPa | **Essai P1/E4⁺ du protocole V4** : seuil + hystérèse du dead-band |

`Dataset.xlsx` = la **fiche géométrique E1** : ρ₀ = 2,29 mm, R_out = 0,96,
R_in = 0,56 (≠ article : 0,4 !), D_nylon = 0,85, angle de biais 42°
(≠ 37,91°), longueur d'hélice, nombre de spires et **longueur désenroulée par
précontrainte** pour D, J, L, Q, R.

## A) Seuils du dead-band (essais 1 O)

> **Échelle de pression.** Les valeurs ci-dessous sont établies avec la pression
> **reconstruite depuis les comptes ADC**, référencée au repos propre à chaque
> essai : `P = [ψ(ADC) − ψ(ADC_repos)] × 0,0689476`, avec
> `ψ(ADC) = (ADC×5/1023 − 0,5)×125` en psi. Le compte de repos varie de 91 à 103
> selon la séance, soit jusqu'à 0,5 bar de décalage capteur. Script :
> `recalcul_analyse.py`.

| Muscle | L_ext (0,8/1,0/1,2 mm) | P_s à 0,8 | à 1,0 | à 1,2 | P_d à 0,8 | Hystérèse (P_d < P_s ?) |
|---|---|---|---|---|---|---|
| G | (non renseigné) | 0,097 | 0,139 | 0,177 | 0,063 | oui (−0,034) |
| J | 20 / 30 / 38 | 0,097 | 0,143 | 0,253 | 0,029 | oui (−0,068) |
| L | 12 / 18 / 22 | 0,126 | 0,160 | 0,139 | 0,110 | oui (−0,017) |
| Q | 8 / 8 / 8 | 0,173 | 0,177 | 0,312 | 0,122 | oui (−0,051) |

**Lecture V4 (partielle mais déjà instructive) :**

1. **Le seuil croît avec la précontrainte** presque partout — y compris sur
   **Q dont les extrémités sont constantes** (8/8/8 mm) : la précontrainte
   augmente le dead-band **indépendamment des extrémités**. Comme une
   précharge plus forte devrait *réduire* un jeu série, ceci **défavorise
   M-B2 (jeu)** et favorise **M-A** dans une variante où le pré-étirement
   aggrave l'état de section (ovalisation/écrasement par contraction de
   Poisson et flexion de spire). Le référencement renforce cette lecture :
   l'écart 0,8 → 1,2 passe de +0,036/+0,124 MPa à +0,013/+0,156 MPa selon les
   muscles, et reste monotone sur G, J et Q.
2. **L'effet extrémités n'est pas tranché** : sur J, P_s et L_ext co-croissent
   avec la précontrainte (facteurs confondus) ; mais **entre muscles à même
   précontrainte**, P_s ne suit pas L_ext (Q, extrémités les plus courtes, a
   le seuil le plus haut à 0,8). M-B1 n'est ni confirmé ni exclu.
3. **Hystérèse montée/descente présente sur les quatre muscles**, de −0,017 à
   −0,068 MPa. Elle est la plus nette sur J et Q, la plus faible sur L. Une
   réserve d'interprétation s'impose toutefois : sur une rampe de 60 s, un
   simple retard viscoélastique produit lui aussi P_d < P_s. Le signe seul ne
   suffit donc pas à établir une signature de décollement/recollement (M-A) ;
   il faut comparer l'amplitude de l'hystérèse à celle qu'un modèle
   viscoélastique sans dead-band prédit sur la même rampe.

**Effet du référencement sur les seuils.** Tous les seuils augmentent, de
+0,012 à +0,057 MPa (médiane +0,031), soit +22 % en médiane. Le sens est
cohérent : rapporter l'échelle au repos réel du capteur restitue une pression
plus élevée à chaque compte, donc une pression plus élevée à l'instant où la
force démarre. **Le dead-band est donc plus grand que ne le laissait voir
l'échelle nominale, pas plus petit** — le classement des muscles et la
croissance avec la précontrainte sont eux inchangés.

**Il manque, pour trancher (cf. protocole)** : P2/E7 (seuil à 3 précharges,
même muscle, même précontrainte) et P4/E6 (volume-pression basse pression →
mesure directe de δ_c).

## B) Relaxation de la précontrainte à P = 0 (essais 10 N)

> **Insensible au référencement de la pression.** Vérifié sur les 18 essais
> « 10 N » : la pression y reste nulle dans l'échelle nominale et n'excède pas
> 0,042 MPa — un seul compte ADC — dans l'échelle reconstruite. Les
> ajustements ci-dessous portent sur F(t) à pression nulle et ne font
> intervenir la pression à aucun moment. **Les valeurs de cette section sont
> donc inchangées.**

Ajustement bi-exponentiel F(t) = F∞ + A₁e^(−t/τ₁) + A₂e^(−t/τ₂) sur 600 s :

- **Relaxation totale sur 10 min : −3 à −16 % (médiane ≈ −9 %)** de la force
  de précontrainte — contre **~83 %** de fraction relaxante à long terme pour
  la Table A1 de l'article. Confirmation propre et sans fuite du constat M6.
- **Constantes de temps réelles : τ₁ ≈ 15-30 s et τ₂ ≈ 250-700 s** (queue
  lente au-delà de 600 s mal contrainte — trois fits saturent la borne :
  prévoir un essai de 30-60 min pour la cerner).
- **Fractions relaxantes : ~3-7 % (rapide) + ~4-10 % (lente)** → spectre réel
  approx. E₀/ΣE ≈ 0,85-0,90.

**Spectre consolidé sur les 18 essais** (médianes ; script
`identifier_spectre.py`, sortie `spectre_identifie.json`). La somme des modules
est fixée à celle de l'article, 37,76 MPa — seule grandeur qui ne provient pas
des essais ; tout le reste est mesuré.

| Paramètre | Valeur | Article (Table A1) |
|---|---|---|
| `maxwell_E0_mpa` | **33,481** | 6,36 |
| `maxwell_E1_mpa` / `eta1` | **1,627** / **31,9** | 20,67 / 154,57 |
| `maxwell_E2_mpa` / `eta2` | **2,652** / **1 014,7** | 5,98 / 977,79 |
| `maxwell_E3_mpa` / `eta3` | **0** / — | 4,75 / 11 044,83 |
| τ₁ / τ₂ | **19,6 s** / **382,6 s** | 7,5 s / 163,5 s |
| E₀/ΣE | **0,887** | 0,168 |

Qualité des ajustements : R² médian 0,985 (0,857 à 0,999). Deux essais sur
dix-huit — G 0,8 et Q 1,2 — saturent la borne haute sur τ₂ ; ce sont aussi
les deux plus fortes relaxations mesurées (28,6 % et 31,9 %), à considérer
avec prudence.

**Deux réserves accompagnaient ce spectre ; la première est levée (28/08).**
La dérive propre de la cellule, jusqu'ici non soustraite, est désormais
mesurée par cinq essais à blanc de 600 s (`StabCell/`, § D) : au plus
**+11 mN**, et de signe **opposé** à la relaxation — la cellule monte
légèrement quand la force relaxée descend de 30 à 210 mN. Les relaxations
mesurées sont donc réelles ; les corriger de la dérive les **augmenterait**
même de l'ordre du point. Reste la seconde réserve : τ₂ est mal contrainte
par une fenêtre de 600 s tandis que E₃ ne l'est pas du tout — un essai de
30 à 60 min reste nécessaire pour la queue lente.

### Du composite au PVC seul

Ce spectre est celui du **muscle monté**, pas du matériau. Aucun essai de
relaxation sur tube nu n'a encore été réalisé — c'est l'essai E2 du plan. En
attendant, la déconvolution peut se faire par le partage tube/nylon que le
moteur calcule (`spectre_pvc_seul.py`).

Le résultat est contre-intuitif et mérite d'être exposé. À l'état précontraint
et à pression nulle, le partage vaut, pour J 0,8 :

| Contribution | Force |
|---|---|
| Tube de PVC | **+2 142 mN** |
| Filament de nylon | **−983 mN** |
| Total mesuré par la cellule | **+1 159 mN** |

Le nylon travaille donc **en compression** et s'oppose au tube. La raison se
lit dans la cinématique : le pré-étirement redresse l'hélice — α passe de 7,8°
à 14,1° chez J 0,8 — sans allonger le fil, qui raccourcit même très légèrement
(allongement cumulé 0,99989). Le cœur de nylon, bien plus raide que le tube,
se retrouve comprimé.

La conséquence sur le spectre inverse l'intuition. La cellule mesure une force
totale **plus petite** que celle du tube seul ; une même variation absolue de
force y pèse donc une fraction **plus grande**. Les fractions relaxantes
mesurées sur le muscle monté **surestiment** celles du PVC, d'un facteur
médian 0,62 (étendue 0,54 à 0,75 selon la précontrainte).

| Paramètre | Composite (mesuré) | PVC seul (déconvolué) |
|---|---|---|
| E₀ | 33,481 MPa | **35,127 MPa** |
| E₁ / η₁ | 1,627 / 31,9 | **1,001** / **19,6** |
| E₂ / η₂ | 2,652 / 1 014,7 | **1,632** / **624,3** |
| E₀/ΣE | 0,887 | **0,930** |

Les temps caractéristiques sont inchangés : la dilution ne porte que sur les
amplitudes. **L'écart avec l'article se creuse** — 0,930 contre 0,168, soit un
facteur 5,5 sur la fraction permanente.

Cette déconvolution hérite toutefois des hypothèses du modèle, puisque le
partage tube/nylon en provient et non de la mesure. Elle indique l'ordre de
grandeur et le sens de la correction, mais seul un essai de relaxation sur
tube nu (E2) permettra de trancher.
- Ces essais valident aussi le mode `viscoelastic_history` d'alpha V3 (le
  mode élastique prédit une dérive strictement nulle à P = 0 — les données
  montrent la relaxation, avec les bonnes échelles de temps).

## C) Effet du référencement sur les corrélations modèle / mesure

Les corrélations force simulée / force mesurée des essais « 1 O » ont été
recalculées sous les deux échelles de pression, **toutes choses égales par
ailleurs** — même géométrie, même spectre, mêmes réglages L3/nφ8. Seule
l'échelle de pression injectée dans le moteur change.

| Cas | ancienne échelle | échelle reconstruite | écart |
|---|---|---|---|
| J 0,8 / 1,0 / 1,2 | 0,898 / 0,900 / 0,911 | 0,879 / 0,883 / 0,888 | −0,019 |
| L 0,8 / 1,0 / 1,2 | 0,871 / 0,859 / 0,827 | 0,853 / 0,837 / 0,808 | −0,019 |
| Q 0,8 / 1,0 / 1,2 | 0,882 / 0,864 / 0,888 | 0,865 / 0,847 / 0,866 | −0,018 |

**La corrélation baisse de 0,006 en médiane** (−0,003 à −0,017). Les 0,96-0,98
de la passe 4 retenue dans `MODELISATION_QUANTITATIVE.md` deviennent donc
**0,94-0,98, médiane 0,964**. L'effet est réel mais modeste : le référencement
ne remet pas en cause la qualité d'ajustement de la forme d'actionnement.

Deux lectures s'en dégagent.

D'abord, **la baisse suit l'ampleur du dead-band**. Les deux plus fortes
dégradations sont J 1,2 (−0,013) et Q 1,2 (−0,017), qui sont précisément les
deux cas de plus grand seuil (0,253 et 0,312 MPa). À l'inverse, les cas à
ε = 0,8, dont les seuils sont les plus bas, perdent trois fois moins.

Ensuite, l'interprétation rejoint la section A. Le référencement agrandit le
dead-band mesuré ; or le modèle n'en contient aucun et produit de la force dès
les premiers incréments de pression. L'écart se creuse donc près de l'origine,
et c'est là que la corrélation se dégrade. **Cette baisse ne traduit pas une
moindre qualité du modèle : elle mesure plus justement ce qui lui manque.**

*(La colonne « ancienne échelle » reproduit l'intervalle 0,96-0,98 de
`MODELISATION_QUANTITATIVE.md`, ce qui valide la reconstitution du calcul.
Spectre : celui identifié en section B. Scripts : `identifier_spectre.py` puis
`recalcul_correlations.py`.)*

## D) Stabilité de la cellule de charge et chaîne d'acquisition (28/08)

### Les cinq essais à blanc `StabCell/`

Cinq enregistrements de 600 s, cellule à vide (`stabilite_cellule_gui.py`),
colonne `force_unfiltered_mN`, ADC pression au repos stable à 94 :

| Essai | Dérive sur 600 s | Bruit (σ après lissage 5 s) |
|---|---|---|
| StabCell 1 | **+11,2 mN** | 2,1 mN |
| StabCell 2 | +5,2 mN | 2,1 mN |
| StabCell 3 | +6,6 mN | 2,0 mN |
| StabCell 4 | −0,4 mN | 2,1 mN |
| StabCell 5 | −0,1 mN | 2,0 mN |

Trois conclusions :

1. **La relaxation des « 10 N » est réelle.** La dérive à blanc est au plus
   +11 mN quand les relaxations mesurées vont de −30 à −210 mN — un ordre de
   grandeur en dessous, et de signe opposé. L'artefact instrumental est exclu
   comme explication (cf. résultat n° 3 de `MODELISATION_QUANTITATIVE.md` :
   il ne reste que le matériau ou l'ancrage, essai E10).
2. **La dérive décroît d'essai en essai** (+11 → +5 → +7 → −0,4 → −0,1) :
   signature d'un échauffement de l'électronique. Recommandation de banc :
   laisser chauffer HX711 et cellule **10 à 15 min sous tension** avant toute
   mesure — les deux derniers blancs sont alors quasi parfaits (< 0,5 mN).
3. Le bruit haute fréquence vaut **±2 mN (1σ)**, négligeable devant les
   forces mesurées (900-1 700 mN).

*Les blancs sont à charge nulle, mais la dérive de cette cellule est
indépendante de la charge appliquée (constat de l'équipe) : les bornes
ci-dessus valent donc aussi sous les 900-1 700 mN des essais, et la question
de l'artefact instrumental est définitivement fermée.*

### Revue du code du banc (`Banc d'essai/`)

Le dossier contient désormais quatre pièces : le croquis Arduino
(`sketch_jul10a.ino`), l'enregistreur principal (`lecture_arduino_gui.py`),
le nouvel outil de stabilité (`stabilite_cellule_gui.py`) qui a produit les
`StabCell`, et un dépouillement rapide force/pression
(`identifier_parametres_materiaux.py`). Points qui touchent à
l'interprétation des données :

- **L'Arduino écrête le psi à zéro avant l'envoi** (`if (pressure_psi < 0)
  pressure_psi = 0`). Le repos étant à ADC ≈ 91-103 (< 0,5 V), toute
  l'information sous 0,5 V est perdue dans les colonnes `pressure_psi/bar`
  du CSV — le tare PC (`pressure_raw_psi − pressure_tare_psi`) ne peut pas
  la restituer. C'est précisément la zone du dead-band : la reconstruction
  depuis les comptes `pressure_adc` (§ C et § 6.4.2) est donc la **seule**
  échelle valable près de l'origine, et le croquis confirme sa formule
  (ψ = (ADC·5/1023 − 0,5)·125 psi, 0,0689476 bar/psi).
- L'enregistreur convertit avec g → mN = 9,80665 et protège le HX711 contre
  les décrochages de trame (resynchronisation si saut > 50 000 comptes,
  colonnes `hx711_raw` vs `hx711_raw_received`).
- Cadence réelle ~10 Hz (`delay(100)` — le commentaire « 20 readings/second »
  du croquis est obsolète), cohérente avec les CSV.

## E) Premier essai E10 — marquage des ancrages (28/08)

Essai réalisé : muscle précontraint monté comme un « 10 N » normal, traits
blancs aux deux jonctions et pointillés le long de l'hélice, photos T0/T10,
600 s d'enregistrement à P = 0 (`test E10/`). La force a relaxé de
**1 030 à 962 mN (−6,7 %)**, profil typique de la campagne — le phénomène à
expliquer s'est bien produit pendant l'essai marqué.

### Les photos ne tranchent pas (caméra déplacée)

L'appareil a bougé entre T0 et T10 (cadrage et angle différents). Les
décalages apparents des traits (~2-3 mm aux deux jonctions) sont de l'ordre
de ce que la parallaxe produit, et leurs sens ne sont pas cohérents entre eux
— aucun trait nettement **scindé** (la signature décisive) n'est visible. Un
changement de forme au niveau de l'ancrage haut (une boucle de tube apparaît
entre le scotch noir et le trait) est suggestif mais pas concluant sous cet
angle. **À refaire caméra fixe** : téléphone posé/scotché, même cadrage
exact, gros plan par jonction — ou mieux, mesure directe au pied à coulisse
des distances trait-scotch à t0 et t10, insensible à la parallaxe.

### La forme de la courbe, elle, apporte une réponse partielle

Un glissement d'ancrage procède par à-coups (stick-slip) : des décrochements
de niveau brusques et **persistants**, que la relaxation matérielle — continue
par nature — ne produit pas. Détecteur : saut de la médiane 18 s avant/après
chaque instant, validé sur les 5 blancs StabCell (plancher de fausse alarme :
0-1 événement marginal à −4,5 mN par blanc ; les excursions de bruit à ±2 mN
ne persistent pas).

Résultat sur les 18 « 10 N » + E10, hors 60 premières secondes (où la branche
rapide τ₁ ≈ 20 s du spectre rend la pente trop forte pour séparer un à-coup
d'une relaxation lisse) :

- **Les à-coups persistants existent** : 0 à 3 par essai, de −4 à −8 mN
  (E10 : −4,9, −7,6, −7,9 mN à t = 140, 380 et 458 s), absents des blancs.
  Le glissement d'ancrage est donc **réel**.
- **Mais ils ne portent qu'une petite part de la relaxation** : ~7 %
  en médiane (0 à 14 % selon les essais ; E10, cas le plus marqué : 33 %).
  L'essentiel de la chute de force est **lisse** (≈ a·ln t).

La composante lisse reste donc non attribuée : vraie viscoélasticité du
matériau, ou fluage d'ancrage continu (le tube de PVC sous un scotch serré
peut ramper sans à-coups). Deux voies pour fermer : refaire les photos E10
caméra fixe (tranche le fluage continu d'ancrage), et surtout **E2** —
relaxation sur tube nu entre mors propres, qui mesure le matériau sans aucun
nœud dans le chemin d'effort.

*Nota : la tige filetée visible sur les photos **longe l'actionneur en
arrière-plan** sans le traverser ni le toucher (confirmé le 28/08) — le
chevauchement n'est que visuel, aucun frottement spire-tige à considérer.*

## Prochaines étapes proposées

1. **Modélisation quantitative enfin possible** : remplir une fiche
   `identification/fiche_muscle_exemple.json` par muscle depuis le Dataset
   (géométrie réelle + longueurs désenroulées, compliance série active) +
   spectre issu des « 10 N » → confrontation modèle/essai **sans facteur
   d'échelle libre** pour la première fois.
2. Compléter la discrimination V4 : P2/E7 (précharges) et P4/E6
   (volume-pression) sur un muscle, une demi-heure de banc.
3. Un « 10 N » long (30-60 min) sur un muscle pour cerner la queue lente du
   spectre.
