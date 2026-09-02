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
> 0,0084 MPa — deux comptes ADC, le quantum valant 0,0042 MPa soit 0,042 bar —
> dans l'échelle reconstruite. Les
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
même de l'ordre du point. La seconde réserve — τ₂ mal contrainte
par une fenêtre de 600 s, E₃ pas du tout — est levée à son tour par les six
essais longs « 30 L » du 01/09 : voir § F, qui révise τ₂ à ≈ 1 200 s et
E₀/ΣE à ≈ 0,83 (borne haute).

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

> **Attention aux versions.** Le croquis a été modifié le 28/08 sur deux
> points qui changent la lecture des CSV. **Tous les fichiers de la campagne
> — `Muscle *`, `StabCell`, `test E10` — ont été acquis avec la version
> antérieure** et doivent être dépouillés comme indiqué ci-dessous.

**Écrêtage du psi — retiré du croquis, mais présent dans les données déjà
acquises.** La version employée pour la campagne appliquait
`if (pressure_psi < 0) pressure_psi = 0` avant l'envoi. Le repos étant à
ADC ≈ 91-103 (< 0,5 V), toute l'information sous 0,5 V est perdue dans les
colonnes `pressure_psi/bar` de ces fichiers, et le tare PC
(`pressure_raw_psi − pressure_tare_psi`) ne peut pas la restituer. C'est
précisément la zone du dead-band : **pour les fichiers existants, la
reconstruction depuis les comptes `pressure_adc` (§ A, § C, § 6.4.2 du
rapport) est la seule échelle valable près de l'origine**. Le croquis actuel
transmet la valeur négative telle quelle et laisse la correction au tare ;
les acquisitions futures n'auront donc plus besoin de cette reconstruction,
mais elle reste exacte dans les deux cas et peut rester la méthode par défaut.
La formule est confirmée par le croquis : ψ = (ADC·5/1023 − 0,5)·125 psi,
0,0689476 bar/psi.

**Cadence et horodatage — également modifiés.** La version de campagne émettait
sur minuterie (`delay(100)`, soit ~10 Hz mesurés ; le commentaire
« 20 readings/second » du croquis était obsolète) et les échantillons étaient
datés par le PC à réception, d'où une gigue de 50 à 133 ms. Le croquis actuel
n'émet que lorsque le HX711 signale une conversion prête — ce qui supprime
aussi les zéros de trame — et joint la date `millis()` relevée juste après la
lecture de la cellule. Les CSV produits désormais portent deux bases de temps,
`time_s` (Arduino) et `time_pc_s`, dont l'écart mesure la latence de liaison.
**Conséquence pour les fichiers de la campagne** : leur colonne `time_s` est
datée PC, donc affectée par la gigue. Aucune incidence sur les relations
force/pression, qui viennent de la même trame ; restait à vérifier l'effet sur
les grandeurs construites sur le temps.

*Quantifié le 28/08 (`effet_gigue.py`).* La gigue réelle des 18 essais
« 10 N » vaut σ(Δt) = **31 ms** pour une période de 127 ms, soit σ = 22 ms sur
chaque date (l'horloge PC ne dérivant pas, les erreurs ne s'accumulent pas :
Var(Δt) = 2·Var(ε)). Propagée par Monte-Carlo — 40 regrilles perturbées par
essai, réajustement complet à chaque fois — elle donne :

| | Écart-type relatif médian | Maximum |
|---|---|---|
| τ₁ | **2,0 × 10⁻⁴** | 3,8 × 10⁻⁴ |
| τ₂ | **1,1 × 10⁻⁴** | 1,7 × 10⁻³ |

Soit **0,02 % sur τ₁ et 0,01 % sur τ₂**, à comparer à la dispersion entre
spécimens : un facteur 2,7 sur τ₁ (12 à 32 s) et un facteur 15 sur τ₂ (108 à
4 000 s). La gigue d'horodatage est donc sans effet mesurable sur le spectre
identifié, et les fichiers datés PC restent pleinement exploitables. La
réserve est levée.

- L'enregistreur convertit avec g → mN = 9,80665 et protège le HX711 contre
  les décrochages de trame (resynchronisation sur 5 valeurs concordantes si
  saut > 50 000 comptes, colonnes `hx711_raw` vs `hx711_raw_received`).
- Le parseur accepte les trames à 5, 6 ou 7 champs : les fichiers des trois
  versions de croquis restent donc lisibles par le même outil.

## E) Premier essai E10 — marquage des ancrages (28/08)

Essai réalisé : muscle précontraint monté comme un « 10 N » normal, traits
blancs aux deux jonctions et pointillés le long de l'hélice, photos T0/T10,
600 s d'enregistrement à P = 0 (`test E10/`). La force a relaxé de
**1 030 à 962 mN (−6,7 %)**, profil typique de la campagne — le phénomène à
expliquer s'est bien produit pendant l'essai marqué. *(Essai identifié
depuis : muscle **G** — fichier `Muscle G E10.csv`, cohérent avec le « G »
inscrit sur le scotch de la photo T0.)*

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

## F) Essais longs « 30 L » — la queue lente enfin contrainte (01/09)

Six relaxations de 1 800 s à P = 0, ε = 1,0 (« L = 30 min N » dans le
Dataset) : muscles A, B, C, I, K, O. Dépouillement par médiane seconde par
seconde, ajustement à deux exponentielles, détecteur d'à-coups persistants
validé sur les blancs (§ E).

| Muscle | F₀ (mN) | à 600 s | à 1 800 s | τ₁ (s) | τ₂ (s) | fractions f₁/f₂ | E₀/ΣE | R² |
|---|---|---|---|---|---|---|---|---|
| C | 1 242 | −12,0 % | −14,8 % | 53 | 1 271 | 0,089 / 0,076 | 0,835 | 0,992 |
| I | 1 333 | −13,1 % | −16,2 % | 9 | 484 | 0,054 / 0,133 | 0,813 | 0,983 |
| K | 1 324 | −8,1 % | −12,8 % | 60 | 1 668 | 0,040 / 0,130 | 0,830 | 0,997 |
| O | 1 372 | −11,3 % | −15,1 % | 43 | 1 157 | 0,073 / 0,100 | 0,827 | 0,997 |
| **médianes** | | | | **48** | **1 214** | **0,064 / 0,115** | **0,829** | |
| A *(écarté)* | 1 534 | −11,6 % | −10,2 % | — | — | remonte après 600 s | — | 0,84 |
| B *(écarté)* | 864 | −3,5 % | **+2,8 %** | — | — | dépasse son niveau initial | — | 0,01 |

Cinq enseignements :

1. **La relaxation ne s'arrête pas à 600 s** : −8 à −13 % à 600 s, −13 à
   −16 % à 1 800 s sur les quatre essais propres.
2. **τ₂ vaut ≈ 1 200 s, pas ≈ 400 s.** L'estimation de la fenêtre 600 s
   (382,6 s) était biaisée vers le bas, mécaniquement : on ne peut pas voir
   une constante de temps plus longue que la fenêtre. τ₁ (10-60 s) est
   confirmée.
3. **E₀/ΣE descend de 0,887 à ≈ 0,83** — et reste une **borne haute** : la
   pente finale en log t (−21 à −60 mN/décade) montre que la queue n'est pas
   close à 30 min ; une branche τ ≳ 5 000 s peut exister (les ajustements à
   trois exponentielles la poursuivent mais extrapolent hors fenêtre, non
   retenus). L'écart à l'article (0,168) demeure : facteur ~5.
4. **Les à-coups d'ancrage restent minoritaires sur 30 min** : −4 à −43 mN
   sur des chutes totales de 150 à 215 mN (médiane ≈ 8 % ; I est le plus
   touché avec un événement isolé de −20 mN à t = 418 s).
5. **Deux essais sur six sont écartés** : A remonte après 600 s et B finit
   **au-dessus** de son niveau initial (+25 mN) — une force qui remonte à
   longueur bloquée et P = 0 n'est pas de la viscoélasticité ; l'ampleur
   (~50 mN crête à crête) dépasse la dérive de chauffe des blancs (+11 mN).
   Suspect n° 1 : la température ambiante sur 30 min (B est aussi un muscle
   « non fonctionnel » à extrémités courtes). Recommandation : noter ou
   logger la température de la pièce pour les essais longs.

**Spectre fenêtre 30 min** (ε = 1,0, muscles C/I/K/O, ΣE = 37,76 MPa
maintenu) : E₀ ≈ **31,3**, E₁ ≈ **2,4** / η₁ ≈ **115** (τ₁ 48 s),
E₂ ≈ **4,3** / η₂ ≈ **5 270** (τ₂ 1 214 s). Pour l'actionnement rien ne
change (branche rapide inchangée) ; pour un maintien ≥ 10 min, c'est ce
spectre qu'il faut mettre dans les champs maxwell_* de l'interface.

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
