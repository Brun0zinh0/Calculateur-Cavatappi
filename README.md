# Interface Cavatappi Alpha V3

Ce dossier contient tout le nécessaire pour lancer l'interface Streamlit du modèle Cavatappi.

Alpha V3 est la copie consolidée d'Alpha V2 (moteur identique,
`2026.07.30-series-compliant-ends-10`) issue de l'audit scientifique complet du
17-18/08/2026. Elle ajoute : le dossier `validation/` (validation quantitative
contre la figure 7 de l'article, avec baseline de non-régression), le dossier
`audit/` (rapport et résultats bruts de l'audit), et `PLAN_DE_CORRECTION.md`
(feuille de route des corrections, phases 0-4). Le présent README a été corrigé
en Phase 1 du plan : les paragraphes marqués « (audit 2026-08) » documentent
les écarts réels entre le logiciel et les articles.

## Contenu du dossier

- `Base.py` : moteur scientifique du modèle.
- `parametres.py` : paramètres, géométrie, pression et cache.
- `pression.py` : lecture des historiques CSV mesurés et conversion des unités.
- `affichage.py` : visualiseur Cavatappi 3D interactif avec perspective, éclairage de profondeur, grille et rotation automatique, ainsi que les graphes Matplotlib.
- `parallel.py` : exécution parallèle des études comparatives.
- `timing.py` : estimation et suivi des temps de calcul.
- `interface.py` : application Streamlit.
- `requirements.txt` : dépendances Python recommandées.
- `installer_dependances.bat` : installateur des dépendances sous Windows.
- `lancer_interface.bat` : lancement de l'interface sous Windows.
- `test_installation.py` : vérification rapide des imports et du modèle.
- `test_scientifique.py` : tests numériques des équilibres, historiques de pression et du mode suspendu, plus le test de non-régression contre la figure 7 (baseline ±2 %).
- `validation/` : validation quantitative contre la figure 7 de l'article (voir `validation/README.md`).
- `audit/` : rapport d'audit 2026-08, constats bruts, baseline de non-régression, scripts de contre-expertise.
- `PLAN_DE_CORRECTION.md` : feuille de route des corrections issues de l'audit.

## Installation sur un PC Windows

1. Copier tout le dossier `alpha V3` sur le PC.
2. Installer Python si nécessaire.
3. Double-cliquer sur `installer_dependances.bat`.
4. Double-cliquer sur `lancer_interface.bat`.

L'interface s'ouvre ensuite dans le navigateur par défaut du PC. Le lanceur
choisit automatiquement un port local libre : plusieurs versions de
l'application peuvent ainsi rester ouvertes sans afficher la mauvaise
interface. L'adresse utilisée est aussi indiquée dans la fenêtre de commande.


## Options d'hystérèse

L'onglet `Hystérèse` peut afficher plusieurs cycles sur le même graphe. Dans
la barre latérale, champ `Cycles à afficher`, entrer par exemple :

```text
1, 2
```

Le même onglet peut aussi comparer :

- plusieurs précontraintes initiales ;
- plusieurs vitesses de pression injectée en `MPa/s`.

Pour les vitesses de pression injectée, l'historique de pression utilisé pour
la comparaison est triangulaire et linéaire afin que la pente corresponde à la
valeur demandée.

## Export des résultats

Chaque onglet de simulation propose un bouton `Exporter les résultats en CSV`.
Le fichier utilise l'encodage UTF-8 et le séparateur `;` pour faciliter son
ouverture dans un tableur en environnement français. L'export d'hystérèse
indique le cas comparé, le numéro du cycle et le temps relatif dans le cycle.

## Affichage d'un essai expérimental

Dans l'onglet des courbes temporelles, le volet `Afficher un essai expérimental
CSV` permet d'importer un relevé contenant le temps, la pression et la force.
Les colonnes usuelles sont reconnues automatiquement, puis peuvent être
sélectionnées manuellement, avec conversion de la pression en MPa et de la
force en mN.

Deux affichages selon la source de pression de la simulation :

- si la pression injectée est un **historique mesuré CSV** (volet `Pression et
  actionnement`) et que l'option `Superposer l'essai expérimental sur le
  graphe principal` (barre latérale, activée par défaut) est cochée,
  simulation et essai partagent la même base de temps : la superposition se
  fait sur **un seul graphe** (force simulée et force mesurée sur l'axe
  principal, pression mesurée injectée sur l'axe secondaire ; le couple n'est
  pas affiché dans cette vue) ;
- sinon, le graphique expérimental est affiché séparément des courbes
  simulées, comme auparavant.

Le volet latéral `Pression et actionnement` permet aussi de masquer le couple
et de superposer la pression à la force sur un second axe vertical.

## Dependances Python

Les versions utilisées pour cette interface sont listées dans `requirements.txt`.
L'interface a été vérifiée avec Python 3.14 sous Windows. Python 3.11 ou une
version plus récente est recommandée, sous réserve de la disponibilité des
versions de NumPy et SciPy indiquées.
Pour une installation manuelle :

```powershell
python -m pip install -r requirements.txt
```

## Portabilite

Le dossier `alpha V2` ne dépend pas du répertoire courant. Les chemins sont résolus à
partir de l'emplacement de `interface.py`.

L'interface, les calculs parallèles et les tests importent directement le
`Base.py` placé dans le dossier `alpha V2`. Le lanceur se place automatiquement
dans ce dossier avant de démarrer l'application.

## Vérification de l'installation

Depuis le dossier `alpha V2`, exécuter :

```powershell
python test_installation.py
python test_scientifique.py
```

Le second test contrôle notamment la convergence de l'équilibre bloqué, la
décomposition force/couple, la durée exacte des profils de pression et
l'équation de déplacement du mode masse suspendue. Il exécute aussi le test de
non-régression contre la figure 7 (baseline ±2 %, ~2-3 minutes ;
`TCPA_SKIP_VALIDATION=1` pour le sauter, `TCPA_SLOW_VALIDATION=1` pour ajouter
la discrétisation de production).

## Extrémités désenroulées

Le paramètre `Longueur totale désenroulée aux extrémités` représente la somme
des portions non hélicoïdales situées aux deux bouts de l'actionneur. La
`Longueur hélicoïdale active initiale` reste la longueur utilisée pour calculer
le nombre de spires et la force produite par la géométrie hélicoïdale.

Les extrémités sont inactives du point de vue de la génération de force, mais
elles sont élastiques. Leur section composite associe le tube et le nylon liés :
leurs rigidités `EA` et `EI` sont calculées directement avec les modules et les
dimensions saisis.

En actionnement bloqué, le solveur impose la compatibilité :

```text
variation de longueur de la spire
+ variation de longueur des extrémités
= 0
```

La déformation des extrémités vaut `ΔF / k_ext`. Elle autorise donc un léger
mouvement de la partie hélicoïdale malgré le blocage global et réduit la force
transmise. Aucun coefficient multiplicatif de force n'est appliqué.

Le mode `Traction et flexion au raccord` considère deux portions de même
longueur quittant la spire suivant sa tangente. Il additionne leur compliance
axiale et leur compliance de flexion de poutre. Le mode `Traction axiale
uniquement` représente des extrémités parfaitement alignées avec l'axe et
produit une raideur beaucoup plus élevée.

Conventions précisées par l'audit 2026-08 (Phase 2) :

- la compliance des extrémités est évaluée à l'angle d'hélice **après
  précontrainte** (α_tk), au moment du verrouillage de la référence série —
  c'est autour de cet état que la compatibilité est linéarisée (item 2.5 ;
  auparavant elle restait figée à α₀, écart ~8 %) ;
- pour les longueurs exportées, les extrémités sont comptées à **pleine
  longueur sur l'axe** (« longueur entre mors »), y compris en mode tangente
  où leur raideur les idéalise quasi transversales — convention d'affichage
  assumée (item 2.6) ;
- la **précontrainte s'applique à la spire seule** : `eps` est la
  pré-déformation de la partie hélicoïdale active, sans compatibilité série
  pendant l'étirement. La correspondance avec le ε_tk expérimental (mesuré
  entre mors) n'est donc exacte que pour des extrémités rigides (item 2.11) ;
- la référence série est verrouillée sur l'état **P = 0** post-précontrainte,
  même si un historique CSV démarre à P(0) > 0 (item 2.3).

## Contraction avec une masse suspendue

La contraction est calculée directement par rapport à la longueur de l'actionneur
au début de l'essai :

```text
contraction(t) = L(t = 0) - L(t)
actionnement(t) = 100 [L(t = 0) - L(t)] / L(t = 0)
```

La valeur reste signée. Une valeur positive indique une contraction et une
valeur négative signifie que l'actionneur est plus long qu'à `t = 0`. Aucun
témoin numérique ni aucune valeur absolue ne sont utilisés.

Par défaut, la masse est appliquée avant le début de l'essai et l'état de Maxwell
est amené à son équilibre à `0 MPa`. L'horloge est ensuite remise à zéro avant la
montée en pression. L'option peut être désactivée pour étudier explicitement la
récupération de la précontrainte immédiatement après la mise en charge, mais ce
transitoire ne doit alors pas être interprété comme une relaxation due à la pression.

(audit 2026-08) **Pour reproduire le protocole d'actionnement libre de
l'article EXP** (figures 11 et 15-17) : mettre la précontrainte `eps` à `0`
(le protocole EXP n'en comporte pas), **désactiver** l'équilibrage préalable
sous la masse (les essais de l'article incluent le fluage sous poids et
n'atteignent jamais l'équilibre asymptotique), et normaliser l'actionnement
par la longueur initiale non chargée `L_T0` (éq. 25 d'EXP) plutôt que par la
`L(0)` chargée/équilibrée utilisée par défaut ici — l'écart entre les deux
normalisations est un facteur ≈ `L(0)/L_T0` (2,11 % contre 3,71 % sur le cas
type de l'audit). **Depuis la Phase 4 (décision D3), la normalisation par
`L_T0` est celle appliquée par défaut** : `free_actuation_percent` divise la
contraction par la longueur naturelle fabriquée de l'actionneur (longueur
hélicoïdale active initiale + extrémités, sans charge ni précontrainte —
c'est le `L_T0` de l'éq. 25 d'EXP, dont le protocole n'a pas de
précontrainte) ; l'ancienne normalisation par la référence chargée reste
exportée sous `free_actuation_percent_loaded_ref`, et la longueur de
référence non chargée est exportée en `reference_unloaded_length_mm`.

## Domaine d'utilisation

- La pression maximale proposée par l'interface est limitée à 1,5 MPa, domaine expérimental de l'article.
- L'intégration exponentielle est le choix recommandé et utilisé par défaut.
- (audit 2026-08) Le profil radial de l'angle de biais par défaut (`paper_linear`, θ linéaire en rayon) suit la phrase descriptive de la section 3.3 de l'article (« varie linéairement ») mais **pas son équation (11)**, qui prescrit `θ_j = arctan((R_j/R_out)·tan θ_f)` — le profil exact d'un tube uniformément torsadé, reprise à l'identique dans l'article compagnon EXP (éq. 7). Ce profil conforme à l'équation est disponible via l'option `uniform_twist`. L'article se contredit entre son texte et son équation ; l'écart mesuré entre les deux profils est d'environ 5-6 % sur le couple crête et 7 % sur l'amplitude de force d'actionnement (décision D1 du plan : défaut conservé, écart documenté).
- (audit 2026-08) Par défaut (`section_update_mode = fixed`), la section du tube — rayons, angles de biais et rigidités par couche — reste **figée** pendant toute la simulation ; seule la géométrie hélicoïdale est mise à jour. C'est un écart au schéma lagrangien réactualisé de l'article, qui re-référence la configuration complète à chaque pas (`r = R + Δu`, éq. 16-17), avec un déplacement radial atteignant ~20 % sous 1,3 MPa. **L'option `updated` est désormais corrigée et validée** (item 3.1, correcteur de point milieu du BVP radial) : l'ancienne dérive de cyclage — un ratchet d'intégration explicite du premier ordre, vérifié O(h), qui contractait la section d'environ −0,3 %/cycle — est réduite d'un facteur ~25 (fermeture de cycle quasi élastique ~0,01 %/cycle, testée par la suite). Sur la figure 7, `updated` corrigé ramène le couple de −33 % à **−6/−12 %** des ancrages expérimentaux et la force à **±1 %** (ε = 0,8) / **−8/−9 %** (ε = 1,0) de la théorie. Le défaut reste `fixed` (référence de la baseline) ; utiliser `updated` pour les études de fidélité. Coût : ~2× sur la partie BVP ; dérive résiduelle O(h) documentée.
- Le module axial par défaut est la somme `E0 + E1 + E2 + E3`, comme indiqué dans l'annexe A.
- (audit 2026-08) La précontrainte construit un état élastique de référence `t_k` **avec des branches de Maxwell vierges** (σᵢ = 0), là où l'article applique la loi viscoélastique dès la phase d'élongation (loi A.4 sur les deux phases ; relaxation post-étirement de la fig. 11 d'EXP). Conséquences mesurées : niveaux absolus de force environ −18 à −22 % en suivant le schéma de l'article ; sur l'atténuation de « training » des 10 premiers cycles (−10,3 % dans l'article), ce schéma n'en reproduit que ~17 % (un schéma viscoélastique en récupère ~2/3) ; le dead-band d'origine fluage et la fig. 11 d'EXP ne sont pas simulables. Les amplitudes pic-vallée d'actionnement restent quasi insensibles (< 5 %). L'option de précontrainte viscoélastique est **disponible** depuis l'item 3.2 (voir « Comparaison avec la figure 7 »).
- L'anisotropie identique de l'article reste le mode par défaut. Le mode `relaxation axiale identifiée` limite les branches Maxwell à la direction caractérisée par l'essai de traction.
- Le profil débit/volume linéaire est utilisé par défaut. Le profil non linéaire est phénoménologique. Depuis la Phase 2 de l'audit 2026-08, ce défaut linéaire vaut aussi pour les fonctions de l'API `Base` (`cyclic_pressure_history`, `ramp_hold_pressure_history`, `run_hold_relaxation`) — auparavant ces points d'entrée imposaient silencieusement le profil non linéaire γ = 3,5 ; les exposants sont centralisés dans `Base.NONLINEAR_GAMMA_LOAD/UNLOAD`.
- Les coefficients du nylon doivent rester égaux à 1 pour reproduire sa loi élastique linéaire complète.
- (audit 2026-08, item 3.4) L'appariement des coefficients de Poisson effectifs suit par défaut la **lettre des articles** (`poisson_pairing = paper_crossed` : ν̄(s→φ) appliqué sur ε_r et ν̄(s→r) sur ε_φ — écriture probablement coquillée, identique dans BLOCKED éq. 18-19 et EXP éq. 19-20). L'option `physical` échange l'appariement ; impact mesuré ≤ 1,5 % sur le couple, < 1 % sur la force. Question à poser aux auteurs pour trancher définitivement.
- Les résultats servent à l'étude et à la comparaison. Un dimensionnement de sécurité nécessite une validation expérimentale du spécimen réel.


## Alpha V4 — mécanismes physiques optionnels (audit → rapport → contre-expertise)

Moteur `2026.09.02-v4-15`. Six améliorations issues de la confrontation aux
essais (chapitre 7 du rapport de stage) ont été proposées sous trois angles,
contre-expertisées de façon adversariale (8 retenues sur 12), puis
implémentées. **Toutes sont inactives par défaut** : avec les valeurs par
défaut, le moteur est bit-identique à l'alpha V3 (vérifié sur les modes
bloqué, série, section réactualisée, suspendu ; suite de 33 tests, baseline
figure 7 intacte). Elles se règlent dans l'expander « Mécanismes Alpha V4 »
des réglages avancés.

| Mécanisme | Réglage (0 = off) | Écart visé (rapport) | Ce que ça change |
|---|---|---|---|
| **V4-1 Pression d'engagement** (reformage de la section ovalisée, § 7.10.3) | `engagement_ovality_e0`, `engagement_ring_factor`, `engagement_unload_ratio` | Seuil de démarrage 1,7 bar vs 0 simulé ; super-linéarité sous 2 bar | Ovalité e comme variable d'état ; P_eff = P(1 − e/e0) — démarrage quadratique P²/P_r0, pleine pression une fois la section ronde (gain conservé). P_r0 = k·E_radius·(t/R_m)³·e0 sur la **section** du tube. r < 1 : hystérésis du seuil. |
| **V4-2 Frottement sec** (élément de Jenkins, § 7.5) | `friction_pressure_coulomb_mpa` | Aire d'hystérésis +62 mN·bar mesurée vs −6 simulée ; seuil de descente négatif (§ 7.5.4) | Élément de Coulomb sur la **transmission de la pression** : P_eff = P − P_f, P_f écrêté à ±P_c. Retard de P_c en charge, avance en décharge : descente au-dessus de la montée, force résiduelle à P = 0, indépendant de la vitesse. Établi par expérience numérique : en mode bloqué, tout patin interne (dw, dv, dκ, ±) est ré-absorbé par l'équilibre géométrique et n'ouvre aucune boucle. Calibration : hystérésis du seuil 0,038 MPa ≈ 2·P_c. |
| **V4-3 Convention de pré-étirement** (item 2.11) | `prestretch_convention = grip_to_grip` | Dérive de F0 avec ε (×1,26 → ×1,76) | ε appliqué à la longueur entre mors ; les extrémités désenroulées s'allongent en série pendant l'étirement (résolution à deux inconnues par incrément). Sans extrémités : identique à `coil_only`. Peut **augmenter** F0 pour des extrémités courtes (physiquement correct). |
| **V4-4 Viscosité d'Eyring** (§ 7.8.7) | `eyring_sigma_star_mpa` | Incompatibilité d'amplitude (relaxation ×0,07) | η_eff = η·(s/σ*)/sinh(s/σ*) par couche et par branche, s = norme de la contrainte de branche. Attention : dans ce modèle les contraintes de branche du tube valent ~0,01-0,05 MPa à la précontrainte (la grande déformation du muscle est géométrique) — σ* doit être de cet ordre. Exige l'intégration exponentielle. |
| **V4-5 Fluage d'ancrage** (§ 7.10.4, essai E10) | `anchor_creep_c_mm`, `anchor_creep_t0_s` | Part d'ancrage de la relaxation (≥ 7 %) | δ(t) = c·ln(1 + t/t0) depuis le blocage, en série dans la longueur bloquée (avec ou sans compliance des extrémités). |
| **V4-6 Identification en boucle fermée** | `identification/identifier_spectre_moteur.py` | Spectre 7.6 dilué ×14 par la chaîne | Ajuste {E_i, τ_i} pour que la force **simulée** en maintien à P = 0 reproduise la mesure (forme normalisée, ΣE maintenu). |

Rejetée à l'implémentation : le « jeu radial tube-nylon » (retenu par la
contre-expertise) — sous pression le rayon interne augmente, le jeu
s'ouvrirait au lieu de se fermer ; mécanisme incohérent avec le sens de la
déformation radiale du moteur.

Première démonstration sur un essai réel (muscle J, rampe lente « 1 O »,
ε = 1,0, géométrie du Dataset, E_nylon 2 065 MPa) : avec e0 = 0,20,
P_c = 0,02 MPa et r = 0,7, la corrélation de forme passe de 0,979 à 0,987
sans perte de gain (0,52 contre 0,44 pour un simple décalage de pression), F0
inchangé à 0,3 % près. Le déficit de gain résiduel (~×0,5) n'est pas
touché : il relève de la pente pression→force du modèle (rapport § 7.9.2),
pas du dead-band.

Calibration attendue : e0(ε) et k par imagerie de la section sous pression
(essai n° 2 du tableau 10.3) ; σ* conjointement avec le spectre sur les essais « 10 N » et l'essai E2 (tube
nu) ; P_c sur l'hystérésis du seuil et l'aire de boucle des rampes lentes ;
c, t0 sur la part à-coups + continue des ancrages (E10, E2).

## Cache des réglages

L'application sauvegarde automatiquement les derniers réglages et les derniers
résultats dans le dossier temporaire Windows. Ce cache n'est pas nécessaire :
s'il n'existe pas sur un autre PC, l'application repart avec les valeurs par
défaut.

Dans l'interface, le bouton `Paramètres par défaut` réinitialise les paramètres
de simulation et efface les derniers résultats mis en cache.

## Transmission de la pression

La pression déforme le tube et modifie la géométrie de l'hélice. Lorsque
l'actionneur est bloqué, cette tendance au mouvement produit la force calculée
par le modèle.

(audit 2026-08) Avec le réglage par défaut `section_update_mode = fixed`, les
directions matérielles anisotropes et la section du tube ne sont **pas**
réorientées pendant la simulation — seule l'option `updated` réoriente les
fibres à chaque pas, comme le prescrit l'article (voir « Domaine
d'utilisation » et l'item 3.1 du plan de correction).

## Comparaison avec la figure 7

Alpha V3 utilise uniquement le modèle de Maxwell généralisé. Par défaut, la
précontrainte forme une base élastique conservée, puis les branches de Maxwell
décrivent les variations de contrainte produites après le début de
l'actionnement. Le mode élastique instantané n'est pas proposé.

(audit 2026-08, item 3.2) La **précontrainte viscoélastique est désormais
proposée en option** (`viscoelastic_history`, sélecteur « Précontrainte ») :
les branches de Maxwell sont actives dès l'élongation à 20 mm/min, comme dans
l'article. Arbitrage mesuré sur la figure 7 : ce mode reproduit le
« training » en ordre de grandeur (−5,0 % d'atténuation des pics sur
11 cycles, contre −10,3 % mesuré par l'article et −0,9 % en mode élastique),
rend la fig. 11 d'EXP simulable (fluage sous poids du bon signe), et amène le
couple à **−2/−1 %** des ancrages expérimentaux (avec la section
réactualisée) ; en revanche il abaisse les niveaux absolus de force
d'environ 22 % (−25/−30 % vs les ancrages), car le spectre de la Table A1
sur-relaxe (constat M6 de l'audit). Les niveaux publiés par l'article — pics
ET vallées — coïncident avec le mode **élastique conservé**, ce qui suggère
que la simulation des auteurs n'a pas non plus relaxé la prétension. Le
défaut reste donc `elastic_tk_reference` ; utiliser `viscoelastic_history`
pour étudier training, dead-band et relaxation de prétension, idéalement
avec un spectre identifié sur le matériau réel (item 3.3 du plan).

La validation distingue maintenant deux cibles : la courbe théorique noire de
l'article et les points expérimentaux colorés. Les erreurs RMSE sont calculées
séparément afin de ne pas confondre reproduction du modèle publié et accord
avec un spécimen expérimental.

(audit 2026-08) La validation quantitative vit dans `validation/`
(`python validation/validation_figure7.py`) avec une baseline de
non-régression exercée par la suite de tests. Écarts résiduels connus au gel
de la baseline (mode `fixed` par défaut), après correction de la calibration
de la numérisation et recalage de la pression sur les pics du texte : force
−4 à −7 % (ε = 0,8) et −11 à −13 % (ε = 1,0) par rapport à la courbe
théorique de l'article, couple −10 % et −25 %. **Depuis l'item 3.1 (correcteur
de point milieu), le mode `updated` referme l'essentiel de ces déficits** :
force −0,7/+0,7 % (ε = 0,8) et −8/−9 % (ε = 1,0), couple −6 à −12 % des
ancrages expérimentaux — mieux que la courbe théorique publiée elle-même,
vraisemblablement affectée du même biais d'intégration. Le résidu de force à
ε = 1,0, croissant avec la précontrainte, pointe vers le schéma de
précontrainte (item 3.2) — pas vers la discrétisation (< 1 %).

## Historique de pression mesuré

Dans `Pression et actionnement`, choisissez `Historique pression/temps mesuré
(CSV)`, puis chargez un fichier contenant une colonne de temps et une colonne
de pression. Les unités `MPa`, `bar`, `kPa` et `psi` sont prises en charge. Les
temps irréguliers sont conservés : aucune rampe de pression n'est reconstruite.
L'option de zéro initial soustrait uniquement l'offset du capteur au premier
échantillon.

## Condition du nylon

Alpha V2 utilise uniquement un nylon linéaire bilatéral lié aux extrémités. Il
participe entièrement à la précontrainte et à l'actionnement, avec des facteurs
de couplage fixés à `1`. Les modes traction seulement et glissement axial ont
été retirés.

## Intégration temporelle

L'intégration temporelle détermine comment la mémoire de contrainte des branches
de Maxwell est mise à jour entre `t` et `t + dt`. Pour chaque branche :

```text
d sigma_i / dt = E_i d epsilon / dt - sigma_i / tau_i
tau_i = eta_i / E_i
```

Le temps `tau_i` contrôle la rapidité de relaxation. Le pas `dt` contrôle la
fréquence à laquelle la pression, la déformation, la géométrie et l'équilibre
mécanique sont recalculés. Cette intégration est donc indispensable pour obtenir
la relaxation, l'hystérésis et la dépendance à la vitesse d'actionnement.

`Euler explicite` approxime la dérivée au début du pas :

```text
sigma_i(t + dt) = sigma_i(t)
                  + E_i Delta epsilon
                  - dt sigma_i(t) / tau_i
```

Cette méthode est d'ordre 1 et impose `dt < 2 tau_min`. Un pas proche de cette
limite peut déjà créer des oscillations numériques ou des contraintes de signe
incorrect. Un petit pas reste nécessaire pour obtenir une courbe précise.

`Intégration exponentielle stable` applique la décroissance exacte de la branche
pendant le pas :

```text
sigma_i(t + dt) = exp(-dt / tau_i) sigma_i(t)
                  + contribution de la déformation du pas
```

Elle évite l'instabilité propre à Euler pour la relaxation de Maxwell et constitue
le choix recommandé. Elle ne rend toutefois pas un grand `dt` précis : il faut
encore échantillonner suffisamment la rampe de pression, les changements de
géométrie, les pics et les boucles d'hystérésis.

## actionnement libre
Rh désigne le rayon de l’hélice, c’est-à-dire la distance entre l’axe central de l’actionneur et la ligne centrale du tube enroulé. Dans l’interface, il est affiché en mm.
beta_h désigne l’angle hélicoïdal instantané de l’actionneur. C’est l’angle formé par la direction locale du tube enroulé par rapport au plan transversal / à la géométrie de l’hélice. Dans l’interface, il est affiché en degrés.
