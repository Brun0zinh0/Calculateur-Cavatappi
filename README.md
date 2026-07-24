# Interface Cavatappi Beta

Ce dossier contient tout le nécessaire pour lancer l'interface Streamlit du modèle Cavatappi

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
- `test_scientifique.py` : tests numériques des équilibres, historiques de pression et du mode suspendu.

## Installation sur un PC Windows

1. Copier tout le dossier `Beta` sur le PC.
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

Le dossier Beta ne dépend pas du répertoire courant. Les chemins sont résolus à
partir de l'emplacement de `interface.py`.

L'interface, les calculs parallèles et les tests importent directement le
`Base.py` placé dans le dossier `Beta`. Le lanceur se place automatiquement
dans ce dossier avant de démarrer l'application.

## Vérification de l'installation

Depuis le dossier `Beta`, exécuter :

```powershell
python test_installation.py
python test_scientifique.py
```

Le second test contrôle notamment la convergence de l'équilibre bloqué, la
décomposition force/couple, la durée exacte des profils de pression et
l'équation de déplacement du mode masse suspendue.

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

## Domaine d'utilisation

- La pression maximale proposée par l'interface est limitée à 1,5 MPa, domaine expérimental de l'article.
- L'intégration exponentielle est le choix recommandé et utilisé par défaut.
- Le profil radial de l'angle de biais est linéaire, conformément à la section 3.3 de l'article.
- La section de référence reste fixe pendant l'intégration et la géométrie hélicoïdale est mise à jour.
- Le module axial par défaut est la somme `E0 + E1 + E2 + E3`, comme indiqué dans l'annexe A.
- La précontrainte construit un état élastique de référence `t_k`. Les branches Maxwell décrivent ensuite les incréments d'actionnement autour de cet état.
- L'anisotropie identique de l'article reste le mode par défaut. Le mode `relaxation axiale identifiée` limite les branches Maxwell à la direction caractérisée par l'essai de traction.
- Le profil débit/volume linéaire est utilisé par défaut. Le profil non linéaire est phénoménologique.
- Les coefficients du nylon doivent rester égaux à 1 pour reproduire sa loi élastique linéaire complète.
- Les résultats servent à l'étude et à la comparaison. Un dimensionnement de sécurité nécessite une validation expérimentale du spécimen réel.


## Cache des réglages

L'application sauvegarde automatiquement les derniers réglages et les derniers
résultats dans le dossier temporaire Windows. Ce cache n'est pas nécessaire :
s'il n'existe pas sur un autre PC, l'application repart avec les valeurs par
défaut.

Dans l'interface, le bouton `Paramètres par défaut` réinitialise les paramètres
de simulation et efface les derniers résultats mis en cache.

## Transmission de la pression

La pression déforme le tube, réoriente ses directions matérielles anisotropes
et modifie la géométrie de l'hélice. Lorsque l'actionneur est bloqué, cette
tendance au mouvement produit la force calculée par le modèle.

## Comparaison avec la figure 7

La Beta utilise uniquement le modèle de Maxwell généralisé. La précontrainte
forme une base élastique conservée, puis les branches de Maxwell décrivent les
variations de contrainte produites après le début de l'actionnement. Le mode
élastique instantané et la précontrainte viscoélastique ne sont pas proposés.

La validation distingue maintenant deux cibles : la courbe théorique noire de
l'article et les points expérimentaux colorés. Les erreurs RMSE sont calculées
séparément afin de ne pas confondre reproduction du modèle publié et accord
avec un spécimen expérimental.

## Historique de pression mesuré

Dans `Pression et actionnement`, choisissez `Historique pression/temps mesuré
(CSV)`, puis chargez un fichier contenant une colonne de temps et une colonne
de pression. Les unités `MPa`, `bar`, `kPa` et `psi` sont prises en charge. Les
temps irréguliers sont conservés : aucune rampe de pression n'est reconstruite.
L'option de zéro initial soustrait uniquement l'offset du capteur au premier
échantillon.

## Condition du nylon

La Beta utilise uniquement un nylon linéaire bilatéral lié aux extrémités. Il
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
