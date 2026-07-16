# Interface Cavatappi Beta

Ce dossier contient tout le nécessaire pour lancer l'interface Streamlit du modèle Cavatappi

## Contenu du dossier

- `Base.py` : moteur scientifique du modèle.
- `parametres.py` : paramètres, géométrie, pression et cache.
- `pression.py` : lecture des historiques CSV mesurés et conversion des unités.
- `affichage.py` : visualiseur Cavatappi et graphes Matplotlib.
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

## option complexe

Dans le modèle principal, la pression ne donne pas directement une force axiale comme un vérin. Elle déforme le tube : le tube gonfle, ses directions matérielles anisotropes se réorientent, l’hélice cherche à changer de géométrie, et comme l’actionneur est bloqué, cette tendance au mouvement devient une force mesurée.

Donc le chemin principal est :
pression -> déformation du tube -> changement de géométrie hélicoïdale empêché -> force bloquée

La force de fond pression, elle, représente un autre chemin possible :
pression -> poussée directe sur une surface fermée -> force axiale

C’est le même principe qu’un vérin pneumatique : si vous avez une pression interne P qui pousse sur une surface A, vous obtenez une force directe :
F = P*A

Dans l'interface, cette contribution utilise uniquement un facteur physique de
zéro (désactivée) ou un (activée). Aucun multiplicateur ajustable n'est utilisé.

## Comparaison avec la figure 7

Le mode Maxwell généralisé reste le modèle temporel principal. Une réponse
`Elastique instantanee` est aussi disponible comme référence sans relaxation :
elle ne constitue pas une calibration et ne doit pas servir à représenter
l'hystérèse. Le mode par défaut conserve plutôt un état élastique `t_k`
séparé puis applique Maxwell aux incréments d'actionnement. Cette formulation
rapproche simultanément les niveaux de force et le couple de la théorie publiée
sans changer les modules du matériau.

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

Trois conditions discrètes sont disponibles, sans coefficient intermédiaire :

- `Linéaire bilatéral` reproduit l'équation (4) de l'article et autorise traction et compression ;
- `Traction seulement` annule la force lorsque le filament devient mou ;
- `Glissant axialement` annule sa force axiale tout en conservant son confinement géométrique en flexion et torsion.

Le paramètre Intégration temporelle correspond à la façon dont le modèle avance dans le temps, pas après pas, pour mettre à jour la partie viscoélastique du tube. À chaque pas de temps dt, le modèle connaît une nouvelle pression, calcule une nouvelle déformation, puis met à jour les contraintes dans les branches de Maxwell.

Euler explicite de l'article met à jour séparément chaque branche :

branche i à t + dt = branche i à t + effet élastique de la déformation - relaxation pendant dt

Cette formule exige un pas inférieur à deux fois le plus petit temps de relaxation. Un pas supérieur au temps de relaxation peut déjà produire des oscillations numériques.

L'intégration exponentielle cohérente utilise la décroissance exacte de chaque branche sur le pas :

contrainte_i(t + dt) = exp(-dt / tau_i) contrainte_i(t) + contribution de la déformation du pas

Le même opérateur tangent exponentiel est utilisé pour la condition de pression radiale et pour la mise à jour des contraintes. C'est le choix par défaut recommandé.

## actionnement libre
Rh désigne le rayon de l’hélice, c’est-à-dire la distance entre l’axe central de l’actionneur et la ligne centrale du tube enroulé. Dans l’interface, il est affiché en mm.
beta_h désigne l’angle hélicoïdal instantané de l’actionneur. C’est l’angle formé par la direction locale du tube enroulé par rapport au plan transversal / à la géométrie de l’hélice. Dans l’interface, il est affiché en degrés.
