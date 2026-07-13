# Interface Cavatappi - livrable portable

Ce dossier contient tout le nécessaire pour lancer l'interface Streamlit du modèle Cavatappi

## Contenu du dossier

- `Base.py` : moteur scientifique du modèle.
- `livrable_base_calcul.py` : chargement  du moteur `Base.py`.
- `livrable_parametres.py` : paramètres, géométrie, pression et cache.
- `livrable_affichage.py` : visualiseur Cavatappi et graphes Matplotlib.
- `livrable_interface.py` : application Streamlit.
- `requirements.txt` : dépendances Python recommandées.
- `installer_dependances.bat` : installateur des dépendances sous Windows.
- `lancer_interface.bat` : lancement de l'interface sous Windows.
- `test_installation.py` : vérification rapide des imports et du modèle.

## Installation sur un PC Windows

1. Copier tout le dossier `livrable` sur le PC.
2. Installer Python si nécessaire.
3. Double-cliquer sur `installer_dependances.bat`.
4. Double-cliquer sur `lancer_interface.bat`.

L'interface s'ouvre ensuite dans le navigateur par défaut du PC


## Options d'hystérèse

L'onglet `Hystérèse` peut afficher plusieurs cycles sur le même graphe. Dans
la barre latérale, champ `Cycles à afficher`, entrer par exemple :

```text 1, 2```

Le même onglet peut aussi comparer :

- plusieurs précontraintes initiales ;
- plusieurs vitesses de pression injectée en `MPa/s`.

Pour les vitesses de pression injectée, l'historique de pression utilisé pour
la comparaison est triangulaire et linéaire afin que la pente corresponde à la
valeur demandée.


## Dependances Python

Les versions utilisees pour ce livrable sont listees dans `requirements.txt`.
Pour une installation manuelle :

```powershell
python -m pip install -r requirements.txt
```

## Portabilite

Le livrable ne dépend pas du répertoire courant. Les chemins sont résolus à
partir de l'emplacement de `livrable_interface.py`.

`livrable_base_calcul.py` cherche `Base.py` dans cet ordre :

1. dans le dossier `livrable` ;
2. dans le dossier parent, en secours.


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

Le paramètre Intégration temporelle correspond à la façon dont le modèle avance dans le temps, pas après pas, pour mettre à jour la partie viscoélastique du tube. À chaque pas de temps dt, le modèle connaît une nouvelle pression, calcule une nouvelle déformation, puis met à jour les contraintes dans les branches de Maxwell.

Incrémentale article : 
nouvelle contrainte = ancienne contrainte + contribution élastique due à la nouvelle déformation - contribution de relaxation pendant dt

Explicite article ( chaque branche est updaté séparement)
branche i à t + dt = branche i à t + effet élastique de la déformation - relaxation de cette branche pendant dt
nouvelle contrainte totale = ressort permanent E0 + branche Maxwell 1 + branche Maxwell 2 + branche Maxwell 3

Exponentielle stable (Cette option utilise une formule exponentielle pour représenter la relaxation :)
contrainte_i(t + dt) ≈ exp(-dt / tau) * contrainte_i(t) + nouvelle contribution élastique

Normalement
Incrémentale article  -> meilleur choix par défaut pour rester proche de l'article
Explicite article     -> utile pour comparer / vérifier la logique branche par branche
Exponentielle stable  -> utile si la simulation devient instable ou si dt est plus grand

## actionnement libre
Rh désigne le rayon de l’hélice, c’est-à-dire la distance entre l’axe central de l’actionneur et la ligne centrale du tube enroulé. Dans l’interface, il est affiché en mm.
beta_h désigne l’angle hélicoïdal instantané de l’actionneur. C’est l’angle formé par la direction locale du tube enroulé par rapport au plan transversal / à la géométrie de l’hélice. Dans l’interface, il est affiché en degrés.
