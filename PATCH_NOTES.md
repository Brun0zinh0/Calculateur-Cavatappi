# Notes de version

## v2.0.0-beta.1 - 16 juillet 2026

Cette version Beta transforme le livrable initial en une application autonome,
modulaire et portable pour l'étude des actionneurs Cavatappi.

### Interface

- Interface entièrement en français avec paramètres classés par catégories
  dépliantes et aides contextuelles.
- Visualiseur géométrique du Cavatappi avec dimensions, rotation à la souris
  et contrôle au clavier.
- Mémorisation des paramètres et des derniers résultats de simulation.
- Bouton de réinitialisation des paramètres par défaut.
- Barre de progression et estimation adaptative du temps de calcul.
- Sélection automatique d'un port local libre au lancement.

### Études disponibles

- Actionnement bloqué avec courbes temporelles de force, couple et pression.
- Étude de l'hystérèse avec sélection de plusieurs cycles, précontraintes et
  vitesses d'injection de pression.
- Étude de la relaxation après actionnement à pression constante.
- Force maximale en fonction de la précontrainte initiale.
- Actionnement libre avec masse suspendue, durée réglable, vitesse
  d'actionnement et maintien de pression.
- Lecture directe d'un historique pression/temps mesuré au format CSV.

### Modèle et calcul

- Séparation de la précontrainte élastique de référence et des incréments
  viscoélastiques d'actionnement.
- Choix de l'anisotropie de Maxwell et des conditions physiques du nylon.
- Choix de l'intégration temporelle et de la réponse élastique instantanée.
- Paramètres géométriques, matériaux et numériques accessibles dans
  l'interface.
- Calcul parallèle multi-cœur pour les études comparatives.
- Contrôles renforcés des entrées et des équilibres numériques.

### Organisation du projet

- Suppression du préfixe `livrable` dans les noms des modules.
- Import direct du moteur local `Base.py`.
- Suppression de la façade redondante `base_calcul.py`.
- Retrait de l'outil optionnel de validation Figure 7 du paquet Beta.
- Sept modules d'exécution et deux fichiers de test.

### Vérifications

- Compilation réussie des neuf fichiers Python.
- Test d'installation et mini-simulation validés.
- Huit tests scientifiques validés.
- Exécution parallèle avec deux processus validée sous Windows.
- Démarrage à froid de l'interface Streamlit validé sans erreur d'import ou
  d'exécution.

### Compatibilité

- Python 3.11 ou version plus récente recommandé.
- Installation automatisée des dépendances sous Windows.
- L'ancienne version reste disponible sur la branche `main` et sous le tag
  `v1.0.0-livrable`.
