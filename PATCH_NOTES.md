Depuis la V1, le modèle du Cavatappi a fait l’objet d’une révision importante. Les changements ne concernent pas uniquement l’interface : la représentation de la précontrainte, la loi viscoélastique, la contribution du nylon, la mécanique de pression et les méthodes numériques ont également été améliorées. L’objectif principal était d’obtenir un modèle plus stable, plus transparent et plus proche de la littérature, sans introduire de coefficient de calibration destiné uniquement à reproduire artificiellement les résultats espérés.

Dans la V1, le tube était déjà représenté par un modèle de Maxwell généralisé composé d’une branche élastique permanente et de trois branches viscoélastiques. Toutefois, la précontrainte initiale et les déformations produites par l’actionnement étaient traitées dans le même historique viscoélastique. Par conséquent, une partie de la contrainte créée lors de la mise en précontrainte pouvait ensuite se relaxer pendant la simulation. Cette formulation pouvait provoquer une diminution excessive de la force de rappel, notamment lors d’un maintien prolongé de la pression.

Dans la nouvelle version, la précontrainte peut être construite comme un état élastique de référence. La force de rappel présente avant l’actionnement est alors conservée comme base mécanique, tandis que les branches de Maxwell décrivent principalement les variations de contrainte produites après le début de l’actionnement. Cette séparation évite que toute la rigidité initiale du Cavatappi soit assimilée à une contrainte viscoélastique destinée à disparaître progressivement. Un second mode conserve néanmoins l’ancienne logique de précontrainte viscoélastique afin de permettre la comparaison entre les deux hypothèses.

L’intégration temporelle a également été revue. La V1 utilisait par défaut une formulation incrémentale explicite, sensible à la valeur du pas de temps et susceptible d’accumuler des erreurs au cours des simulations longues. La nouvelle version utilise par défaut une intégration exponentielle stable pour chaque branche de Maxwell. La décroissance visqueuse sur un pas de temps est ainsi calculée à partir du facteur exponentiel associé au temps de relaxation de la branche. Cette méthode reste stable lorsque le pas de temps devient relativement grand. L’intégration d’Euler explicite est toujours disponible à titre de comparaison, mais le programme vérifie désormais que le pas de temps choisi respecte ses conditions de stabilité.

La définition de la rigidité axiale a aussi été clarifiée. Dans la V1, le module axial était essentiellement donné par une valeur saisie manuellement de 31,24 MPa. La nouvelle version permet soit de conserver cette convention, soit d’utiliser la somme des modules du modèle de Maxwell :
Avec les paramètres actuels, cette somme vaut 37,76 MPa. Ce choix, utilisé par défaut, augmente physiquement la rigidité instantanée sans ajouter de facteur correctif arbitraire. Il contribue notamment à obtenir une force de rappel initiale plus cohérente avec la rigidité totale décrite par les branches de Maxwell.

La gestion de l’anisotropie viscoélastique est devenue configurable. Le premier mode applique une relaxation proportionnelle dans les différentes directions matérielles, conformément à l’hypothèse d’anisotropie identique. Le second limite la relaxation identifiée expérimentalement à la direction axiale, tandis que les autres directions conservent leur comportement principalement élastique. Cette option permet de ne pas déduire automatiquement une relaxation en cisaillement ou dans la direction radiale à partir d’un essai de traction axial. Elle ne remplace cependant pas une identification expérimentale indépendante des propriétés de cisaillement.
La contribution mécanique du nylon a été profondément clarifiée. Dans la V1, son couplage avec la déformation axiale produite par l’actionnement était nul par défaut. Le nylon contribuait donc à la force associée à la précontrainte, mais ne réagissait pas complètement aux changements de géométrie provoqués par la pression. Dans la version actuelle, ce couplage vaut physiquement un par défaut : le nylon se comporte comme un élément élastique dont la force dépend de sa déformation réelle. Trois conditions sont proposées :
nylon lié aux extrémités avec comportement linéaire bilatéral ;
nylon non collé travaillant uniquement en traction ;
nylon glissant axialement mais restant confiné dans le tube.
Ces conditions sont discrètes et physiquement interprétables. Aucun coefficient intermédiaire n’est utilisé pour ajuster artificiellement la contribution du nylon.

La mécanique du tube et de la pression a aussi été rendue plus explicite. Le modèle permet de choisir si la section radiale reste celle de référence ou si elle évolue avec la déformation. Le profil radial de l’angle de biais du tube peut suivre une variation linéaire avec le rayon ou une loi de torsion uniforme. Après réévaluation de ses conditions physiques, l’ancienne option ajoutant une poussée directe sur un fond pressurisé a été entièrement supprimée. La force totale calculée est désormais strictement la somme des contributions mécaniques du tube et du nylon, aussi bien en actionnement bloqué que dans le mode avec masse suspendue.

La génération de la pression a été améliorée. Le modèle peut construire les cycles à partir du débit et du volume injecté, produire une rampe suivie d’un maintien à pression constante ou utiliser directement un historique expérimental pression-temps provenant d’un fichier CSV. Dans ce dernier cas, les instants de mesure irréguliers sont conservés et aucune rampe artificielle n’est reconstruite. Cela permet de comparer la simulation à un essai réel en lui appliquant exactement la pression mesurée.

Le mode d’actionnement libre avec masse suspendue a également été complété. Il repose sur l’équilibre entre les efforts du tube, ceux du nylon et la charge suspendue. La durée, la vitesse d’augmentation de la pression et le maintien éventuel de la pression sont maintenant réglables. Il devient ainsi possible d’étudier aussi bien le déplacement pendant le gonflage que son évolution viscoélastique pendant une phase de maintien.
Enfin, la fiabilité numérique a été renforcée par des contrôles sur la géométrie, les propriétés matérielles, le pas de temps et les résidus d’équilibre. Les études comparatives peuvent utiliser plusieurs cœurs du processeur. La version actuelle est accompagnée d’un test d’installation et de neuf tests scientifiques portant notamment sur les profils de pression, l’équilibre bloqué, le mode masse suspendue, les conditions du nylon, l’état de précontrainte, la lecture des mesures CSV et l’export des résultats.

Pour finir la possibilité d'exporter les résultats en fichier CSV pour les traiter en dehors de l'environement du logiciel à était ajouté.

## Mise Ã  jour Beta - interface fiabilisÃ©e (24 juillet 2026)

- Persistance durable des paramÃ¨tres et des rÃ©sultats, avec Ã©critures atomiques et migration des anciens caches.
- Import et export des paramÃ¨tres en JSON ou CSV depuis la barre latÃ©rale.
- Les CSV de rÃ©sultats incluent dÃ©sormais les paramÃ¨tres d'entrÃ©e et la version du modÃ¨le.
- Interface responsive avec action principale accessible immÃ©diatement, mÃ©triques compactes et onglets adaptÃ©s au mobile.
- Tous les rÃ©glages avancÃ©s sont directement disponibles dans la barre latÃ©rale.
- Visualiseur 3D et graphiques harmonisÃ©s avec le thÃ¨me sombre.
- LisibilitÃ© corrigÃ©e dans la coupe du tube, notamment l'annotation du diamÃ¨tre du nylon.
- Progression, estimation de durÃ©e, rÃ©initialisation et messages d'Ã©tat amÃ©liorÃ©s.
- Conservation des modes constitutifs fixes de la Beta et renforcement des validations d'entrÃ©e.

## Mise Ã  jour visuelle du visualiseur 3D (28 juillet 2026)

- Le tube PVC utilise dÃ©sormais un rendu blanc cassÃ© semi-transparent.
- Les ombres, reflets et contours ont Ã©tÃ© retravaillÃ©s pour amÃ©liorer la perception de la gÃ©omÃ©trie.
- Un liserÃ© discret distingue les Ã©tats fabriquÃ© et prÃ©contraint.
- Le rendu a Ã©tÃ© vÃ©rifiÃ© sur ordinateur et mobile, ainsi que pendant la rotation interactive.


## Mise a jour Audit 2026-08 - alpha V3 (24 aout 2026)

Cette version est issue d'un audit scientifique complet du moteur contre les
trois articles de reference (equation par equation) et les donnees
experimentales de l'equipe, suivi de l'execution integrale du plan de
correction (PLAN_DE_CORRECTION.md, rapport complet dans audit/rapport_audit.html).
Moteur 2026.08.21-audit-phase4-14, suite de tests passee de 14 a 26 tests.

- Chaine de validation quantitative : validation/validation_figure7.py
  (calibration de la numerisation corrigee, baseline de non-regression +/-2 %
  exercee par la suite de tests), validation/validation_fig11_exp.py
  (mode suspendu valide contre la fig. 11 d'EXP : fluage a ~10 %),
  validation/etude_convergence.py + CONVERGENCE.md (reglages recommandes :
  interactif 3 couches/phi 8/dt 0,5 ; production 6 couches/phi 8/dt 0,1).
- Mode section reactualisee (updated) corrige : le ratchet de cyclage etait
  un biais O(h) de l'integrateur geometrique explicite, elimine par un
  correcteur de point milieu du BVP radial (fermeture de cycle x28). Sur la
  figure 7, le couple passe de -33 % a -6/-12 % des ancrages experimentaux.
- Precontrainte viscoelastique en option (selecteur « Precontrainte ») :
  branches de Maxwell actives des l'elongation, training et fig. 11
  simulables ; le defaut reste la reference elastique (fidele aux niveaux
  publies).
- Outil d'identification materiau (identification/identifier_materiau.py) :
  paliers -> constantes de temps et fractions de relaxation, facteur
  d'echelle en forme fermee, verdict « spectre vs structure » sur la
  courbure F(P). Constat chiffre : le spectre de l'article n'est pas
  transferable aux muscles de l'equipe (relaxation 4-9 % mesuree vs ~83 %).
- Interface : superposition d'un essai experimental sur UN SEUL graphe quand
  la pression injectee est l'historique mesure (option activable/desactivable
  dans la barre laterale), option d'appariement des coefficients de Poisson,
  actionnement en % normalise par la longueur non chargee L_T0 (ancienne
  normalisation exportee en _loaded_ref).
- Corrections diverses de l'audit : profils de pression lineaires par defaut
  sur tous les points d'entree, verrou de reference serie a P=0, garde
  theta=0, inference d'unites par mots entiers, rampe de mise en charge de la
  masse, garde de residu d'equilibre, compliance des extremites evaluee a
  l'angle post-precontrainte, documentation honnete des ecarts a l'article
  (profil d'angle de biais, section figee) directement dans le README.
