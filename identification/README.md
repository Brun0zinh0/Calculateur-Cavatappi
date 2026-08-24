# Identification matériau — alpha V3

Outil de l'item 3.3 du plan de correction (audit 2026-08) : confronter le
spectre de Maxwell de l'article aux essais des muscles de l'équipe, et
identifier ce qui est identifiable avec les données disponibles.

## Usage

```powershell
python identifier_materiau.py --essai "..\..\expérimentale\Bruno\essai_force_pression_20260730_153353.csv"
python identifier_materiau.py --essai <essai> --fiche fiche_muscle_exemple.json
python identifier_materiau.py --essai <essai> --sans-simulation   # volet paliers seul
```

Rapport JSON et console ; sorties dans `out/`.

## Ce que l'outil identifie — et ce qu'il n'identifie pas

| Volet | Identifié | Condition |
|---|---|---|
| A. Paliers | Constantes de temps τ_k et fractions de relaxation (sans géométrie) | Paliers de pression **sans fuite** (essai E3) ; un palier où la force monte est signalé, pas forcé |
| B. Échelle | Un facteur d'échelle global de force + corrélation de forme | Fiche géométrique (E1) absente → les modules absolus ne sont **pas** identifiés |
| C. Courbure | Verdict « spectre vs structure » sur la courbure F(P) de la première montée | Rampe couvrant ≥ 0,1 MPa |

## Dépendances au plan d'essais (E1-E4)

- **E1 fiche géométrique** : remplir `fiche_muscle_exemple.json` pour chaque
  muscle (rayons du tube, hélice, angle de biais, longueurs, prédéformation
  réellement appliquée) — sans elle, l'échelle absorbe tout.
- **E2 relaxation uniaxiale / DMA sur le tube seul** : la vraie identification
  du spectre {E_k, η_k} se fait là, pas sur l'actionneur complet.
- **E3 maintien sans fuite** : l'essai « maintien 20 min » actuel mesure la
  fuite du circuit (corr(F,P) = 0,997) — inutilisable pour les τ_k.
- **E4 rampe quasi statique lente** : sépare la courbure F(P) de la
  viscoélasticité.

## Résultats du dégrossissage sur les données existantes (21/08/2026)

Exécuté sans fiche géométrique (E1 manquant), mode `elastic_tk_reference` :

| Essai | Paliers | Relaxation identifiée | Échelle | Corr. | Courbure F(P) (b/a, /MPa) |
|---|---|---|---|---|---|
| Bruno 153353 (0-0,40 MPa, 200 s) | 2 | 5,3 % (τ~313 s) sur un palier ; l'autre à force **croissante** | 6,2 | 0,47 | mesuré **+100** vs simulé −2,8 → **structurelle** |
| Bruno 155135 (0-0,58 MPa, 101 s) | 1 | aucune (force croissante) | 5,8 | 0,69 | mesuré **+118** vs simulé −1,3 → **structurelle** |
| Sacha muscle B (cadence 1 s supposée) | 2 | 4,0 % et 9,2 % (τ dispersés) | 5,6 | 0,87 | −0,8 vs −1,1 → compatible |

Conclusions :

- Les muscles de l'équipe relaxent de **4 à 9 %** aux paliers, contre **~83 %**
  à long terme pour le spectre Table A1 de l'article : le spectre publié est
  confirmé **non transférable** (constat M6 de l'audit, chiffré). Les données
  actuelles pointent vers un spectre quasi élastique (E0/ΣE ≈ 0,9 contre 0,17
  dans l'article) — à confirmer par les essais E2/E3.
- La courbure F(P) des montées de Bruno est **structurelle** (30 à 90 fois la
  courbure maximale du modèle, quel que soit le facteur d'échelle) : limite du
  modèle de l'article, documentée — à ne pas compenser par un fit de spectre.
- Les facteurs d'échelle ~5,6-6,2 absorbent la géométrie inconnue (essai E1
  indispensable pour des modules absolus).
