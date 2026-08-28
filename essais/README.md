# Campagne d'essais — notes d'analyse

Ce dossier contient les notes d'analyse et les scripts de dépouillement de la
campagne d'essais d'août 2026 (muscles D/G/J/L/Q/R, banc HX711 + capteur de
pression). **Les données brutes (CSV, Dataset.xlsx, essais StabCell) ne sont
pas versionnées** : les scripts sont fournis pour documenter la méthodologie
et se relancent en plaçant le dossier de données à côté.

| Fichier | Contenu |
|---|---|
| `ANALYSE_PRELIMINAIRE.md` | Dead-band mesuré, spectre de Maxwell identifié sur les « 10 N », déconvolution PVC seul, effet du référencement pression, stabilité de la cellule (StabCell) et revue de la chaîne d'acquisition |
| `MODELISATION_QUANTITATIVE.md` | Confrontation modèle/essai sans facteur d'échelle libre : niveaux absolus ±13 %, corrélations 0,94-0,98, incompatibilité d'amplitude, configuration recommandée par usage |
| `identifier_spectre.py` | Ajustement du spectre de Maxwell sur les relaxations « 10 N » |
| `recalcul_analyse.py` | Reconstruction de l'échelle de pression depuis les comptes ADC |
| `recalcul_correlations.py` | Corrélations modèle/mesure sous les deux échelles de pression |
| `spectre_pvc_seul.py` | Déconvolution tube/nylon du spectre composite |
| `spectre_identifie.json` | Spectre consolidé (muscle monté) — utilisable dans l'interface |
| `spectre_pvc_seul.json` | Spectre déconvolué (PVC seul) |

Les fiches de géométrie mesurée par muscle sont dans
`identification/fiches/` (utilisables directement par l'interface).
