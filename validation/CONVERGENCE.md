# Étude de convergence du moteur Base.py (item 4.2, plan de correction — audit 2026-08)

**Date :** 2026-08-24 · **Moteur :** `Base.py` version **2026.08.21-audit-phase4-14** · **Script :** `validation/etude_convergence.py` · **Résultats bruts :** `validation/out/convergence_results.json` et `validation/out/convergence_tables.md` (arbre d'exécution `C:\Users\b.pereiraazevedo\AppData\Local\Temp\tcpa_v3\alpha V3\validation\out\`)

> **Note de version.** Le plan de correction citait la version `2026.08.21-audit-phase32-13`. `Base.py` a été mis à jour vers `phase4-14` le 2026-08-24 à 08:44 (dans l'espace de travail et dans l'arbre court, synchronisés à l'identique) **avant** le lancement de la grille : les 30 runs de cette étude ont donc tous tourné sur la même version `phase4-14`, l'étude est cohérente.

## 1. Protocole

- **Bloqué (balayages A–E, G)** : pression cyclique **linéaire** (défaut projet depuis l'audit phase 2), Pmax = 1,4 MPa, pré-étirement eps = 0,8, 3 cycles ; débit 10 mL/min, volume 1,50 mL → période 18 s, durée 54 s. Historique de pression **généré explicitement** par `cyclic_pressure_history(dt=…)` et passé via `pressure_time`/`pressure_MPa`.
- **Suspendu (balayage F)** : `run_suspended_actuation`, charge 1 N, rampe simple 0 → 1,4 MPa en 9 s puis maintien 60 s, dt = 0,5 s.
- **Métriques** (dernier cycle) : pic de force (`force_total_mN`, max), vallée de force (min), pic de couple d'actionnement (max |`torque_act_microNm`|) ; en suspendu : contraction max et finale (%). Temps de calcul mur par run (machine locale, Python 3.14, NumPy 2.4.4).
- **Erreur relative** : par rapport au réglage **le plus raffiné de chaque balayage** (colonne « réf. »).

### Rappels critiques (acquis de l'audit + confirmés ici)

1. **Le `dt` de la configuration est IGNORÉ dès qu'une histoire de pression est fournie** : le pas de temps effectif du solveur est la grille de pression. Pour contrôler dt, générer l'historique avec le dt voulu (ce que fait ce script).
2. **`n_phi` est sans effet en actionnement bloqué (intégration axisymétrique)** — précision apportée par ce travail : l'effet est strictement nul **dès n_phi = 8** (écart ≤ 0,002 % vs n_phi = 32) ; **n_phi = 4 n'est PAS suffisant** (+1,94 % sur le couple).
3. **Garde d'Euler explicite** (`paper_explicit`) : `Base._validate_time_step` refuse dt ≥ 2·τ_min = 2·(154,57/20,67) ≈ **14,96 s**. Tous les dt testés (≤ 2 s) sont admissibles.
4. **Coût du mode `section_update_mode='updated'`** : l'audit annonçait ~2× ; **mesuré ici : 1,14–1,26×** seulement sur ce protocole (voir § 3.5 et anomalies § 5).

## 2. Tableaux de convergence

### 2.1 Balayage A — n_layers (bloqué, n_phi = 16, dt = 0,5 s, exponential, section fixe)

| n_layers | pic force (mN) | err. (%) | pic couple (µN·m) | err. (%) | vallée force (mN) | err. (%) | temps (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1  | 1550,36 | **−6,84** | 431,02 | −2,31 | 1371,62 | −6,36 | 2,7 |
| 2  | 1636,30 | −1,67 | 438,81 | −0,55 | 1442,06 | −1,55 | 5,2 |
| 3  | 1652,40 | −0,71 | 440,20 | −0,23 | 1455,21 | −0,66 | 7,6 |
| 4  | 1658,05 | −0,37 | 440,69 | −0,12 | 1459,83 | −0,34 | 10,0 |
| 6  | 1662,09 | −0,12 | 441,04 | −0,04 | 1463,13 | −0,12 | 15,1 |
| 8  | 1663,50 | −0,04 | 441,16 | −0,01 | 1464,28 | −0,04 | 19,9 |
| 10 | 1664,16 | réf. | 441,22 | réf. | 1464,82 | réf. | 24,6 |

**Genou : n_layers = 3–4.** L'erreur est à peu près divisée par 2 à chaque couche ajoutée. n_layers = 2 tient < 2 %, n_layers = 4 tient < 0,5 % sur les trois métriques, n_layers = 6 tient < 0,15 %. n_layers = 1 est inutilisable (−6,8 % de force).

### 2.2 Balayage B — n_phi (bloqué, n_layers = 4, dt = 0,5 s)

| n_phi | pic force (mN) | err. (%) | pic couple (µN·m) | err. (%) | vallée force (mN) | err. (%) | temps (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4  | 1661,04 | +0,18 | 449,25 | **+1,94** | 1464,90 | +0,35 | 3,7 |
| 8  | 1658,05 | 0,000 | 440,70 | +0,002 | 1459,84 | +0,001 | 5,9 |
| 16 | 1658,05 | 0,000 | 440,69 | 0,000 | 1459,83 | 0,000 | 10,4 |
| 32 | 1658,05 | réf. | 440,69 | réf. | 1459,83 | réf. | 18,7 |

**Genou : n_phi = 8** (convergence exacte, la quadrature angulaire résout exactement les harmoniques du problème axisymétrique dès 8 points). Le temps de calcul croît quasi linéairement avec n_phi : passer de 16 à 8 divise le coût par ~1,8 **sans aucune perte**. n_phi = 4 est à proscrire.

### 2.3 Balayage C — dt (bloqué, n_layers = 4, n_phi = 16, exponential)

| dt (s) | pic force (mN) | err. (%) | pic couple (µN·m) | err. (%) | vallée force (mN) | err. (%) | temps (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2    | 1666,48 | +0,66 | 470,74 | **+9,05** | 1455,79 | −0,35 | 3,6 |
| 1    | 1661,20 | +0,34 | 451,91 | **+4,69** | 1458,44 | −0,17 | 5,6 |
| 0,5  | 1658,05 | +0,15 | 440,69 | **+2,09** | 1459,83 | −0,08 | 10,3 |
| 0,25 | 1656,46 | +0,06 | 435,06 | +0,78 | 1460,53 | −0,03 | 19,1 |
| 0,1  | 1655,51 | réf. | 431,67 | réf. | 1460,95 | réf. | 47,2 |

**Le pic de couple est la métrique limitante : convergence du 1er ordre en dt** (l'erreur est divisée par 2 quand dt est divisé par 2). Extrapolation de Richardson (ajustement T(dt) = T₀ + c·dt sur dt = 0,5/0,25, vérifié à ±0,01 µN·m sur dt = 0,1) : **limite dt→0 ≈ 429,4 µN·m**, pente c ≈ 22,5 µN·m/s. Même dt = 0,1 s garde donc **+0,52 % de biais résiduel sur le couple** ; la force, elle, est convergée < 0,7 % dès dt = 2 s et < 0,2 % dès dt = 0,5 s.

### 2.4 Balayage D — intégrateurs (bloqué, n_layers = 4, n_phi = 16)

| dt (s) | intégrateur | pic force (mN) | pic couple (µN·m) | vallée force (mN) | temps (s) |
|---:|---|---:|---:|---:|---:|
| 0,5 | exponential    | 1658,05 | 440,69 | 1459,83 | 10,3 |
| 0,5 | paper_explicit | 1656,65 | 435,77 | 1459,46 | 9,4 |
| 2   | exponential    | 1666,48 | 470,74 | 1455,79 | 3,6 |
| 2   | paper_explicit | 1662,37 | 456,30 | 1453,79 | 3,3 |

Écart `paper_explicit` vs `exponential` à dt identique : **−0,08 % force / −1,11 % couple** à dt = 0,5 s ; −0,25 % / −3,07 % à dt = 2 s. Les deux intégrateurs convergent vers la même limite (429,4 µN·m sur le couple) ; à dt fini, `paper_explicit` est ici légèrement *moins* biaisé sur le couple, mais il est conditionnellement stable (garde dt < 14,96 s) et à peine plus rapide. **`exponential` reste le choix par défaut** (inconditionnellement stable, même coût) ; `paper_explicit` ne sert qu'à la fidélité stricte à l'article.

### 2.5 Balayage E — section réactualisée (`section_update_mode='updated'`, n_phi = 16, exponential)

| n_layers | dt (s) | pic force (mN) | err. (%)¹ | pic couple (µN·m) | err. (%)¹ | vallée force (mN) | err. (%)¹ | temps (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1   | 1706,55 | −1,20 | 577,43 | +2,78 | 1408,15 | −1,62 | 3,6 |
| 2 | 0,5 | 1700,56 | −1,54 | 556,43 | −0,96 | 1410,90 | −1,43 | 6,5 |
| 4 | 1   | 1729,20 | +0,12 | 582,21 | +3,63 | 1425,33 | −0,42 | 6,6 |
| 4 | 0,5 | 1723,01 | −0,24 | 560,96 | −0,15 | 1428,14 | −0,22 | 11,7 |
| 6 | 1   | 1733,42 | +0,36 | 583,11 | +3,79 | 1428,52 | −0,20 | 9,6 |
| 6 | 0,5 | 1727,19 | réf. | 561,81 | réf. | 1431,35 | réf. | 17,3 |

¹ vs le plus raffiné du balayage (L6, dt = 0,5).

Le comportement de convergence est **le même qu'en mode `fixed`** : mêmes genoux (n_layers = 4, couple dominé par dt). Écart physique fixed → updated à réglages identiques (dt = 0,5) : **+3,9 % sur le pic de force, +27 % sur le pic de couple** — c'est un effet de modèle, pas de discrétisation. Surcoût mesuré : **×1,14–1,26** seulement (voir anomalie § 5.2).

### 2.6 Balayage F — mode suspendu (charge 1 N, rampe 9 s → 1,4 MPa + maintien 60 s, dt = 0,5, n_phi = 16)

| n_layers | contraction max (%) | err. (%)² | contraction finale (%) | err. (%)² | longueur finale (mm) | temps (s) |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 5,739 | −0,26 | 2,579 | **+2,71** | 46,34 | 10,4 |
| 4 | 5,751 | −0,04 | 2,521 | +0,40 | 45,94 | 18,9 |
| 6 | 5,754 | réf. | 2,511 | réf. | 45,86 | 28,0 |

² vs n_layers = 6.

Sanity : la force axiale totale reste exactement égale à la charge (1000 mN, résidu ~5·10⁻¹⁵ N·mm) sur tous les runs. Le **pic** de contraction converge très vite (n_layers = 2 suffit à 0,3 %), mais la **contraction finale après relaxation** est plus sensible : n_layers = 2 dévie de +2,7 %, il faut **n_layers = 4** pour tenir < 0,5 %.

### 2.7 Balayage G — validation croisée des réglages recommandés (erreurs combinées)

Référence : n_layers = 10, n_phi = 16, dt = 0,1 s (le run le plus raffiné de l'étude, 118,5 s).

| Réglage | n_layers | n_phi | dt (s) | pic force err. (%) | pic couple err. (%) | vallée err. (%) | temps (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Interactif**  | 3 | 8 | 0,5 | −0,55 | +1,87 | −0,73 | **4,5** |
| **Production** | 6 | 8 | 0,1 | −0,12 | −0,04 | −0,12 | **39,8** |
| Référence      | 10 | 16 | 0,1 | réf. (1662 mN) | réf. (432,1 µN·m) | réf. (1466 mN) | 118,5 |

Les erreurs combinées mesurées confirment que les erreurs par paramètre s'additionnent approximativement et que les deux réglages tiennent leurs objectifs.

## 3. Réglages de production recommandés

| Usage | n_layers | n_phi | dt (s) | intégrateur | section | erreur max mesurée³ | temps (protocole 3 cycles / 54 s) |
|---|---:|---:|---:|---|---|---|---:|
| **Interface interactive** | **3** | **8** | **0,5** | exponential | fixed | **< 2 %** (−0,55 % force, +1,87 % couple, −0,73 % vallée) | **≈ 4,5 s** |
| **Production / publication** | **6** | **8** | **0,1** | exponential | fixed | **< 0,15 %** (−0,12 % force, −0,04 % couple) | **≈ 40 s** |

³ Erreur combinée mesurée (balayage G) vs la référence n_layers = 10 / n_phi = 16 / dt = 0,1.

**Notes d'application :**
- Ces dt s'entendent comme **pas de la grille de pression** (le dt de config est ignoré dès qu'un historique est fourni). Pour un CSV de pression expérimental, rééchantillonner à ≤ 0,5 s (interactif) / ≤ 0,1 s (publication) si le couple est exploité.
- Si **seule la force** compte, dt = 0,5 s suffit même en publication (erreur force < 0,2 %) : c'est le couple qui impose dt = 0,1 s.
- Toute valeur de couple publiée à dt = 0,1 s porte encore **+0,5 % de biais de discrétisation temporelle** (limite Richardson 429,4 µN·m sur ce protocole) ; pour un résultat de couple critique, faire deux dt (0,2 et 0,1) et extrapoler linéairement à dt → 0.
- **Mode `updated`** : mêmes réglages (genoux identiques), surcoût ~×1,2.
- **Mode suspendu** : n_layers = 4 minimum si l'état relaxé final est exploité (n_layers = 2 acceptable pour le seul pic de contraction en interactif).
- **`paper_explicit`** : uniquement pour comparaison à l'article, avec dt < 14,96 s (garde automatique du moteur) ; à dt ≤ 0,5 s l'écart avec `exponential` est ≤ 1,1 % (couple).
- Les défauts actuels de `default_simulation_config` (n_layers = 4, n_phi = 24, dt = 0,25) sont sûrs mais sous-optimaux : n_phi = 24 gaspille ~×3 de coût angulaire sans gain (8 suffit), et dt = 0,25 laisse +0,8 % sur le couple.

## 4. Analyse — genoux de convergence

- **n_layers** : convergence géométrique régulière (erreur ÷ ~2 par couche) ; genou à **3–4 couches**. La force et le couple convergent par le bas.
- **n_phi** : convergence **exacte dès 8** (quadrature angulaire du problème axisymétrique) ; aucun raffinement au-delà n'a d'effet, conformément à l'audit (qui l'avait mesuré sur 16 → 32) — l'étude étend le constat à 8.
- **dt** : la force converge vite (< 0,2 % dès dt = 0,5 s), le **pic de couple converge au 1er ordre en dt** et domine le budget d'erreur ; c'est le paramètre coûteux (temps ∝ 1/dt).
- **Intégrateur** : effet du second ordre devant celui de dt ; `exponential` par défaut.
- Le classement des sensibilités sur ce protocole : **dt (couple) > n_layers > intégrateur ≫ n_phi (nul dès 8)**.

## 5. Anomalies détectées

1. **Changement de version du moteur en cours de mission** : `Base.py` est passé de `2026.08.21-audit-phase32-13` (citée par le plan) à `2026.08.21-audit-phase4-14` le 2026-08-24 08:44, de façon synchronisée dans les deux arbres. Tous les runs de l'étude ont tourné sur `phase4-14` (version enregistrée dans le champ `meta` du JSON) ; l'étude est auto-cohérente mais **n'est pas opposable à `phase32-13`**. À re-vérifier si les recommandations doivent couvrir une autre version.
2. **Surcoût du mode `updated` bien plus faible qu'annoncé** : l'audit indiquait ~2×, la mesure donne **×1,14–1,26** sur ce protocole (3,6–17,3 s vs 2,7–15,1 s à réglages égaux). Soit le coût du recalcul de section est devenu marginal dans `phase4-14`, soit le facteur ~2× de l'audit provenait d'un autre protocole — la pénalité n'est en tout cas plus un argument contre le mode `updated`.
3. **Biais de couple non résorbé à dt = 0,1 s** (+0,52 % vs la limite dt → 0) : le pic de couple est une métrique intrinsèquement du 1er ordre en dt dans ce schéma ; documenté ci-dessus avec la procédure d'extrapolation.
4. **n_phi = 4 casse la neutralité angulaire** (+1,9 % sur le couple) : la règle « n_phi sans effet » ne vaut qu'à partir de 8.
5. Aucune anomalie de robustesse : 30/30 runs convergés, résidus ≤ 2,7·10⁻⁷ N·mm (bloqué) et ≤ 5,2·10⁻¹⁵ N·mm (suspendu), équilibre de la charge suspendue exact à la précision machine.

## 6. Reproduction

```bash
cd "alpha V3/validation"
python etude_convergence.py --list          # état de la grille (30 runs)
python etude_convergence.py --max-seconds 480   # exécute les runs manquants (relançable, flush incrémental)
python etude_convergence.py --report        # régénère out/convergence_tables.md
```

Les résultats bruts de la présente étude sont dans l'arbre d'exécution :
`C:\Users\b.pereiraazevedo\AppData\Local\Temp\tcpa_v3\alpha V3\validation\out\convergence_results.json` (et `convergence_tables.md`).
