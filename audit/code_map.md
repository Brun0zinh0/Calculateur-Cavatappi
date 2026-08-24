# Cartographie du moteur scientifique `Base.py` (alpha V2) — rapport d'audit

Version du modèle : `MODEL_VERSION = "2026.07.30-series-compliant-ends-10"` (Base.py l.28). Fichiers analysés : `alpha V2/Base.py` (2071 lignes), `alpha V2/parametres.py` (891 l.), `alpha V2/pression.py` (167 l.), `alpha V2/README.md`, et `alpha V2/interface.py` (survol ciblé). Chemin racine : `C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/alpha V2/`.

---

## 1. Structures de données et paramètres par défaut

Le moteur duplique les structures : `Base.py` fournit des **fabriques de `SimpleNamespace`** (autonomie du module), `parametres.py` fournit des **dataclasses équivalentes** utilisées par l'interface. Les deux jeux de défauts ne coïncident pas partout (voir §5).

### 1.1 `default_maxwell_tensile_params` (Base.py l.36-47) / `MaxwellTensileParams` (parametres.py l.200-220)
| Champ | Défaut | Unité (implicite) |
|---|---|---|
| `E0` | 6.36 | MPa (ressort permanent) |
| `E1, E2, E3` | 20.67, 5.98, 4.75 | MPa |
| `eta1, eta2, eta3` | 154.57, 977.79, 11044.83 | MPa·s |

Les unités ne sont **pas déclarées** dans Base.py ; elles le sont via les suffixes des clés de `DEFAULT_SETTINGS` (`maxwell_E0_mpa`, `maxwell_eta1_mpa_s`, parametres.py l.179-185). Temps de relaxation τᵢ = ηᵢ/Eᵢ ≈ 7.5 s, 163 s, 2325 s. Accès génériques : `_maxwell_E/_eta/_E_total/_rates/_branch_count` (Base.py l.119-142), compatibles namespace (attributs E1…) et dataclass (propriétés `E`, `eta`, `E_total`).

### 1.2 `default_material_params` (Base.py l.50-66) / `MaterialParams` (parametres.py l.223-236)
| Champ | Défaut Base.py | Défaut parametres.py | Unité |
|---|---|---|---|
| `E_axial` | **37.76** | **31.24** | MPa |
| `E_radius` | 8.82 | 8.82 | MPa |
| `G12` | 7.24 | 7.24 | MPa |
| `nu12`, `nu23` | 0.205, 0.422 | idem | — |
| `E_nylon`, `G_nylon` | 3690, 790 | idem | MPa |
| `maxwell_anisotropy_mode` | `"paper_equal"` | idem | (`"paper_equal"` ou `"axial_test_only"`) |
| `nylon_condition_mode` | `"bonded_linear"` (seul accepté, l.333) | idem | — |
| `nylon_axial_prestrain/actuation_coupling` | 1.0 (verrouillés à 1.0, l.341-346) | idem | — |

37.76 = E0+E1+E2+E3 (mode `maxwell_sum`) ; 31.24 = valeur « paper_table ». `build_config` (parametres.py l.789-794) écrase `E_axial` par `maxwell.E_total` si `axial_modulus_mode == "maxwell_sum"` (défaut).

### 1.3 `default_geometry_params` (Base.py l.69-84) / `GeometryParams` (parametres.py l.239-251)
| Champ | Défaut | Unité |
|---|---|---|
| `Rout`, `Rin` | 1.0, 0.4 | mm |
| `r_nylon` | 0.77/2 | mm |
| `rho0` (rayon d'hélice) | 2.16 | mm |
| `alpha0_deg` (angle d'hélice initial) | 10.53 | **degrés** |
| `theta_f_deg` (angle de biais en surface) | 37.91 | **degrés** |
| `initial_length` | 32.45 | mm (longueur axiale active) |
| `uncoiled_length` | 0.0 | mm (extrémités non enroulées, somme des deux bouts) |
| `uncoiled_compliance_mode` | `"tangent_beam"` | (`"axial_rod"` possible) |
| `bias_angle_profile` | `"paper_linear"` | (`"uniform_twist"` possible) |
| `section_update_mode` | `"fixed"` | (`"updated"` possible) |

### 1.4 `default_discretization` (Base.py l.87-95)
`n_layers=18`, `n_phi=72`, `pre_steps=120`, `dw_bracket=(-0.08, 0.08)` (sans unité, borne de la déformation incrémentale dw). **Jamais utilisés par les runners** : `run_blocked_actuation` (l.1671-1676) et `run_suspended_actuation` (l.1755-1760) reconstruisent la discrétisation depuis la config avec `dw_bracket=(-0.05, 0.05)`.

### 1.5 `default_simulation_config` (Base.py l.98-116) / `SimulationParams` (parametres.py l.254-273)
Divergences notables : Base.py → `n_cycles=3`, `Pmax=1.3` MPa, `dt=0.25` s, `n_layers=4`, `nonlinear_pressure=False` ; parametres.py → `n_cycles=11`, `Pmax=1.5` (`P_MAX_MPA`), `dt=0.5`, `n_layers=8`, `nonlinear_pressure=True`. `flow_rate_mL_min=10.0` (mL/min), `volume_mL=1.5` (mL) identiques. `DEFAULT_SETTINGS` de l'interface (parametres.py l.122-197) : `dt=0.5`, `n_layers=4`, `n_phi=16`, `nonlinear_pressure=False`.

### 1.6 États et résultats
- `HelixState` (Base.py l.145-151) : `rho` [mm], `alpha` [rad], `h` [mm/rad — hauteur par radian, pas = 2πh], `pressure` [MPa], `time` [s].
- `StepResult` (l.154-174) : `dw` [-], `dv` [rad/mm, incrément de torsion linéique], `dkappa` [1/mm], `Ft` [N], `Tt` [N·mm], `residual` [N·mm], `Ftube/Fnylon` [N, projetés sur l'axe], `Mtube/Mnylon/Ttube/Tnylon` [N·mm], `axial_stretch` [-], `Rin/Rout` [mm].
- Champs de contrainte : tableaux `(n_layers, n_phi, 6)` en composantes de Voigt **[s, φ, r, rφ?, sr?, sφ]** — seules les positions 0,1,2,5 sont peuplées (l.807) : `sigma_reference` (base élastique de précontrainte), `sigma0` (branche E0), `sigma_i` (branches de Maxwell, dim. supplémentaire n_maxwell), `sigma_total` (l.386-390), en MPa.

### 1.7 `DEFAULT_SETTINGS` et validation (parametres.py)
Clés suffixées par l'unité (`rout_mm`, `p_max_mpa`, `maxwell_eta1_mpa_s`, `suspended_mass_g`, `flow_rate_mL_min`, `alpha0_deg`…). Bornes dures dans `normalize_settings` (l.362-410), p.ex. `p_max_mpa ∈ [0, 1.5]`. Verrouillages Alpha V2 (l.354-359) : mode constitutif, référence de précontrainte, nylon, couplages = 1. Validations physiques : `geometry_error` l.525, `material_error` l.549 (reconstruit la matrice C et vérifie sa positivité — duplication du code de `ti_stiffness_from_paper`), `numerical_error` l.621 (dont stabilité d'Euler explicite l.651-675). `derived_geometry` (l.685-761) : grandeurs dérivées en mm/mm²/mm⁴, volume interne en mL via `× 1.0e-3` (l.729, mm³→mL).

---

## 2. Flux du solveur `TCPAMaxwellBlockedModel` (Base.py l.306-1511)

### 2.1 Construction (`__init__`, l.309-405)
- Conversion degrés→radians : `alpha0 = deg2rad(alpha0_deg)`, `theta_f = deg2rad(theta_f_deg)` (l.358-359).
- Géométrie : `h0 = rho0·tan(alpha0)` [mm/rad], `turns = initial_length/(2π h0)` (l.360-361).
- Maillage : couches radiales `R_edges = linspace(Rin, Rout, n_layers+1)` (l.365-367), secteurs angulaires `phi` (l.368-369).
- Raideur locale `C_local_total` (l.371) ; angles de biais par couche `_bias_angles` (l.378, cf. l.530-534 : linéaire `r/Rout·θf` ou `arctan(r/Rout·tanθf)`).
- Décomposition Maxwell + rotation par couche : `_rebuild_section_properties` (l.384, corps l.505-528).
- Compliance série des extrémités : `_uncoiled_series_compliance` (l.400, corps l.468-482) [mm/N].
- Validation : `_validate_inputs` (l.407-455), inclut la positivité spectrale de C (l.453-455).

### 2.2 Déroulé d'un pas bloqué — `step(pressure_new, dt, h_target)` (l.1135-1176)
1. **Validation du pas** : `_validate_time_step` (l.536-549) — pour Euler explicite, exige `dt < 2·τ_min`.
2. **Incrément de pression** : `dP = pressure_new − helix.pressure` (l.1141), en MPa. La pression n'entre que comme **condition limite radiale incrémentale** sur la paroi interne (`b_vec[0] = −dP − k_sig`, l.686) ; il n'y a **aucun terme de poussée de fond** (P·π·Rin²) dans l'équilibre axial.
3. **Recherche de l'équilibre** : `dw = _find_dw(dP, dt, h_target)` (l.1142 → corps l.973-1006). Balayage de 65 valeurs de dw dans `dw_bracket`, `brentq` sur chaque changement de signe du résidu de moment (l.988-990, xtol=1e-10), choix de la racine de plus petit |dw| (l.992) ; repli `minimize_scalar` borné avec rejet si résidu > 1e-6 N·mm ou solution en butée (l.994-1006).
4. **État d'essai** : `_trial_state(dw, dP, dt, h_target)` (l.781-874), qui enchaîne :
   - **Cinématique hélicoïdale** (l.782 → `_kinematic_increments` l.569-574 + `_new_geometry_from_dw_and_h` l.560-567) : la fibre s'allonge de (1+dw) à h imposé → `rho_new = √(l² − h²)`, `alpha_new = atan2(h, rho_new)` ; incréments de torsion `dv = sin2α/(2ρ)|new − old` et de courbure `dkappa = cos²α/ρ|new − old`.
   - **Tangente algorithmique viscoélastique** : `_algorithmic_data(dt)` (l.783 → l.616-638). Euler explicite : C_alg = C_total, historique = −dt·(Eᵢ/ηᵢ)·σᵢ ; exponentielle : facteurs `(1−e^(−dt/τ))/(dt/τ)` sur les Cᵢ et historique `(e^(−dt/τ)−1)·σᵢ` (l.627-637). Pendant la précontrainte (`_building_reference_state`), tangente purement élastique (l.618-619).
   - **Problème radial de Lekhnitskii** : `_solve_radial_constants` (l.784 → l.661-721). Système 2n×2n pour les constantes (C1, C2) du déplacement radial `u = C1·R^μ + C2·R^(−μ) + A·dv·R² + B·dw·R` (`_u_base` l.586-592, `_layer_AB_mu` l.576-584 avec `μ = √(C̄22/C̄33)`). CL : σ_r incrémental = −dP en R_in (l.679-687), continuité de u (l.689-700) et de σ_r (l.702-711) aux interfaces, σ_r = 0 en R_out (l.713-720). L'historique visqueux moyenné en φ (`CSmean`, l.673) entre dans la CL radiale.
   - **Mise à jour de section** : `_trial_section_geometry` (l.785 → l.738-779). Mode `fixed` : géométrie de section inchangée (l.745-751) ; mode `updated` : rayons + angles de biais recalculés (l.752-779), garde-fous l.761, 777.
   - **Champ de déformation et contraintes** (l.796-818) : pour chaque (couche j, angle φ) — terme de courbure `curv = (dκ·R·cosφ + u·K_old·cosφ)/(1 + K_old·R·cosφ)` avec `K_old = cos²α/ρ` (l.794) ; corrections par coefficients de Poisson effectifs `vbar` (issus de `effective_poissons` l.264-269) : `ε_r = du − ν̄12·curv`, `ε_φ = u/R − ν̄13·curv`, `ε_s = dw + curv`, `γ_sφ = dv·R/denom − ν̄14·curv` (l.803-806). Puis mise à jour des contraintes : base de référence pendant la précontrainte (l.808-811) ou `σ0 += C0·dε` + branches de Maxwell `_update_maxwell_branches` (l.640-659 : schéma explicite l.644-646 ou exponentiel exact l.647-654).
   - **Résultantes de section** (l.820-829) : `Ftube = ∮∫σ_s R dR dφ` [N], `Mtube = ∮∫σ_s·(R cosφ) dA` [N·mm], `Ttube = ∮∫σ_sφ·R dA` [N·mm]. NB : contraintes évaluées aux anciens rayons `self.R_centers` (l.796) mais intégrées avec les poids des rayons **nouveaux** (l.823) — sans effet en mode `fixed`.
   - **Nylon lié** (l.831-835) : `Fny += π·E_nylon·r² · dw` [N], `Mny += (π/4)·E_nylon·r⁴ · dκ`, `Tny += (π/2)·G_nylon·r⁴ · dv` (EA, EI, GJ d'une tige circulaire).
   - **Équilibre bloqué** (l.837-848) : `Ft = (Ftube+Fny)/sinα` ; `Tt = (M + Ft·ρ·sinα)/cosα` ; **résidu de fermeture** = `Tt·sinα + Ft·ρ·cosα − (Ttube+Tny)` [N·mm] — c'est ce résidu qu'annule `_find_dw`.
5. **Commit** : `_commit_trial` (l.1120-1133) copie les états de contrainte, forces nylon, géométrie de section, puis `_rebuild_section_properties`. Mise à jour de `helix` (l.1146-1152) et enregistrement du `StepResult` (l.1154-1176).

### 2.3 Variantes de pas
- **`step_blocked_series`** (l.1178-1215) : blocage avec extrémités compliantes. `_find_blocked_series_trial` (l.1008-1118) résout par `least_squares` 2×2 en (dw, h_target) les résidus normalisés [moment, compatibilité de longueur] avec `_series_length_residual` (l.496-503 : Δlongueur active + compliance·ΔF = 0). Redémarrages alternatifs l.1073-1102, tolérances finales 1e-5 N·mm / 1e-6 mm (l.1111).
- **`step_suspended`** (l.1316-1359) : actionnement libre sous charge. `_find_suspended_state` (l.1233-1314) : `least_squares` 3×3 en (dw, ρ, α), résidus `_suspended_residuals` (l.1217-1231) = les trois équilibres de l'article (F = F_load·sinβ, M = −F_load·Rh·sinβ, T = F_load·Rh·cosβ), via `_trial_state_from_geometry` (l.876-971, duplication quasi complète de `_trial_state` avec ρ, α imposés). Bornes : dw∈[−0.25, 0.25], ρ∈[1.001·Rout, max(3ρ0, 2Rout)], α∈[0.5°, 85°] (l.1236-1237). Après convergence, `h = ρ·tanα` (l.1328).
- **`prestretch_to(eps_tk, strain_rate_mm_min=20.0)`** (l.1361-1379) : rampe de h de h0 à (1+eps)·h0 en `pre_steps` pas, durée `60·eps·L0/vitesse` [s], à pression nulle, avec `_building_reference_state=True` → construit la base élastique `sigma_reference` (précontrainte « élastique conservée »).

### 2.4 Boucles temporelles
`run_pressure_history` (l.1381-1392, choisit step ou step_blocked_series selon `series_compliance_enabled`) ; `run_pressure_history_suspended` (l.1394-1402). Les deux exigent un temps strictement croissant et intègrent de proche en proche (dt variable par échantillon).

---

## 3. Localisation des composants demandés

| Composant | Localisation | Notes |
|---|---|---|
| Raideur transversalement isotrope | `ti_stiffness_from_paper` Base.py **l.182-213** | Résout un système 4×4 pour (C11, C12, C22, C23) à partir de (E_axial, E_radius, ν12, ν23), C44 = (C22−C23)/2, C55=C66=G12. Ordre Voigt « engineering shear ». Dupliqué dans `material_error` (parametres.py l.594-613). |
| Rotation de Voigt corrigée | `rotate_stiffness_bias` Base.py **l.239-261**, avec `voigt_to_tensor` l.219-228 / `tensor_to_voigt` l.231-236, `VOIGT_PAIRS` l.216 | Rotation tensorielle complète (einsum 4 indices) autour de l'axe r ; le remplissage symétrique sans facteur 2 est cohérent avec la convention engineering (commentaire l.220 « pas de facteur shear supplementaire » = correction V6). Symétrisation finale l.261. Contrôle d'invariance isotrope : `stiffness_sanity_report` l.1994-2010. |
| Poissons effectifs | `effective_poissons` l.264-269 | Inversion du sous-bloc [0,1,2,5]. |
| Intégrateurs temporels Maxwell | `_algorithmic_data` **l.616-638** (tangente + historique) et `_update_maxwell_branches` **l.640-659** (mise à jour des branches) ; noms normalisés par `_normalize_integration_name` l.272-282 ; garde de stabilité `_validate_time_step` l.536-549 | Deux schémas : `paper_explicit` (Euler avant) et `exponential` (exact par branche). |
| Solveur d'équilibre bloqué | `_find_dw` **l.973-1006** : balayage 65 points + **`brentq`** (l.990) + repli `minimize_scalar` (l.998) | Racine de plus petit |dw| retenue si plusieurs (l.992). |
| Solveur série (extrémités) | `_find_blocked_series_trial` **l.1008-1118** : **`least_squares`** 2D | Verrou du référentiel : `lock_blocked_series_reference` l.488-494. |
| Solveur suspendu | `_find_suspended_state` **l.1233-1314** : **`least_squares`** 3D multi-départs | |
| Extrémités compliantes | rigidités `_uncoiled_section_rigidities` l.457-466 ; compliance `_uncoiled_series_compliance` **l.468-482** (mode `tangent_beam` : traction projetée sin²α + flexion de poutre cos²α, deux demi-longueurs ; mode `axial_rod` : L/EA) ; résidu de compatibilité `_series_length_residual` l.496-503 ; sorties l.1432-1471 | Équivalent interface : `derived_geometry` parametres.py l.713-728. |
| Conventions de sortie | `history_arrays` **l.1404-1511** (conversions N→mN, N·mm→µN·m, longueurs) puis `add_corrected_output_conventions` **l.1593-1614** (grandeurs « actuation » = valeur − valeur au 1er échantillon post-précontrainte) | |
| `run_blocked_actuation` | **l.1658-1704** | prestretch → pas initial dt=0 à P(0) (l.1696) → verrou série (l.1697) → historique → découpe à partir de `i_act0` et re-zéro du temps (l.1701-1702). |
| `run_suspended_actuation` | **l.1737-1810** | Équilibrage préalable de la masse `_equilibrate_suspended_load` l.1707-1734 (pas de relaxation τ₁…τ₃, 5τmax, 20τmax en intégration exponentielle forcée, puis horloge/historique remis à zéro l.1732-1733) ; longueur de référence l.1785-1791 ; contraction signée l.1799-1806. |
| `run_hold_relaxation` | **l.1813-1843** | Rampe+maintien via `ramp_hold_pressure_history` l.1555-1590 ; indices de maintien et courbes de relaxation l.1835-1842. |
| Historiques de pression | `cyclic_pressure_history` l.1519-1552, `ramp_hold_pressure_history` l.1555-1590, `_time_grid_with_events` l.285-298, `_prepare_actuation_history` l.1617-1639 | Équivalent interface : `make_pressure_history` parametres.py l.858-891. |
| Graphes | `plot_response` l.1861, `plot_decomposition` l.1888, `plot_hysteresis` l.1935, `plot_all` l.1974 | |
| Validation rapide | `quick_validation` l.2013-2051, `summary` l.2054-2062 | Vérifie notamment l'additivité tube+nylon. |

---

## 4. Conventions d'unités

**Système cohérent du moteur : mm – N – MPa – s – rad** (MPa·mm² = N ; MPa·mm³ = N·mm). Aucune unité n'est déclarée dans Base.py hors suffixes des clés de sortie ; les entrées portent leurs unités via les noms de `DEFAULT_SETTINGS`.

- **Contraintes/modules** : MPa partout (E, G, η en MPa·s). Pression en MPa de bout en bout (`pressure_MPa`).
- **Longueurs** : mm (Rout, Rin, rho0, initial_length, uncoiled_length, u, R). `h` en **mm/rad** (pas géométrique = 2πh, cf. clé `h_mm_per_rad` l.1454 et `h0_mm_per_rad` parametres.py l.736).
- **Angles** : entrées utilisateur en **degrés** (`alpha0_deg`, `theta_f_deg`), converties une seule fois en radians à l'init (l.358-359). Interne tout en radians ; sorties doubles `alpha_rad`/`alpha_deg` (l.1410-1411).
- **Forces** : internes en N (`Ft`, `Ftube`…) ; sorties `force_N` et `force_mN = 1000·Ft` (l.1415-1416). Décomposition tube/nylon en mN (l.1503-1504).
- **Couples** : internes en N·mm ; sorties `torque_Nmm` et `torque_microNm = 1000·Tt` (l.1417-1418) — correct car 1 N·mm = 1000 µN·m. Le résidu reste en N·mm (`max_abs_residual_Nmm`, l.2061).
- **dv** : rad/mm (torsion linéique) ; **dkappa** : 1/mm ; **dw** : adimensionnel.
- **Temps** : s. Débit `flow_rate_mL_min` en mL/min, converti par `60·volume/flow` (demi-période en s, l.1533 ; parametres.py l.850) ; vitesse de précontrainte en mm/min convertie l.1371 (`60·eps·L0/rate`).
- **Volume** : mL ; mm³→mL par ×1.0e-3 (parametres.py l.729).
- **Masse suspendue** : g côté interface, convertie en N par `masse_g × 1.0e-3 × 9.80665` (interface.py **l.536 et l.2104**). Le moteur ne reçoit que `load_N`.
- **Conversions de pression CSV** (pression.py l.106/142) : vers MPa — MPa×1, bar×0.1, kPa×0.001, psi×0.006894757293168361 ; plafond 1.5 MPa vérifié (pression.py l.149-150). `PSI_TO_MPA` aussi dans parametres.py l.17.
- **Conversions de force CSV** (pression.py l.107) : vers mN — mN×1, N×1000, g×9.80665, kg×9806.65 (poids). Temps : ms×1e-3 (l.105). Inférence d'unités depuis les en-têtes : `infer_pressure_unit`/`infer_force_unit`/`infer_time_unit` (pression.py l.45-73) — attention, `infer_force_unit` détecte « mn » dans n'importe quel en-tête minuscule.
- **Transmission interface → modèle** : les champs Streamlit sont saisis directement en mm/MPa/deg (interface.py l.826-849), rassemblés dans `settings` (l.1244-1278), validés/bornés par `normalize_settings`, puis convertis en dataclasses par `build_config` (parametres.py l.777-844) **sans aucune conversion d'unité** (mm reste mm, MPa reste MPa ; seuls le diamètre nylon → rayon l.812 et le choix de E_axial l.789-794 transforment les valeurs). Les runs passent par `run_model` (interface.py l.456-470), `run_relaxation_model` (l.478-497), `run_suspended_model` (l.533-556) ; vitesses de pression en MPa/s (`make_pressure_rate_history` l.371-391 : demi-période = Pmax/rate ; rampe suspendue l.508, l.2105).

---

## 5. Bizarreries et points d'attention relevés

1. **Défauts `E_axial` incohérents entre modules** : 37.76 MPa (Base.py l.53, = ΣE Maxwell) contre 31.24 MPa (parametres.py l.174 et l.225, valeur « table de l'article »). Sans passer par `build_config`, un appel direct `Base.run_blocked_actuation(...)` et une config interface en mode `paper_table` ne simulent pas le même matériau. `stiffness_sanity_report` (l.2004) utilise en dur 31.24 alors que le défaut du module est 37.76.
2. **Doubles jeux de défauts divergents** : `default_simulation_config` (Base.py l.98 : n_cycles=3, Pmax=1.3, dt=0.25, n_layers=4, nonlinear_pressure=False) vs `SimulationParams` (parametres.py l.254 : n_cycles=11, Pmax=1.5, dt=0.5, n_layers=8, nonlinear_pressure=True) vs `DEFAULT_SETTINGS` (n_layers=4, n_phi=16, nonlinear_pressure=False). Trois « défauts » différents selon le point d'entrée.
3. **`default_discretization` en grande partie mort** : n_layers=18/n_phi=72/pre_steps=120 et `dw_bracket=(-0.08,0.08)` ne servent que si on instancie `TCPAMaxwellBlockedModel` sans runner ; les runners forcent `dw_bracket=(-0.05,0.05)` (l.1675, 1759).
4. **Constantes magiques du profil de pression non linéaire** : `gamma_load=3.5`, `gamma_unload=2.8` (Base.py l.1526-1527 et re-dupliquées en dur dans parametres.py l.877-878) — phénoménologiques, non tracées vers l'article (le README les qualifie seulement de « phénoménologique »).
5. **Période d'hystérèse en dur** : `plot_hysteresis` (l.1945) suppose `period = 2·60·1.50/10.0 = 18 s` si non fournie — valable uniquement pour le débit/volume par défaut.
6. **Absence de poussée de fond hydraulique** : la pression n'entre que par la CL radiale interne (l.686) ; aucun terme axial P·π·Rin². À confronter à l'article BLOCKED (choix de modélisation, pas nécessairement un bug).
7. **Duplication quasi intégrale `_trial_state` / `_trial_state_from_geometry`** (l.781-874 vs 876-971) : ~90 lignes clonées ; tout correctif de physique doit être appliqué deux fois (risque de divergence).
8. **Incohérence potentielle en mode `updated`** : les contraintes sont évaluées aux anciens rayons (`self.R_centers`, l.796) mais intégrées avec les rayons/épaisseurs *nouveaux* (l.823) ; les angles de biais par couche sont recalculés puis les raideurs reconstruites au commit (l.1133) — cohérent seulement au premier ordre. Sans effet en mode `fixed` (défaut).
9. **Dégénérescence non traitée du problème radial** : dans `_layer_AB_mu` (l.576-584), `A = (C26−2C36)/(4C33−C22)` et `B = (C12−C13)/(C33−C22)` divisent par des différences qui tendent vers 0 pour un angle de biais faible (couches internes : C̄22→C̄33, μ→1, et R^μ devient colinéaire à la solution particulière B·dw·R). Cas limite μ=1 de Lekhnitskii non géré ; seul garde-fou : `mu = sqrt(max(C22/C33, 1e-14))` (l.581).
10. **Vitesse de précontrainte dupliquée en dur** : 20 mm/min est le défaut de `prestretch_to` (l.1361) et réapparaît sous forme du « /20.0 » dans le test de stabilité d'Euler pendant la précontrainte (parametres.py l.664-670). Si l'un change, l'autre ment.
11. **`run_hold_relaxation` (Base.py) ≠ `run_relaxation_model` (interface)** : la version Base utilise `nonlinear_ramp=True` par défaut (l.1560) alors que l'interface transmet `config.nonlinear_pressure` (interface.py l.487) — mêmes noms de sorties, rampes différentes.
12. **Attribut dynamique sur dataclass** : `run_model` écrit `config.measured_cycle_period_s` (interface.py l.464) hors définition de `SimulationParams` — fonctionne mais fragile (pickle/replace).
13. **`_find_dw` : politique multi-racines** : si plusieurs équilibres existent dans le bracket, la racine de plus petit |dw| est choisie silencieusement (l.992) ; aucune trace de multi-stabilité n'est remontée.
14. **Historique inclut la précontrainte** : `model.history` contient les pas de prestretch ; le re-zéro se fait par découpe `i_act0` dans les runners (l.1698-1702). Un usage direct de `history_arrays()` sur le modèle donne des séries incluant la précontrainte avec un temps non re-basé.
15. **`_equilibrate_suspended_load` force l'intégration exponentielle** (l.1721-1728) même si l'utilisateur a demandé `paper_explicit`, puis restaure — documenté nulle part hors du code.
16. **Résidus « erreur » sentinelles 1e6/1e100** dans les fonctions de coût (l.996, 1056, 1266) : standard mais peut piéger `least_squares` sur des plateaux ; compensé par les multi-départs.
17. **`infer_force_unit`** (pression.py l.58-66) : la détection `"mn" in header` matche aussi des mots contenant « mn » ; l'utilisateur peut corriger dans l'interface, risque faible mais réel de mauvaise unité par défaut.
18. **Commentaire l.2 du docstring** (« module autonome final ») et README cohérents avec le code ; pas de contradiction majeure commentaire/code détectée hormis le point 11.

**Chaîne d'appel de référence (interface)** : Streamlit → `settings` → `normalize_settings`/`settings_error` (parametres.py l.344/681) → `build_config` (l.777) → `Base.run_blocked_actuation` / `run_suspended_actuation` → `TCPAMaxwellBlockedModel.prestretch_to` → `step`/`step_blocked_series`/`step_suspended` → `history_arrays` → `add_corrected_output_conventions` → affichage/export (mN, µN·m, mm, MPa).