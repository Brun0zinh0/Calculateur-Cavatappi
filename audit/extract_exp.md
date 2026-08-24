# Rapport de référence — Article EXP
**« Twisted and coiled tube actuators driven by hydraulic pressure: Experiment and theory »**
Lei Liu, Zhiya Zhang, Dabiao Liu — *Thin-Walled Structures* 209 (2025) 112957, doi:10.1016/j.tws.2025.112957 — 14 pages (aucune annexe : le corps se termine p. 13, les références p. 14).
Mots-clés : pression hydraulique, actionneur torsadé-enroulé (TCPA), hystérésis, modèle de Maxwell généralisé.

**Note de lecture** : cet article est le compagnon « expérience + théorie » des deux autres PDF. Son modèle (Maxwell généralisé + description lagrangienne réactualisée + décomposition en couches concentriques) est structurellement identique à celui de l'article BLOCKED ; ici l'actionnement est **libre** (charge morte constante, extrémité basse libre), pas bloqué.

---

## 1. Résumé du modèle et de la démarche

Les auteurs fabriquent des actionneurs hydrauliques à partir d'un tube PVC étiré à froid (rapport 3:1), renforcé par un fil nylon 6 traversant, puis torsadé (TPA = actionneur de torsion) et enfin enroulé en hélice sur mandrin (TCPA = actionneur de contraction), avec recuit à 90 °C pendant 180 min et 10–15 cycles d'entraînement pour éliminer la déformation plastique. Deux bancs d'essai sont montés : (i) torsion libre d'un TPA vertical avec une palette (poids 0,02 N, inertie 1,126 kg·mm²) filmée par caméra, (ii) contraction d'un TCPA sous charge morte de 1 N mesurée par capteur laser ; l'eau est injectée par pousse-seringue (volume imposé, 1 à 20 mL/min) et la pression est mesurée par capteur. Les essais montrent un comportement viscoélastique marqué : fluage (la course de torsion continue de croître après arrêt de l'injection alors que la pression chute), hystérésis asymétrique sous cycles, et dérive avec le nombre de cycles. Le tube PVC est modélisé comme un composite biphasé (microfibrilles cristallines = ressorts, matrice amorphe = amortisseurs) par un modèle de Maxwell généralisé : un ressort C₀ en parallèle avec m = 3 branches de Maxwell, chaque composant étant transversalement isotrope avec **le même rapport d'anisotropie** (hypothèse simplificatrice assumée). Les paramètres (E₀…E₃, η₁…η₃) sont identifiés par moindres carrés sur des essais de relaxation à trois vitesses de déformation. Comme la déformation atteint ~20 %, la théorie des petites déformations est jugée insuffisante : le modèle est écrit en **incréments** (description lagrangienne réactualisée), la configuration courante servant de référence à chaque pas Δt. Le TPA est traité comme un cylindre axisymétrique multicouche (n couches concentriques, angle de biais variant « linéairement » via arctan), avec la solution de déplacement radial de Pipes & Hubert, conditions aux limites en pression et continuité inter-couches ; l'équilibre force/moment avec l'inertie de la palette donne la course de torsion γ = ∫v₀dt. Le TCPA est traité comme une hélice (courbure K = cos²β_h/R_h) : les incréments de déformation incluent la courbure initiale, les efforts du tube et du fil nylon sont sommés et équilibrés avec la charge dans le repère de Frenet–Serret, ce qui donne le déplacement δ puis la déformation d'actionnement ε_AF. Les prédictions reproduisent bien le fluage et les cycles (Figs. 10–11). Les études paramétriques montrent : la performance décroît quand la vitesse de pressurisation ou la charge augmente ; elle croît fortement avec l'anisotropie E_axial/E_radial ; l'angle de biais optimal de la couche intermédiaire reste ≈ 41° quelle que soit la géométrie ; la capacité de charge et la course se règlent par l'angle et le rayon d'hélice.

---

## 2. Inventaire exhaustif des équations

Notation : indices E/S = élastique/visqueux ; barre supérieure notée `X̄` = grandeur transformée dans le repère du cylindre (ou effective) ; gras = matrices/vecteurs (6×1 ou 6×6 en notation de Voigt, implicite).

### Modèle de Maxwell généralisé (§4.1)

**Éq. (1), p. 4** — Relations constitutives des éléments :
- Ressort : `σ₀ = C₀ ε`
- Branches de Maxwell : `σᵢ + (ηᵢ/Cᵢ) σ̇ᵢ = ηᵢ ε̇`, i = [1, m]

ε : matrice de déformation (–) ; m : nombre de branches de Maxwell (= 3 ici) ; σ₀, C₀ : matrice de contrainte (MPa) et matrice de raideur (MPa) du ressort ; σᵢ, Cᵢ, ηᵢ : contrainte (MPa), raideur (MPa) et viscosité (MPa·s) de la i-ème branche. Hypothèse : viscoélasticité linéaire par branche ; l'écriture matricielle (ηᵢ/Cᵢ apparaissant comme rapport) suppose implicitement que ηᵢ et Cᵢ sont proportionnelles (même anisotropie — cf. §4.2).

**Éq. (2), p. 4** — Contrainte totale : `σ = σ₀ + Σᵢ₌₁ᵐ σᵢ` (assemblage en parallèle, MPa).

**Éq. (3), p. 5** — Forme incrémentale de (1) :
`Δσ₀ = C₀ Δε` ; `σᵢ Δt + (ηᵢ/Cᵢ) Δσᵢ = ηᵢ Δε`, i = [1, m]
Δt : incrément de temps (s). C'est la discrétisation explicite (σᵢ évaluée au pas courant).

**Éq. (4), p. 5** — Relation constitutive incrémentale globale :
`Δσ = Δσ₀ + Σᵢ Δσᵢ = C_E Δε + C_S`
avec `C_E = C₀ + Σᵢ₌₁ᵐ Cᵢ` (raideur instantanée, MPa) et `C_S = −Δt Σᵢ₌₁ᵐ Cᵢ σᵢ / ηᵢ` (terme visqueux, MPa — dépend de l'état de contrainte courant). Deux contributions : partie élastique C_EΔε, partie visqueuse C_S. C'est l'équation-pivot du schéma incrémental.

### Cinématique du TPA (§4.3)

**Éq. (5), p. 6** — Ansatz de déplacement (torsion axisymétrique, incrément Δt) :
`u = u(r)`, `v = v₀ r z`, `w = w₀ z`
u : déplacement radial (mm) ; v : tangentiel ; w : axial ; v₀ : incrément de torsion par unité de longueur (rad/mm) ; w₀ : incrément de déformation axiale (–) ; r, z : coordonnées cylindriques dans la configuration de référence courante.

**Éq. (6), p. 6** — Incréments de déformation associés :
`Δε_r = du/dr`, `Δε_θ = u/r`, `Δε_z = w₀`, `Δε_zθ = v₀ r`, `Δε_rθ = Δε_rz = 0`
(sans dimension ; Δε_zθ en notation ingénieur).

**Éq. (7), p. 7** — Angle de biais de la couche j :
`α_j = arctan( (R_j / R_TO) · tan α_f )`
R_j : rayon de la j-ème couche concentrique (mm) ; R_TO : rayon extérieur courant du tube (mm) ; α_f : angle de biais à la surface extérieure (°). Le texte dit « varie linéairement dans la direction radiale » mais la formule exacte est en arctan (c'est tan α qui varie linéairement en r). NB : le texte imprime « radius of the *i*-th concentric layer » (coquille j/i).

**Éq. (8), p. 7** — Transformation de la matrice de raideur U (repère matériau, axes le long des fibres) vers Ū (repère du cylindre), rotation d'angle α_j, avec `p = cos α_j`, `q = sin α_j` :
```
Ū₁₁ = U₁₁p⁴ + U₂₂q⁴ + 2(U₁₂ + 2U₆₆)p²q²
Ū₁₂ = U₁₂(p⁴ + q⁴) + (U₁₁ + U₂₂ − 4U₆₆)p²q²
Ū₁₃ = U₁₃p² + U₂₃q²
Ū₁₆ = (U₁₂ − U₂₂ + 2U₆₆)pq³ + (U₁₁ − U₁₂ − 2U₆₆)p³q
Ū₂₂ = U₁₁p⁴ + U₂₂q⁴ + 2(U₁₂ + 2U₆₆)p²q²   [tel qu'imprimé — voir §6 : coquille quasi certaine]
Ū₂₃ = U₂₃p² + U₁₃q²
Ū₂₆ = (U₁₂ − U₂₂ + 2U₆₆)p³q + (U₁₁ − U₁₂ − 2U₆₆)pq³
Ū₃₃ = U₃₃
Ū₃₆ = (U₁₃ − U₂₃)pq
Ū₆₆ = U₆₆(p² − q²)² + (U₁₁ + U₂₂ − 2U₁₂)p²q²
```
(MPa). La construction de U à partir des constantes ingénieur (E_axial, E_radial, G₁₂, ν₁₂, ν₂₃) n'est pas donnée : renvoi à la réf. [25] (Yang & Li 2016).

**Éq. (9), p. 7** — Relation constitutive incrémentale dans le repère cylindre :
`Δσ = C̄_E Δε + C̄_S`
avec `C̄_E = C̄₀ + Σᵢ₌₁ᵐ C̄ᵢ`, `C̄_S = −Δt Σᵢ₌₁ᵐ C̄ᵢ σᵢ / η̄ᵢ` (mêmes formes que (4), matrices transformées par (8)).

**Éq. (10), p. 7** — Solution du déplacement radial dans la couche j (d'après Pipes & Hubert [35]) :
`u = C₁ r^μ + C₂ r^(−μ) + [(C̄_E26 − 2C̄_E36) / (4C̄_E33 − C̄_E22)] v₀ r² + [(C̄_E12 − C̄_E13) / (C̄_E33 − C̄_E22)] w₀ r`
avec `μ = sqrt(C̄_E22 / C̄_E33)` (–). C₁, C₂ : constantes d'intégration par couche (déterminées par (11)–(12)). Les C̄_Eab sont les composantes ab de C̄_E (MPa). Attention à l'indexation des axes sous-jacente (1 = direction fibre, 2/3 = transverses) héritée de [35].

**Éq. (11), p. 7** — Conditions aux limites en pression :
`Δσ_r |_(R=R_TI) = −ΔP`, `Δσ_r |_(R=R_TO) = 0`
ΔP : incrément de pression hydraulique interne (MPa) ; R_TI, R_TO : rayons intérieur/extérieur courants.

**Éq. (12), p. 8** — Continuité inter-couches :
`u^j |_(r=R_j) = u^(j+1) |_(r=R_j)` et `Δσ_r^j |_(R=R_j) = Δσ_r^(j+1) |_(R=R_j)`, j = [1, n−1]
n : nombre de couches (valeur numérique jamais donnée).

### Équilibre du TPA et course de torsion

**Éq. (13), p. 8** — Efforts résultants du tube (TPA, axisymétrique) :
`F_tube = 2π ∫₀^t1 ∫_RTI^RTO σ_z r dr dt` et `T_tube = 2π ∫₀^t1 ∫_RTI^RTO σ_θz r² dr dt`
F_tube : tension (N) ; T_tube : couple (N·mm). ⚠ L'intégrale sur le temps rend les unités littéralement incohérentes (N·s) : il faut comprendre l'accumulation temporelle des incréments de contrainte (σ dans l'intégrande = dσ accumulé sur l'historique) — voir §6. Le texte parle de « torque **H**_tube » mais l'équation écrit T_tube.

**Éq. (14), p. 8** — Équilibre des forces : `F = F_tag`
(F_tag = poids de la palette = 0,02 N ; la longueur du TPA est supposée quasi constante ⇒ inertie axiale négligée).

**Éq. (15), p. 8** — Équilibre des moments : `T = J_tag v₀ / Δt`
J_tag : inertie de rotation de la palette (kg·mm²) ; v₀/Δt ≈ accélération angulaire par unité de longueur — l'inertie n'est retenue que dans le bilan de moment (extrémité basse libre en rotation).

**Éq. (16), p. 8** — Course de torsion : `γ = ∫₀^t1 v₀ dt` (rad/mm).
La résolution combinée des Éqs. (6)–(15) donne u, v₀, w₀ à chaque pas.

### Modèle du TCPA (§4.4)

**Éq. (17), p. 8** — Courbure de l'hélice : `K(t) = cos²β_h(t) / R_h(t)`
R_h : rayon d'hélice (mm) ; β_h : angle d'hélice (°/rad). Configuration intermédiaire à l'instant t = référence.

**Éq. (18), p. 8** — Incréments de torsion et de courbure entre t et t+Δt :
`v₀ = sin2β_h(t+Δt) / (2R_h(t+Δt)) − sin2β_h(t) / (2R_h(t))`
`κ₀ = cos²β_h(t+Δt) / R_h(t+Δt) − cos²β_h(t) / R_h(t)`
(unités : 1/mm). Ce sont les formules classiques torsion/courbure d'une hélice (τ = sin2β/2R, κ = cos²β/R).

**Éq. (19), p. 9** — Incréments de déformation du tube torsadé avec courbure initiale (d'après [36]) :
```
Δε_r  = ∂u/∂r − ν̄₁₂ (κ₀ r cosθ + u K cosθ)/(1 + K r cosθ)
Δε_θ  = u/r  − ν̄₁₃ (κ₀ r cosθ + u K cosθ)/(1 + K r cosθ)
Δε_z  = w₀ + (κ₀ r cosθ + u K cosθ)/(1 + K r cosθ)
Δε_θz = v₀ r/(1 + K r cosθ) − ν̄₁₄ (κ₀ r cosθ + u K cosθ)/(1 + K r cosθ)
Δε_rθ = Δε_rz = 0
```
θ : coordonnée angulaire dans la section ; le dénominateur (1 + K r cosθ) traduit la métrique d'une poutre courbe. ν̄₁₂, ν̄₁₃, ν̄₁₄ : coefficients de Poisson effectifs.

**Éq. (20), p. 9** — Détermination des Poisson effectifs :
`C̄⁻¹ [Ē₁ 0 0 0 0 0]ᵀ = [1 −ν̄₁₂ −ν̄₁₃ −ν̄₁₄ 0 0]ᵀ`
Ē₁ : module d'Young effectif (MPa). La matrice C̄ utilisée ici (raideur 6×6 dans le repère courant) n'est pas explicitement reliée à C̄_E — laissé implicite.

**Éq. (21), p. 9** — Efforts résultants du tube du TCPA (intégration complète sur la section, sans axisymétrie) :
`F_tube = ∫₀^t1 ∫₀^2π ∫_RTI^RTO σ_z r dr dθ dt`
`M_tube = ∫₀^t1 ∫₀^2π ∫_RTI^RTO σ_z r² cosθ dr dθ dt`
`T_tube = ∫₀^t1 ∫₀^2π ∫_RTI^RTO σ_θz r² dr dθ dt`
F : tension (N), M : moment de flexion (N·mm), T : couple (N·mm). Même remarque d'unités que (13) sur ∫dt.

**Éq. (22), p. 9** — Contributions du fil nylon 6 (élastique linéaire, rayon r_nylon supposé constant pendant l'actionnement) :
`F_nylon = π E_nylon r²_nylon w₀` ; `M_nylon = (1/4) π E_nylon r⁴_nylon k₀` ; `T_nylon = (1/2) π G_nylon r⁴_nylon v₀`
(k₀ imprimé — vraisemblablement κ₀ de l'Éq. (18)). Formules poutre circulaire : EA·ε, EI·κ, GJ·τ.

**Éq. (23), p. 10** — Équilibre dans le repère de Frenet–Serret {n, b, t} (inertie négligée : accélération ~10⁻⁴ m/s² ≪ g) :
`F_tube + F_nylon = F_load sinβ_h`
`M_tube + M_nylon = −F_load R_h sinβ_h`
`T_tube + T_nylon = F_load R_h cosβ_h`
F_load : charge morte (N). La combinaison (21)–(23) donne w₀, θ_h (sic — probablement β_h) et R_h au cours de l'actionnement.

**Éq. (24), p. 10** — Déplacement axial du TCPA :
`δ = (1 + w₀) L_T0 · sinβ_h / sinβ_h0 − L_T0`
L_T0 : longueur initiale du TCPA (mm) ; β_h0 : angle d'hélice initial. (δ > 0 = allongement.)

**Éq. (25), p. 10** — Déformation d'actionnement :
`ε_AF = (δ(P=0) − δ(P)) / L_T0`
(positive en contraction sous pression ; négative si allongement par fluage).

### Relations importantes non numérotées

| Page | Relation | Rôle |
|---|---|---|
| p. 2 | `γ = φ / L_P0` | course de torsion mesurée : φ = angle de la palette (rad), L_P0 = longueur initiale du TPA (mm) → rad/mm |
| p. 3 | `ε_AF = (L_T1(P) − L_T1(P=0)) / L_T0` | définition **expérimentale** de la déformation d'actionnement (L_T1 = longueur courante) — signe opposé à l'Éq. (25), cf. §6 |
| p. 6 | `E_axial = E₀ + E₁ + E₂ + E₃` | lien module macroscopique / branches de Maxwell (réponse instantanée) |
| p. 7 | `p = cos α_j`, `q = sin α_j` | pour l'Éq. (8) |
| p. 7 | `μ = sqrt(C̄_E22/C̄_E33)` | exposant de l'Éq. (10) |

---

## 3. Tables de paramètres

### Table 1 (p. 2) — Géométrie des actionneurs fabriqués

| Type | Paramètre | Valeur |
|---|---|---|
| TPA | Diamètre extérieur initial du tube PVC 2R_TO0 | 1,87 ± 0,02 mm |
| TPA | Diamètre intérieur initial 2R_TI0 | 0,80 ± 0,02 mm |
| TPA | Angle de biais initial α_f0 | 16,42 ± 0,01° |
| TCPA | Diamètre du fil nylon 6, 2r_nylon | 0,77 ± 0,02 mm |
| TCPA | Diamètre extérieur initial 2R_TO0 | 2,00 ± 0,02 mm |
| TCPA | Diamètre intérieur initial 2R_TI0 | 0,80 ± 0,02 mm |
| TCPA | Angle de biais initial α_f0 | 35,69 ± 0,01° |
| TCPA | Rayon d'hélice initial R_h0 | 2,66 ± 0,02 mm |
| TCPA | Angle d'hélice initial β_h0 | 6,87 ± 0,01° |

### Table 2 (p. 6) — Paramètres de traction du modèle de Maxwell (identifiés par moindres carrés sur la relaxation, m = 3)

| Paramètre | Valeur |
|---|---|
| Module du ressort parallèle E₀ | 6,36 MPa |
| Module de la branche 1, E₁ | 20,67 MPa |
| Module de la branche 2, E₂ | 5,98 MPa |
| Module de la branche 3, E₃ | 4,75 MPa |
| Viscosité de la branche 1, η₁ | 154,57 MPa·s |
| Viscosité de la branche 2, η₂ | 977,79 MPa·s |
| Viscosité de la branche 3, η₃ | 11 044,83 MPa·s |

(Temps de relaxation implicites τᵢ = ηᵢ/Eᵢ ≈ 7,5 s ; 163,5 s ; 2325 s.)

### Table 3 (p. 6) — Paramètres macroscopiques matériau (PVC étiré + recuit, nylon 6)

| Paramètre | Valeur |
|---|---|
| Module axial du tube PVC E_axial | 37,8 MPa |
| Module radial du tube PVC E_radial | 9,1 MPa |
| Module de cisaillement du PVC G₁₂ | 7,2 MPa |
| Poisson du PVC ν₁₂ | 0,205 (repris de [27]) |
| Poisson du PVC ν₂₃ | 0,422 (repris de [27]) |
| Module du fil nylon 6 E_nylon | 3694 MPa |
| Module de cisaillement du nylon G_nylon | 790 MPa |

Vérification : E₀+E₁+E₂+E₃ = 37,76 ≈ E_axial = 37,8 MPa. Anisotropie E_axial/E_radial ≈ 4,15.

### Autres valeurs dispersées dans le texte

| Grandeur | Valeur | Page |
|---|---|---|
| Tube PVC brut : OD / ID / longueur | 3,20 ± 0,02 mm / 1,60 ± 0,02 mm / ~100 mm | 2 |
| Tubes acier d'extrémité : OD / ID | 1,60 / 0,80 ± 0,02 mm | 2 |
| Fil nylon 6 : longueur / diamètre | 300 mm / 0,77 ± 0,02 mm | 2 |
| Étirage à froid | jusqu'à ~300 mm, rapport 3:1 | 2 |
| Recuit | 90 °C, 180 min (tube aplani : 90 °C, 3 h) | 2, 6 |
| Entraînement | 10 à 15 cycles | 2 |
| Longueur initiale TPA L_P0 | 51,88 ± 0,10 mm | 2 |
| Longueur initiale TCPA L_T0 | 25,07 ± 0,10 mm | 3 |
| Poids palette F_tag / inertie J_tag | 0,02 N / 1,126 kg·mm² (1,13 dans Fig. 12b) | 2 |
| Charge TCPA F_load | 1 N (0,2–1,0 N en simulation) | 2, 12 |
| Pompe / capteur pression / caméra / laser | LEAD FLUID TYD02-01 / NanJing Aire AE-T / GP-460H / Keyence LK-G80 | 2–3 |
| Essai de relaxation : étirement | 30 → 31,2 mm (4 %) aux vitesses 2×10⁻⁴, 1×10⁻³, 5×10⁻³ /s | 6 |
| Contraintes au pic de relaxation | 0,61 / 0,81 / 1,19 MPa → convergence ~0,3 MPa après ~60 min | 6 |
| Testeur de traction | CARE IBTC-300SL | 5 |
| Rayon d'hélice pour Fig. 11(c) | R_h0 = 2,15 mm | 10 |
| Angle de biais alternatif testé | α_f0 = 40,02° (Fig. 10(c)) | 8, 10 |

---

## 4. Hypothèses et domaine de validité déclarés

1. **Matériau** : PVC = composite biphasé (microfibrilles cristallines ↔ élasticité, matrice amorphe ↔ viscosité) ; Maxwell généralisé à m = 3 branches suffit à décrire la relaxation.
2. **Même anisotropie pour tous les composants** (ressorts et amortisseurs transversalement isotropes avec le même rapport) — hypothèse explicitement identifiée comme cause de l'écart théorie/expérience : en réalité la phase amorphe (visqueuse) est moins anisotrope, d'où course de torsion sous-estimée en pressurisation et surestimée en fluage, et boucle d'hystérésis théorique plus petite que l'expérimentale (p. 10).
3. **Grandes déformations traitées par incréments** : description lagrangienne réactualisée ; Δt suffisamment petit pour que la déformation par pas soit « mineure » ; déformation totale jusqu'à ~20 %.
4. **TPA** : torsion axisymétrique ; ansatz (5) ; n couches concentriques ; tan(angle de biais) linéaire en r, maximum α_f en surface ; longueur ≈ constante pendant la torsion (inertie axiale négligée) ; inertie uniquement dans le bilan de moment (palette) ; rotation libre en bas.
5. **TCPA** : hélice à courbure K = cos²β_h/R_h ; effets de la courbure initiale via (19) ; rayon du fil nylon constant pendant l'actionnement ; nylon élastique linéaire ; inertie négligeable (accélération ~10⁻⁴ m/s² ≪ g) ; équilibre quasi-statique dans le repère de Frenet–Serret.
6. **Auto-contact** : la raideur en compression du TCPA auto-contacté ≫ raideur en traction ⇒ la déformation d'actionnement max est atteinte à l'auto-contact ; à β_h0 = 10° l'espacement inter-spires empêche l'auto-contact à 1,2 MPa sans charge.
7. **Domaine exploré** : pression 0–1,2 MPa ; vitesses de pressurisation 1–20 mL/min (essais) et 0,01–1 MPa/s (simulations) ; charges 0,2–1 N ; angles de biais 0–70° ; anisotropie E_axial/E_radial 3–8. Les actionneurs doivent être « entraînés » (10–15 cycles) : le modèle décrit le régime stabilisé, pas le premier cycle vierge.
8. L'entrée du modèle est la **courbe pression–temps mesurée** (« substituting the pressurization-time curve ») : la relation volume injecté → pression n'est PAS modélisée.

---

## 5. Figures clés

- **Fig. 1 (p. 2)** : photos de fabrication — (a) tube PVC + tubes acier + fil nylon ; (b) étirage ; (c) torsadage ; (d) enroulement. Échelles 10–30 mm.
- **Fig. 2 (p. 3)** : schémas des deux bancs — (a) torsion : tube précurseur pendu, palette (« Tag »), caméra en dessous, pression par le haut ; (b) contraction : TCP muscle + charge + capteur laser.
- **Fig. 3 (p. 3)** : **fluage du TPA**. 3 sous-graphes vs temps (0–300 s) : volume d'eau (0→2 mL, rampe à 5 mL/min puis palier), pression (MPa, pic **1,238 MPa** à 2 mL puis décroissance par fluage), course de torsion (rad/mm, **0,0447** au pic de pression puis croissance continue : **0,0734 / 0,0815 / 0,0868 / 0,0915** rad/mm à 1/2/3/4 min). Point crucial : la pression baisse pendant que la torsion continue de monter (signature viscoélastique).
- **Fig. 4 (p. 4)** : **cycles de torsion du TPA** (1 mL/min, 0↔2 mL, 3 cycles, 0–750 s). Volume et pression périodiques symétriques (pic ~1,2 MPa) ; course asymétrique : **0,0117** rad/mm à 0,5 MPa en montée vs **0,0325** en descente (hystérésis) ; résidu **0,00134** rad/mm à P = 0 ; pics à 2 mL : **0,0511 / 0,0529 / 0,0555** rad/mm (dérive de fluage).
- **Fig. 5 (p. 5)** : **cycles de contraction du TCPA** (5 mL/min, 0↔3 mL, 10 cycles, 0–720 s). Axe « Tensile actuation » de −0,4 à 0,1 (sans unité) : pendant la pressurisation (0–36 s) l'actionnement **diminue d'abord puis augmente** (compétition fluage sous poids / contraction sous pression, minimum −0,236 à 18 s) ; en dépressurisation (36–72 s) il diminue. Longueur moyenne dérive avec les cycles (fluage). (b) photos de la charge (règle 0–16 cm).
- **Fig. 6 (p. 6)** : schéma du modèle de Maxwell généralisé (C₀ ∥ (η₁,C₁) ∥ (η₂,C₂) ∥ … ∥ (η_m,C_m)).
- **Fig. 7 (p. 6)** : relaxation de contrainte, contrainte (MPa) vs temps (0–3600 s) aux 3 vitesses de déformation avec ajustements (pics 1,19/0,81/0,61 MPa → ~0,3 MPa) ; valide m = 3.
- **Fig. 8 (p. 7)** : (a) angle de biais α_f sur le tube (coordonnées r, θ, z) ; (b) section multicouche : R_TI, R₁, …, R_{n−1}, R_TO, déplacements u et v₀z.
- **Fig. 9 (p. 7)** : géométrie du TCPA : R_h, β_h, 2R_TO, F_load aux deux extrémités, repère (n, b, t) ; triangle de l'hélice reliant L_T0, L_T1, β_h0, β_h et δ (δ = variation de hauteur ; L/tanβ = développement).
- **Fig. 10 (p. 8)** : **validation TPA**. (a) fluage : théorie vs essai, course (0–0,10 rad/mm) vs temps (0–300 s) — bon accord, théorie légèrement au-dessus après arrêt d'injection ; (b) cycles à 1 mL/min (pics ~0,05–0,056 rad/mm) — bon accord ; (c) cycles à 5 mL/min pour α_f0 = 16,42° et 40,02° : pic **0,0437** rad/mm à 5 mL/min contre **0,0515** à 1 mL/min (16,42°), et **0,0603** rad/mm à 40,02° (courbes rouges plus hautes) → la course croît avec l'angle de biais et décroît avec la vitesse de pressurisation.
- **Fig. 11 (p. 9)** : **validation TCPA**. (a) fluage sous 1 N : déformation → −0,34 à 600 s (élongation 0,262 à 72 s) ; (b) 10 cycles (0–720 s), pics ~0 / creux ~−0,35, accord théorie/essai ; l'actionnement à 3 mL aux cycles 2/5/10 : **0,28 / 0,287 / 0,284** (≈ constant) ; (c) hystérésis déformation–pression (0–1,2 MPa) à 5 et 20 mL/min (R_h0 = 2,15 mm) : boucle théorique plus petite que l'expérimentale (anisotropie visqueuse) ; (d) cycles pour R_C0 = 2,15 vs 2,66 mm : la déformation croît avec le rayon d'hélice (pics ~0,27 vs ~0,17).
- **Fig. 12 (p. 10)** : simulations TPA. (a) course vs pression pour Ṗ = 1 / 0,1 / 0,01 MPa/s : courses **0,0266 / 0,0403 / 0,0782** rad/mm ; boucle d'hystérésis et course croissent quand Ṗ décroît ; (b) effet de l'inertie J_tag = 1,13 / 10 / 100 kg·mm² : l'inertie freine en pressurisation, amplifie en dépressurisation, aire de boucle croît avec J_tag.
- **Fig. 13 (p. 11)** : course vs angle de biais (0–70°, P = 1,2 MPa). (a) Ṗ = 0,01…1 MPa/s : maximum ~0,15 rad/mm à 0,01 MPa/s ; l'optimum (~48–50°) varie peu avec Ṗ ; (b) anisotropie E_axial/E_radial = 3…8 (à 0,1 MPa/s) : course max de ~0,042 à ~0,155 rad/mm ; l'angle optimal **décroît** quand l'anisotropie croît.
- **Fig. 14 (p. 11)** : géométrie du TPA vs angle de biais. (a) R_TI0 constant, R_TO0 = 0,735…1,535 mm : épaissir le tube réduit la course, l'angle optimal augmente ; (b) épaisseur constante (R_TO0−R_TI0), diamètre croissant : course augmente, angle optimal diminue. Conclusion clé : **l'angle de biais optimal de la couche intermédiaire ≈ 41°**, indépendant de Ṗ et de la géométrie (p. 12).
- **Fig. 15 (p. 12)** : TCPA, déformation vs pression (0–1,2 MPa). (a) Ṗ = 1/0,1/0,01 MPa/s sous 1 N : à Ṗ décroissant, déformation et boucle croissent (à 0,01 MPa/s : ε = −0,037 à P = 0 par fluage) ; à 1 MPa/s la boucle est minime (≈ actionnement linéaire) ; (b) F_load = 0,2/0,6/1,0 N à 0,1 MPa/s : la déformation décroît avec la charge ; boucle du TCPA nettement plus petite que celle du TPA (élasticité linéaire du nylon).
- **Fig. 16 (p. 12)** : TCPA, déformation vs angle de biais. (a) l'angle optimal **augmente** avec Ṗ (max ~0,38 à 0,01 MPa/s vers ~40°) ; (b) quasi insensible à la charge (max ~0,29–0,37 vers ~45°).
- **Fig. 17 (p. 13)** : TCPA, déformation vs charge (0–1 N, 0,1 MPa/s, 1,2 MPa). (a) R_C0 = 2,00/2,66/3,00 mm : montée jusqu'à l'auto-contact puis lente décroissance ; pour R_h0 = 2,66 mm, maximum ε_AF^max à **F_load = 0,21 N** ; ε_AF ∈ [0,9·ε^max, ε^max] pour F_load ∈ [0,18 ; 0,65] N ; rayon ↑ ⇒ déformation ↑ mais capacité de charge ↓ ; (b) θ_h0 = 6,87/8/10° : angle d'hélice ↑ ⇒ capacité de charge ↑, déformation ↓ ; à 10° pas d'auto-contact à 1,2 MPa sans charge.

---

## 6. Ambiguïtés et points laissés implicites (à l'attention des auditeurs)

1. **Coquille quasi certaine dans l'Éq. (8)** : Ū₂₂ est imprimée **identique à Ū₁₁** (`U₁₁p⁴ + U₂₂q⁴ + 2(U₁₂+2U₆₆)p²q²`). La transformation orthotrope standard donne `Ū₂₂ = U₁₁q⁴ + U₂₂p⁴ + 2(U₁₂+2U₆₆)p²q²` (échange p↔q). Vérifier ce que fait le code : s'il recopie la formule imprimée, il reproduit la coquille.
2. **Intégrales temporelles des Éqs. (13) et (21)** : telles qu'écrites, F = 2π∫₀^t1∫σ_z r dr **dt** a la dimension d'une force×temps. Il faut comprendre σ comme le cumul des incréments Δσ sur l'historique (l'intégrale en temps matérialise l'accumulation incrémentale), c'est-à-dire F(t₁) = 2π∫ σ_z(t₁) r dr avec σ_z(t₁) = Σ Δσ_z. Le schéma exact d'accumulation n'est pas écrit.
3. **Deux définitions de ε_AF de signes opposés** : p. 3 (expérimental) `ε_AF = (L_T1(P) − L_T1(P=0))/L_T0` vs Éq. (25) `ε_AF = (δ(P=0) − δ(P))/L_T0`. Les figures 5 et 11(a,b) tracent des valeurs négatives (élongation), les figures 11(c,d), 15, 16, 17 des valeurs positives (contraction). La convention effective de chaque figure doit être déduite du contexte.
4. **Matrice U non explicitée** : la construction de la raideur transversalement isotrope U à partir de (E_axial, E_radial, G₁₂, ν₁₂, ν₂₃) est renvoyée à [25] ; de même la répartition de l'anisotropie sur chaque branche (C̄ᵢ, η̄ᵢ) est seulement décrite verbalement (« même anisotropie pour chaque composant »).
5. **Éq. (20)** : la matrice C̄ (6×6) utilisée pour extraire ν̄₁₂, ν̄₁₃, ν̄₁₄ et Ē₁ n'est pas définie explicitement (C̄_E ? raideur totale courante ?) ; l'ordre des composantes du vecteur (r, θ, z, θz, …) est implicite.
6. **Paramètres numériques du schéma jamais donnés** : nombre de couches n, pas de temps Δt, critère de convergence de l'itération sur la configuration de référence. Aucune indication d'implémentation.
7. **Entrée pression** : le modèle consomme P(t) mesuré ; la relation volume injecté ↔ pression (compressibilité, fuite, conformité du circuit) n'est pas modélisée. Les simulations des §5.2–5.3 utilisent des rampes en MPa/s (0,01–1), pas des volumes.
8. **Incohérences de notation** : couple noté H_tube dans le texte vs T_tube dans (13) ; k₀ dans (22) vs κ₀ dans (18) ; R_C0 (Figs. 11d, 17a) vs R_h0 (texte, Table 1) pour le rayon d'hélice ; θ_h0 (Fig. 17b) vs β_h0 pour l'angle d'hélice ; « R_j radius of the i-th layer » (Éq. 7) ; l'Éq. (23) est suivie de « we obtain w₀, θ_h, R_h » (θ_h non défini, lire β_h).
9. **Légendes de figures erronées** : Fig. 5(a) légende « torsional stroke » alors que l'axe est « Tensile actuation » ; Fig. 10 : deux sous-figures étiquetées « (b) » (la seconde est (c), à 5 mL/min).
10. **Angle de biais initial du TPA** : Table 1 donne α_f0 = 16,42° (TPA) et 35,69° (TCPA) ; Fig. 10(c) compare 16,42° et 40,02° — l'échantillon à 40,02° n'apparaît pas dans la Table 1.
11. **J_tag** : 1,126 kg·mm² (texte p. 2) vs 1,13 kg·mm² (Fig. 12b) — arrondi.
12. **« Variation linéaire » de l'angle de biais** (p. 6) contredite par la forme arctan de l'Éq. (7) (c'est tanα qui est linéaire en r) — le code doit choisir l'une des deux.
13. **Auto-contact** : le traitement mécanique de l'auto-contact (raideur de compression « significativement plus grande ») n'est pas formalisé par une équation — la réf. [36] (Liu et al., Smart Mater. Struct. 33 (2024) 065022) est invoquée ; le critère d'arrêt de contraction du modèle est implicite.
14. **État initial après entraînement** : le modèle part d'une configuration « stabilisée » (après 10–15 cycles) mais l'état de contrainte interne initial (σᵢ(0) des branches de Maxwell, précontrainte inter-spires évoquée p. 4) n'est pas spécifié.
15. **Fig. 11(c)** : R_h0 = 2,15 mm pour ces essais (différent de la Table 1 : 2,66 mm) — deux spécimens de TCPA distincts coexistent donc dans l'article sans que la géométrie complète du second (angle d'hélice, biais) soit tabulée.

**Fichiers utiles** : pages rendues en PNG dans `C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad\pages\p01.png … p14.png` (zooms d'équations : `eq8_zoom.png`, `eq10_zoom.png`, `eq19_zoom.png` ; texte brut : `exp_fulltext.txt`).