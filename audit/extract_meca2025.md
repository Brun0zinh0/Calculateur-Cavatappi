# Rapport de référence — Article « MECA2025 »

**Fichier source** : `C:/Users/b.pereiraazevedo/OneDrive - House Of HR NV/Documents/Stage muscle artificièle/modèle/modèle chinois/Espace de travail/Article support/2025--Mechanicsoftwistedandcoiledtube-basedartificialmusclesdrivenbyhydraulicpressure.pdf` (11 pages PDF ; p. 1 = couverture ResearchGate, p. 2–11 = pages d'article 1420203:1 à 1420203:10).

**Référence bibliographique** : Liu H., Zhang Z.Y., Liu L., Liu D.* — *Mechanics of twisted and coiled tube-based artificial muscles driven by hydraulic pressure*, Science China Technological Sciences, avril 2025, Vol. 68, Iss. 4, art. 1420203, DOI 10.1007/s11431-024-2871-y. Huazhong University of Science and Technology, Wuhan. Reçu 04/07/2024, accepté 02/01/2025, publié en ligne 12/03/2025.

**Convention de pagination du présent rapport** : « p. N » = page du fichier PDF ; la page d'article correspondante est N−1 (ex. p. 6 PDF = page 1420203:5).

**Important pour l'audit** : cet article traite de l'**actionnement libre** (rotation libre du tube torsadé, contraction libre du muscle TCP sous charge morte constante) avec un modèle **élastique linéaire transversalement isotrope, incrémental**. Il ne contient **aucun modèle viscoélastique de Maxwell** (celui-ci est dans l'article BLOCKED). La « Supporting Information » (détails de fabrication et d'essais) est un document en ligne séparé, **absent de ce PDF** — les auditeurs ne peuvent pas s'y référer via ce fichier.

---

## 1. Résumé du modèle et de la démarche

L'article étudie, par expérience, théorie et simulation éléments finis (Abaqus), la mécanique et l'actionnement de muscles artificiels fabriqués à partir de tubes en PVC torsadés (actionneurs de torsion) puis enroulés en hélice (actionneurs de traction, dits TCP ou « cavatappi »), actionnés par pression hydraulique interne. Les tubes PVC (renforcés intérieurement par un monofilament de nylon 6) sont étirés à froid à 3× leur longueur initiale, ce qui aligne les chaînes moléculaires et rend le matériau transversalement isotrope (rapport E1/E2 passant de 1 à 3,64), puis recuits à 90 °C pendant 90 min après chaque étape (étirage, torsion, enroulement) pour éliminer les contraintes résiduelles. La pression interne produit une contrainte de cisaillement σ_sφ dans la paroi du tube torsadé, qui provoque un dé-torsadage local : ce dé-torsadage se manifeste comme une rotation pour le tube torsadé et comme une contraction du pas d'hélice pour le muscle TCP. Le modèle phénoménologique décrit le tube comme une poutre à courbure initiale K en coordonnées curvilignes hélicoïdales, avec un tenseur de déformation issu du tenseur de Cauchy-Green (formulation de la réf. [27] des auteurs). La paroi du tube est discrétisée en n couches cylindriques concentriques homogénéisées, chacune transversalement isotrope avec la direction matérielle inclinée de l'angle de biais local β (maximal β_f = 40° en surface) ; la solution radiale de chaque couche suit la forme analytique de Pipes & Hubert [31]. La résolution est incrémentale (lagrangienne actualisée) : la pression croît linéairement de 0 à P sur un temps t, et à chaque pas Δt la configuration courante devient la nouvelle configuration de référence. Les conditions aux limites sont la continuité du déplacement radial et de la contrainte radiale aux 2n−2 interfaces, −ΔP sur la paroi interne et 0 sur la paroi externe. Pour le tube torsadé (K = 0), les conditions T = 0 et H = 0 (extrémité basse libre en rotation, cible de masse négligeable) donnent la déformation de rotation ε_Rot = Θ·R_out/L0, linéaire en pression au-delà de 0,5 MPa. Pour le muscle TCP, la géométrie hélicoïdale (angle α, rayon ρ) évolue à nombre de spires constant ; l'équilibre en effort/moment/couple sur la section (avec la contribution du mandrin de nylon, non négligeable en traction) et la cinématique δ = L·sinα − L0·sinα0 donnent la déformation d'actionnement ε_Ten = δ/(L0·sinα0), qui croît non linéairement avec P. La théorie prédit un maximum d'actionnement en torsion à un angle de biais critique β_f ≈ 47° sous 1,4 MPa, et un actionnement en traction croissant avec le diamètre d'enroulement (paramètre C = (R_core+R_out)/R_out), au prix d'une capacité de charge réduite. Les écarts théorie/expérience (< 6 % en torsion, < 20 % en traction, maximaux vers 0,5 MPa) sont attribués au fluage viscoélastique du PVC et à la non-linéarité de la montée en pression réelle, le modèle étant strictement élastique linéaire en petites déformations.

---

## 2. Inventaire exhaustif des équations

Notation : `Φ, φ` angles polaires de section (référence/courant) ; `K, κ` courbures (référence/courant) ; `κ̄` incrément de courbure ; barre sur C et ν = grandeurs effectives ; indices T = tension/torsion, B = flexion.

### 3.1 — Cadre théorique du tube torsadé (p. 5 PDF = page article 4)

**Éq. (1)** — Coordonnées curvilignes de référence :
`X^1 = R,  X^2 = Φ,  X^3 = S`
- R : coordonnée radiale dans la section [m] ; Φ : angle polaire de section [rad] ; S : abscisse curviligne le long de l'axe du tube [m]. Configuration de référence (hélice de courbure K).

**Éq. (2)** — Vecteurs de base covariants de référence :
```
G1 = −sinΦ e1 + cosΦ cos(KS) e2 + cosΦ sin(KS) e3
G2 = −R cosΦ e1 − R sinΦ cos(KS) e2 − R sinΦ sin(KS) e3
G3 = −(1 + K R cosΦ) sin(KS) e2 + (1 + K R cosΦ) cos(KS) e3
```
- {e1, e2, e3} : repère cartésien global ; K : courbure initiale du tube [1/m]. Hypothèse : tube modélisé comme tore/hélice de courbure K (torsion géométrique de l'hélice non explicitée dans la métrique — voir §6).

**Éq. (3)** — Tenseur métrique de référence :
`[G_ij] = diag[1, R², (1+KR cosΦ)²]` ; `[G^ij] = diag[1, R⁻², (1+KR cosΦ)⁻²]`
- Sans unité (composantes mixtes) ; la base est orthogonale.

**Éq. (4)** — Coordonnées curvilignes courantes : `x^1 = r, x^2 = φ, x^3 = s`.

**Éq. (5)** — Vecteurs de base covariants courants (mêmes formes que (2) avec r, φ, κ, s) :
```
g1 = −sinφ e1 + cosφ cos(κs) e2 + cosφ sin(κs) e3
g2 = −r cosφ e1 − r sinφ cos(κs) e2 − r sinφ sin(κs) e3
g3 = −(1 + κ r cosφ) sin(κs) e2 + (1 + κ r cosφ) cos(κs) e3
```
- κ : courbure du tube en configuration courante [1/m].

**Éq. (6)** — Métrique courante : `[g_ij] = diag[1, r², (1+κ r cosφ)²]`.

**Éq. (7)** — Élément d'arc de référence :
`dS² = G_ij dX^i dX^j = c_ij dx^i dx^j`, avec `c_ij = G_αβ X^α,_i X^β,_j`, `X^α,_i = ∂X^α/∂x^i`, `x^α,_i = ∂x^α/∂X^i`.
- c_ij : métrique de référence tirée en coordonnées courantes (pull-back).

**Éq. (8)** — Tenseur de déformation de Cauchy-Green (droit, tel que nommé par les auteurs) :
`(c⁻¹)^i_j = (c⁻¹)^{ik} g_kj = G^{αβ} x^i,_α x^k,_β g_kj`

**Éq. (9)** (p. 6 PDF) — Tenseur des déformations :
`ε^i_j = ½ [ (c⁻¹)^i_j − δ^i_j ]`
- δ^i_j : Kronecker. **Attention** : forme de type Almansi mais avec le signe (c⁻¹ − I)/2, et non (I − c⁻¹)/2 ; convention à vérifier vis-à-vis du code (voir §6). Sans unité.

**Schéma incrémental (texte, p. 6 PDF, non numéroté)** : la pression croît linéairement de 0 à P sur un temps t ; à chaque pas Δt, l'état à t_i sert de configuration initiale pour résoudre t_i + Δt (lagrangien actualisé). Les variables de position (R, Φ, S, K) de la configuration intermédiaire à t_i sont supposées connues.

**Éq. (10)** — Relation référence → courant sur un incrément :
`r = R + u,  φ = Φ + vS,  s = (1+w)S,  κ = K + κ̄`
- u : déplacement radial [m] ; v : angle tangentiel par unité de longueur (incrément de taux de torsion) [rad/m] ; w : déformation axiale (incrément) [–] ; κ̄ : variation de courbure sur Δt [1/m].

**Décomposition (texte, non numérotée)** : `u = u_T(R) + u_B(R, Φ)` — u_T dû à la tension/torsion (axisymétrique), u_B dû à la flexion. Hypothèse : influence de u_B sur ε_s ignorée ; u_B n'est jamais explicité dans l'article.

**Éq. (11)** — Incréments de déformations (d'après la réf. [27] des auteurs) :
```
Δε_r  = ∂u_T/∂R − ν̄12 · (κ̄ R cosΦ + u_T K cosΦ)/(1 + K R cosΦ)
Δε_φ  = u_T/R  − ν̄13 · (κ̄ R cosΦ + u_T K cosΦ)/(1 + K R cosΦ)
Δε_s  = w + (κ̄ R cosΦ + u_T K cosΦ)/(1 + K R cosΦ)
Δε_φs = vR/(1 + K R cosΦ) − ν̄14 · (κ̄ R cosΦ + u_T K cosΦ)/(1 + K R cosΦ)
Δε_rs = Δε_rφ = 0
```
- ν̄12, ν̄13, ν̄14 : coefficients de Poisson effectifs (couplage axial→radial, axial→circonférentiel, axial→cisaillement). Le terme commun (κ̄RcosΦ + u_T K cosΦ)/(1+KRcosΦ) est la part de déformation axiale due à la courbure. Sans unité.

**Éq. (12)** — Définition des Poisson effectifs :
`C̄⁻¹ [Ē1, 0, 0, 0]^T = [1, −ν̄12, −ν̄13, −ν̄14]^T`
- Ē1 : module d'Young axial effectif [Pa] ; C̄ : matrice de rigidité effective 4×4 (formulation de la réf. [25], Yang & Li). Hypothèse : chaque couche homogénéisée transversalement isotrope, direction matérielle inclinée de l'angle de biais local β(R).

**Éq. (13)** — Loi de comportement incrémentale par couche (ordre des composantes s, φ, r, sφ) :
```
[Δσ_s ]   [C̄11 C̄12 C̄13 C̄16] [Δε_s ]
[Δσ_φ ] = [C̄12 C̄22 C̄23 C̄26] [Δε_φ ]
[Δσ_r ]   [C̄13 C̄23 C̄33 C̄36] [Δε_r ]
[Δσ_sφ]   [C̄16 C̄26 C̄36 C̄66] [Δε_sφ]
```
- Δσ en [Pa]. Matrice symétrique 4×4 ; les C̄ dépendent de la couche (angle β_i). Le tube est découpé en n couches cylindriques concentriques (Fig. 4(d)) ; n **n'est jamais spécifié**.

**Éq. (14)** — Déplacement radial de la couche i (solution de Pipes & Hubert [31]) :
`u_T^i = C1^i R^μ + C2^i R^(−μ) + [(C̄26 − 2C̄36)/(4C̄33 − C̄22)] v R² + [(C̄12 − C̄13)/(C̄33 − C̄22)] w R`
avec `μ = sqrt(C̄22/C̄33)`.
- C1^i, C2^i : constantes d'intégration de la couche i (2n inconnues au total, plus v et w) [unités assurant u en m]. Hypothèse : solution d'équilibre radial axisymétrique par couche.

**Éq. (15)** — Continuité inter-couches (2n−2 conditions) :
`u_T^i |_(R=R_i) = u_T^(i+1) |_(R=R_i)` ; `Δσ_r^i |_(R=R_i) = Δσ_r^(i+1) |_(R=R_i)` ; `i = [1, n−1]`
- R_i : rayon de l'interface i.

**Conditions aux limites radiales (texte, non numérotées, p. 6 et 7 PDF)** :
`Δσ_r^1 |_(R=R_in) = −ΔP` et `σ_r^n |_(R=R_out) = 0` (p. 6, tube torsadé) ; réécrit `Δσ_r^n |_(R=R_out) = 0` p. 7 (TCP). ΔP : incrément de pression interne [Pa].

**Éq. (16)** — Efforts résultants du tube multicouche :
```
T = ∫₀ᵗ ∫₀^{2π} ∫_{R_in}^{R_out} σ_s  R dR dΦ dt
G = ∫₀ᵗ ∫₀^{2π} ∫_{R_in}^{R_out} σ_s R cosΦ · R dR dΦ dt
H = ∫₀ᵗ ∫₀^{2π} ∫_{R_in}^{R_out} σ_sφ R² dR dΦ dt
```
- T : tension [N] ; G : moment de flexion [N·m] ; H : couple de torsion [N·m]. **Attention notation** : l'intégrale sur dt avec σ (et non σ̇ ou Δσ) est dimensionnellement incohérente telle qu'imprimée ; l'interprétation cohérente avec le schéma incrémental est le cumul temporel des incréments de contrainte (voir §6).

### 3.2.1 — Actionnement du tube torsadé (p. 6–7 PDF)

**Hypothèses** : K = 0 et κ̄ = 0 (tube droit) ; influence du mandrin nylon négligée (tube et nylon droits) ; extrémité haute fixe, masse de la cible négligeable.

**Éq. (17)** — Conditions de charge : `T = 0` et `H = 0`.

**Éq. (18)** — Déformation d'actionnement en torsion :
`ε_Rot = (Θ/L0) · R_out`, avec `Θ = L0 ∫₀ᵗ v dt`
- Θ : angle de rotation total de l'extrémité libre [rad] ; L0 : longueur initiale du tube torsadé [m] ; R_out : rayon externe [m]. Dans les figures et le texte, cette grandeur est notée **ε_Tor** (incohérence de notation ε_Rot/ε_Tor). v est ici traité comme un taux [rad/(m·s)] — notation incrémentale relâchée.

### 3.2.2 — Actionnement du muscle TCP (p. 7 PDF)

**Configuration** : charge F appliquée au muscle ; rotation bloquée à l'extrémité basse ⇒ couple de réaction M. Indice « 0 » = configuration initiale à t_i.

**Éq. (19)** — Courbure initiale du muscle TCP :
`K = cos²α0 / ρ0`
- α0 : angle d'hélice initial [rad ou °] ; ρ0 : rayon d'hélice initial [m]. (Courbure de Frenet d'une hélice.)

**Éq. (20)** — Incréments cinématiques hélicoïdaux sur Δt :
`v = sinα·cosα/ρ − sinα0·cosα0/ρ0` ; `κ̄ = cos²α/ρ − cos²α0/ρ0`
- α, ρ : angle et rayon d'hélice courants. v = variation de la torsion géométrique (tortuosité) de l'hélice [1/m] ; κ̄ = variation de courbure [1/m].

**Éq. (21)** — Conservation du nombre de spires pendant l'actionnement :
`ρ/cosα = (1+w)·ρ0/cosα0`
- Hypothèse forte : le nombre de spires reste constant.

**Contribution du nylon (texte, non numérotée)** :
`ΔH_nylon = G_n · π R_n⁴ · v / 2` ; `ΔG_nylon = E_n · π R_n⁴ · κ̄ / 4`
- Rigidités de torsion (G·J, J = πR⁴/2) et de flexion (E·I, I = πR⁴/4) du monofilament ; R_n = 0,19 mm ; E_n = 2,8 GPa ; G_n = 0,63 GPa (mesurés). Hypothèses : dimensions du nylon inchangées ; couple/flexion du nylon déterminés directement par la géométrie hélicoïdale ; nylon négligé en torsion pure (3.2.1) mais **non négligeable** en traction TCP (module ≫ PVC). Repère de Frenet-Serret {n, b, t} utilisé pour l'hélice.

**Éq. (22)** — Équilibre du brin hélicoïdal (efforts projetés sur la section du tube) :
```
T = F sinα
G + G_nylon = M cosα − F ρ sinα
H + H_nylon = M sinα + F ρ cosα
```
- F : charge axiale appliquée [N] ; M : couple de réaction (rotation bloquée en bas) [N·m] ; T, G, H : cf. éq. (16). Fermeture du système avec (15), (20), (21) et les CL radiales (−ΔP intérieur, 0 extérieur).

**Éq. (23)** — Cinématique déplacement/angle d'hélice :
`δ = L sinα − L0 sinα0`
- δ : déplacement d'actionnement axial [m] ; L, L0 : longueurs courante et initiale du **tube torsadé** (longueur curviligne du brin) [m]. L0·sinα0 = longueur axiale initiale du ressort.

**Éq. (24)** — Déformation d'actionnement en traction :
`ε_Ten = δ / (L0 sinα0)`
- Normalisation par la longueur axiale initiale du muscle. Les figures tracent des valeurs positives pour la contraction (convention de signe implicite, voir §6).

### Relations non numérotées supplémentaires

- **Paramètre d'enroulement** (p. 10 PDF) : rayon d'hélice = `R_core + R_out` (R_core : rayon du mandrin rigide) ; `C = (R_core + R_out)/R_out`. C ↑ ⇒ course ↑ mais capacité de charge ↓.
- **FE (p. 7–8 PDF)** : constitutif élastique linéaire transversalement isotrope ; éléments C3D8R ; tube torsadé de longueur 40 mm ; muscle TCP réduit à **une seule spire** avec nylon ; axes matériaux locaux : 1 = direction hélicoïdale (fibre), 2 = direction tangentielle de la face d'extrémité, 3 = normale à la paroi interne ; angle de déflexion de fibre 40° appliqué autour de l'axe 3 ; pressions et CL imposées via points de référence.

---

## 3. Tables de paramètres

### Table 1 (p. 4 PDF = page article 3) — Matériau et géométrie

| Paramètre | Symbole | Valeur | Unité | Source |
|---|---|---|---|---|
| Module axial du PVC (étiré, λ=3) | E1 | 31,24 ± 1,9 | MPa | mesuré (IBTC-300SL) |
| Module radial du PVC | E2 = E3 | 8,82 ± 0,46 | MPa | mesuré |
| Module de cisaillement du PVC | G12 = G13 | 7,24 ± 0,32 | MPa | mesuré (machine de torsion maison [30]) |
| Poisson PVC | ν23 | 0,422 | – | réf. [29] (non mesuré) |
| Poisson PVC | ν12 = ν13 | 0,205 | – | réf. [29] (non mesuré) |
| Rayon externe du tube PVC | R_out | 0,48 ± 0,02 | mm | – |
| Rayon interne du tube PVC | R_in | 0,20 ± 0,02 | mm | – |
| Angle de biais maximal (surface) | β_f | 40,0° ± 0,01° | ° | – |
| Angle d'hélice du muscle TCP | α0 | 13,7° ± 0,01° | ° | – |
| Rayon d'hélice du muscle TCP | ρ0 | 2,72 ± 0,02 | mm | – |
| Rayon du nylon 6 | R_n | 0,19 ± 0,02 | mm | – |

### Autres valeurs dispersées dans le texte

| Grandeur | Valeur | Page PDF |
|---|---|---|
| Module élastique du nylon 6, E_n | 2,8 GPa | 7 |
| Module de cisaillement du nylon 6, G_n | 0,63 GPa | 7 |
| Élongation d'étirage à froid λ (fabrication) | 3 (rupture si > 3) | 3–4 |
| Rapport E1/E2 : λ=1 → 1 ; λ=3 → 3,64 | – | 4 |
| Recuit (après chaque étape) | 90 °C, 90 min | 3–4 |
| Cycles d'entraînement avant essai | 10 | 5 |
| Débit de pressurisation | 1 mL/min | 5 |
| Charge morte pour l'essai TCP (Fig. 10) | F = 0,98 N | 9 |
| Pression des cartographies de contraintes (Fig. 7) | P = 1,4 MPa | 8 |
| Seuil de linéarité torsion | P ≈ 0,5 MPa | 9 |
| Angle de biais critique (max de ε_Tor à 1,4 MPa) | β_f ≈ 47° | 9 (et « β = 47° » en conclusion, p. 10) |
| Longueur du tube torsadé en FE | 40 mm | 7 |
| Valeurs de C testées (Fig. 11) | 1,63 ; 2,58 ; 3,13 | 10 |
| Écart max théorie/exp (torsion) | < 6 % | 9 |
| Écart max théorie/exp (traction, vers 0,5 MPa) | < 20 % | 10 |

---

## 4. Hypothèses et domaine de validité déclarés

1. **Élasticité linéaire, petites déformations** — théorie et FE ; les auteurs reconnaissent que le PVC est viscoélastique (fluage sous charge morte et sous pression croissante) et que cela explique les écarts (p. 9–10).
2. **Isotropie transverse** de chaque couche, direction matérielle = direction de fibre inclinée de β(R) ; homogénéisation par couche cylindrique (matrice C̄ de la réf. [25]).
3. **Schéma incrémental** : pression supposée croître **linéairement** dans le temps ; en pratique le débit imposé (1 mL/min) produit une montée en pression non linéaire — source d'écart identifiée (p. 9).
4. **Tube torsadé** : K = 0, κ̄ = 0 ; nylon négligé ; masse de la cible négligée ; T = 0, H = 0.
5. **Muscle TCP** : nombre de spires constant (éq. 21) ; influence de u_B sur ε_s négligée ; dimensions du nylon inchangées ; rotation bloquée à l'extrémité basse (couple M) ; nylon pris en compte via ΔH_nylon, ΔG_nylon.
6. **Fabrication** : λ ≤ 3 (rupture au-delà) ; β_f ≤ ~40° (auto-enroulement ou rupture au-delà) ; recuits supposés éliminer toutes les contraintes résiduelles.
7. **Limites explicitement admises (conclusion, p. 10)** : effets dépendant de la vitesse non modélisés ; à haute pression, grandes déformations ⇒ non-linéarité géométrique nécessaire ; le modèle reste « qualitativement » applicable pour guider la conception.
8. Validité démontrée sur : P ∈ [0 ; ~1,4–1,6] MPa, β_f ∈ [~17° ; 40°] (expérimental), C ∈ [1,63 ; 3,13].

---

## 5. Figures clés

- **Figure 1 (p. 4 PDF)** — Chaîne de fabrication : (a) tube PVC avec nylon 6 à l'intérieur et axes matériaux 1/2/3 ; (b) étirage + recuit ; (c) torsion sous charge axiale + recuit (définit β_f) ; (d) enroulement sur mandrin + recuit (définit α0, ρ0).
- **Figure 2 (p. 4 PDF)** — Modules du PVC vs élongation d'étirage λ (0,8–3,2). Axe gauche étiqueté « Elastic modulus (MPa) », gradué **0–0,4** ; axe droit E1/E2 (0–4). Points : à λ=1, E1=E2≈0,05 (unité de l'axe), E1/E2=1 ; à λ≈2,3, E1≈0,10, E1/E2≈1,65 ; à λ=3, E1≈0,32, E2≈0,09, E1/E2≈3,64. **Incohérence d'unités avec la Table 1** (voir §6).
- **Figure 3 (p. 5 PDF)** — Bancs d'essai : (a) tube torsadé vertical, pompe d'injection en haut, cible légère en bas, rotation Θ mesurée par caméra CCD ; (b) muscle TCP avec charge morte F en bas, déplacement δ mesuré par capteur laser. Pression P à 1 mL/min.
- **Figure 4 (p. 5 PDF)** — Géométrie du modèle : (a) configuration initiale (courbure 1/K, base G_i) ; (b) configuration courante (1/κ, base g_i) ; (c) angle de biais β_f sur le tube ; (d) section multicouche, couches i = 1…n entre R_in et R_out.
- **Figure 5 (p. 7 PDF)** — Géométrie du muscle TCP : rayon ρ, efforts F et M aux deux extrémités, angle α, repère de Frenet {n, b, t} ; schéma cinématique reliant L0, L, α0, α et δ (δ = décalage axial).
- **Figure 6 (p. 8 PDF)** — Maillages FE : (a) tube torsadé (section avec coordonnées R, θ, z ; nylon jaune au centre) ; (b) une spire de TCP avec champ d'orientations matérielles locales (axes 1/2/3).
- **Figure 7 (p. 8 PDF)** — Cartographies des contraintes dans la section du tube torsadé à P = 1,4 MPa, théorie (gauche) vs FEM (droite), en MPa : (a) σ_r compressive, min ≈ −1,31 à la paroi interne → ≈ −0,012 à l'externe ; (b) σ_φ tractive, ≈ +2,13 interne → +0,72 externe ; (c) σ_s : +0,124 (traction, interne) → −0,073 (compression, externe) — résultante axiale compressive ; (d) σ_sφ : +0,539 (interne) → −0,215 (externe), décroissant radialement — moteur du dé-torsadage. Bon accord théorie/FEM.
- **Figure 8 (p. 9 PDF)** — ε_Tor (0–0,08) vs P (0–1,6 MPa), tube torsadé β_f = 40° : théorie quasi linéaire au-delà de 0,5 MPa (≈ 0,062 à 1,4 MPa), FE (≈ 0,056), expérience non linéaire (≈ 0,011 à 0,45 MPa ; ≈ 0,028 à 0,9 ; ≈ 0,059 à 1,37). Écart max < 6 % (attribué au fluage du PVC).
- **Figure 9 (p. 9 PDF)** — ε_Tor vs angle de biais β_f (0–90°) à P = 1,4 MPa : courbe théorique en cloche, **maximum ≈ 0,064 à β_f ≈ 47°**, retombant vers 0 à 0° et 90° ; points expérimentaux ≈ (17°, 0,037), (28°, 0,054), (40°, 0,059) — en dessous de la théorie près de 40° (dégradation matériau par sur-torsion).
- **Figure 10 (p. 9 PDF)** — ε_Ten (0–0,30) vs P (0–1,6 MPa), muscle TCP, F = 0,98 N : croissance non linéaire ; expérience ≈ 0,05 à 0,6 MPa, ≈ 0,13 à 0,9, ≈ 0,27 à 1,3 ; théorie légèrement au-dessus à basse pression, FE proche ; écart max < 20 % vers 0,5 MPa.
- **Figure 11 (p. 10 PDF)** — ε_Ten (0–0,6) vs P (0–1,6 MPa) pour C = 1,63, 2,58, 3,13 (théorie + expérience) : à ≈ 1,3 MPa, ε_Ten ≈ 0,20 (C=1,63), ≈ 0,27 (C=2,58), ≈ 0,48 (C=3,13). Course ↑ avec C ; capacité de charge ↓ avec C.

---

## 6. Ambiguïtés et points laissés implicites (à l'attention des auditeurs)

1. **Unités de la Figure 2 incohérentes avec la Table 1** : l'axe indique « MPa » gradué 0–0,4, avec E1 ≈ 0,32 à λ=3, alors que la Table 1 donne E1 = 31,24 MPa. Facteur ≈ 100 entre figure et table (l'axe devrait probablement être ×10² MPa). Le rapport E1/E2 (axe droit) est, lui, cohérent (3,64 ≈ 31,24/8,82 = 3,54). Retenir les valeurs de la **Table 1** comme référence.
2. **Nombre de couches n jamais spécifié** : la discrétisation multicouche (éq. 13–15) est un choix numérique laissé libre ; le profil β(R) par couche (vraisemblablement tanβ(R) = R·tanβ_f/R_out, classique) n'est pas donné.
3. **Notation incrémentale relâchée** : éq. (16) intègre σ_s, σ_sφ sur dt sans point de dérivée (dimensionnellement faux tel qu'imprimé) ; de même Θ = L0∫v dt (éq. 18) avec v défini comme un incrément. Interprétation cohérente : cumul des incréments Δσ, Δv sur les pas de temps.
4. **Double notation ε_Rot / ε_Tor** pour la même grandeur (éq. 18 vs figures 8–9 et texte).
5. **Signe du tenseur de déformation (éq. 9)** : ε = (c⁻¹ − I)/2, opposé de la convention Euler-Almansi usuelle e = (I − c⁻¹)/2 ; vérifier quelle convention le code implémente réellement et la cohérence des signes en aval (u_B jamais explicité non plus).
6. **Conditions aux limites radiales** écrites une fois `σ_r^n|R_out = 0` (p. 6) et une fois `Δσ_r^n|R_out = 0` (p. 7) — la version incrémentale (Δ) est la cohérente.
7. **Convention de signe de δ et ε_Ten** : δ = L sinα − L0 sinα0 devrait être négatif pour une contraction, mais les figures 10–11 tracent des valeurs positives ; l'article ne précise jamais si ε_Ten est une amplitude ou si un signe est absorbé quelque part.
8. **Incohérence apparente ρ0 / C** : la Table 1 donne ρ0 = 2,72 mm, soit C = ρ0/R_out ≈ 5,7, hors de la plage C = 1,63–3,13 de la Figure 11 ; l'article ne précise pas quel C correspond au muscle des Figures 3(b)/10 ni si α0 = 13,7° vaut pour toutes les valeurs de C. Les muscles de Fig. 10 et Fig. 11 semblent être des échantillons différents sans que ce soit dit.
9. **L0 (longueur du tube torsadé du muscle TCP) et le nombre de spires jamais donnés** ; la FE ne modélise qu'une spire. La masse morte (0,98 N ≈ 100 g) n'est donnée que pour la Fig. 10.
10. **Relation torsion insérée ↔ β_f non explicitée** (combien de tours/m pour obtenir β_f = 40°) — renvoyée à la Supporting Information, absente du PDF.
11. **G23 et E3 ne sont pas utilisés/donnés séparément** : E2 = E3, G12 = G13 sont donnés ; G23 (=E2/2(1+ν23) si isotropie transverse) est implicite. Les ν sont repris de la réf. [29] (cavatappi), pas mesurés sur ce PVC.
12. **Éq. (14)** : les coefficients (C̄26−2C̄36)/(4C̄33−C̄22) et (C̄12−C̄13)/(C̄33−C̄22) proviennent de Pipes & Hubert [31] ; toute vérification du code doit contrôler ces dénominateurs exacts (risque élevé de coquille de recopie ; l'article BLOCKED, qui reprend ce cadre, doit être comparé terme à terme).
13. **La métrique (3) est celle d'un tore** (courbure K seule) : la torsion géométrique de l'hélice n'apparaît pas dans la métrique de référence, elle n'entre que par la cinématique incrémentale v (éq. 20). C'est un choix de modélisation non discuté.
14. **β = 47° vs β_f ≈ 47°** : la conclusion (p. 10) écrit « β = 47° » sans l'indice f ; il s'agit bien de l'angle de biais de surface.
15. **Le modèle est purement élastique** : aucun paramètre viscoélastique (temps de relaxation, modules de Maxwell) dans cet article ; le fluage n'est invoqué que qualitativement pour expliquer les écarts. Toute fonctionnalité viscoélastique du code audité doit être tracée vers l'article BLOCKED, pas vers MECA2025.

---

**Fichiers de travail produits (scratchpad)** : texte brut extrait `meca2025_text.txt`, rendus pleine page `meca_p01.png`–`meca_p11.png`, zooms `fig2_zoom.png`, `eq16_zoom.png`, `eq18_zoom.png` dans `C:/Users/B3DCB~1.PER/AppData/Local/Temp/claude/C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail/c5852e10-6d38-4703-b218-a1649a5b12c9/scratchpad`.