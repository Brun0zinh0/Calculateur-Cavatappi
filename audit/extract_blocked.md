# Rapport de référence — Article « BLOCKED »
**Hu J., Liu L., Liu H., Teng J., Liu D., « Blocked actuation of twisted and coiled tube-based polymer actuators driven by hydraulic pressure », International Journal of Smart and Nano Materials, 17 sept. 2025, DOI 10.1080/19475411.2025.2558660.**

Convention de pagination : le PDF comporte 20 pages ; la page 1 du PDF est la page de garde Taylor & Francis. La page imprimée *n* de l'article correspond à la page *n+1* du PDF. Les pages citées ci-dessous sont les **pages imprimées de l'article** (p. 1–19).

---

## 1. Résumé du modèle et de la démarche

L'article étudie l'actionnement **bloqué** (deux extrémités encastrées, rotation et translation empêchées) de muscles artificiels TCPA fabriqués à partir d'un tube PVC (TYGON ND 100-65) étiré à froid (rapport 3:1), torsadé (angle de biais visé 37°), enroulé sur mandrin, renforcé par un monofilament nylon 6 interne servant de mandrin anti-collapse, puis recuit à 90 °C pendant 90 min. Sous pression hydraulique interne (eau injectée par pousse-seringue), le tube torsadé tend à se détordre ; comme les extrémités sont bloquées, il génère un couple de réaction T et une force axiale F mesurés simultanément sur un banc tension-torsion micrométrique maison, avec capteur de pression et carte NI PCI-6280. Le protocole : pré-étirement à une pré-déformation ε_tk donnée (0,5 à 1,2 expérimentalement), puis 100 cycles injection/retrait à 10 mL/min, pression crête limitée à 1,5 MPa (1,2–1,3 MPa en régime long), les 10 premiers cycles servant de « training » pour stabiliser l'hystérésis (acquisition exploitée à partir du 11e cycle). Le modèle théorique est incrémental, en **description lagrangienne réactualisée** (la configuration à t sert de référence pour t+Δt), car les déformations atteignent ~20 %. Le PVC est décrit par un **modèle de Maxwell généralisé** : un ressort C₀ en parallèle avec trois branches de Maxwell (Cᵢ, ηᵢ), i = 1..3, écrit sous forme matricielle anisotrope (isotropie transverse) ; le nylon 6 est élastique linéaire. Le tube torsadé est discrétisé en *n* couches cylindriques concentriques ; l'angle de biais varie linéairement avec le rayon (θⱼ = arctan(Rⱼ/R_out·tan θ_f)) ; chaque couche suit la solution de déplacement radial de Pipes-Hubert avec continuité de Δu et Δσ_r aux interfaces, et conditions aux limites en pression sur les faces interne/externe. La cinématique de l'hélice est fermée par trois relations : conservation du nombre de tours N (extrémités encastrées), relation ρ/cos α liée à l'allongement du fil, et condition de compatibilité bloquée ρ·tan α = constante pendant l'actionnement. Les efforts internes du tube (F_tube, M_tube, T_tube) s'obtiennent par intégration des incréments de contrainte sur la section, ceux du nylon par des formules de poutre élastique ; l'équilibre de l'hélice (théorie des ressorts) relie ces efforts internes aux sorties mesurables F_t et T_t. Le modèle reproduit bien les figures 7–9 (validation quantitative après 10 cycles de training) mais **sous-estime l'aire des boucles d'hystérésis**, ce que les auteurs attribuent à l'hypothèse (ii) d'anisotropie identique des éléments élastiques et visqueux. Les études paramétriques théoriques montrent : couple normalisé maximal à un indice de ressort intermédiaire (pic en C), force normalisée croissante quand C diminue et quand ε_tk augmente, et réduction de l'hystérésis quand la vitesse de pressurisation augmente (0,1 → 1 MPa/s). Le training cyclique réduit la pression de seuil (« dead-band ») de 0,38 à 0,23 MPa et améliore la répétabilité (<3 % sur 100 cycles).

---

## 2. Inventaire exhaustif des équations

### Corps de l'article

**Éq. (1) — p. 6** — Lois incrémentales du modèle de Maxwell généralisé :
```
Δσ₀ = C₀ Δε
σᵢ Δt + ηᵢ Cᵢ⁻¹ Δσᵢ = ηᵢ Δε ,   i = [1,3]
```
- σ₀ [MPa] : matrice (vecteur de Voigt) des contraintes du ressort parallèle ; C₀ [MPa] : sa matrice de rigidité.
- σᵢ, Cᵢ, ηᵢ : contrainte [MPa], rigidité [MPa] et viscosité [MPa·s] (matricielles) de la i-ème branche de Maxwell.
- Δε [–] : incrément de déformation (commun à toutes les branches, montage en parallèle) ; Δt [s] : pas de temps.
- Hypothèses : viscoélasticité linéaire par incrément ; branches élastiques = microfibrilles cristallines, dashpots = matrice amorphe ; même anisotropie pour toutes les matrices (hyp. ii).
- Conséquence physique invoquée p. 10 : pendant la pressurisation σ et σ̇ ont le même signe, signes opposés en dépressurisation → asymétrie charge/décharge (hystérésis).

**Éq. (2) — p. 6** — Incrément de contrainte totale :
```
Δσ = Δσ₀ + Σᵢ₌₁³ Δσᵢ = C_E Δε + C_S
```
avec, **relation non numérotée immédiatement sous (2)** :
```
C_E = C₀ + Σᵢ₌₁³ Cᵢ          (rigidité effective, MPa)
C_S = − Δt Σᵢ₌₁³ Cᵢ σᵢ ηᵢ⁻¹   (terme « visqueux », homogène à une contrainte, MPa)
```
- Δσ se décompose en partie élastique C_E Δε et partie visqueuse C_S (qui dépend de l'état de contrainte courant σᵢ de chaque branche — mémoire du matériau).

**Éq. (3) — p. 7** — Équilibre de l'hélice (extrémités bloquées) :
```
F_tube + F_nylon = F_t sin α_t
M_tube + M_nylon = T_t cos α_t − F_t ρ_t sin α_t
T_tube + T_nylon = T_t sin α_t + F_t ρ_t cos α_t
```
- F_t [N] : force axiale au bout du TCPA ; T_t [N·m] : couple au bout du TCPA (les sorties mesurées/prédites).
- F_tube, F_nylon [N] : efforts normaux le long de l'axe centroïdal du tube/du filament ; M [N·m] : moment de fléchissement ; T_tube, T_nylon [N·m] : couples de torsion autour de l'axe du fil.
- α_t [rad] : angle d'hélice à l'instant t ; ρ_t [m] : rayon d'hélice à t.
- Hypothèse : théorie des ressorts hélicoïdaux ; le tube et le nylon travaillent en parallèle (mêmes déformations généralisées de la fibre centrale).

**Éq. (4) — p. 8** — Incréments d'effort interne du filament nylon 6 (élastique linéaire, section circulaire invariante — hyp. i et iii) :
```
ΔF_nylon = π E_nylon r_nylon² Δw
ΔM_nylon = (1/4) π E_nylon r_nylon⁴ Δκ
ΔT_nylon = (1/2) π G_nylon r_nylon⁴ Δv
```
- E_nylon [Pa] : module d'Young du nylon 6 ; G_nylon [Pa] : module de cisaillement ; r_nylon [m] : rayon du filament.
- Δw [–] : incrément de déformation axiale (le long du fil) ; Δv [rad/m] : incrément d'angle tangentiel par unité de longueur (taux de torsion) ; Δκ [1/m] : incrément de courbure, pendant Δt.

**Éq. (5) — p. 8** — Lien cinématique entre (Δv, Δκ) et la géométrie de l'hélice :
```
Δv = sin(2α_{t+Δt}) / (2ρ_{t+Δt}) − sin(2α_t) / (2ρ_t)
Δκ = cos²α_{t+Δt} / ρ_{t+Δt} − cos²α_t / ρ_t
```
- Formules classiques de la torsion géométrique τ = sin(2α)/(2ρ) et de la courbure κ = cos²α/ρ d'une hélice.

**Éq. (6) — p. 8** — Conservation du nombre de tours (extrémités encastrées ; réf. [21]) :
```
ρ_{t+Δt} / cos α_{t+Δt} = (1 + Δw) ρ_t / cos α_t
```
- Le rapport ρ/cos α est proportionnel à la longueur de fil par tour ; il ne varie que par l'étirement du fil Δw.

**Éq. (7) — p. 8** — Déplacement de l'extrémité pendant la phase de pré-étirement (t ∈ [0, t_k], sans pression) :
```
δ_t = 2πN ( ρ_t tan α_t − ρ_{t₀} tan α_{t₀} )
```
- N [–] : nombre de tours d'hélice ; δ_t [m] : déplacement d'extrémité ; ρ tan α = pas/(2π) de l'hélice.
- Indices : t₀ = état initial auto-contact ; t_k = fin d'étirement/début de pressurisation.

**Éq. (8) — p. 8** — Définition de la pré-déformation :
```
ε_tk = ( ρ_tk tan α_tk − ρ_{t₀} tan α_{t₀} ) / ( ρ_{t₀} tan α_{t₀} )
```
- ε_tk [–] : pré-déformation appliquée avant actionnement (0,5–1,2 en essai ; 0,5–2,0 en simulation).

**Éq. (9) — p. 9** — Condition de compatibilité bloquée pendant l'actionnement (t > t_k) :
```
ρ_{tk+Δt} tan α_{tk+Δt} = ρ_{tk} tan α_{tk}
```
- La longueur axiale du TCPA est fixe (encastrement) : le pas de l'hélice est constant. Avec (6), ceci détermine α et ρ à chaque pas connaissant Δw.

**Éq. (10) — p. 9** — Efforts internes du tube par intégration des contraintes sur la section (réf. [20]) :
```
F_tube = ∫₀ᵗ ∫₀²π ∫_{R_in}^{R_out} σ_s R dR dΦ dt
M_tube = ∫₀ᵗ ∫₀²π ∫_{R_in}^{R_out} σ_s R² cos Φ dR dΦ dt
T_tube = ∫₀ᵗ ∫₀²π ∫_{R_in}^{R_out} σ_φs R² dR dΦ dt
```
- σ_s [MPa] : contrainte axiale (direction s du fil) ; σ_φs [MPa] : cisaillement circonférentiel-axial ; R [m], Φ [rad], S : coordonnées curvilignes locales de la section (cf. Fig. 6 et Éq. (16)).
- L'intégrale sur t exprime le cumul des incréments de contrainte issus de (A.4) — notation à interpréter comme une somme des contributions incrémentales (voir §6, ambiguïtés).

**Éq. (11) — p. 9** — Angle de biais de la couche j (hyp. iv, variation linéaire radiale) :
```
θⱼ = arctan( (Rⱼ / R_out) tan θ_f )
```
- θ_f [°] : angle de biais maximal (surface externe) ; Rⱼ [m] : rayon de la couche j ; θⱼ : angle des microfibres de la couche j.

**Éq. (12) — p. 9** — Déplacement radial d'une couche (solution de Pipes & Hubert [34]) :
```
Δu = C₁ R^μ + C₂ R^(−μ)
     + [ (C̄_E26 − 2 C̄_E36) / (4 C̄_E33 − C̄_E22) ] Δv R²
     + [ (C̄_E12 − C̄_E13) / (C̄_E33 − C̄_E22) ] Δw R
```
avec, **relation non numérotée** : `μ = sqrt( C̄_E22 / C̄_E33 )`
- Δu [m] : déplacement radial incrémental de la couche j ; C₁, C₂ : constantes d'intégration propres à chaque couche (déterminées par (13)–(15)) ; C̄_Eij : composantes de la matrice de rigidité effective transformée (Éq. (A.3)).

**Éq. (13) — p. 9** — Continuité aux interfaces des n couches :
```
Δuʲ|_{R=Rⱼ} = Δuʲ⁺¹|_{R=Rⱼ} ,  Δσ_rʲ|_{R=Rⱼ} = Δσ_rʲ⁺¹|_{R=Rⱼ} ,  j = [1, n−1]
```

**Éq. (14) — p. 9** — Conditions aux limites radiales, **phase d'élongation** (sans pression) :
```
Δσ_r|_{R=R_out} = Δσ_r|_{R=R_in} = 0
```

**Éq. (15) — p. 9** — Conditions aux limites radiales, **phase d'actionnement** :
```
Δσ_r|_{R=R_out} = 0 ,   Δσ_r|_{R=R_in} = −ΔP
```
- ΔP [MPa] : incrément de pression hydraulique interne. Hyp. v : aucune pression externe.

**Éq. (16) — p. 10** — Coordonnées curvilignes des configurations de référence et courante :
```
X¹ = R, X² = Φ, X³ = S ;   x¹ = r, x² = φ, x³ = s
```

**Éq. (17) — p. 10** — Passage référence → courante :
```
r = R + Δu ,  φ = Φ + Δv·S ,  s = (1 + Δw) S ,  κ = K + Δκ
```
avec, **relation non numérotée** : `K = cos²α_{t₀} / ρ_{t₀}` (courbure initiale du tube torsadé, [1/m]).

**Éq. (18) — p. 10** — Incréments de déformation avec courbure initiale (réf. [35]) :
```
ε_φr = ε_rs = 0
ε_r  = ∂u/∂R − ν̄₁₂ (κR cos Φ + uK cos Φ) / (1 + KR cos Φ)
ε_φ  = u/R  − ν̄₁₃ (κR cos Φ + uK cos Φ) / (1 + KR cos Φ)
ε_s  = ω + (κR cos Φ + uK cos Φ) / (1 + KR cos Φ)
ε_φs = νR / (1 + KR cos Φ) − ν̄₁₄ (κR cos Φ + uK cos Φ) / (1 + KR cos Φ)
```
- u : déplacement radial ; ω (= w) : déformation axiale ; ν (= v) : torsion par unité de longueur ; κ : courbure ; K : courbure initiale. Le terme (1 + KR cos Φ) est le facteur métrique dû à la courbure initiale du fil (position Φ dans la section).
- ν̄₁₂, ν̄₁₃, ν̄₁₄ [–] : coefficients de Poisson effectifs (Éq. (19)).
- Nota : notation « totale » (sans Δ) alors que le cadre est incrémental — voir §6.

**Éq. (19) — p. 10** — Détermination des coefficients de Poisson effectifs (réf. [21]) :
```
C̄_E⁻¹ [ Ē_L  0  0  0 ]ᵀ = [ 1  −ν̄₁₂  −ν̄₁₃  −ν̄₁₄ ]ᵀ
```
- Ē_L [MPa] : module d'Young effectif du tube torsadé ; C̄_E : matrice de rigidité effective [17] (Éq. (A.3)) réduite ici à 4 composantes (r, φ, s, φs) — dimension implicite, voir §6.

**Éq. (20) — p. 13** — Couple d'actionnement normalisé :
```
T* = T cos α_tk / ( 2π E_nylon (R_out)⁴ )
```
- T* a la dimension d'une longueur⁻¹ ; l'échelle de couleur de la Fig. 10(a) le donne en **rad/mm** (0,22–0,52).

**Relation non numérotée — p. 13** — Indice de ressort :
```
C = ρ / R_out
```

**Éq. (21) — p. 13** — Force d'actionnement normalisée (réf. [35]) :
```
F* = (1 + ε_tk) F sin α_tk / ( π R_out² )
```
- F* [MPa] (échelle Fig. 10(b) : 0,15–1,25 MPa).

**Éq. (22) — p. 13** — Forme équivalente via l'équilibre (3) :
```
F* = (1 + ε_tk) ( F_tube + F_nylon ) / ( π R_out² )
```
- F_tube s'obtient en intégrant σ_s (qui décroît quand C croît [35]) ; F_nylon dépend de w, qui croît avec ε_tk — mécanisme expliquant les tendances de la Fig. 10(b).

### Annexe A (p. 17–18) — Relation constitutive viscoélastique incrémentale du PVC

**Éq. (A.1) — p. 17** — Matrice de rigidité d'une couche dans les axes matériaux (isotropie transverse, réf. [17]) :
```
      | C₁₁ C₁₂ C₁₂    0        0    0  |
      | C₁₂ C₂₂ C₂₃    0        0    0  |
C_E = | C₁₂ C₂₃ C₂₂    0        0    0  |
      |  0   0   0  (C₂₂−C₂₃)/2  0    0  |
      |  0   0   0     0       C₅₅   0  |
      |  0   0   0     0        0   C₆₆ |
```
- Axe 1 = direction axiale du matériau (microfibres), axes 2, 3 = plan transverse isotrope. **C₅₅ n'est défini nulle part** (voir §6).

**Éq. (A.2) — p. 17** — Système donnant C₁₁, C₁₂, C₂₃, C₂₂, C₆₆ :
```
C₁₁ = E_axial + 2 ν₁₂ C₁₂
C₁₂ = ν₁₂ (C₂₂ + C₂₃)
C₁₂ = C₁₁ ν₂₁ + C₁₂ ν₂₃
C₂₂ = C₁₂ ν₂₁ + C₂₃ ν₂₃ + E_radius
ν₂₁ = ν₁₂ E_radius / E_axial
ν₃₁ = ν₁₂ E_radius / E_axial
ν₃₂ = ν₂₃ E_radius / E_axial
C₆₆ = G₁₂
```
- E_axial, E_radius [MPa] : modules axial et radial du PVC étiré ; G₁₂ [MPa] : module de cisaillement ; ν₁₂, ν₂₃ [–] : coefficients de Poisson (Table B1). Les 3e et 4e lignes, combinées aux deux premières, ferment le système linéaire en (C₁₁, C₁₂, C₂₂, C₂₃) ; ν₃₁, ν₃₂ sont donnés mais inutilisés ailleurs.

**Éq. (A.3) — p. 17** — Rotation de C_E vers les coordonnées « colonne » (cylindriques), pour la couche j :
```
C̄₁₁ = C₁₁ m⁴ + 2(C₁₂ + C₁₆) m²n² + C₂₂ n⁴
C̄₁₂ = (C₁₁ + C₂₂ − 4C₆₆) m²n² + C₁₂ (m⁴ + n⁴)
C̄₁₃ = C₁₃ m² + C₂₃ n²
C̄₁₆ = −C₂₂ m n³ + C₁₁ m³ n − (C₁₂ + 2C₆₆) m n (m² − n²)
C̄₂₂ = C₁₁ n⁴ + 2(C₁₂ + 2C₄₄) m²n² + C₂₂ m⁴
C̄₂₃ = C₁₃ n² + C₂₃ m²
C̄₂₆ = −C₂₂ m³ n + C₁₁ m n³ + (C₁₂ + 2C₆₆) m n (m² − n²)
C̄₃₃ = C₃₃
C̄₃₆ = (C₁₃ − C₂₃) m n
C̄₆₆ = (C₁₁ + C₂₂ − 2C₁₂) m²n² + C₆₆ (m² − n²)²
```
avec `m = cos θⱼ , n = sin θⱼ`.
- Ici C₁₃, C₃₃, C₄₄, C₁₆ désignent des **entrées de la matrice (A.1)** : C₁₃ = C₁₂, C₃₃ = C₂₂, C₄₄ = (C₂₂−C₂₃)/2, et C₁₆ = 0 dans (A.1) — la présence de C₁₆ dans C̄₁₁ est très probablement une coquille pour 2C₆₆ (voir §6). Rotation autour de l'axe radial, d'angle θⱼ (angle de biais de la couche).

**Éq. (A.4) — p. 17** — Relation constitutive incrémentale finale (combinaison de (2) et (A.3)) :
```
Δσ = C̄_E Δε + C̄_S
```

**Éq. (A.5) — p. 17** :
```
C̄_E = C̄₀ + Σᵢ₌₁³ C̄ᵢ ,     C̄_S = −Δt Σᵢ₌₁³ C̄ᵢ σᵢ / η̄ᵢ
```
avec, **relation non numérotée cruciale (répartition des branches de Maxwell)** :
```
C̄_k = [ E_k / (E₀+E₁+E₂+E₃) ] · C̄_E ,   k = [0,3]
η̄ᵢ  = [ ηᵢ  / (E₀+E₁+E₂+E₃) ] · C̄_E ,   i = [1,3]
```
- Chaque branche k porte une fraction E_k/ΣE de la rigidité anisotrope totale C̄_E ; la « viscosité matricielle » η̄ᵢ est la même matrice pondérée par ηᵢ/ΣE (d'où l'hypothèse ii : anisotropie identique pour élasticité et viscosité). Attention à la circularité apparente de l'écriture (C̄_k défini via C̄_E, lui-même somme des C̄_k) — voir §6.
- **Relation non numérotée (p. 17–18)** : `E_axial = E₀ + E₁ + E₂ + E₃` (module axial du PVC pour le modèle de Maxwell en parallèle).

---

## 3. Tables de paramètres

### Table A1 (p. 17) — Paramètres du modèle de Maxwell généralisé (traction, PVC)
| Paramètre | Symbole | Valeur |
|---|---|---|
| Module du ressort parallèle | E₀ | 6,36 MPa |
| Module de l'élément 1 | E₁ | 20,67 MPa |
| Viscosité de l'élément 1 | η₁ | 154,57 MPa·s |
| Module de l'élément 2 | E₂ | 5,98 MPa |
| Viscosité de l'élément 2 | η₂ | 977,79 MPa·s |
| Module de l'élément 3 | E₃ | 4,75 MPa |
| Viscosité de l'élément 3 | η₃ | 11 044,83 MPa·s |

Identification (p. 17–18) : traction de 30 → 31,2 mm (ε = 4 %) à trois vitesses de déformation (2×10⁻⁴, 1×10⁻³, 5×10⁻³ s⁻¹), puis relaxation ~60 min ; ajustement par moindres carrés (machine CARE IBTC-300SL). Temps caractéristiques induits : τ₁ = η₁/E₁ ≈ 7,5 s ; τ₂ ≈ 163,5 s ; τ₃ ≈ 2325 s. **Nota : E₀+E₁+E₂+E₃ = 37,76 MPa ≠ E_axial = 31,24 MPa de la Table B1 (voir §6).**

### Table B1 (p. 18) — Paramètres matériau et géométrie des TCPA
| Paramètre | Symbole | Valeur |
|---|---|---|
| Module axial du PVC | E_axial | 31,24 ± 1,90 MPa |
| Module radial du PVC | E_radius | 8,82 ± 0,46 MPa |
| Module de cisaillement du PVC | G₁₂ | 7,24 ± 0,32 MPa |
| Poisson PVC | ν₁₂ | 0,205 (repris de [15]) |
| Poisson PVC | ν₂₃ | 0,422 (repris de [15]) |
| Module axial nylon 6 | E_nylon | 3,69 ± 0,25 GPa |
| Module de cisaillement nylon 6 | G_nylon | 0,79 ± 0,08 GPa |
| Diamètre du nylon 6 | 2r_nylon | 0,77 mm |
| Diamètre extérieur du tube PVC (après étirage) | 2R_out0 | 2,0 ± 0,02 mm |
| Diamètre intérieur du tube PVC | 2R_in0 | 0,8 ± 0,02 mm |
| Angle de biais initial des microfibres | θ_f0 | 37,91° ± 0,01° |
| Rayon d'hélice initial | ρ₀ | 2,16 ± 0,02 mm |
| Angle d'hélice initial | α₀ | 10,53° ± 0,01° |

Mesures : microscope PZ-CS3500H (géométrie) ; testeur tension-torsion micrométrique maison [32] + IBTC-300SL (modules). PVC et nylon modélisés transversalement isotropes.

### Fabrication (p. 3–4)
- Tube PVC TYGON ND 100-65 : longueur initiale ~100 mm, Ø ext. 3,2 mm, Ø int. 1,6 mm.
- Monofilament nylon 6 (Klinskaya New FISHLINE) : longueur ~300 mm, Ø 0,75 mm (0,77 mm en Table B1), lubrifié à l'huile alimentaire.
- Étirage à froid à 300 mm, rapport 3:1 ; angle de biais cible 37° ; enroulement homochiral sur mandrin (auto-contact) ; 3 torsions ajoutées toutes les 3 spires pour maintenir l'angle de biais ; recuit 90 °C / 90 min [20].

### Protocole expérimental (p. 4–6)
- Longueur initiale du TCPA testé : 32,45 ± 0,02 mm ; pré-étirement à ε = 0,5 à 20 mm/min ; injection cyclique 1,22 mL/cycle à 10 mL/min, 100 cycles ; pression crête limitée 1,5 MPa (plasticité), 1,2–1,3 MPa en usage prolongé.
- Matériel : pousse-seringue Leapet TYD02-01 ; capteur de pression AE-T-G (Nanjing Aire) ; DAQ NI PCI-6280.
- Dégradations mesurées sur les 10 premiers cycles : force crête 1041,74 → 934,37 mN ; couple crête 1208,95 → 1187,40 µNm ; pression crête 1,20 → 1,11 MPa ; après training : P ≈ 1,16 MPa, T ≈ 1185,07 µNm.
- Dead-band : 0,38 MPa (1er cycle) → 0,23 MPa (11e cycle) (« réduction 65,8 % » selon le texte — voir §6).

### Conditions des figures de validation et simulations
| Cas | Mandrin | Débit | Volume/cycle | Pré-déformation | Page |
|---|---|---|---|---|---|
| Fig. 7 | 2,0 mm | 10 mL/min | 1,50 mL | 0,8 et 1,0 | p. 10–11 |
| Fig. 8 | 3,0 mm | 10 mL/min | 1,76 mL | 0,5 et 1,2 | p. 11 |
| Fig. 9 | 3,0 mm | — | — | 0,5 / 1,0 / 1,2 (11e cycle) | p. 12 |
| Fig. 10 (simulation) | — | 10 mL/min | 1,5 mL, P crête 1,4 MPa | 0,5–2,0 ; C = 2–3 ; « R_out = 1 mm » ; « α_tk = 37° » | p. 13 |
| Fig. 11 (simulation) | C = 2,25 | Ṗ = 0,1 / 0,5 / 1 MPa/s | — | 2,0 ; relaxation 60 s après étirement | p. 14 |

### Table C1 (p. 19) — Benchmark des performances
Notre travail (PVC) : puissance spécifique moyenne 0,8 W/g ; travail spécifique 0,11–0,38 J/g ; contrainte d'actionnement max 0,70 MPa ; déformation max 50 % [20,21] ; efficacité contractile 45 % (annoncée 22,5× le record antérieur [37]) ; fréquence > 1 Hz (* limitée par le régulateur hydraulique). Comparaisons : Haines [1] nylon 6,6 (27 W/g ; 2,48 J/g ; 100 MPa ; 46 % ; 1 % ; 1 Hz) ; Chen [9] nylon 6 (1,1 ; 2,10 ; 0,30 ; 50 % ; 1 % ; 0,02 Hz) ; Park [38] PET (~0,21 ; 0,39–0,85 ; – ; 12,1 %) ; Hiraoka [37] PE (1,6 ; 1,7 ; – ; ~23 % ; 1,2–2,0 %) ; Dong [39] CNT (18 % de déformation).

---

## 4. Hypothèses et domaine de validité déclarés

**Les cinq approximations du modèle (p. 6, §3) :**
1. (i) Le tube PVC est viscoélastique ; le nylon est élastique linéaire.
2. (ii) Tous les éléments du modèle de Maxwell ont **le même degré d'anisotropie** (cause déclarée de la sous-estimation des aires d'hystérésis, p. 12–13).
3. (iii) La géométrie du filament nylon 6 est invariante pendant la déformation.
4. (iv) L'angle de biais varie linéairement dans la direction radiale.
5. (v) Aucune pression externe appliquée sur le tube pendant l'actionnement.

**Autres hypothèses/conditions structurantes :**
- Description lagrangienne réactualisée obligatoire : déformations globales jusqu'à ~20 % (la théorie des petites déformations est déclarée insuffisante, p. 7).
- Extrémités encastrées → N constant (Éq. 6) et longueur bloquée (Éq. 9) ; la configuration à t_k sert de référence pour l'actionnement.
- Domaine exploré : pression ≤ 1,5 MPa (au-delà : déformation plastique) ; pré-déformation 0,5–2,0 (« pour prévenir la déformation plastique ») ; indice de ressort 2–3.
- Auto-contact des spires avant essai ; intervalle pré-étirement → pompage minimisé (éviter la relaxation) ; 10 cycles de préconditionnement avant exploitation des données.
- Limites reconnues (p. 14–15, conclusions) : effets de température non modélisés (proposition future : TTSP + cycles thermiques) ; anisotropies élastique/visqueuse non découplées ; dépendance à une pompe externe.

---

## 5. Figures clés

- **Fig. 1 (p. 3)** : étapes de fabrication — (a) tube PVC initial, (b) étirage, (c) torsion, (d) enroulement + recuit.
- **Fig. 2 (p. 4)** : banc d'essai cyclique — (a) TCPA fixé aux deux extrémités (état initial), (b) actionnement cyclique après pré-étirement ; pousse-seringue, capteur de pression, mesure synchrone T/F/P.
- **Fig. 3 (p. 5)** : réponse type sur 100 cycles à ε = 0,5. (a) trois panneaux vs temps 0–1500 s : Couple (µNm, 0–5000, encart 0–1500 sur 150–400 s), Force (mN, 600–1800, encart 600–1000 ; dérive descendante marquée de la force par fluage), Pression (MPa, 0–3,0, encart 0–1,5). (b) boucles d'hystérésis 1er vs 11e cycle : couple vs pression (≈ −150–1350 µNm ; 0–1,2 MPa) et force vs pression (600–1050 mN) ; lignes pointillées bleues = dead-band 0,38 MPa (1er) et 0,23 MPa (11e) ; recouvrement non monotone en décharge (remontée transitoire puis décroissance).
- **Fig. 4 (p. 5)** : schéma du modèle de Maxwell généralisé — ressort C₀ en parallèle avec 3 branches (ressort Cᵢ en série avec amortisseur ηᵢ), contrainte σ appliquée, déformation commune ε.
- **Fig. 5 (p. 7)** : géométrie du TCPA — (a) état initial t₀ (ρ_t₀, α_t₀), (b) état étiré t_k (ρ_tk, α_tk, δ_tk, F_tk, T_tk), (c) état actionné t sous pression P (ρ_t, α_t, F_t, T_t), (d) triangle de l'hélice déroulée (2πρ_t₀, 2πρ_tk, δ_tk).
- **Fig. 6 (p. 8)** : tube torsadé à courbure initiale — (a) configuration initiale (S, 1/K, vecteurs e₁..e₃, G₁..G₃, angle Φ), (b) configuration courante (1/κ, g₁..g₃), (c) angle de biais θ_f, R_out, R_in, (d) section multicouche (Rⱼ, j = 1..n, Δu, ΔvS).
- **Fig. 7 (p. 11) — CIBLE DE VALIDATION** : comparaison théorie/expérience, actionnement bloqué, mandrin 2,0 mm, 10 mL/min, 1,50 mL/cycle. Chaque sous-figure = 3 panneaux empilés vs temps **180–380 s** (≈ 9 cycles, période ≈ 22 s) : Force (mN), Couple (µNm), Pression (MPa, 0–1,5, crêtes ≈ 1,4, vallées ≈ 0). (a) ε_tk = 0,8 : force sur axe 1200–2200 mN, pics ≈ 1800 mN décroissant légèrement, vallées expérimentales ≈ 1300–1400 mN (la théorie surestime un peu les vallées ≈ 1450 mN) ; couple sur axe 0–1650 µNm (graduations 550/1100), pics ≈ 800–900 µNm, vallées ≈ 0. (b) ε_tk = 1,0 : force 2000–2800 mN, pics ≈ 2400–2460 mN ; couple axe gradué 480/960/1440 µNm, pics ≈ 870 µNm. **Ancrages chiffrés du texte (11e cycle)** : P crête 1,41–1,43 MPa ; F crête 1824,60 mN (0,8) → 2424,67 mN (1,0), +32,8 % ; T crête 877,01 → 873,41 µNm, −0,4 %.
- **Fig. 8 (p. 11)** : idem, mandrin 3,0 mm, 1,76 mL/cycle. (a) ε_tk = 0,5 : force 600–1200 mN, pics ≈ 950–1000 mN ; couple 0–960 µNm (320/640), pics ≈ 700 µNm ; temps 180–380 s. (b) ε_tk = 1,2 : temps **180–500 s** (~19 cycles) ; force 1300–1800 mN, pics ≈ 1650–1700 mN avec dérive descendante ; **couple sur axe −250 à 1250 µNm : les vallées passent sous zéro** (couple négatif en fin de décharge) ; pics ≈ 750 µNm. Ancrages : P crête 1,35 (0,5) / 1,32 (1,2) MPa ; F 946,56 → 1671,67 mN (+76,6 %) ; T 729,73 → 674,94 µNm (−7,5 %).
- **Fig. 9 (p. 12)** : boucles du 11e cycle, mandrin 3,0 mm, légende **λ_k** (= pré-déformation). (a) couple vs pression : axes 0,05–1,30 MPa / −150–900 µNm ; λ_k = 0,5 et 1,2 ; croissance quasi linéaire en pressurisation, boucle étroite, remontée transitoire puis décroissance en décharge ; l'aire de boucle change peu de 0,5 à 1,2. (b) force vs pression : 0–1,5 MPa / 500–2100 mN ; λ_k = 0,5 (≈ 700 mN à P≈0), 1,0 (≈ 1130 mN), 1,2 (≈ 1440 mN) ; plateau initial (dead-band ≈ 0–0,25 MPa, attribué au fluage : l'injection est compensée par l'expansion du diamètre interne) puis montée dominée par la contraction.
- **Fig. 10 (p. 13)** : cartes théoriques (11 cycles dont 10 de training ; 10 mL/min, 1,5 mL, 1,4 MPa) en fonction de C (2,0–3,0) et ε_tk (0,5–2,0). (a) T* en **rad/mm**, échelle 0,22–0,52 : pic de T* à C intermédiaire (position du pic indépendante de ε_tk), croissance monotone avec ε_tk. (b) F* en **MPa**, échelle 0,15–1,25 : F* croît quand C diminue et quand ε_tk croît ; C ≈ 2 + ε_tk > 1,5 = zone de force maximale.
- **Fig. 11 (p. 14)** : effet de la vitesse de pressurisation (C = 2,25, ε_tk = 2,0, relaxation 60 s post-étirement). (a) couple vs pression, 600–2600 µNm / 0–1,4 MPa : à 1,3 MPa, T = 2454,91 (1 MPa/s), 2352,72 (0,5), 2033,75 (0,1) µNm ; à 0,1 MPa/s, retour à 0,01 MPa : T = 701,62 µNm (−27,7 %). (b) force vs pression, 3,55–3,95 N : 3,88 / 3,93 / 3,94 N à 1,3 MPa ; retour 3,51 N (−3,1 %). Boucle quasi linéaire à 1 MPa/s ; le couple est plus sensible à Ṗ que la force.
- **Fig. B1 (p. 18)** : photo microscope du TCPA (échelle 100 µm) annotée : 2(ρ₀ + R_out0), 2R_out0, α₀, θ_f0, microfibres visibles.

---

## 6. Ambiguïtés et points laissés implicites (à l'attention des auditeurs)

1. **Incohérence E_axial** : Table B1 donne E_axial = 31,24 ± 1,90 MPa, mais l'annexe A affirme E_axial = E₀+E₁+E₂+E₃ = **37,76 MPa** (Table A1). Écart ~21 %, hors barre d'erreur. Le papier ne dit pas laquelle est utilisée dans les simulations.
2. **Circularité de la définition C̄_k** (p. 17) : C̄_k = [E_k/ΣE]·C̄_E alors que C̄_E = C̄₀ + ΣC̄ᵢ. Lecture cohérente : C̄_E est la matrice effective complète issue de (A.1)–(A.3) (avec E_axial = ΣE_k), et chaque branche en porte la fraction E_k/ΣE — mais ce n'est jamais dit explicitement.
3. **C₅₅ jamais défini** dans (A.1)/(A.2) (par isotropie transverse on attendrait C₅₅ = C₆₆ = G₁₂, non écrit).
4. **Coquilles probables dans (A.3)** : C̄₁₁ contient « C₁₆ », entrée nulle de (A.1) (attendu : 2C₆₆ dans la formule standard de rotation) ; C̄₂₂ contient C₄₄ (= (C₂₂−C₂₃)/2 en tant qu'entrée de (A.1)) là où C̄₁₁ utilise C₆₆ — asymétrie suspecte. C₁₃ et C₃₃ doivent se lire comme entrées de (A.1) (C₁₃ = C₁₂, C₃₃ = C₂₂). Un code fidèle au papier et un code fidèle à la mécanique standard peuvent différer ici.
5. **« Réduction de 65,8 % » du dead-band** (p. 6) : 0,38 → 0,23 MPa correspond arithmétiquement à −39,5 % (0,23/0,38 = 60,5 %). Le chiffre 65,8 % est inexpliqué.
6. **(A.2) sur-écrit** : deux lignes distinctes définissent « C₁₂ » ; C₂₃ n'a pas de ligne propre (il se déduit du système) ; ν₃₁, ν₃₂ définis mais inutilisés.
7. **Notation incrémentale incohérente** : Éq. (10) intègre « σ_s dt » sur le temps (cumul d'incréments non formalisé) ; Éq. (18) écrit u, ω, ν, κ sans Δ alors que le cadre est incrémental ; le texte renvoie plusieurs fois à « Eq. (A.4) » pour des objets différents (constitutive, incréments de déformation — en réalité Éq. 18 —, efforts internes — en réalité Éq. 10).
8. **Éq. (19) sous-déterminée en dimensions** : C̄_E y est implicitement une matrice 4×4 réduite (composantes r, φ, s, φs) ; la formule de Ē_L (module effectif du tube torsadé) n'est pas donnée (renvoi à [17]).
9. **Renvois de figures erronés** : §3.3 renvoie à « Figure 5(c),(d) » pour la décomposition en couches (c'est la Fig. 6) ; §3.2 renvoie à « Figure 6(a),(b) » pour les états t₀/t_k (c'est la Fig. 5) et à « Figure 6(d) » pour la phase de traction ; « Figure 4(b) » est cité (p. 7) alors que la Fig. 4 n'a pas de panneaux. Paragraphe corrompu p. 14 (« C%MathType!End!... ») : résidu d'éditeur d'équations recouvrant la phrase « minimiser l'indice de ressort (p. ex. C = 2) ».
10. **Notations concurrentes pour la pré-déformation** : ε_tk dans le texte/Éq. (8), **λ_k** dans la légende de la Fig. 9.
11. **Fig. 10, paramètres douteux** : « The outer diameter of the PVC tube R_out = 1 mm » (R_out est un rayon ; 2R_out0 = 2,0 mm en Table B1 — cohérent si l'on lit « rayon extérieur = 1 mm ») et « bias angle α_tk = 37° » (α désigne partout ailleurs l'angle d'hélice, 10,53° ; 37° est l'angle de biais θ_f).
12. **Enveloppes de pression incohérentes** : résumé « 0–1,2 MPa » ; §2.3 « crête 1,5 MPa, long terme 1,2–1,3 MPa » ; Fig. 7 « 1,41–1,43 MPa » ; Fig. 8 « 1,32–1,35 MPa » ; simulations Fig. 10 « 1,4 MPa ».
13. **Volumes injectés variables selon les cas** : 1,22 mL (§2.3), 1,50 mL (Fig. 7), 1,76 mL (Fig. 8), 1,5 mL (Fig. 10) — sans explication du choix.
14. **Entrée pression du modèle non spécifiée** : le modèle prend P(t) (via ΔP dans Éq. 15) mais l'article ne dit pas si P(t) simulé provient du signal expérimental ou d'un modèle volume→pression (aucun couplage volume injecté/pression n'est formulé) ; à la Fig. 11, Ṗ est idéalisé constant.
15. **Paramètres numériques d'implémentation absents** : nombre de couches n, pas de temps Δt, nombre de tours N, longueur du fil — jamais chiffrés ; N n'est pas dans la Table B1 (l'initial length 32,45 mm et ρ₀, α₀ permettraient de le reconstituer, non fait dans le papier).
16. **ν₁₂ = 0,205 et ν₂₃ = 0,422 ne sont pas mesurés** : repris de la réf. [15] (Cavatappi, autre matériau/procédé).
17. **Le dead-band n'est pas une équation du modèle** : il est présenté comme émergent (fluage compensant l'injection) ; aucun critère explicite n'est donné.
18. **Sous-estimation déclarée de l'hystérésis** : conséquence directe de l'hyp. (ii) ; les auteurs indiquent que découpler les anisotropies élastique/visqueuse corrigerait ce biais (piste, non implémentée).
19. **Nylon transversalement isotrope dans B1 mais isotrope dans (4)** : seuls E_nylon et G_nylon apparaissent dans le modèle (poutre) ; le caractère transverse isotrope annoncé n'est pas exploité.
20. **Fig. 8(b)** : couples négatifs en fin de décharge (axe jusqu'à −250 µNm) — comportement que le modèle reproduit ; un code qui borne le couple à 0 serait infidèle au papier.

Fichiers de travail (scratchpad) : texte brut extrait `blocked_text.txt` ; rendus pleine page `blocked_p01.png` … `blocked_p20.png` ; zooms `p18_A2_zoom.png`, `p18_Ck_zoom.png` dans `C:\Users\B3DCB~1.PER\AppData\Local\Temp\claude\C--Users-b-pereiraazevedo-OneDrive---House-Of-HR-NV-Documents-Stage-muscle-artifici-le-mod-le-mod-le-chinois-Espace-de-travail\c5852e10-6d38-4703-b218-a1649a5b12c9\scratchpad`.