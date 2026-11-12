# Phenol–water: 3D water-rotation results (S0)

Surface: ωB97X-D/def2-TZVP relaxed (θ, α, τ) grid. 200 3D points plus the 142-point τ = 0 grid; all nodes from data.

MP2 correction: CP-MP2/aug-cc-pVTZ − ωB97X-D at 33 points, as a smooth fit:
- `ext` model: 16 parameters, rms 9.9 cm⁻¹
- `base` model: 12 parameters, rms 15.4 cm⁻¹

Method: hindered rigid rotor in a Wigner basis with the exact G(Ω). It is converged with basis size (jmax 16/20/24 agree to 1e-4 cm⁻¹).

| Surface | H2O splitting (cm⁻¹) | D2O splitting (cm⁻¹) | ZPE H2O / D2O (cm⁻¹) |
|---|---|---|---|
| ωB97X-D | 0.464 | 0.0160 | 209.5 / 156.6 |
| + MP2 corr. (ext fit), **recommended** | **0.487** | **0.0186** | 190.1 / 141.7 |
| + MP2 corr. (base fit) | 0.536 | 0.0225 | 185.7 / 139.3 |
| 2D model (τ frozen), ωB97X-D | 0.44 | 0.013 | 121 / – |

First excited rotor pair (H2O, ext): 262.42 (even) / 258.29 (odd). Its splitting is −4.1 cm⁻¹ (inverted; it was −5.0 on ωB97X-D).

## Size of the MP2 correction (cm⁻¹, relative to the minimum)
- Shoulder side (α ≥ 135, τ = 0): ωB97X-D too low by +22 to +199. The error is largest at θ = 110 and decreases toward the collinear geometry.
- In-plane water (τ = 90, α = 0): ωB97X-D too high by −106 (θ150), −227 (θ130), −366 (θ110).
- Minimum region and α ≤ 90 at θ ≥ 150: within ±25.

## Experiment
Berden et al. 1996 measure a 0.85 cm⁻¹ S1−S0 difference. Their rigid 1D model implies about 0.7 cm⁻¹ in S0. The S1 calculation is still to do.
