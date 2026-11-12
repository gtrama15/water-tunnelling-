# 3D water-rotation calculation (phenol–water)

`dvr3d.py`: full rotation of a rigid water about O_w on the (θ, α, τ) surface.

- **Basis.** Wigner functions D^j_mk, j ≤ jmax. This has no coordinate singularity, so the collinear geometry θ = 180 (where α and τ describe the same motion) is handled exactly.
- **Kinetic energy.** Exact orientation-dependent G for rotation about the water axes, with the complex's translation and rotation projected out (same recipe as `dvr2b.py`).
- **H exchange.** Equal to (−1)^k. Even k gives the symmetric states, odd k the antisymmetric ones. The splitting is E0(odd) − E0(even).
- **Potential.** V0(θ, α) from `grid_pt.dat` (τ = 0), plus dV(θ, α, τ) from the step-3 grid (τ = 45, 90, 135; 4-point Fourier in τ, period 180°), plus a wall K(90 − θ)² where there is no data, plus the optional MP2 correction.

```
python3 dvr3d.py --jmax 16 20 --iso H2O D2O                 # wB97X-D surface
python3 dvr3d.py --jmax 20 --mp2                            # + CP-MP2 correction (fit_mp2corr.py, mp2corr.dat)
python3 dvr3d.py --grid3d grid_3d.dat --nofill ...          # final: refuse to run if grid nodes are missing
python3 fit_mp2corr.py [mp2corr.dat]                        # fit report: point-by-point MP2, wB97X-D, diff, fit
```

## Checks

| Check | Result |
|---|---|
| (θ, α, τ) → geometry → (θ, α, τ) round trip | 1e-13° |
| Basis orthonormality on the quadrature grid | 2e-14 |
| Hermiticity | 3e-13 |
| Constant-G free rotor vs. exact asymmetric-top levels | exact to 1e-6 cm⁻¹ (all levels and degeneracies) |
| jmax 12 → 16 → 20 | splitting 0.4800 → 0.4814 → 0.4808 cm⁻¹ (ZPE 213.40 → 213.43 → 213.42) |
| Wall K = 1 / 5 / 20 cm⁻¹ deg⁻² | splitting 0.4805 / 0.4814 / 0.4830 cm⁻¹ |

Note: ARPACK cannot find an *exactly* zero eigenvalue (the free-rotor j = 0 state). This is irrelevant for any real potential.
