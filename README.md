# Water H-exchange tunnelling in aromatic–water complexes

Full 3D water-rotation (hindered rigid rotor, Wigner basis, exact G matrix) tunnelling splittings for
phenol–H2O, p-difluorobenzene–H2O and anisole–H2O, on relaxed wB97X-D/def2-TZVP (theta, alpha, tau) grids
with CP-MP2/aug-cc-pVTZ energies, and S1 from TD-wB97X-D excitation energies scaled to EOM-CCSD.

| folder | contents |
| --- | --- |
| `rotor/` | 3D rotor code: `dvr3d.py` (phenol–water), `dvr3d_pdfb.py`, `dvr3d_anis.py`; `fit_mp2corr.py` (phenol MP2 correction fit) |
| `grids/` | Z-matrix input generators, collect and restart scripts, SLURM scripts for the relaxed grids of each complex |
| `mp2/` | CP-MP2 single-point workflows (phenol: correction points; pDFB and anisole: every grid point) |
| `s1/` | TD-DFT vertical excitation workflows and the EOM-CCSD collector |
| `results/` | collected grids (`.dat`) and excitation energies used for the published numbers |

Gaussian logs and checkpoint files are not included. Note: some scripts contain absolute paths of the original
working directories (`/home/rama/tunnel/...`); adjust them to your layout before rerunning.

Example (phenol–water, MP2-corrected S0):
    cd rotor && MP2MODEL=ext python3 dvr3d.py --grid2d ../results/phenol_water/grid_pt.dat \
        --grid3d ../results/phenol_water/grid_3d.dat ../results/phenol_water/grid_3dfine.dat --nofill --jmax 16 --mp2
