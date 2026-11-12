anisole-water: CP-MP2/aug-cc-pVTZ at every converged grid geometry (wB97X-D relaxed) + the wB97X-D minimum (reference).
Unzip next to tunnel_anis3d:
  cd ~/Rama/Non-Covalent/mp2anis
  python3 make_mp2grid.py ../tunnel_anis3d/grid_anis3d.xyz        -> N pending
  sbatch --array=1-N%3 slurm_mp2grid.sh                             (16 cores, 40 GB each; ~20 min per point)
Later (restarted grid points): rerun make_mp2grid.py with the new .xyz ONLY when no array of this folder is running.
Collect:  python3 collect_mp2grid.py ../tunnel_anis3d/grid_anis3d.dat  -> grid_anis3d_mp2.dat  (send it)
