anisole-water: relaxed 3D grid of the water orientation (theta, alpha, tau), wB97X-D/def2-TZVP; same design as
tunnel_pdfb3d (see anis_geom.py). u = O_w -> O(ether), n = normal of the (O_w, O_ether, C_ipso) plane.
594 constrained optimisations, theta 2..120 (15), alpha 0..180 (22.5), tau 0..157.5 (22.5); job list theta 15-90 first.
The inputs already use CalcFC and MaxCycles=150 (what the pDFB grid needed in its restarts).
Step-1 points in these coordinates: minimum theta 48 / alpha -65 / tau 93; flip TSs (288 cm-1) theta 62 / alpha 180 / tau 90;
bifurcated TS (343 cm-1) theta 3 / alpha 180 / tau 178.

  sbatch -c 4 --mem=12G --array=1-594%6 slurm_anis.sh
  python3 collect.py                    -> grid_anis3d.dat, grid_anis3d.xyz ; lists unconverged
  restarts as for pDFB:  python3 make_restart.py <folder with the failed .log/.gjf>   ->  <folder>_restart
Checks done: all 594 Z-matrices rebuild theta/alpha/tau to < 1e-8 deg; closest water-H contact in the start
geometries 1.90 A; collect.py and make_restart.py tested on synthetic logs.
