p-difluorobenzene-water: relaxed 3D grid of the water orientation (theta, alpha, tau), wB97X-D/def2-TZVP.
Same design as the phenol-water step 3 (see pdfb_geom.py for the coordinate definitions).
594 constrained optimisations: theta = 2..120 (15-degree steps), alpha = 0..180 (22.5), tau = 0..157.5 (22.5).
theta, alpha, tau frozen; O_w position, H-bond geometry, water and ring relax. Job list ordered theta 15-90 first.

  sbatch -c 4 --mem=12G --array=1-594%12 slurm_pdfb.sh        (12 x 4 = 48 cores; lower %12 to share the cluster)
  python3 collect.py                                         -> grid_pdfb3d.dat, grid_pdfb3d.xyz ; lists unconverged
  python3 make_restart.py inputs
  sbatch -c 4 --mem=12G --array=1-N%12 slurm_pdfb.sh restart
  python3 collect.py inputs inputs_restart
Send grid_pdfb3d.dat and grid_pdfb3d.xyz.

Checks done: all 594 Z-matrices rebuild the target theta/alpha/tau to < 1e-8 deg; closest water-H contact in the
start geometries 1.90 A (O_w moved back along O_w->F where needed, which leaves theta/alpha/tau unchanged);
collect.py and make_restart.py tested on synthetic logs (converged and unconverged).
