CP-MP2/aug-cc-pVTZ at every pDFB-water grid geometry (the MP2 - wB97X-D difference is too irregular to fit:
leave-one-out error >= 60 cm-1 on the 40 test points). Unzip next to tunnel_pdfb3d and mp2pdfb:
  cd ~/Rama/Non-Covalent/mp2grid
  python3 make_mp2grid.py ../tunnel_pdfb3d/grid_pdfb3d.xyz          (reuses the 36 points already in ../mp2pdfb)
  sbatch --array=1-N%3 slurm_mp2grid.sh                              (N printed by make_mp2grid.py; 16 cores, 40 GB each)
When the restart round finishes: rerun make_mp2grid.py with the new grid .xyz (adds only new points) and submit again.
Collect:  python3 collect_mp2grid.py ../tunnel_pdfb3d/grid_pdfb3d.dat   -> grid_pdfb3d_mp2.dat   (send it)
