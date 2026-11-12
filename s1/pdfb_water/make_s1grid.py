#!/usr/bin/env python3
"""TD-wB97X-D/def2-TZVP vertical excitations at every pDFB-water grid geometry.
  python3 make_s1grid.py ../tunnel_pdfb3d/grid_pdfb3d.xyz
Writes inputs/s1_pd_thTTT_aAAAA_tTTTT.gjf for geometries without a finished log, and joblist.txt (pending ones;
the reference minimum s1_pd_sp_min comes first). Safe to rerun when more grid points converge."""
import sys, os, re
os.makedirs('inputs', exist_ok=True); jobs = []
ROUTE = '#p TD(NStates=4) wB97XD/Def2TZVP SCF=Tight NoSymm Int=UltraFine'
done = lambda n: os.path.exists(f'inputs/{n}.log') and 'Normal termination' in open(f'inputs/{n}.log').read()[-3000:]
if not done('s1_pd_sp_min'): jobs.append('s1_pd_sp_min')
for xyz in sys.argv[1:] or ['../tunnel_pdfb3d/grid_pdfb3d.xyz']:
    L = open(xyz).read().split('\n'); i = 0
    while i < len(L) and L[i].strip():
        n = int(L[i]); t = L[i + 1]; at = L[i + 2:i + 2 + n]; i += 2 + n
        th, al, ta = [float(x) for x in re.search(r'theta=([\d.]+) alpha=([\d.]+) tau=([\d.]+)', t).groups()]
        nm = f's1_pd_th{th:03.0f}_a{al*10:04.0f}_t{ta*10:04.0f}'
        if done(nm): continue
        open(f'inputs/{nm}.gjf', 'w').write(f'%nprocshared=4\n%mem=8GB\n%chk={nm}.chk\n{ROUTE}\n\nS1 vertical theta={th} alpha={al} tau={ta} '
                                            f'(pDFB-water S0 wB97XD relaxed geometry)\n\n0 1\n' + '\n'.join(at) + '\n\n')
        jobs.append(nm)
open('joblist.txt', 'w').write(''.join(j + '\n' for j in jobs))
print(f'{len(jobs)} pending in joblist.txt;  submit:  sbatch -c 4 --mem=12G --array=1-{len(jobs)}%4 slurm_s1.sh')
