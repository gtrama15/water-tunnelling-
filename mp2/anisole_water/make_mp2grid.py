#!/usr/bin/env python3
"""CP-MP2/aug-cc-pVTZ single point at every converged anisole-water grid geometry (wB97X-D relaxed).
  python3 make_mp2grid.py ../tunnel_anis3d/grid_anis3d.xyz [more .xyz]
Writes inputs/m_thTTT_aAAAA_tTTTT.gjf for geometries without a finished log, and joblist.txt (pending points; the
reference minimum m_sp_min first). Safe to rerun later, but ONLY when no array of this folder is running."""
import sys, os, re
import numpy as np
done = lambda n: os.path.exists(f'inputs/{n}.log') and 'Normal termination' in open(f'inputs/{n}.log').read()[-3000:]
os.makedirs('inputs', exist_ok=True); jobs = [] if done('m_sp_min') else [(-1, 0, 0, 0, 'm_sp_min')]; new = 0
for xyz in sys.argv[1:] or ['../tunnel_anis3d/grid_anis3d.xyz']:
    L = open(xyz).read().split('\n'); i = 0
    while i < len(L) and L[i].strip():
        n = int(L[i]); t = L[i + 1]; at = [l.split() for l in L[i + 2:i + 2 + n]]; i += 2 + n
        th, al, ta = [float(x) for x in re.search(r'theta=([\d.]+) alpha=([\d.]+) tau=([\d.]+)', t).groups()]
        nm = f'm_th{th:03.0f}_a{al*10:04.0f}_t{ta*10:04.0f}'
        if done(nm): continue
        el = [a[0] for a in at]; X = np.array([[float(v) for v in a[1:4]] for a in at])
        D = np.linalg.norm(X[:, None] - X[None], axis=2)
        Ow = [k for k, e in enumerate(el) if e == 'O' and sum(el[j] == 'H' and D[k, j] < 1.15 for j in range(n)) == 2]
        assert len(Ow) == 1, nm
        W = {Ow[0]} | {j for j in range(n) if el[j] == 'H' and D[Ow[0], j] < 1.15}
        body = ''.join(f'{el[k]}(Fragment={2 if k in W else 1}) {X[k,0]:14.8f} {X[k,1]:14.8f} {X[k,2]:14.8f}\n' for k in range(n))
        if not os.path.exists(f'inputs/{nm}.gjf'): new += 1
        open(f'inputs/{nm}.gjf', 'w').write(f'%nprocshared=16\n%mem=32GB\n%chk={nm}.chk\n#p MP2/aug-cc-pVTZ Counterpoise=2 SCF=Tight NoSymm MaxDisk=170GB\n\n'
                                             f'anisole-water grid theta={th} alpha={al} tau={ta}: CP-MP2/aug-cc-pVTZ at the wB97XD geometry\n\n0 1 0 1 0 1\n{body}\n')
        jobs.append(((0 if 15 <= th <= 90 else 1), th, al, ta, nm))
jobs.sort(); open('joblist.txt', 'w').write(''.join(j[4] + '\n' for j in jobs))
print(f'{new} new inputs; {len(jobs)} pending in joblist.txt\nsubmit:  sbatch --array=1-{len(jobs)}%3 slurm_mp2grid.sh')
