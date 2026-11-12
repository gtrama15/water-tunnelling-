#!/usr/bin/env python3
"""CP-MP2/aug-cc-pVTZ single point at every converged pDFB-water grid geometry (wB97X-D relaxed).
  python3 make_mp2grid.py ../tunnel_pdfb3d/grid_pdfb3d.xyz [more .xyz]
Writes inputs/m_thTTT_aAAAA_tTTTT.gjf for every geometry that has no finished log yet and is not one of the
36 points already done in ../mp2pdfb (those are reused), and joblist.txt (all pending points).
Safe to rerun after more grid points converge: only new geometries are added."""
import sys, os, re, glob
import numpy as np
done_old = {}
for f in glob.glob('../mp2pdfb/inputs/mp2cp_g_*.log'):
    if 'Normal termination' in open(f).read()[-3000:]:
        th, al, ta = [float(x) for x in re.findall(r'th(\d+)_a(\d+)_t(\d+)', f)[0]]; done_old[(th, al, ta)] = f
os.makedirs('inputs', exist_ok=True); new = 0; jobs = []
for xyz in sys.argv[1:] or ['../tunnel_pdfb3d/grid_pdfb3d.xyz']:
    L = open(xyz).read().split('\n'); i = 0
    while i < len(L) and L[i].strip():
        n = int(L[i]); t = L[i + 1]; at = [l.split() for l in L[i + 2:i + 2 + n]]; i += 2 + n
        th, al, ta = [float(x) for x in re.search(r'theta=([\d.]+) alpha=([\d.]+) tau=([\d.]+)', t).groups()]
        if (th, al, ta) in done_old: continue
        nm = f'm_th{th:03.0f}_a{al*10:04.0f}_t{ta*10:04.0f}'
        log = f'inputs/{nm}.log'
        if os.path.exists(log) and 'Normal termination' in open(log).read()[-3000:]: continue
        el = [a[0] for a in at]; X = np.array([[float(v) for v in a[1:4]] for a in at])
        O = el.index('O'); W = {O} | {k for k, e in enumerate(el) if e == 'H' and np.linalg.norm(X[k] - X[O]) < 1.15}
        assert len(W) == 3, nm
        body = ''.join(f'{el[k]}(Fragment={2 if k in W else 1}) {X[k,0]:14.8f} {X[k,1]:14.8f} {X[k,2]:14.8f}\n' for k in range(len(el)))
        if not os.path.exists(f'inputs/{nm}.gjf'): new += 1
        open(f'inputs/{nm}.gjf', 'w').write(f'%nprocshared=16\n%mem=32GB\n%chk={nm}.chk\n#p MP2/aug-cc-pVTZ Counterpoise=2 SCF=Tight NoSymm MaxDisk=170GB\n\n'
                                             f'pDFB-water grid theta={th} alpha={al} tau={ta}: CP-MP2/aug-cc-pVTZ at the wB97XD geometry\n\n0 1 0 1 0 1\n{body}\n')
        jobs.append(((0 if 15 <= th <= 90 else 1), th, al, ta, nm))
jobs.sort(); open('joblist.txt', 'w').write(''.join(j[4] + '\n' for j in jobs))
print(f'{len(done_old)} points reused from ../mp2pdfb; {new} new inputs; {len(jobs)} pending in joblist.txt')
print(f'submit:  sbatch --array=1-{len(jobs)}%3 slurm_mp2grid.sh')
