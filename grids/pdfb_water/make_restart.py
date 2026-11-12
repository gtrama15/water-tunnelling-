#!/usr/bin/env python3
"""Restart inputs for unconverged pDFB-water grid points:  python3 make_restart.py [inputs]
Each restart starts from the LOWEST-energy step of the failed run, keeps the frozen constants exactly, and adds
CalcFC. The variables are recomputed from that step's geometry using the input's own Z-matrix definition.
Writes inputs_restart/ and joblist_restart.txt."""
import sys, os, glob, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect import parse
import pdfb_geom as G
d = (sys.argv[1:] or ['inputs'])[0].rstrip('/'); out = d + '_restart'; os.makedirs(out, exist_ok=True); names = []
for log in sorted(glob.glob(os.path.join(d, 'pd_th*_a*_t*.log'))):
    name = os.path.basename(log)[:-4]; p = parse(log)
    if p is not None and p['ok']: continue
    src = open(os.path.join(d, name + '.gjf')).read()
    if p is None: open(os.path.join(out, name + '.gjf'), 'w').write(src); names.append(name); print(name, ': no steps, original input'); continue
    E = np.array([s['E'] for s in p['steps']]); k = int(E.argmin()); P = p['steps'][k]['Xfull']
    head, rest = src.split('\n\n', 1); title, body = rest.split('\n\n', 1)
    zlines = body.split(' Variables:')[0].strip().split('\n')[1:]          # after '0 1'
    if len(P) != len(zlines): print(f'  {name}: {len(P)} centres in log vs {len(zlines)} z-matrix lines; skipped'); continue
    consts = body.split(' Constants:')[1]; cn = {l.split()[0] for l in consts.strip().split('\n') if l.strip()}
    var = {}
    for i, l in enumerate(zlines):
        t = l.split()
        if len(t) >= 3 and t[2] not in cn: var.setdefault(t[2], G.bond(P[i], P[int(t[1]) - 1]))
        if len(t) >= 5 and t[4] not in cn: var.setdefault(t[4], G.ang(P[i], P[int(t[1]) - 1], P[int(t[3]) - 1]))
        if len(t) >= 7 and t[6] not in cn: var.setdefault(t[6], G.dih(P[i], P[int(t[1]) - 1], P[int(t[3]) - 1], P[int(t[5]) - 1]))
    head = head.replace('Opt=(Z-Matrix,', 'Opt=(Z-Matrix,CalcFC,') if 'CalcFC' not in head else head
    head = head.replace(name + '.chk', name + '_rs.chk')
    txt = head + '\n\n' + title + ' RESTART\n\n0 1\n' + '\n'.join(zlines) + '\n Variables:\n' + \
          ''.join(f'  {k:6s} {v:14.8f}\n' for k, v in var.items()) + ' Constants:' + consts
    open(os.path.join(out, name + '.gjf'), 'w').write(txt); names.append(name)
    print(f'{name}: restart from step {k + 1}/{len(E)}')
open('joblist_restart.txt', 'w').write(''.join(n + '\n' for n in names))
print(f'{len(names)} restart inputs in {out}/ ; joblist_restart.txt\nsubmit:  sbatch -c 4 --mem=12G --array=1-{len(names)}%12 slurm_pdfb.sh restart')
