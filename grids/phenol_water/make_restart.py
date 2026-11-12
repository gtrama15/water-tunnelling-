#!/usr/bin/env python3
"""Write restart inputs for scan points that did not converge (and were not rescued by collect.py).
The restart starts from the LOWEST-energy geometry of the failed run, keeps the frozen r1/tau constants
exactly as in the original input, and computes the force constants once at the start (CalcFC), which
usually makes a floppy optimisation converge within ~20-40 steps.
Usage:  python3 make_restart.py [inputs] [--only-unrescued] [--tol 1.0] [--maxcyc 100] [--nocalcfc]
        default: restart EVERY unconverged run (rescued values stay as fallback in collect.py)
        --only-unrescued : skip the runs collect.py accepts as RESCUED
Writes  inputs_restart/<name>.gjf  and  joblist_restart.txt
Run them:  sbatch --array=1-N%5 slurm_array.sh restart   (the script prints the exact command)
"""
import sys, os, glob, re, argparse
import numpy as np
from collect import parse, rescue
from make_scan2 import zmatrix

ap=argparse.ArgumentParser(); ap.add_argument('dir',nargs='?',default='inputs'); ap.add_argument('--prefix',default='ps')
ap.add_argument('--only-unrescued',action='store_true',help='skip runs that collect.py would accept as RESCUED'); ap.add_argument('--tol',type=float,default=1.0)
ap.add_argument('--maxcyc',type=int,default=100); ap.add_argument('--nocalcfc',action='store_true')
a=ap.parse_args()
out=a.dir.rstrip('/')+'_restart'; os.makedirs(out,exist_ok=True)
names=[]
for log in sorted(glob.glob(os.path.join(a.dir,a.prefix+'_th*_a*.log'))):
    name=os.path.basename(log)[:-4]
    p=parse(log)
    if p is not None and p['ok']: continue
    gjf=os.path.join(a.dir,name+'.gjf')
    if not os.path.exists(gjf): print('  missing',gjf); continue
    src=open(gjf).read()
    if p is None:                                   # crashed before any step: just resubmit the original
        open(os.path.join(out,name+'.gjf'),'w').write(src); names.append(name); print(f'{name}: no steps, original input'); continue
    k,still,prc=rescue(p,a.tol)
    if k is not None and a.only_unrescued: continue
    E=np.array([s['E'] for s in p['steps']]); kmin=int(E.argmin())
    P=list(p['steps'][kmin]['Xfull'])
    if len(P)!=18: print(f'  {name}: expected 18 centres incl. dummies, found {len(P)}; skipped'); continue
    _,var,_=zmatrix(P)
    # header (Link0 + route) and z-matrix body + constants from the original input; only the variables change
    head,rest=src.split('\n\n',1)
    lines=head.split('\n')
    route=[i for i,l in enumerate(lines) if l.startswith('#')][0]
    r=lines[route]
    m=re.search(r'Opt=\(([^)]*)\)',r,re.I)
    opts=[o for o in m.group(1).split(',') if not o.lower().startswith('maxcycles')]
    if not a.nocalcfc and not any(o.lower()=='calcfc' for o in opts): opts.append('CalcFC')
    opts.append(f'MaxCycles={a.maxcyc}')
    if not any(o.lower().startswith('maxstep') for o in opts): opts.append('MaxStep=10')
    lines[route]=r[:m.start()]+'Opt=('+','.join(opts)+')'+r[m.end():]
    lines=[l.replace(name+'.chk',name+'_rs.chk') for l in lines]
    body=rest.split(' Variables:')[0]
    consts=rest.split(' Constants:')[1]
    txt='\n'.join(lines)+'\n\n'+body.replace('relaxed scan','relaxed scan RESTART',1)+' Variables:\n'
    cnames={l.split()[0] for l in consts.strip().split('\n') if l.strip()}
    txt+=''.join(f'  {k:6s} {v:14.8f}\n' for k,v in var.items() if k not in cnames)+' Constants:'+consts
    open(os.path.join(out,name+'.gjf'),'w').write(txt); names.append(name)
    print(f'{name}: restart from lowest-energy step {kmin+1}/{len(E)}  (E={E[kmin]:.8f})')
JL=('joblist_restart.txt' if a.dir.rstrip('/')=='inputs' else 'joblist_'+os.path.basename(a.dir.rstrip('/'))+'_restart.txt'); MODE={'inputs':'restart','inputs_tau':'taurestart','inputs_3d':'3drestart'}.get(a.dir.rstrip('/'),'restart')
open(JL,'w').write(''.join(f'{n} 1\n' for n in names))
print(f'{len(names)} restart inputs in {out}/ ; job list {JL}')
print(f'submit on the cluster:   sbatch --array=1-{len(names)}%5 slurm_array.sh {MODE}')
