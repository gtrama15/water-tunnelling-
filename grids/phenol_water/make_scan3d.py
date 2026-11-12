#!/usr/bin/env python3
"""Step 3: relaxed 3D grid of phenol-water in the water orientation (theta, alpha, tau).
theta, alpha as in step 2; tau = rotation of the water about its own bisector, from the co-rotating
orientation (tau = 0 is the step-2b slice, already computed).  theta, alpha AND tau are frozen;
everything else (O...O distance, H-bond angles, water and phenol geometry) relaxes.
Symmetry: H-label exchange tau -> tau+180, and the phenol-plane mirror (alpha,tau) -> (-alpha,-tau),
so alpha = 0..180 with tau = 0..180 covers everything; on the in-plane lines alpha = 0, 180 also tau -> -tau.
Grid: tau = 45, 90, 135 (tau=135 skipped on alpha = 0, 180 where it equals tau=45)
      theta = 110..170 step 10 and 179;  alpha = 0..180 step 22.5  ->  200 constrained optimisations.
"""
import numpy as np, os, argparse, sys
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import make_scan2 as M
THETAS=[110,120,130,140,150,160,170,179]
ALPHAS=[0,22.5,45,67.5,90,112.5,135,157.5,180]
TAUS=[45,90,135]
def priority(th,al,ta):
    if ta==90: return 1                       # in-plane torsion route (water in the phenol plane at alpha=0)
    if 120<=th<=160: return 2
    return 3
if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--route',default='#p wB97XD/def2TZVP Opt=(Z-Matrix,MaxCycles=100,MaxStep=10) NoSymm Int=UltraFine')
    p.add_argument('--nproc',default='10'); p.add_argument('--mem',default='4GB'); p.add_argument('--out',default='inputs_3d')
    a=p.parse_args(); os.makedirs(a.out,exist_ok=True); M.FROZEN.add('D17'); jobs=[]
    for ta in TAUS:
        for th in THETAS:
            for al in ALPHAS:
                if ta==135 and al in (0,180): continue
                name=f'p3_th{th:03d}_a{int(round(al*10)):04d}_t{ta:03d}'
                M.write_input(f'{a.out}/{name}.gjf',f'phenol-water relaxed 3D grid theta={th} alpha={al} tau={ta} (all three frozen)',
                              M.geometry(th,al,ta),a.route,a.nproc,a.mem)
                jobs.append((priority(th,al,ta),th,al,ta,name))
    jobs.sort()
    # split: cluster gets every 5th job from the priority-ordered list... simpler: alternate 2 of 5 to cluster
    cl=[j for k,j in enumerate(jobs) if k%5 in (1,3)]; ws=[j for k,j in enumerate(jobs) if k%5 not in (1,3)]
    open('joblist_3d.txt','w').write(''.join(f'{j[4]} {j[0]}\n' for j in jobs))
    open('joblist_3d_cluster.txt','w').write(''.join(f'{j[4]} {j[0]}\n' for j in cl))
    open('joblist_3d_ws.txt','w').write(''.join(f'{j[4]} {j[0]}\n' for j in ws))
    print(len(jobs),'inputs in',a.out,'| cluster',len(cl),'| workstation',len(ws),'| priorities',{k:sum(j[0]==k for j in jobs) for k in (1,2,3)})
