#!/usr/bin/env python3
"""Step 3b: tau refinement of the relaxed 3D grid. tau = 22.5, 67.5, 112.5, 157.5 (between the step-3 nodes),
theta = 110..170 step 10 and 179, alpha = 0..180 step 22.5.  On alpha = 0, 180 the mirror gives
tau=157.5 == 22.5 and 112.5 == 67.5, so those are skipped -> 256 constrained optimisations (theta, alpha, tau frozen).
Names: p3_thTTT_aAAAA_tTTTT  (alpha and tau x10)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_scan2 as M
THETAS=[110,120,130,140,150,160,170,179]; ALPHAS=[0,22.5,45,67.5,90,112.5,135,157.5,180]
TAUS=[22.5,67.5,112.5,157.5]
ROUTE='#p wB97XD/def2TZVP Opt=(Z-Matrix,MaxCycles=100,MaxStep=10) NoSymm Int=UltraFine'
os.makedirs('inputs_3dfine',exist_ok=True); M.FROZEN.add('D17'); jobs=[]
for th in THETAS:
    for al in ALPHAS:
        for ta in TAUS:
            if al in (0,180) and ta>90: continue
            n=f'p3_th{th:03d}_a{int(round(al*10)):04d}_t{int(round(ta*10)):04d}'
            M.write_input(f'inputs_3dfine/{n}.gjf',f'phenol-water relaxed 3D grid theta={th} alpha={al} tau={ta} (all three frozen)',
                          M.geometry(th,al,ta),ROUTE,'4','8GB')
            jobs.append(((0 if 130<=th<=170 else 1),th,al,ta,n))
jobs.sort()                                   # theta 130-170 first (most of the ground-state density)
open('joblist_3dfine.txt','w').write(''.join(j[4]+'\n' for j in jobs))
print(len(jobs),'inputs in inputs_3dfine/')
