#!/usr/bin/env python3
"""Relaxed 3D grid of anisole-water in (theta, alpha, tau), all three frozen (see anis_geom.py).
theta = 2, 15 ... 120 (step 15), alpha = 0 ... 180 (22.5), tau = 0 ... 157.5 (22.5); on alpha = 0 and 180 the ring
mirror gives tau == -tau, so only tau <= 90 there  ->  9 x (7x8 + 2x5) = 594 constrained optimisations.
Names an_thTTT_aAAAA_tTTTT (alpha, tau x10). Job list ordered: theta 15-90 first (minimum, both saddles)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anis_geom as G
THETAS = [2, 15, 30, 45, 60, 75, 90, 105, 120]
ALPHAS = [i * 22.5 for i in range(9)]
TAUS = [i * 22.5 for i in range(8)]
ROUTE = '#p wB97XD/def2TZVP Opt=(Z-Matrix,CalcFC,MaxCycles=150,MaxStep=10) NoSymm Int=UltraFine'
if __name__ == '__main__':
    os.makedirs('inputs', exist_ok=True); jobs = []
    for th in THETAS:
        for al in ALPHAS:
            for ta in TAUS:
                if al in (0, 180) and ta > 90: continue
                n = f'an_th{th:03d}_a{round(al*10):04d}_t{round(ta*10):04d}'
                G.write_input(f'inputs/{n}.gjf', f'anisole-water relaxed 3D grid theta={th} alpha={al} tau={ta} (all three frozen)',
                              th, al, ta, ROUTE)
                jobs.append((0 if 15 <= th <= 90 else 1, th, al, ta, n))
    jobs.sort()
    open('joblist.txt', 'w').write(''.join(j[4] + '\n' for j in jobs))
    print(len(jobs), 'inputs in inputs/ ; joblist.txt')
