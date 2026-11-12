#!/usr/bin/env python3
"""Build the MP2//wB97X-D grid:  python3 collect_mp2grid.py ../tunnel_pdfb3d/grid_pdfb3d.dat
Writes grid_pdfb3d_mp2.dat in the same format as grid_pdfb3d.dat, with E = E_wB(min) + [E_MP2 - E_MP2(min)]
(so energies stay relative to the same reference), for every grid point that has a finished MP2 log."""
import sys, os, re, glob
H = 219474.6354
def mp2(f):
    if not os.path.exists(f): return None
    t = open(f, errors='replace').read()
    if 'Normal termination' not in t: return None
    m = re.findall(r'[Cc]ounterpoise corrected energy\s*=\s*(-?\d+\.\d+)', t)
    return float(m[-1]) if m else None
Emp2_min = mp2('../mp2pdfb/inputs/mp2cp_sp_min.log'); E_wb_min = -507.2002529730
assert Emp2_min, 'need ../mp2pdfb/inputs/mp2cp_sp_min.log'
out = []; have = miss = 0
for l in open(sys.argv[1] if len(sys.argv) > 1 else '../tunnel_pdfb3d/grid_pdfb3d.dat'):
    if l.startswith('#'): out.append(l); continue
    s = l.split(); th, al, ta = map(float, s[:3])
    e = mp2(f'inputs/m_th{th:03.0f}_a{al*10:04.0f}_t{ta*10:04.0f}.log')
    if e is None and al == round(al) and ta == round(ta): e = mp2(f'../mp2pdfb/inputs/mp2cp_g_th{th:03.0f}_a{al:03.0f}_t{ta:03.0f}.log')
    if e is None: miss += 1; continue
    have += 1; E = E_wb_min + (e - Emp2_min); s[3] = f'{E:.10f}'; s[4] = f'{(E - E_wb_min) * H:.2f}'
    out.append('  '.join(s) + '\n')
open('grid_pdfb3d_mp2.dat', 'w').writelines(out)
print(f'{have} grid points with MP2 -> grid_pdfb3d_mp2.dat ; {miss} without MP2 yet')
