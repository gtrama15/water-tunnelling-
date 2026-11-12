#!/usr/bin/env python3
"""python3 collect_eom.py  -> eom.dat: job, ground-state energy, then each excited state: E(eV) f main-transition"""
import glob, re, os
rows = []
for f in sorted(glob.glob('inputs/*.log')):
    t = open(f, errors='replace').read(); n = os.path.basename(f)[:-4]
    if 'Normal termination' not in t: rows.append(f'{n:34s} NOT FINISHED'); continue
    E0 = re.findall(r'(?:E\(CORR\)=\s*|SCF Done:\s+E\(\S+\)\s+=\s+)(-?\d+\.\d+)', t)
    st = re.findall(r'Excited State\s+\d+:\s+\S+\s+(-?\d+\.\d+) eV\s+\S+ nm\s+f=(-?\d+\.\d+)[^\n]*\n((?:\s+\d+[AB]?\s*->\s*\d+[AB]?\s+-?\d\.\d+\n)+)', t)
    s = []
    for e, fo, tr in st:
        m = max(re.findall(r'(\d+)[AB]?\s*->\s*(\d+)[AB]?\s+(-?\d\.\d+)', tr), key=lambda x: abs(float(x[2])))
        s.append(f'{float(e):8.4f} {float(fo):7.4f} {m[0]}->{m[1]}')
    rows.append(f'{n:34s} {E0[-1] if E0 else "?":>18s}   ' + '   '.join(s))
open('eom.dat', 'w').write('\n'.join(rows) + '\n'); print(len(rows), 'jobs -> eom.dat')
