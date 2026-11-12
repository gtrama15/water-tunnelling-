#!/usr/bin/env python3
"""Collect CP-MP2 energies of the correction points:  python3 collect_mp2corr.py  -> mp2corr.dat"""
import re,glob,os
rows=[]
for f in sorted(glob.glob('inputs/mp2cp_*.log')):
    t=open(f).read(); cp=re.findall(r'[Cc]ounterpoise corrected energy\s*=\s*(-?\d+\.\d+)',t)
    if cp: rows.append((os.path.basename(f)[6:-4],float(cp[-1])))
open('mp2corr.dat','w').write(''.join(f'{n:24s} {E:.10f}\n' for n,E in rows))
print(len(rows),'finished -> mp2corr.dat')
