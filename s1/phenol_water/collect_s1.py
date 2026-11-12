#!/usr/bin/env python3
"""python3 collect_s1.py  -> s1_vert.dat : name  E0(hartree)  then for states 1-4: E(eV) f ; plus the main
orbital transition of each state (to check which one is the phenol pi-pi*)"""
import glob,re,os
rows=[];bad=[]
for f in sorted(glob.glob('inputs/s1_*.log')):
    t=open(f).read(); nm=os.path.basename(f)[:-4]
    if 'Normal termination' not in t: bad.append(nm); continue
    E0=float(re.findall(r'SCF Done:\s+E\(\S+\)\s+=\s+(-?\d+\.\d+)',t)[-1])
    st=re.findall(r'Excited State\s+(\d+):\s+\S+\s+(-?\d+\.\d+) eV\s+\S+ nm\s+f=(-?\d+\.\d+)[^\n]*\n((?:\s+\d+\s*->\s*\d+\s+-?\d\.\d+\n)+)',t)
    s=[]
    for k,e,fo,tr in st[:4]:
        m=max(re.findall(r'(\d+)\s*->\s*(\d+)\s+(-?\d\.\d+)',tr),key=lambda x:abs(float(x[2])))
        s.append(f'{float(e):9.5f} {float(fo):7.4f} {m[0]}->{m[1]}')
    rows.append(f'{nm:26s} {E0:.10f}   '+'   '.join(s))
open('s1_vert.dat','w').write('\n'.join(rows)+'\n')
print(len(rows),'finished -> s1_vert.dat ;',len(bad),'not finished')
