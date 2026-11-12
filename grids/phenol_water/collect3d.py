#!/usr/bin/env python3
"""Collect the 3D-grid logs:  python3 collect3d.py [dirs ...]   (default inputs_3d inputs_3d_restart)
Writes grid_3d.dat  (theta alpha tau E(hartree) dE(cm-1, vs lowest in file) checks status) and geoms_3d.xyz"""
import numpy as np, glob, re, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from collect import parse, rescue, HART
from check_scan2 import active as _active, frame_from
def tau_of(X,al):
    Ow=X[13]; u,n,w=frame_from(X); b=(X[14]+X[15])/2-Ow; b/=np.linalg.norm(b); a=np.radians(al)
    e0=np.cos(a)*n-np.sin(a)*w; e=X[14]-X[15]; e-=(e@b)*b; e/=np.linalg.norm(e)
    return np.degrees(np.arctan2(np.cross(e0,e)@b,e0@e))%180
dirs=sys.argv[1:] or ['inputs_3d','inputs_3d_restart']
files={}
for d in dirs:
    for f in glob.glob(os.path.join(d,'p3_th*_a*_t*.log')): files.setdefault(os.path.basename(f),[]).append(f)
rows=[];bad=[]
fx=open('geoms_3d.xyz','w')
for base in sorted(files):
    th,al,ta=re.findall(r'th(\d+)_a(\d+)_t(\d+)',base)[0]; th,al,ta=float(th),float(al)/10,float(ta)
    best=None
    for f in files[base]:
        q=parse(f)
        if q is None: continue
        st='OK' if q['ok'] else None
        if not q['ok']:
            k,still,prc=rescue(q,1.0)
            if k is not None: s=q['steps'][k]; q.update(E=s['E'],X=s['X'],lab=s['lab']); st='RESCUED'
        if st and (best is None or q['E']<best[0]['E']): best=(q,st)
    if best is None: bad.append(base); continue
    q,st=best; tc,ac=_active(q['X']); tt=tau_of(q['X'],al)
    rows.append((th,al,ta,q['E'],tc,abs(ac),tt,st))
    fx.write('16\ntheta=%.1f alpha=%.1f tau=%.1f E=%.10f\n'%(th,al,ta,q['E']))
    for l,x in zip(q['lab'],q['X']): fx.write('%-2s %14.8f %14.8f %14.8f\n'%(l,*x))
fx.close()
if not rows: sys.exit('no converged points yet')
Emin=min(r[3] for r in rows)
with open('grid_3d.dat','w') as f:
    f.write('#  theta   alpha     tau        E(hartree)       dE(cm-1)  th_chk  a_chk  tau_chk status\n')
    for r in sorted(rows,key=lambda r:(r[2],r[1],r[0])):
        f.write('%7.1f %7.1f %7.1f %18.10f %10.2f %7.2f %7.2f %7.2f  %s\n'%(r[0],r[1],r[2],r[3],(r[3]-Emin)*HART,r[4],r[5],r[6],r[7]))
print(f'{len(rows)} usable points -> grid_3d.dat, geoms_3d.xyz ;  {len(bad)} still to (re)run')
dr=[(abs(r[4]-r[0]),abs(r[5]-r[1]),min(abs(r[6]-r[2]),180-abs(r[6]-r[2]))) for r in rows if r[0]<179]
print('max drift theta %.3f alpha %.3f tau %.3f deg'%tuple(np.max(np.array(dr),axis=0)))
if bad: print('missing/unconverged:',' '.join(b[:-4] for b in bad[:30]),'...' if len(bad)>30 else '')
