#!/usr/bin/env python3
"""Collect the step-2 relaxed-scan logs (theta, alpha) into a grid file + geometries.
Usage:  python3 collect.py [inputs_dir ...] [--rescue 1.0]
Writes  grid_scan2.dat   theta alpha E(hartree) dE(cm-1) theta_check alpha_check tau status
        geoms_scan2.xyz  final geometry of every converged point (16 atoms, z-matrix order)
Checks that the frozen coordinates really stayed at their target values.
"""
import numpy as np, glob, re, sys, os
HART=219474.6354
Z2S={1:'H',6:'C',8:'O'}
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))

def parse(fn):
    """Return final data plus every optimisation step (energy, geometry without dummies, geometry with dummies)."""
    t=open(fn,errors='replace').read()
    ok_term='Normal termination' in t
    conv=('Optimization completed' in t) or ('Stationary point found' in t)
    orient=[(m.start(),m.group(1)) for m in re.finditer(r'(?:Input|Standard|Z-Matrix) orientation:.*?-{20,}\n.*?-{20,}\n(.*?)\n\s*-{20,}',t,re.S)]
    scf=[(m.start(),float(m.group(1))) for m in re.finditer(r'SCF Done:\s+E\(\S+\)\s+=\s+(-?\d+\.\d+)',t)]
    pred=[(m.start(),float(m.group(1).replace('D','E'))) for m in re.finditer(r'Predicted change in Energy=\s*(-?[\d.]+D[-+]\d+)',t)]
    if not orient or not scf: return None
    steps=[]
    for i,(pos,E) in enumerate(scf):                    # pair each energy with the geometry printed before it
        nxt=scf[i+1][0] if i+1<len(scf) else len(t)
        dE=[v for p,v in pred if pos<p<nxt]               # Gaussian's predicted energy gain for this step
        prev=[b for p,b in orient if p<pos]
        if not prev: continue
        rows=[l.split() for l in prev[-1].strip().split('\n')]
        full=np.array([[float(v) for v in r[3:6]] for r in rows])
        real=[i for i,r in enumerate(rows) if int(r[1])>0]          # drop dummy atoms (atomic number -1)
        lab=[Z2S.get(int(rows[i][1]),'X') for i in real]
        steps.append(dict(E=E,X=full[real],Xfull=full,lab=lab,pred=dE[0] if dE else None))
    if not steps: return None
    last=steps[-1]
    return dict(E=last['E'],lab=last['lab'],X=last['X'],ok=ok_term and conv,term=ok_term,conv=conv,steps=steps)

def rescue(p,tol,window=20):
    """Unconverged constrained optimisation. The frozen coordinates are exact at every step, so every step lies on
    the constrained surface; the lowest-energy step is accepted only if
      (a) the energy gained during the last `window` steps is < tol (cm-1), and
      (b) Gaussian's predicted remaining energy change at that step is < tol (cm-1)
          (guards against runs that hit their minimum early and then wandered uphill).
    Returns (step index or None, energy still gained over the last window, predicted change at best step)."""
    E=np.array([s['E'] for s in p['steps']])
    if len(E)<window+5: return None,None,None
    k=int(E.argmin()); still=(E[:-window].min()-E.min())*HART
    pr=p['steps'][k]['pred']; prc=abs(pr)*HART if pr is not None else None
    ok=still<tol and (prc is None or prc<tol)
    return (k if ok else None), still, prc

from check_scan2 import active as _active
def active(X):
    """theta, |alpha| and the water rotation tau about its own bisector (0 = H-H as in the minimum)."""
    th,al=_active(X)
    from check_scan2 import frame_from
    u,n,w=frame_from(X); a=np.radians(abs(al)); Ow=X[13]
    b=(X[14]+X[15])/2-Ow; b/=np.linalg.norm(b); e0=np.cos(a)*n-np.sin(a)*w
    e=X[14]-X[15]; e-=(e@b)*b; e/=np.linalg.norm(e)
    tau=np.degrees(np.arccos(min(1,abs(e@e0))))
    return th,abs(al),tau

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('dirs',nargs='*',default=['inputs'],help='e.g.  inputs inputs_restart  (for each point the best result is used)')
    ap.add_argument('--rescue',type=float,default=1.0,help='accept unconverged runs whose energy changed less than this (cm-1) over their last 20 steps; 0 = off')
    ap.add_argument('--prefix',default='ps',help='ps = tau relaxed (step 2), pt = tau frozen (step 2b)')
    a=ap.parse_args(); d=' '.join(a.dirs)
    OUT={'ps':('grid_scan2.dat','geoms_scan2.xyz')}.get(a.prefix,(f'grid_{a.prefix}.dat',f'geoms_{a.prefix}.xyz'))
    rows=[];bad=[]
    fx=open(OUT[1],'w')
    files={}
    for dd in a.dirs:
        for fn in glob.glob(os.path.join(dd,a.prefix+'_th*_a*.log')): files.setdefault(os.path.basename(fn),[]).append(fn)
    for base in sorted(files):
        fn=files[base][-1]
        r1,t=map(float,re.findall(r'_th(\d+)_a(\d+)',base)[0])
        cands=[]
        for f in files[base]:
            q=parse(f)
            if q is None: continue
            st='OK' if q['ok'] else 'NOTCONV'; info=None if q['ok'] else 'not converged'
            if not q['ok'] and a.rescue>0:
                k,still,prc=rescue(q,a.rescue)
                kb=int(np.argmin([s_['E'] for s_ in q['steps']]))
                pinfo=(', predicted gain %.2f cm-1'%prc) if prc is not None else ''
                if k is not None:
                    s_=q['steps'][k]; q.update(E=s_['E'],X=s_['X'],lab=s_['lab']); st='RESCUED'
                    info='RESCUED: lowest of %d steps (step %d)%s'%(len(q['steps']),k+1,pinfo)
                else:
                    info=('not converged' if q['term'] else 'stopped')+(': lowest at step %d/%d%s -> restart'%(kb+1,len(q['steps']),pinfo) if still is not None else ': only %d steps -> restart'%len(q['steps']))
            cands.append(({'OK':0,'RESCUED':1,'NOTCONV':2}[st],q['E'],f,q,st,info))
        if not cands: bad.append((fn,'no energy/geometry -> resubmit')); continue
        _,_,fn,p,status,info=min(cands,key=lambda c:(c[0],c[1]))
        if info: bad.append((fn,info))
        rc,ac,tc=active(p['X'])
        rows.append((r1,t,p['E'],rc,ac,tc,status))
        if status in ('OK','RESCUED'):
            fx.write('16\ntheta=%.1f alpha=%.1f E=%.10f\n'%(r1,t,p['E']))
            for l,x in zip(p['lab'],p['X']): fx.write('%-2s %14.8f %14.8f %14.8f\n'%(l,*x))
    fx.close()
    if not rows:
        for fn,why in bad: print('  ',os.path.basename(fn),':',why)
        sys.exit('no usable logs found in '+d)
    Emin=min(r[2] for r in rows if r[6] in ('OK','RESCUED'))
    with open(OUT[0],'w') as f:
        f.write('#  theta   alpha        E(hartree)        dE(cm-1) th_check  a_check    tau  status\n')
        for r in sorted(rows,key=lambda r:(r[1],r[0])):
            f.write('%7.2f %7.2f %18.10f %10.2f %8.2f %8.2f %6.1f  %s\n'%(r[0],r[1],r[2],(r[2]-Emin)*HART,r[3],r[4],r[5],r[6]))
    nok=sum(r[6]=='OK' for r in rows); nres=sum(r[6]=='RESCUED' for r in rows)
    print(f'{nok} converged, {nres} rescued, {len(rows)-nok-nres} need restart  -> {OUT[0]}, {OUT[1]}')
    good=[r for r in rows if r[6] in ('OK','RESCUED')]
    if good: print('max drift of frozen coords: theta %.3f deg, alpha %.3f deg (should be ~0)'%(max(abs(r[3]-r[0]) for r in good if r[0]<179.5),max(abs(r[4]-r[1]) for r in good if r[0]<179.5)))
    for fn,why in bad: print('  ',os.path.basename(fn),':',why)
    for al in (0.0,105.0,180.0):
        sel=[r for r in good if r[1]==al]
        if sel: m=min(sel,key=lambda r:r[2]); print('  alpha=%5.0f: lowest point theta=%.0f  dE=%.1f cm-1'%(al,m[0],(m[2]-Emin)*HART))
