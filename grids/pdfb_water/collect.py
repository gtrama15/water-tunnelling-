#!/usr/bin/env python3
"""Collect the pDFB-water 3D-grid logs:  python3 collect.py [dirs ...]   (default: inputs inputs_restart)
Writes grid_pdfb3d.dat (theta alpha tau E(hartree) dE(cm-1) theta/alpha/tau as found  status) and grid_pdfb3d.xyz.
Unconverged runs are accepted as RESCUED only if their energy is flat (< 1 cm-1 over the last 20 steps)."""
import numpy as np, glob, re, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdfb_geom as G
HART = 219474.6354
Z2S = {1: 'H', 6: 'C', 8: 'O', 9: 'F'}
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


if __name__ == '__main__':
    dirs = sys.argv[1:] or ['inputs', 'inputs_restart']
    files = {}
    for d in dirs:
        for f in glob.glob(os.path.join(d, 'pd_th*_a*_t*.log')): files.setdefault(os.path.basename(f), []).append(f)
    rows, bad, dr = [], [], []
    fx = open('grid_pdfb3d.xyz', 'w')
    for base in sorted(files):
        th, al, ta = [float(x) for x in re.findall(r'th(\d+)_a(\d+)_t(\d+)', base)[0]]; al /= 10; ta /= 10
        best = None
        for f in files[base]:
            q = parse(f)
            if q is None: continue
            st = 'OK' if q['ok'] else None
            if not q['ok']:
                k, still, prc = rescue(q, 1.0)
                if k is not None: s = q['steps'][k]; q.update(E=s['E'], X=s['X'], lab=s['lab']); st = 'RESCUED'
            if st and (best is None or q['E'] < best[0]['E']): best = (q, st)
        if best is None: bad.append(base[:-4]); continue
        q, st = best; tc, ac, tt = G.angles(q['X'])
        d = [abs(tc - th), (0 if th < 5 else abs((ac - al + 180) % 360 - 180)), abs((tt - ta + 90) % 180 - 90)]; dr.append(d)
        rows.append((th, al, ta, q['E'], tc, ac, tt, st))
        fx.write('15\ntheta=%.1f alpha=%.1f tau=%.1f E=%.10f\n' % (th, al, ta, q['E']))
        for l, x in zip(q['lab'], q['X']): fx.write('%-2s %14.8f %14.8f %14.8f\n' % (l, *x))
    fx.close()
    if not rows: sys.exit('no converged points yet')
    Emin = min(r[3] for r in rows)
    with open('grid_pdfb3d.dat', 'w') as f:
        f.write('#  theta   alpha     tau        E(hartree)       dE(cm-1)  th_chk  a_chk  tau_chk status\n')
        for r in sorted(rows, key=lambda r: (r[2], r[1], r[0])):
            f.write('%7.1f %7.1f %7.1f %18.10f %10.2f %7.2f %7.2f %7.2f  %s\n' % (r[0], r[1], r[2], r[3], (r[3] - Emin) * HART, *r[4:]))
    print(f'{len(rows)} usable points -> grid_pdfb3d.dat, grid_pdfb3d.xyz ;  {len(bad)} still to (re)run')
    print('max drift theta %.3f alpha %.3f tau %.3f deg' % tuple(np.max(np.array(dr), axis=0)))
    if bad: print('missing/unconverged:', ' '.join(bad))
