#!/usr/bin/env python3
"""
Fit the CP-MP2/aug-cc-pVTZ - wB97X-D correction surface dV(theta, alpha, tau) from the single points in
mp2corr.dat, and provide it as a callable for dvr3d.py (--mp2).

  dV = sum_{l=0..2} (a_l + b_l x + c_l x^2) cos(l alpha)  +  (d0 + d1 x + d2 cos alpha) (1 - cos 2tau) / 2
  x  = (theta - 140) / 40,  theta clamped to [110, 179] (the range of the data)

Energies are relative to the minimum (tf_th130_a000) at both levels, so dV = 0 there.
(theta, alpha, tau) of each point are recomputed from its own coordinates.
"""
import numpy as np, re, os, sys, glob
sys.path.insert(0, '/home/claude/tunnel/step2'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_scan2 as M
HC = 219474.6354
HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(HERE, '../mp2corr/inputs')
GPT = os.path.join(HERE, 'grid_pt.dat')
GBIS = '/home/claude/pw_results/relaxed/grid_relaxed.dat'

_u, _n, _w = M.frame(np.array(M.REF, float)); SW = np.sign(_w @ np.cross(_n, _u))

def coords(name):
    X = []
    for l in open(os.path.join(INP, f'mp2cp_{name}.gjf')):
        s = l.split()
        if len(s) == 4 and '(Fragment' in s[0]: X.append([float(v) for v in s[1:]])
    return np.array(X)

def angles(X):
    u, n, _ = M.frame(X); Ow = X[13]
    w = SW * np.cross(n, u)                     # sign of w fixed by the reference minimum, not by X
    h1, h2 = X[14] - Ow, X[15] - Ow
    b = h1 / np.linalg.norm(h1) + h2 / np.linalg.norm(h2); b /= np.linalg.norm(b)
    e = h1 - h2; e -= (e @ b) * b; e /= np.linalg.norm(e)
    th = np.degrees(np.arccos(b @ u)); al = np.degrees(np.arctan2(b @ n, b @ w)); a = np.radians(al)
    e0 = np.cos(a) * n - np.sin(a) * w
    tau = np.degrees(np.arctan2(e @ np.cross(b, e0), e @ e0))
    return th, al, tau

def edft():
    E = {}
    for l in open(GPT):
        if l.startswith('#'): continue
        s = l.split(); E[f'tf_th{float(s[0]):03.0f}_a{float(s[1]):03.0f}'] = float(s[2])
    for l in open(GBIS):
        if l.startswith('#'): continue
        s = l.split(); E[f'bis_r{float(s[0]):03.0f}_t{float(s[1]) + 90:02.0f}'] = float(s[2])
    return E

MODEL = os.environ.get('MP2MODEL', 'ext')     # 'ext' (16 parameters) or 'base' (12)

def basis(th, al, tau):
    th, al, tau = map(np.atleast_1d, (th, al, tau))
    x = (np.clip(th, 110, 179) - 140) / 40; a = np.radians(al); t = (1 - np.cos(2 * np.radians(tau))) / 2
    cols = [x ** p * np.cos(l * a) for l in range(3) for p in range(3)]
    t2 = (1 - np.cos(4 * np.radians(tau))) / 2
    cols += [t, t * x, t * np.cos(a)]
    if MODEL == 'ext':
        cols += [t * x ** 2, t * x * np.cos(a), t2, t2 * np.cos(a)]
    return np.stack(cols, -1)

class Corr:
    def __init__(self, f=os.path.join(HERE, '../mp2corr/mp2corr.dat'), verbose=True):
        Em = {l.split()[0]: float(l.split()[1]) for l in open(f) if l.strip()}
        Ed = edft(); ref = 'tf_th130_a000'
        rows = []
        for nme, e in Em.items():
            if nme not in Ed: print('  no wB97X-D energy for', nme); continue
            th, al, tau = angles(coords(nme))
            d = ((e - Em[ref]) - (Ed[nme] - Ed[ref])) * HC
            rows.append((nme, th, al, tau, (e - Em[ref]) * HC, (Ed[nme] - Ed[ref]) * HC, d))
        self.rows = rows
        A = basis(*np.array([r[1:4] for r in rows]).T); y = np.array([r[6] for r in rows])
        has_tau = np.any(np.abs(np.sin(np.radians(2 * A[:, 0] * 0 + 2 * np.array([r[3] for r in rows])))) > 0.1)
        self.ncol = A.shape[1] if has_tau else 9
        ncol = min(self.ncol, len(rows) - 1)
        self.c = np.zeros(A.shape[1])
        self.c[:ncol] = np.linalg.lstsq(A[:, :ncol], y, rcond=1e-3)[0]   # drop directions the data cannot fix
        res = A @ self.c - y
        if verbose:
            print(f'MP2 correction fit: {len(rows)} points, {ncol} parameters, rms residual {np.sqrt(np.mean(res**2)):.1f} cm-1')
            print(f'  {"point":15s} {"theta":>6s} {"alpha":>6s} {"tau":>6s} {"MP2":>8s} {"wB97XD":>8s} {"diff":>7s} {"fit":>7s}')
            tg, ag, ug = np.meshgrid(np.linspace(110, 179, 24), np.linspace(-180, 180, 37), np.linspace(0, 180, 13))
            cg = self(tg, ag, ug)
            print(f'  correction over theta 110-179, all alpha, tau: {cg.min():.0f} to {cg.max():.0f} cm-1')
            for r, f_ in zip(rows, A @ self.c):
                print(f'  {r[0]:15s} {r[1]:6.1f} {r[2]:6.1f} {r[3]:6.1f} {r[4]:8.1f} {r[5]:8.1f} {r[6]:7.1f} {f_:7.1f}')

    def __call__(self, th, al, tau):
        sh = np.shape(th)
        return (basis(th, al, tau) @ self.c).reshape(sh)

if __name__ == '__main__':
    Corr(*(sys.argv[1:2]))
