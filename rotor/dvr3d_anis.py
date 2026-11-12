#!/usr/bin/env python3
"""
3D hindered-rotor calculation for the water in anisole-water, on the relaxed (theta, alpha, tau) grid of
tunnel_pdfb3d (wB97X-D/def2-TZVP). Same method as dvr3d.py (Wigner basis, exact orientation-dependent G, H exchange =
(-1)^k), with the pDFB frame:  u = O_w -> O7 (ether O), n = normal of the (O_w, O7, C1) plane, w = n x u.

Potential: for each tau node (0, 22.5, ..., 157.5) a bicubic spline in (theta, alpha) over the data
(theta 2..120, alpha -180..180 via the ring mirror V(th,-a,tau) = V(th,a,-tau)), then periodic monotone Hermite
(or Fourier) interpolation in tau (period 180: H exchange).  theta < 2 is clamped to 2 (next to the pole);
above theta = 120 (no data, all points there > 300 cm-1) a wall K (theta-120)^2.
"""
import numpy as np, os, sys, argparse
from scipy.interpolate import RectBivariateSpline, griddata
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '../anis/tunnel_anis3d'))
import dvr3d as D
import anis_geom as PG

E_REF = -423.2274270740          # anisole-water wB97X-D minimum (step 1)
HC = D.HC

# ---------------------------------------------------------------- frame, geometry and masses for pDFB-water
REF = PG.REF.copy()                                   # 15 atoms: C1-C6 F7 F8 H9-H12 Ow13 H H
U, N_, W = PG.frame(REF)
OW = REF[-3]
_r = np.mean([np.linalg.norm(REF[-2] - OW), np.linalg.norm(REF[-1] - OW)])
_g = np.arccos(np.dot(REF[-2] - OW, REF[-1] - OW) / (np.linalg.norm(REF[-2] - OW) * np.linalg.norm(REF[-1] - OW))) / 2
LAB = PG.LAB
MASS0 = {'C': 12.0, 'O': 15.994915, 'H': 1.007825, 'F': 18.998403}
YL = np.cross(U, W)
D.U, D.N_, D.W, D.OW, D.LABF = U, N_, W, OW, np.array([W, YL, U])      # used by D.internal / D.body_axes
NA = len(LAB)

def Gbody(e, y, b, mH):
    m = np.array([MASS0[l] for l in LAB]); m[-2:] = mH; sm = np.sqrt(m)
    n = e.shape[0]
    X = np.broadcast_to(REF, (n, NA, 3)).copy()
    h1 = _r * (np.cos(_g) * b + np.sin(_g) * e); h2 = _r * (np.cos(_g) * b - np.sin(_g) * e)
    X[:, NA - 2] = OW + h1; X[:, NA - 1] = OW + h2
    J = np.zeros((n, 3, NA, 3))
    for a, ax in enumerate((e, y, b)):
        J[:, a, NA - 2] = np.cross(ax, h1) * sm[NA - 2]; J[:, a, NA - 1] = np.cross(ax, h2) * sm[NA - 1]
    J = J.reshape(n, 3, 3 * NA)
    com = np.einsum('i,nik->nk', m, X) / m.sum(); Rr = X - com[:, None]
    P = []
    for i in range(3):
        v = np.zeros((n, NA, 3)); v[:, :, i] = sm[None, :]; P.append(v.reshape(n, 3 * NA))
        w = np.zeros(3); w[i] = 1
        P.append((np.cross(w, Rr) * sm[None, :, None]).reshape(n, 3 * NA))
    Q, _ = np.linalg.qr(np.stack(P, -1))
    Mm = np.einsum('nai,nbi->nab', J, J) - np.einsum('nai,nik,nbj,njk->nab', J, Q, J, Q)
    return np.linalg.inv(Mm)
D.Gbody = Gbody

# ---------------------------------------------------------------- potential
class PotPDFB:
    def __init__(self, files, wallK=5.0, fill=True, tauinterp='pchip', corr=None, verbose=True):
        self.wallK, self.tauinterp, self.corr = wallK, tauinterp, corr
        rows = []
        for f in files:
            for l in open(f):
                s = l.split()
                if l.startswith('#') or len(s) < 9 or s[8] not in ('OK', 'RESCUED'): continue
                rows.append((float(s[0]), float(s[1]), float(s[2]), (float(s[3]) - E_REF) * HC))
        self.TH = np.array(sorted({r[0] for r in rows})); self.AL = np.arange(0, 180.01, 22.5)
        self.NT = 8; self.HT = 22.5; TK = self.HT * np.arange(self.NT)
        kidx = lambda t: int(np.argmin(abs(TK - (t % 180.0))))
        V = np.full((self.NT, len(self.TH), len(self.AL)), np.nan)
        for t, a, u, e in rows:
            i, j, k = np.argmin(abs(self.TH - t)), np.argmin(abs(self.AL - a)), kidx(u)
            if np.isnan(V[k, i, j]) or e < V[k, i, j]: V[k, i, j] = e
        mir = {k: kidx(-TK[k]) for k in range(self.NT)}
        for j in (0, len(self.AL) - 1):                      # on alpha = 0, 180: V(tau) = V(-tau)
            for k in range(self.NT):
                m = np.isnan(V[k, :, j]); V[k, m, j] = V[mir[k], m, j]
        self.nmiss = int(np.isnan(V).sum()); self.nnode = V.size
        if self.nmiss:
            if not fill: raise SystemExit(f'{self.nmiss} grid nodes missing')
            T, A = np.meshgrid(self.TH, self.AL, indexing='ij')
            for k in range(self.NT):
                ok = ~np.isnan(V[k]); pts = np.c_[T[ok], A[ok]]
                lin = griddata(pts, V[k][ok], (T, A), method='linear'); nea = griddata(pts, V[k][ok], (T, A), method='nearest')
                V[k] = np.where(ok, V[k], np.where(np.isnan(lin), nea, lin))
        self.V = V
        ALf = np.concatenate([-self.AL[::-1][:-1], self.AL]); ALe = np.concatenate([ALf[:-1] - 360, ALf[:-1], ALf + 360])
        self.sp = []
        for k in range(self.NT):
            full = np.hstack([V[mir[k]][:, ::-1][:, :-1], V[k]])
            self.sp.append(RectBivariateSpline(self.TH, ALe, np.hstack([full[:, :-1], full[:, :-1], full]), kx=3, ky=3))
        if verbose:
            print(f'anisole potential: {len(rows)} points, {self.nnode - self.nmiss}/{self.nnode} (theta,alpha,tau) nodes from data'
                  + (f', {self.nmiss} FILLED BY INTERPOLATION' if self.nmiss else '') + f'; {tauinterp} in tau')

    def __call__(self, th, al, tau):
        th, al, tau = map(np.asarray, (th, al, tau))
        thc = np.clip(th, self.TH[0], self.TH[-1])
        Y = np.stack([s.ev(thc, al) for s in self.sp])
        N = self.NT
        if self.tauinterp == 'fourier':
            c = np.fft.fft(Y, axis=0) / N; ph = 2 * np.pi * np.mod(tau, 180.0) / 180.0
            v = np.real(c[0]) + 0 * ph
            for m in range(1, (N + 1) // 2):
                v = v + 2 * (np.real(c[m]) * np.cos(m * ph) - np.imag(c[m]) * np.sin(m * ph))
            if N % 2 == 0: v = v + np.real(c[N // 2]) * np.cos(N // 2 * ph)
        else:
            u = np.mod(tau, 180.0) / self.HT; k = np.minimum(np.floor(u).astype(int), N - 1); sx = u - k
            g = lambda j: np.take_along_axis(Y, (np.mod(j, N))[None], 0)[0]
            def slope(j):
                a, b = g(j) - g(j - 1), g(j + 1) - g(j)
                with np.errstate(divide='ignore', invalid='ignore'):
                    return np.where(a * b > 0, 2 * a * b / (a + b), 0.0)
            y0, y1, m0, m1 = g(k), g(k + 1), slope(k), slope(k + 1)
            v = ((2 * sx**3 - 3 * sx**2 + 1) * y0 + (sx**3 - 2 * sx**2 + sx) * m0 + (-2 * sx**3 + 3 * sx**2) * y1 + (sx**3 - sx**2) * m1)
        v = v + self.wallK * np.clip(th - self.TH[-1], 0, None) ** 2
        if self.corr is not None: v = v + self.corr(th, al, tau)
        return v

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid', nargs='+', default=[os.path.join(HERE, '../anis/res/grid_anis3d.dat')])
    ap.add_argument('--jmax', type=int, nargs='+', default=[16])
    ap.add_argument('--iso', nargs='+', default=['H2O', 'D2O'])
    ap.add_argument('--wall', type=float, default=5.0)
    ap.add_argument('--tauinterp', default='pchip', choices=['pchip', 'fourier'])
    ap.add_argument('--nofill', action='store_true')
    a = ap.parse_args()
    pot = PotPDFB(a.grid, wallK=a.wall, fill=not a.nofill, tauinterp=a.tauinterp)
    for iso in a.iso:
        mH = {'H2O': 1.007825, 'D2O': 2.014102}[iso]
        for J in a.jmax: D.run(pot, mH, J, f'anisole-{iso} jmax={J}')
