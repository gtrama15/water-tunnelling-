#!/usr/bin/env python3
"""
3D hindered-rotor calculation for the water in phenol-water: full rotation of a rigid water about O_w
(theta, alpha, tau all active) on the step-3 surface, with the exact orientation-dependent G matrix.

Method
  * basis: Wigner functions D^j_{mk}(phi,beta,chi) of the water orientation in the phenol frame, j <= jmax.
    No coordinate singularity anywhere (the collinear geometry theta=180, where alpha and tau become the same
    motion, is an ordinary point).
  * kinetic energy  T = C * sum_ab L_a G_ab(Omega) L_b, L_a = body-fixed angular momentum of the water,
    G = (J'^T J')^-1 with J' the mass-weighted Cartesian velocities for unit rotation about each water axis,
    after projection of overall translation/rotation of the complex (same recipe as the 2D code dvr2b.py;
    phenol and O_w fixed, rigid water at the reference geometry).
  * potential: quadrature grid (Gauss-Legendre in cos beta, uniform in phi and chi), V(theta,alpha,tau) from
      V = V0(theta,alpha)                     2D tau=0 surface, grid_pt.dat (theta 90-179)
        + dV(theta,alpha,tau)                 3D grid (tau = 45, 90, 135), 4-point Fourier in tau (period 180)
        + wall below theta = 90 (no data; K (90-theta)^2)            [+ optional MP2 correction]
  * H-exchange (H1<->H2) = rotation of the water by 180 deg about its bisector = (-1)^k.
    Even k: exchange-symmetric states (para for H2O ... spin weight 1); odd k: antisymmetric (weight 3).
    The tunnelling splitting is E0(odd k) - E0(even k).
"""
import numpy as np, sys, os, argparse, time
from scipy.interpolate import RectBivariateSpline, griddata
from scipy.sparse.linalg import LinearOperator, eigsh
sys.path.insert(0, '/home/claude/tunnel/step2')
import make_scan2 as M

C = 16.857629            # hbar^2/(2 amu A^2) in cm-1
HC = 219474.6354
E_REF = -383.9305807380  # global minimum of the tau=0 surface (theta130, alpha0)
HERE = os.path.dirname(os.path.abspath(__file__))

# ======================================================================================== potential
class Pot3D:
    def __init__(self, f2d, f3d, wallK=5.0, fill=True, corr=None, verbose=True, tauinterp='pchip'):
        self.wallK, self.corr, self.tauinterp = wallK, corr, tauinterp
        # ---- V0(theta, alpha), tau = 0 (as dvr2b.py)
        d = np.loadtxt(f2d, usecols=(0, 1, 2))
        V = {(a, b): (e - E_REF) * HC for a, b, e in d}
        if (90.0, 90.0) not in V: V[(90.0, 90.0)] = 0.5 * (V[(90.0, 75.0)] + V[(90.0, 105.0)])
        TH = sorted(set(d[:, 0])); AL = sorted(set(d[:, 1]))
        Z = np.array([[V[(t, a)] for a in AL] for t in TH])
        Z = np.vstack([Z, 0.5 * (Z[-1] + Z[-1][::-1])]); THg = np.array(TH + [180.0])
        self.th0min = TH[0]
        self.s0 = RectBivariateSpline(THg, *self._ext(np.array(AL), Z), kx=3, ky=3)
        # ---- dV(theta, alpha, tau_k), tau_k = 45, 90, 135
        rows = []
        for f in f3d:
            for l in open(f):
                s = l.split()
                if l.startswith('#') or len(s) < 9 or s[8] not in ('OK', 'RESCUED'): continue
                rows.append((float(s[0]), float(s[1]), float(s[2]), (float(s[3]) - E_REF) * HC))
        self.TH3 = np.array([110., 120., 130., 140., 150., 160., 170., 179.])
        self.AL3 = np.arange(0, 180.01, 22.5)
        tk = sorted({r[2] for r in rows})                     # tau nodes in the data (uniform, excl. 0)
        self.NT = len(tk) + 1; self.HT = 180.0 / self.NT
        if not np.allclose(tk, self.HT * np.arange(1, self.NT)):
            raise SystemExit(f'tau nodes {tk} are not a uniform grid 0..180')
        self.TK = np.array(tk); nk = len(tk)
        kidx = lambda t: int(np.argmin(abs(self.TK - (t % 180.0))))
        D = np.full((nk, len(self.TH3), len(self.AL3)), np.nan)
        for t, a, tau, e in rows:          # keep lowest if duplicated (restart + original)
            i, j, k = np.argmin(abs(self.TH3 - t)), np.argmin(abs(self.AL3 - a)), kidx(tau)
            v = e - self.v0(t, a)
            if np.isnan(D[k, i, j]) or v < D[k, i, j]: D[k, i, j] = v
        # mirror: dV(th,-a,tau) = dV(th,a,-tau); on a = 0, 180 this gives tau == 180 - tau
        mir = {k: kidx(180.0 - self.TK[k]) for k in range(nk)}
        for j in (0, len(self.AL3) - 1):
            for k in range(nk):
                m = np.isnan(D[k, :, j]); D[k, m, j] = D[mir[k], m, j]
        self.nmiss = int(np.isnan(D).sum()); self.nnode = D.size
        if self.nmiss:
            if not fill: raise SystemExit(f'{self.nmiss} grid nodes missing')
            D = self._fill(D)
        self.D = D
        self.sd = []
        for k in range(nk):
            full = np.hstack([D[mir[k]][:, ::-1][:, :-1], D[k]])          # alpha -180..180
            ALf = np.concatenate([-self.AL3[::-1][:-1], self.AL3])
            ALe = np.concatenate([ALf[:-1] - 360, ALf[:-1], ALf + 360])
            Ze = np.hstack([full[:, :-1], full[:, :-1], full])
            self.sd.append(RectBivariateSpline(self.TH3, ALe, Ze, kx=3, ky=3))
        if verbose:
            print(f'potential: {len(rows)} 3D points read, {self.nnode - self.nmiss}/{self.nnode} '
                  f'(theta,alpha,tau) nodes from data, tau step {self.HT:g} deg, {self.tauinterp} in tau' + (f', {self.nmiss} FILLED BY INTERPOLATION' if self.nmiss else ''))

    @staticmethod
    def _ext(AL, Z):    # alpha 0..180 -> -540..540 by mirror + periodicity
        ALf = np.concatenate([-AL[::-1][:-1], AL]); Zf = np.hstack([Z[:, ::-1][:, :-1], Z])
        return (np.concatenate([ALf[:-1] - 360, ALf[:-1], ALf + 360]),
                np.hstack([Zf[:, :-1], Zf[:, :-1], Zf]))

    def v0(self, th, al):
        return self.s0.ev(np.clip(th, self.th0min, 180.), al)

    def _fill(self, D):
        """provisional: fill missing nodes of each tau slice by linear interpolation in (theta, alpha)"""
        T, A = np.meshgrid(self.TH3, self.AL3, indexing='ij')
        for k in range(len(D)):
            ok = ~np.isnan(D[k])
            pts = np.c_[T[ok], A[ok]]
            lin = griddata(pts, D[k][ok], (T, A), method='linear')
            nea = griddata(pts, D[k][ok], (T, A), method='nearest')
            lin = np.where(np.isnan(lin), nea, lin)
            D[k] = np.where(ok, D[k], lin)
        return D

    def __call__(self, th, al, tau):
        """theta in [0,180], alpha in (-180,180], tau in degrees (any)"""
        th, al, tau = map(np.asarray, (th, al, tau))
        v = self.v0(th, al)
        thc = np.clip(th, self.TH3[0], self.TH3[-1])
        F = [s.ev(thc, al) for s in self.sd]
        Y = np.stack(np.broadcast_arrays(np.zeros_like(F[0]), *F))          # tau = 0, h, 2h, ... (f(0) = 0)
        N = self.NT
        if self.tauinterp == 'fourier':      # N-point trig interpolation (rings next to steep tau walls)
            c = np.fft.fft(Y, axis=0) / N; ph = 2 * np.pi * np.mod(tau, 180.0) / 180.0
            f = np.real(c[0]) + 0 * ph
            for m in range(1, (N + 1) // 2):
                f = f + 2 * (np.real(c[m]) * np.cos(m * ph) - np.imag(c[m]) * np.sin(m * ph))
            if N % 2 == 0: f = f + np.real(c[N // 2]) * np.cos(N // 2 * ph)
            v = v + f
        else:                                # periodic monotone (Fritsch-Carlson) Hermite: no overshoot
            u = np.mod(tau, 180.0) / self.HT; k = np.minimum(np.floor(u).astype(int), N - 1); sx = u - k
            g = lambda j: np.take_along_axis(Y, (np.mod(j, N))[None], 0)[0]
            def slope(j):
                a, b = g(j) - g(j - 1), g(j + 1) - g(j)
                with np.errstate(divide='ignore', invalid='ignore'):
                    return np.where(a * b > 0, 2 * a * b / (a + b), 0.0)
            y0, y1, m0, m1 = g(k), g(k + 1), slope(k), slope(k + 1)
            v = v + ((2 * sx**3 - 3 * sx**2 + 1) * y0 + (sx**3 - 2 * sx**2 + sx) * m0
                     + (-2 * sx**3 + 3 * sx**2) * y1 + (sx**3 - sx**2) * m1)
        v = v + self.wallK * np.clip(self.th0min - th, 0, None) ** 2
        if self.corr is not None: v = v + self.corr(th, al, tau)
        return v

# ======================================================================================== geometry
REF = np.array(M.REF, dtype=float)
U, N_, W = M.frame(REF)
OW = REF[13]
_r = np.mean([np.linalg.norm(REF[14] - OW), np.linalg.norm(REF[15] - OW)])
_g = np.arccos(np.dot(REF[14] - OW, REF[15] - OW) / (np.linalg.norm(REF[14] - OW) * np.linalg.norm(REF[15] - OW))) / 2
LAB = ['C'] * 6 + ['O'] + ['H'] * 6 + ['O', 'H', 'H']
MASS0 = {'C': 12.0, 'O': 15.994915, 'H': 1.007825}
YL = np.cross(U, W)                                       # lab frame: X = w, Y = u x w, Z = u
LABF = np.array([W, YL, U])                               # rows: lab axes in Cartesian

def rot(phi, beta, chi):
    """R = Rz(phi) Ry(beta) Rz(chi), arrays -> (...,3,3)"""
    def Rz(a):
        c, s = np.cos(a), np.sin(a); o, z = np.ones_like(a), np.zeros_like(a)
        return np.stack([np.stack([c, -s, z], -1), np.stack([s, c, z], -1), np.stack([z, z, o], -1)], -2)
    def Ry(a):
        c, s = np.cos(a), np.sin(a); o, z = np.ones_like(a), np.zeros_like(a)
        return np.stack([np.stack([c, z, s], -1), np.stack([z, o, z], -1), np.stack([-s, z, c], -1)], -2)
    return Rz(phi) @ Ry(beta) @ Rz(chi)

def body_axes(R):
    """water axes in Cartesian (phenol) frame: e (H2->H1 direction), y = b x e, b (bisector)"""
    A = np.einsum('...ij,ik->...jk', R, LABF)   # row j: body axis j in Cartesian
    return A[..., 0, :], A[..., 1, :], A[..., 2, :]

def internal(e, b):
    """(theta, alpha, tau) in degrees, as in make_scan2.geometry"""
    th = np.degrees(np.arccos(np.clip(b @ U, -1, 1)))
    al = np.degrees(np.arctan2(b @ N_, b @ W))
    a = np.radians(al)
    e0 = np.cos(a)[..., None] * N_ - np.sin(a)[..., None] * W
    tau = np.degrees(np.arctan2(np.einsum('...i,...i', e, np.cross(b, e0)), np.einsum('...i,...i', e, e0)))
    return th, al, tau

def Gbody(e, y, b, mH):
    """3x3 G (amu^-1 A^-2) for rotations about the water body axes (x=e, y, z=b), projected"""
    m = np.array([MASS0[l] for l in LAB]); m[-2:] = mH; sm = np.sqrt(m)
    n = e.shape[0]
    X = np.broadcast_to(REF[:16], (n, 16, 3)).copy()
    h1 = _r * (np.cos(_g) * b + np.sin(_g) * e); h2 = _r * (np.cos(_g) * b - np.sin(_g) * e)
    X[:, 14] = OW + h1; X[:, 15] = OW + h2
    J = np.zeros((n, 3, 16, 3))
    for a, ax in enumerate((e, y, b)):
        J[:, a, 14] = np.cross(ax, h1) * sm[14]; J[:, a, 15] = np.cross(ax, h2) * sm[15]
    J = J.reshape(n, 3, 48)
    com = np.einsum('i,nik->nk', m, X) / m.sum(); Rr = X - com[:, None]
    P = []
    for i in range(3):
        v = np.zeros((n, 16, 3)); v[:, :, i] = sm[None, :]; P.append(v.reshape(n, 48))
        w = np.zeros(3); w[i] = 1
        P.append((np.cross(w, Rr) * sm[None, :, None]).reshape(n, 48))
    Q, _ = np.linalg.qr(np.stack(P, -1))                   # (n,48,6)
    Mm = np.einsum('nai,nbi->nab', J, J) - np.einsum('nai,nik,nbj,njk->nab', J, Q, J, Q)
    return np.linalg.inv(Mm)

# ======================================================================================== rotor basis
def jmats(j):
    k = np.arange(-j, j + 1, dtype=float)
    jp = np.diag(np.sqrt(j * (j + 1) - k[:-1] * (k[:-1] + 1)), -1)   # <k+1|J+|k>
    jx = (jp + jp.T) / 2; jy = (jp - jp.T) / (2j)
    return jx.astype(complex), jy, np.diag(k).astype(complex)

class Rotor:
    def __init__(self, jmax, nb_extra=4):
        self.jmax = J = jmax; L = 2 * J + 1
        self.Nf = 2 * J + 2; self.Nb = J + 1 + nb_extra
        x, wb = np.polynomial.legendre.leggauss(self.Nb)
        self.beta = np.arccos(x); self.wb = wb
        self.phi = 2 * np.pi * np.arange(self.Nf) / self.Nf
        # d^j_{mk}(beta) padded to (J+1, Nb, L, L), with normalisation sqrt((2j+1)/8pi^2)
        self.d = np.zeros((J + 1, self.Nb, L, L))
        self.Jm = np.zeros((3, J + 1, L, L), complex)
        self.mask = np.zeros((J + 1, L, L), bool)
        for j in range(J + 1):
            jx, jy, jz = jmats(j); s = slice(J - j, J + j + 1)
            lam, vec = np.linalg.eigh(jy)
            for ib, bt in enumerate(self.beta):
                dj = (vec * np.exp(-1j * bt * lam)) @ vec.conj().T
                self.d[j, ib, s, s] = dj.real * np.sqrt((2 * j + 1) / (8 * np.pi ** 2))
            for a, mat in enumerate((jx, jy, jz)): self.Jm[a, j, s, s] = mat
            self.mask[j, s, s] = True
        self.kk = np.arange(-J, J + 1)
        self.midx = np.arange(-J, J + 1) % self.Nf        # m -> FFT index
        self.w3 = (2 * np.pi / self.Nf) ** 2

    def grid(self):
        P, B, X = np.meshgrid(self.phi, self.beta, self.phi, indexing='ij')
        return P, B, X                                     # (Nf, Nb, Nf): phi, beta, chi

    def fwd(self, Cf):           # coefficients (J+1,L,L) -> grid (Nf, Nb, Nf)
        A = np.einsum('jbmk,jmk->bmk', self.d, Cf)
        G = np.zeros((self.Nb, self.Nf, self.Nf), complex)
        G[:, self.midx[:, None], self.midx[None, :]] = A
        G = np.fft.fft(np.fft.fft(G, axis=1), axis=2)       # sum_m e^{-i m phi} sum_k e^{-i k chi}
        return G.transpose(1, 0, 2)

    def adj(self, g):            # grid -> coefficients (weighted adjoint of fwd)
        g = g.transpose(1, 0, 2) * self.wb[:, None, None]
        G = np.fft.ifft(np.fft.ifft(g, axis=1), axis=2) * self.Nf ** 2 * self.w3
        A = G[:, self.midx[:, None], self.midx[None, :]]
        return np.einsum('jbmk,bmk->jmk', self.d, A)

    def L(self, a, Cf):          # body-fixed angular momentum on the k index
        return np.einsum('jmk,jlk->jml', Cf, self.Jm[a])

# ======================================================================================== Hamiltonian
class Ham:
    def __init__(self, pot, jmax, mH, verbose=True):
        t0 = time.time()
        self.r = r = Rotor(jmax)
        P, B, X = r.grid(); shp = P.shape
        R = rot(P.ravel(), B.ravel(), X.ravel())
        e, y, b = body_axes(R)
        self.th, self.al, self.tau = internal(e, b)
        self.V = pot(self.th, self.al, self.tau).reshape(shp) if pot is not None else np.zeros(shp)
        G = np.concatenate([Gbody(e[i:i + 20000], y[i:i + 20000], b[i:i + 20000], mH)
                            for i in range(0, len(e), 20000)])
        self.G = G.reshape(shp + (3, 3))
        self.nmv = 0
        if verbose:
            n = r.mask.sum()
            print(f'  jmax={jmax}: {n} basis functions, grid {shp} = {np.prod(shp)} points, '
                  f'setup {time.time() - t0:.1f} s')

    def apply(self, Cf):
        r = self.r; self.nmv += 1
        out = r.adj(self.V * r.fwd(Cf))
        psi = [r.fwd(r.L(b, Cf)) for b in range(3)]
        for a in range(3):
            ga = sum(self.G[..., a, b] * psi[b] for b in range(3))
            out = out + C * r.L(a, r.adj(ga))
        return out

    def solve(self, parity, nev=3, tol=1e-11):
        r = self.r
        sel = r.mask & ((r.kk[None, None, :] % 2 == 0) if parity == 'even' else (r.kk[None, None, :] % 2 == 1)) \
            if parity != 'all' else r.mask
        n = sel.sum()
        def mv(v):                       # real-symmetric form [[Re,-Im],[Im,Re]] (each level twice)
            Cf = np.zeros(sel.shape, complex); Cf[sel] = v[:n] + 1j * v[n:]
            w = self.apply(Cf)[sel]
            return np.concatenate([w.real, w.imag])
        op = LinearOperator((2 * n, 2 * n), matvec=mv, dtype=float)
        v0 = np.random.default_rng(1).standard_normal(2 * n)
        e, v = eigsh(op, k=2 * nev, which='SA', tol=tol, v0=v0, ncv=max(6 * nev, 40), maxiter=20000)
        e, v = e[::2], v[:, ::2]          # degenerate pairs -> one copy
        o = np.argsort(e)
        return e[o], v[:, o], sel

# ======================================================================================== run
def run(pot, mH, jmax, label, nev=3):
    H = Ham(pot, jmax, mH)
    print(f'  V on grid: min {H.V.min():.1f} cm-1 (data minimum = 0)')
    res = {}
    for par in ('even', 'odd'):
        t0 = time.time(); H.nmv = 0
        e, v, sel = H.solve(par, nev)
        res[par] = e
        print(f'  {par:4s} k ({"exchange-symmetric" if par=="even" else "antisymmetric"}): '
              + '  '.join(f'{x:10.4f}' for x in e) + f'   [{H.nmv} H*v, {time.time()-t0:.0f} s]')
    print(f'{label}: ZPE = {res["even"][0]:.2f} cm-1   tunnelling splitting = {res["odd"][0]-res["even"][0]:.5f} cm-1')
    return res

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid2d', default=os.path.join(HERE, 'grid_pt.dat'))
    ap.add_argument('--grid3d', nargs='+', default=[os.path.join(HERE, '../res3d/grid_3d_part1.dat'),
                                                    os.path.join(HERE, '../res3d/grid_3d_restart.dat')])
    ap.add_argument('--jmax', type=int, nargs='+', default=[20])
    ap.add_argument('--wall', type=float, default=5.0)
    ap.add_argument('--iso', nargs='+', default=['H2O', 'D2O'])
    ap.add_argument('--nofill', action='store_true')
    ap.add_argument('--tauinterp', default='pchip', choices=['pchip', 'fourier'])
    ap.add_argument('--mp2', nargs='?', const=os.path.join(HERE, '../mp2corr/mp2corr.dat'),
                    help='add the fitted CP-MP2 - wB97X-D correction (fit_mp2corr.py)')
    a = ap.parse_args()
    corr = None
    if a.mp2:
        from fit_mp2corr import Corr
        corr = Corr(a.mp2)
    pot = Pot3D(a.grid2d, a.grid3d, wallK=a.wall, fill=not a.nofill, corr=corr, tauinterp=a.tauinterp)
    for iso in a.iso:
        mH = {'H2O': 1.007825, 'D2O': 2.014102}[iso]
        for J in a.jmax:
            run(pot, mH, J, f'{iso} jmax={J}')
