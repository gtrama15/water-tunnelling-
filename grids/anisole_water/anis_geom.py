#!/usr/bin/env python3
"""
anisole-water: water-orientation coordinates and Z-matrix (same design as pdfb_geom.py / phenol-water step 3).
  u = O_w -> O7 (the accepting ether O), n = normal of the (O_w, O7, C1) plane, w = n x u; theta, alpha, tau as below.
Atom order (19): C1-C6 ring (C1 ipso, C2 ortho on the water side), O7, C8 (CH3), H9-H13 ring, H14-H16 CH3, Ow17, H, H.
Z-matrix: ... Ow17, X18 (normal), X19 (bisector), H20, H21; frozen R18 A18 D18 R19 A19 D19 D20 D21.
----- original pDFB description follows (indices there: F7 -> O7 here, Ow13 -> Ow17, X14/X15 -> X18/X19) -----

  u = O_w -> F7 (the accepting F), n = normal of the (O_w, F7, C1) plane (direction of dummy X14), w = n x u
  theta = angle between the water bisector and O_w->F7
  alpha = rotation of the bisector about the O_w...F7 axis (alpha = 0: in the (O_w,F7,C1) plane, w side;
          90: along n; 180: in plane, -w side)
  tau   = rotation of the water about its own bisector, from e0 = cos(a) n - sin(a) w
Reference structures (wB97X-D/def2-TZVP):  minimum theta 61, |alpha| 131, tau 63.5 (mirror: 116.5)
  free-H flip TS (60 cm-1) theta 79, alpha 180, tau 90;  bifurcated H-exchange TS (63 cm-1) theta 30, alpha ~0, tau ~0/180
Symmetry: ring mirror (alpha,tau) -> (-alpha,-tau); H exchange tau -> tau+180.

Atom order (15 atoms): C1..C6, F7 (on C1), F8 (on C4), H9 (C2) H10 (C3) H11 (C5) H12 (C6), Ow13, H16, H17
Z-matrix centres: ... Ow13, X14 (on Ow, perpendicular to the Ow-F7-C1 plane), X15 (on Ow, along the bisector),
                  H16 (angle to X15, dihedral about Ow-X15 = tau-type), H17 (same angle, dihedral 180 from H16)
Frozen: R14 A14 D14 R15 A15 D15 D16 D17  ->  theta, alpha, tau fixed; O_w position, H-bond geometry, water and ring relax.
Where the bisector is nearly along n (theta~90, alpha~90) the H16 dihedral is referred to F7 instead of X14.
"""
import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))

def bond(a, b): return np.linalg.norm(a - b)
def ang(a, b, c):
    u, v = a - b, c - b; return np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
def dih(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c; b1n = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1n) * b1n; w = b2 - (b2 @ b1n) * b1n
    return np.degrees(np.arctan2(np.cross(b1n, v) @ w, v @ w))

def read_ref(fn=os.path.join(HERE, 'ref_min.xyz')):
    L = open(fn).read().split('\n'); n = int(L[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([l.split()[1:4] for l in L[2:2 + n]], float)
LAB, REF = read_ref()

def frame(X):
    """X: 15 real atoms. Returns u, n, w."""
    C1, F7, Ow = X[0], X[6], X[16]
    u = (F7 - Ow) / bond(F7, Ow)
    n = np.cross(F7 - Ow, C1 - Ow); n /= np.linalg.norm(n)
    return u, n, np.cross(n, u)

def angles(X):
    """(theta, alpha, tau) in degrees from 15 real atoms; tau in [0,180)."""
    u, n, w = frame(X); Ow = X[16]; h1, h2 = X[17] - Ow, X[18] - Ow
    b = h1 / np.linalg.norm(h1) + h2 / np.linalg.norm(h2); b /= np.linalg.norm(b)
    th = np.degrees(np.arccos(np.clip(b @ u, -1, 1))); al = np.degrees(np.arctan2(b @ n, b @ w)); a = np.radians(al)
    e0 = np.cos(a) * n - np.sin(a) * w; e = h1 - h2; e -= (e @ b) * b; e /= np.linalg.norm(e)
    return th, al, np.degrees(np.arctan2(e @ np.cross(b, e0), e @ e0)) % 180

def geometry(theta, alpha, tau, dmin=1.9):
    """17 centres: the reference minimum with the water re-oriented (rigid water, O_w fixed) + dummies X14, X15."""
    X = REF.copy(); Ow = X[16]; u, n, w = frame(REF)
    r = np.mean([bond(X[17], Ow), bond(X[18], Ow)]); g = np.radians(ang(X[17], Ow, X[18])) / 2
    t, a, tt = np.radians(theta), np.radians(alpha), np.radians(tau)
    b = np.cos(t) * u + np.sin(t) * (np.cos(a) * w + np.sin(a) * n)
    e0 = np.cos(a) * n - np.sin(a) * w; e = np.cos(tt) * e0 + np.sin(tt) * np.cross(b, e0)
    h1, h2 = r * (np.cos(g) * b + np.sin(g) * e), r * (np.cos(g) * b - np.sin(g) * e)
    # start geometry only: if a water H comes closer than dmin to the pDFB atoms, move O_w back along -u.
    # This keeps u and n (O_w stays on the O_w-F7 line, in the same O_w-F7-C1 plane), so theta, alpha, tau are unchanged;
    # O_w's position is optimised anyway.
    for _ in range(40):
        if min(bond(Ow + h, X[j]) for h in (h1, h2) for j in range(16)) >= dmin: break
        Ow = Ow - 0.05 * u
    X[16] = Ow
    return np.array([*X[:17], Ow + n, Ow + b, Ow + h1, Ow + h2])

def rows(theta, alpha):
    ref20 = 7 if np.sin(np.radians(theta)) * abs(np.sin(np.radians(alpha))) > 0.94 else 18
    return [('C', None, None, None), ('C', 1, None, None), ('C', 2, 1, None), ('C', 3, 2, 1), ('C', 4, 3, 2), ('C', 5, 4, 3),
            ('O', 1, 2, 3), ('C', 7, 1, 2), ('H', 2, 1, 6), ('H', 3, 2, 1), ('H', 4, 3, 2), ('H', 5, 4, 3), ('H', 6, 5, 4),
            ('H', 8, 7, 1), ('H', 8, 7, 14), ('H', 8, 7, 14),
            ('O', 7, 1, 2), ('X', 17, 7, 1), ('X', 17, 7, 18), ('H', 17, 19, ref20), ('H', 17, 19, 20)]
FROZEN = {'R18', 'A18', 'D18', 'R19', 'A19', 'D19', 'D20', 'D21'}

def zmatrix(P, R):
    lines, var, const = [], {}, {}
    for i, (lab, a, b, c) in enumerate(R):
        k = i + 1; s = lab
        if a:
            key = f'R{k}'; s += f' {a} {key}'; (const if key in FROZEN else var)[key] = 1.0 if lab == 'X' else bond(P[i], P[a - 1])
        if b:
            if k == 21: s += f' {b} A20'
            else:
                key = f'A{k}'; s += f' {b} {key}'; v = 90.0 if key == 'A18' else ang(P[i], P[a - 1], P[b - 1])
                (const if key in FROZEN else var)[key] = v
        if c:
            key = f'D{k}'; s += f' {c} {key}'; v = 180.0 if key == 'D21' else dih(P[i], P[a - 1], P[b - 1], P[c - 1])
            (const if key in FROZEN else var)[key] = v
        lines.append(s)
    return lines, var, const

def write_input(fn, title, theta, alpha, tau, route, nproc='4', mem='8GB'):
    P = geometry(theta, alpha, tau); L, var, const = zmatrix(P, rows(theta, alpha))
    base = os.path.splitext(os.path.basename(fn))[0]
    with open(fn, 'w') as f:
        f.write(f'%nprocshared={nproc}\n%mem={mem}\n%chk={base}.chk\n{route}\n\n{title}\n\n0 1\n')
        f.write('\n'.join(L) + '\n Variables:\n' + ''.join(f'  {k:6s} {v:14.8f}\n' for k, v in var.items()))
        f.write(' Constants:\n' + ''.join(f'  {k:6s} {v:14.8f}\n' for k, v in const.items()) + '\n')
