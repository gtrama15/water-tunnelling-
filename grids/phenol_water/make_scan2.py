#!/usr/bin/env python3
"""Step 2: RELAXED 2D scan of phenol-water in the coordinates of the water's H-exchange path.

  theta = angle between the water bisector and the O_w->O_ph axis (cone angle, deg)
          (theta ~131 at the minimum; 180 = bisector collinear with O...O, pointing away from phenol)
  alpha = rotation of the bisector about the O...O axis (deg)
          alpha =   0 : bisector in the phenol plane, on the side of the minimum
          alpha =  90 : bisector perpendicular to the phenol plane
          alpha = 180 : bisector in the phenol plane, on the other side (the ~110 cm-1 shoulder)
  The H-exchange saddle (228 cm-1, two independent TS searches) lies at theta=133.6, alpha=109.5.

theta and alpha are FROZEN; everything else relaxes, including the rotation of the water about
its own bisector (tau), the O...O distance and all phenol/water internal coordinates.
Mirror symmetry through the phenol plane maps alpha -> -alpha, so alpha = 0..180 covers the full surface.
This is the relaxed version of the old PCCP grid (r1, r2):  theta=r1 for r1<=180, alpha = r2+90.

Z-matrix:  X15 = dummy on O_w, perpendicular to the (O_w, O_ph, C1) plane
           Xb  = dummy on O_w along the bisector:  angle(Xb,Ow,O_ph)=theta, dihedral(Xb,Ow,O_ph,X15)=D16
           D16 = Dmin*(1 - alpha/90), where Dmin (= +-90) is the in-plane direction of the minimum.
           H17: bond, angle to Xb (= half HOH), dihedral about Xb = tau (VARIABLE);  H18: same angle, dihedral 180.
"""
import numpy as np, os, argparse
def bond(a,b): return np.linalg.norm(a-b)
def ang(a,b,c):
    u,v=a-b,c-b; return np.degrees(np.arccos(np.clip(u@v/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))
def dih(a,b,c,d):
    b0,b1,b2=a-b,c-b,d-c; b1n=b1/np.linalg.norm(b1)
    v=b0-(b0@b1n)*b1n; w=b2-(b2@b1n)*b1n
    return np.degrees(np.arctan2(np.cross(b1n,v)@w, v@w))
# order: C1..C6 O7 H8..H12 (ring H) H13 (phenol OH) Ow14 X15 Xb16 H17 H18
ROWS=[('C',None,None,None),('C',1,None,None),('C',2,1,None),('C',3,2,1),('C',4,3,2),('C',5,4,3),
      ('O',1,2,3),('H',3,2,1),('H',4,3,2),('H',5,4,3),('H',6,5,4),('H',2,1,6),('H',7,1,2),
      ('O',7,1,2),('X',14,7,1),('X',14,7,15),('H',14,16,15),('H',14,16,17)]
FROZEN={'R15','A15','D15','R16','A16','D16','D18'}   # + 'D17' with --freeze-tau

def read_ref(fn):
    L=open(fn).read().split('\n'); return np.array([l.split()[1:4] for l in L[2:18]],float)
REF=read_ref(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ref_min.xyz'))

def frame(X):
    C1,O7,Ow=X[0],X[6],X[13]
    u=(O7-Ow)/bond(O7,Ow)
    n=np.cross(O7-Ow,C1-Ow); n/=np.linalg.norm(n)          # = direction of X15
    w=np.cross(n,u)                                        # in-plane, perpendicular to O...O
    b=(X[14]+X[15])/2-Ow
    if b@w<0: w=-w                                         # alpha=0 on the side of the minimum bisector
    return u,n,w

def geometry(theta,alpha,tau=0.0):
    """Start geometry: reference minimum with the water re-oriented (rigid water, O_w fixed).
    tau = rotation of the water about its own bisector, measured from the co-rotating orientation e0
    (tau=0: H-H along e0 = cos(a) n - sin(a) w;  tau=90 with alpha=0: water lies in the phenol plane)."""
    X=REF.copy(); Ow=X[13]; u,n,w=frame(REF)
    r=np.mean([bond(X[14],Ow),bond(X[15],Ow)]); g=np.radians(ang(X[14],Ow,X[15]))/2
    t,a=np.radians(theta),np.radians(alpha)
    b=np.cos(t)*u+np.sin(t)*(np.cos(a)*w+np.sin(a)*n)
    e0=np.cos(a)*n-np.sin(a)*w                             # H-H direction for tau=0 (as in the minimum)
    tt=np.radians(tau); e=np.cos(tt)*e0+np.sin(tt)*np.cross(b,e0)
    X[14]=Ow+r*(np.cos(g)*b+np.sin(g)*e); X[15]=Ow+r*(np.cos(g)*b-np.sin(g)*e)
    X15=Ow+n; Xb=Ow+b
    return [*X[:14],X15,Xb,X[14],X[15]]

def zmatrix(P):
    lines=[];var={};const={}
    for i,(lab,a,b,c) in enumerate(ROWS):
        k=i+1; s=lab
        if a:
            key=f'R{k}'; s+=f' {a} {key}'; (const if key in FROZEN else var)[key]=1.0 if lab=='X' else bond(P[i],P[a-1])
        if b:
            key=f'A{k}'
            if k==18: s+=f' {b} A17'
            else:
                s+=f' {b} {key}'; v=ang(P[i],P[a-1],P[b-1])
                if key=='A15': v=90.0
                (const if key in FROZEN else var)[key]=v
        if c:
            key=f'D{k}'; s+=f' {c} {key}'; v=dih(P[i],P[a-1],P[b-1],P[c-1])
            if key=='D18': v=180.0
            (const if key in FROZEN else var)[key]=v
        lines.append(s)
    return lines,var,const

def write_input(fn,title,P,route,nproc,mem):
    L,var,const=zmatrix(P); base=os.path.splitext(os.path.basename(fn))[0]
    with open(fn,'w') as f:
        f.write(f'%nprocshared={nproc}\n%mem={mem}\n%chk={base}.chk\n{route}\n\n{title}\n\n0 1\n')
        f.write('\n'.join(L)+'\n Variables:\n'+''.join(f'  {k:6s} {v:14.8f}\n' for k,v in var.items()))
        f.write(' Constants:\n'+''.join(f'  {k:6s} {v:14.8f}\n' for k,v in const.items())+'\n')
    return var,const

THETAS=[90,100,110,120,130,140,150,160,170,175,179]
ALPHAS=list(range(0,181,15))
def priority(th,al):
    if al in (0,180): return 1                  # in-plane lines: minimum, inversion, shoulder
    if al in (90,105,120) and 120<=th<=150: return 1   # around the exchange saddle (theta 134, alpha 110)
    if al in (75,135,150): return 2
    return 3

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--route',default='#p wB97XD/def2TZVP Opt=(Z-Matrix,MaxCycles=100) NoSymm Int=UltraFine')
    p.add_argument('--nproc',default='8'); p.add_argument('--mem',default='32GB'); p.add_argument('--out',default='inputs')
    p.add_argument('--freeze-tau',action='store_true',help='also freeze D17: the water co-rotates rigidly with alpha (old PCCP coordinates)')
    p.add_argument('--prefix',default='ps'); p.add_argument('--joblist',default='joblist.txt')
    a=p.parse_args(); os.makedirs(a.out,exist_ok=True); jobs=[]
    if a.freeze_tau: FROZEN.add('D17')
    for th in THETAS:
        for al in ALPHAS:
            if a.freeze_tau and th==90 and al==90: continue      # bisector along dummy X15: tau undefined
            name=f'{a.prefix}_th{th:03d}_a{al:03d}'
            what='tau frozen (co-rotating water)' if a.freeze_tau else 'tau relaxed'
            write_input(f'{a.out}/{name}.gjf',f'phenol-water relaxed scan theta={th} alpha={al}, {what}',geometry(th,al),a.route,a.nproc,a.mem)
            pr=priority(th,al) if not a.freeze_tau else (1 if th>=150 else 2 if al in (90,105,120,135) else 3)
            jobs.append((pr,th,al,name))
    jobs.sort()
    open(a.joblist,'w').write(''.join(f'{n} {p}\n' for p,_,_,n in jobs))
    print(len(jobs),'inputs in',a.out,'; priorities:',{k:sum(j[0]==k for j in jobs) for k in (1,2,3)})
