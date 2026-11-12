#!/usr/bin/env python3
"""Rebuild every z-matrix and check it reproduces the intended geometry and (theta, alpha)."""
import numpy as np, glob, re, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from make_scan2 import geometry, frame, REF
def place(a,b,c,r,th,ph):
    th,ph=np.radians(th),np.radians(ph); bc=a-b; bc/=np.linalg.norm(bc)
    n=np.cross(b-c,bc); n/=np.linalg.norm(n); m=np.cross(n,bc)
    d=np.array([-r*np.cos(th),r*np.sin(th)*np.cos(ph),r*np.sin(th)*np.sin(ph)])
    return a+d[0]*bc+d[1]*m+d[2]*n
def read(fn):
    t=open(fn).read().split('\n'); i=t.index('0 1'); z=[];k=i+1
    while not t[k].strip().startswith('Variables'): z.append(t[k].split()); k+=1
    v={}
    for l in t[k+1:]:
        p=l.split()
        if len(p)==2: v[p[0]]=float(p[1])
    return z,v
def build(z,v):
    P=[]
    for row in z:
        if len(row)==1: P.append(np.zeros(3))
        elif len(row)==3: P.append(np.array([v[row[2]],0,0]))
        elif len(row)==5:
            a,b=int(row[1])-1,int(row[3])-1; r,th=v[row[2]],np.radians(v[row[4]])
            d=P[b]-P[a]; d/=np.linalg.norm(d); P.append(P[a]+r*np.array([d[0]*np.cos(th)-d[1]*np.sin(th),d[0]*np.sin(th)+d[1]*np.cos(th),0]))
        else:
            a,b,c=int(row[1])-1,int(row[3])-1,int(row[5])-1
            P.append(place(P[a],P[b],P[c],v[row[2]],v[row[4]],v[row[6]]))
    return np.array(P)
def kabsch(A,B):
    A=A-A.mean(0);B=B-B.mean(0);U,S,Vt=np.linalg.svd(A.T@B);d=np.sign(np.linalg.det(U@Vt))
    R=U@np.diag([1,1,d])@Vt; return np.sqrt(((A@R-B)**2).sum(1).mean())
def active(X):
    """theta, alpha from a 16-atom geometry (no dummies)."""
    Ow=X[13]; u,n,w=frame_from(X)
    b=(X[14]+X[15])/2-Ow; b/=np.linalg.norm(b)
    th=np.degrees(np.arccos(np.clip(b@u,-1,1)))
    al=np.degrees(np.arctan2(b@n,b@w))
    return th,al
_u,_n,_w=frame(REF)
def frame_from(X):
    C1,O7,Ow=X[0],X[6],X[13]; u=(O7-Ow)/np.linalg.norm(O7-Ow)
    n=np.cross(O7-Ow,C1-Ow); n/=np.linalg.norm(n); w=np.cross(n,u)
    # same handedness as the reference frame (w = n x u, with the reference sign choice)
    if (np.cross(_n,_u)@_w)<0: w=-w
    return u,n,w
if __name__=='__main__':
    worst=0; worst_a=0
    for fn in sorted(glob.glob('inputs/ps_th*_a*.gjf')):
        th,al=map(float,re.findall(r'th(\d+)_a(\d+)',fn)[0])
        z,v=read(fn); P=build(z,v); ref=np.array(geometry(th,al)); worst=max(worst,kabsch(P,ref))
        X=np.delete(P,[14,15],axis=0); t2,a2=active(X)
        worst_a=max(worst_a,abs(t2-th),abs(((a2-al)+180)%360-180) if th<179.5 else 0)
    print('max RMSD rebuilt vs target = %.1e A ; max error in theta/alpha = %.1e deg ; files %d'%(worst,worst_a,len(glob.glob('inputs/ps_*.gjf'))))
