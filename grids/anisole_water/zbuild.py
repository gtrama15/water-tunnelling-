#!/usr/bin/env python3
"""Rebuild Cartesian centres from a Gaussian Z-matrix input (read, build) and Kabsch RMSD."""
import numpy as np
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
