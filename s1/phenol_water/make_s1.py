import re,os
out='tunnel_s1'
ROUTE='#p TD(NStates=4) wB97XD/Def2TZVP SCF=Tight NoSymm Int=UltraFine'
def frames(f):
    L=open(f).read().split('\n'); i=0
    while i<len(L) and L[i].strip():
        n=int(L[i]); t=L[i+1]; yield t,L[i+2:i+2+n]; i+=2+n
names=[]
for f,kind in (('../dvr/geoms_pt.xyz','pt'),('../res3d/final/geoms_3d.xyz','3d')):
    for t,at in frames(f):
        th=float(re.search(r'theta=([\d.]+)',t).group(1)); al=float(re.search(r'alpha=([\d.]+)',t).group(1))
        if kind=='pt': nm=f's1_pt_th{th:03.0f}_a{al:03.0f}'; tau=0.0
        else:
            tau=float(re.search(r'tau=([\d.]+)',t).group(1)); nm=f's1_3d_th{th:03.0f}_a{al*10:04.0f}_t{tau:03.0f}'
        with open(f'{out}/inputs/{nm}.gjf','w') as g:
            g.write(f'%nprocshared=8\n%mem=8GB\n%chk={nm}.chk\n{ROUTE}\n\n'
                    f'S1 vertical theta={th} alpha={al} tau={tau} (S0 wB97XD relaxed geometry)\n\n0 1\n'+'\n'.join(at)+'\n\n')
        names.append(nm)
test='s1_pt_th130_a000'
names.remove(test); names=[test]+names
half=(len(names)+1)//2
open(f'{out}/joblist_cluster.txt','w').write('\n'.join(names[:half])+'\n')
open(f'{out}/joblist_karplus.txt','w').write('\n'.join(names[half:])+'\n')
print(len(names),'inputs;',half,'cluster,',len(names)-half,'Karplus')
