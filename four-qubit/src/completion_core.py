"""Exact monomial arithmetic and binary lemmas used by the compact checks."""
from pathlib import Path
P=Path(__file__).resolve().parents[1]
def mul(a,b):
    pa,fa=a;pb,fb=b
    return tuple(pa[pb[x]] for x in range(8)),tuple((fb[x]+fa[pb[x]])%8 for x in range(8))
def inv(a):
    p,f=a;q=[0]*8;g=[0]*8
    for x in range(8):q[p[x]]=x;g[p[x]]=-f[x]%8
    return tuple(q),tuple(g)
def pauli(a,b):
    return tuple(x^a for x in range(8)),tuple((2*(a&b).bit_count()+4*(b&x).bit_count())%8 for x in range(8))
ID=pauli(0,0)
def conj(u,v):return mul(mul(u,v),inv(u))
def label(u):
    perm,phase=u;a=perm[0]
    if any(perm[x]!=(x^a) for x in range(8)):return None
    b=0
    for j in range(3):
        d=(phase[1<<j]-phase[0])%8
        if d not in (0,4):return None
        b|=(d//4)<<j
    if any((phase[x]-phase[0])%8 != 4*((b&x).bit_count()%2) for x in range(8)):return None
    return a+8*b
def symp(u):
    c=[label(conj(u,pauli(1<<j,0) if j<3 else pauli(0,1<<(j-3)))) for j in range(6)]
    return None if None in c else tuple(c)
def core(toffoli=0,cs=0,cubic=0,mark=0):
    perm=[];phase=[]
    for n in range(8):
        x,y,z=(n>>j&1 for j in range(3))
        perm.append(n^(4*x*y if toffoli else 0))
        phase.append((2*(n&cs==cs if cs else 0)+4*cubic*x*y*z+(n&mark).bit_count()%2)%8)
    return tuple(perm),tuple(phase)
def group(u):
    return sorted({s for a in range(8) for b in range(8) if (s:=symp(conj(u,pauli(a,b)))) is not None})

