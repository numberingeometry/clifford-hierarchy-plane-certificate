"""Exact small monomials and chart relations for the compact four-qubit checks."""
if not __debug__: raise RuntimeError("Run without Python -O: assertions are certificate checks")
from pathlib import Path
import importlib.util,itertools
P=Path(__file__).resolve().parents[2]
sp=importlib.util.spec_from_file_location('mc',P/'src/completion_core.py')
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
mul,inv,conj,pauli,symp=m.mul,m.inv,m.conj,m.pauli,m.symp
I=tuple(1<<i for i in range(6))
def app(a,v):
 r=0
 for j in range(6):
  if v>>j&1:r^=a[j]
 return r
def prod(a,b):return tuple(app(a,x) for x in b)
def inverse(a):return tuple(next(x for x in range(64) if app(a,x)==1<<j) for j in range(6))
def om(a,b):return (((a&7)&(b>>3)).bit_count()+((b&7)&(a>>3)).bit_count())%2
def span(bs):
 s={0}
 for b in bs:s|={x^b for x in list(s)}
 return frozenset(s)
LAGS=sorted({span((a,b,c)) for a in range(1,64) for b in range(a+1,64) if not om(a,b)
 for c in range(b+1,64) if c!=a^b and not om(a,c) and not om(b,c)},key=lambda x:sorted(x))
assert len(LAGS)==135
def basis(l):
 out=[]
 for v in sorted(l):
  if v not in span(out):out.append(v)
 return tuple(out)
BASES=[basis(l) for l in LAGS]
LMASK=[sum(1<<v for v in l) for l in LAGS]
PAULIS=[pauli(v&7,v>>3) for v in range(64)]
from functools import lru_cache
@lru_cache(None)
def coeffs(u):
 p,f=u;out={}
 for v,(pp,ff) in enumerate(PAULIS):
  a=[0]*4
  for x in range(8):
   if p[x]==pp[x]:
    e=(f[x]-ff[x])%8;a[e%4]+=1 if e<4 else -1
  if any(a):out[v]=tuple(a)
 return out
@lru_cache(None)
def rows(u):return tuple(coeffs(conj(u,p)) for p in PAULIS)
@lru_cache(None)
def charts(u):
 rr=rows(u);masks=[sum(1<<v for v in r) for r in rr];out=[]
 for i,bs in enumerate(BASES):
  s=masks[bs[0]]|masks[bs[1]]|masks[bs[2]]
  choices=[j for j,t in enumerate(LMASK) if s&~t==0]
  assert len(choices)<=1
  if choices:out.append((i,choices[0]))
 return tuple(out)
def mono(t,a,b,c,d,p):
 return tuple(n^(4*(n&1)*((n>>1)&1) if t else 0) for n in range(8)),tuple(
 (2*(a*(n&1)*((n>>1)&1)+b*(n&1)*((n>>2)&1)+c*((n>>1)&1)*((n>>2)&1)+d*(n&1)*((n>>1)&1)*((n>>2)&1))+(n&p).bit_count()%2)%8 for n in range(8))
def affine(p):
 return all(p[n]==p[0] ^ ((p[1]^p[0]) if n&1 else 0) ^ ((p[2]^p[0]) if n&2 else 0) ^ ((p[4]^p[0]) if n&4 else 0) for n in range(8))
def signature(f):
 a=f[3]-f[1]-f[2]+f[0];b=f[5]-f[1]-f[4]+f[0];c=f[6]-f[2]-f[4]+f[0]
 d=f[7]-f[3]-f[5]-f[6]+f[1]+f[2]+f[4]-f[0]
 return a%4,b%4,c%4,d%8
LINEARS=[(p,(0,)*8) for p in itertools.permutations(range(8)) if p[0]==0 and affine(p)]
assert len(LINEARS)==168
DIAG={}
for base in (mono(0,0,0,0,0,0),mono(0,1,0,0,0,0),mono(0,0,0,0,2,0),mono(0,1,0,0,2,0)):
 for c in LINEARS:
  d=conj(c,base);DIAG.setdefault(signature(d[1]),d)
assert len(DIAG)==16
# One linear change for each nonlinear radical direction.
DIR_COORDS=[]
for k in range(1,8):
 g=next(g for g in LINEARS if g[0][k]==4)
 DIR_COORDS.append(g)
POINTS=[mono(0,0,0,0,0,0)[0]]+[conj(inv(g),mono(1,0,0,0,0,0))[0] for g in DIR_COORDS]
assert len(set(POINTS))==8
NONLIN=[conj(inv(g),mono(1,a,0,0,d,0)) for g in DIR_COORDS for a,d in ((0,0),(1,0),(1,2),(2,2))]
STARS=list(dict.fromkeys([*DIAG.values(),*NONLIN]))
STAR_CHARTS={d:charts(d) for d in STARS}
assert all(all(i==j for i,j in cs) for cs in STAR_CHARTS.values())
STAR_E={d:m.group(d) for d in STARS}
assert all(len(g) in (1,4,8) for g in STAR_E.values())
@lru_cache(None)
def decompose(u):
 p,f=u
 if affine(p):
  d=DIAG.get(signature(f))
  if d is None:return None
  g=symp(mul(u,inv(d)))
  return (g,d) if g is not None else None
 for d in NONLIN:
  v=mul(u,inv(d))
  if not affine(v[0]):continue
  g=symp(v)
  if g is not None:return g,d
 return None
def is_c3(u):return all(symp(conj(u,PAULIS[v])) is not None for v in (1,2,4,8,16,32))
def one_label(r):
 if len(r)!=1:return None
 return next(iter(r))
def balanced(r):
 # Numerator 8/sqrt(2)=4sqrt(2)=4(zeta-zeta^3), up to signs.
 return len(r)==2 and all(a in ((0,4,0,-4),(0,-4,0,4)) for a in r.values())
def witnesses(a,b):
 ar,br=rows(a),rows(b);ret=[]
 for i in range(1,64):
  x=one_label(ar[i]);y=one_label(br[i])
  if x is not None and y is not None:ret.append((x,y))
 aa,bb=rows(inv(a)),rows(inv(b))
 bal=[(i,j) for i in range(1,64) if balanced(aa[i]) for j in range(1,64)
  if balanced(bb[j]) and aa[i].keys()==bb[j].keys()]
 ca,cb=dict(charts(a)),dict(charts(b))
 edges=[(ca[l],cb[l]) for l in ca.keys()&cb.keys()]
 return edges,ret,bal
