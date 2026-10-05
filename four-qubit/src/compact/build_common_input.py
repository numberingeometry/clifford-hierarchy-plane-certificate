"""Generate exact common-input cases; no classification is inferred from counts."""
if not __debug__: raise RuntimeError("Run without Python -O: assertions are certificate checks")
from pathlib import Path
import itertools,json,struct,time
import compact_core as q
t=time.perf_counter();out=Path(__file__).resolve().parents[2]/'evidence/compact'
alist=[]
for sign in (1,3):
 alist.extend([(f'odd_aff_{sign}',q.mono(0,0,0,0,sign,4),4),
 (f'odd_tan_{sign}',q.mono(1,0,0,0,sign,1),1),
 (f'odd_trans_{sign}',q.mono(1,0,0,0,sign,4),4)])
for eta in (0,1):
 alist.extend([(f'dim5_{mark}_{cross}_{eta}',q.mono(1,0,0,cross,2*eta,mark),mark)
 for mark,cross in ((1,1),(2,1),(4,0),(4,1))])
for sign in (1,3):
 alist.extend([(f'unmarked_{perm}_{sign}',q.mono(perm,0,0,0,sign,0),0) for perm in (0,1)])
alist.extend([(f'unmarked_dim5_{eta}',q.mono(1,0,0,1,2*eta,0),0) for eta in (0,1)])
assert len(alist)==20
cases=[];excluded=[];cache={};derivs=0
for name,a,mark in alist:
 for perm,x,y,z,d in itertools.product(range(8),range(2),range(2),range(2),range(4)):
  b=(q.POINTS[perm],q.mono(0,x,y,z,d,mark)[1]);bn=f'{perm}{x}{y}{z}{d}'
  ds=[q.mul(q.mul(b,p),q.inv(a)) for p in q.PAULIS]
  dec=[];reason=None
  for v,u in enumerate(ds):
   if u not in cache:
    dec_u=q.decompose(u)
    # No non-C3 matrix may be silently rejected merely by failed decomposition.
    c3=q.is_c3(u)
    if c3 and dec_u is None:raise RuntimeError(('uncovered C3 monomial',name,bn,v,u))
    if not c3 and dec_u is not None:raise RuntimeError('invalid C3 decomposition')
    cache[u]=dec_u;derivs+=1
   if cache[u] is None:reason=v;break
   dec.append(cache[u])
  if reason is not None:excluded.append([name,bn,reason]);continue
  g0i=q.inverse(dec[0][0])
  tests=set();station=set()
  for g,star in dec:
   gg=q.prod(g0i,g)
   tests.update(q.prod(gg,e) for e in q.STAR_E[star])
   station.add((gg,tuple(i for i,j in q.STAR_CHARTS[star])))
  assert q.I in tests
  edges,ret,bal=q.witnesses(a,b)
  edges2=[(q.LMASK[i],tuple(q.app(g0i,v) for v in q.BASES[j])) for i,j in edges]
  ret2=sorted(set((i,q.app(g0i,j)) for i,j in ret))
  bal2=sorted(set((i,q.app(g0i,j)) for i,j in bal))
  tests=sorted(tests);tests.remove(q.I);tests.insert(0,q.I)
  cases.append({'name':name+'_'+bn,'tests':tests,'station':sorted(station),'edges':edges2,'returns':ret2,'balanced':bal2})
 print(name,'cumulative',len(cases),len(excluded),round(time.perf_counter()-t,1),flush=True)
raw=bytearray()
def u32(v):raw.extend(struct.pack('<I',v))
u32(135)
for mask,bs in zip(q.LMASK,q.BASES):raw.extend(struct.pack('<Q',mask));raw.extend(bs)
u32(len(cases))
for c in cases:
 for k in ('tests','station','edges','returns','balanced'):u32(len(c[k]))
 for mat in c['tests']:raw.extend(mat)
 for mat,chartids in c['station']:
  raw.extend(mat);u32(len(chartids))
  for i in chartids:u32(i)
 for mask,bs in c['edges']:raw.extend(struct.pack('<Q',mask));raw.extend(bs)
 for k in ('returns','balanced'):
  for a,b in c[k]:raw.extend((a,b))
(out/'common-input.bin').write_bytes(raw)
(out/'common-input.json').write_text(json.dumps({'cases':cases,'excluded':excluded,'elapsed':time.perf_counter()-t,'decompositions':derivs},indent=2)+'\n',encoding='utf-8')
print('DONE',len(cases),len(excluded),len(raw),time.perf_counter()-t,flush=True)
