"""Generate all Clifford conjugates of the Gaussian third-level core tests."""
if not __debug__: raise RuntimeError("Run without Python -O: assertions are certificate checks")
from pathlib import Path
import json,struct,time
import compact_core as q
from functools import lru_cache
t=time.perf_counter();out=Path(__file__).resolve().parents[2]/'evidence/compact'
def lift_lag(bs):
 xs=[]
 for i in range(3):
  xs.append(next(v for v in range(1,64) if all(q.om(v,z)==(i==j) for j,z in enumerate(bs))
   and all(not q.om(v,x) for x in xs)))
 c=tuple(xs)+tuple(bs)
 assert all(q.om(c[i],c[j])==q.om(1<<i,1<<j) for i in range(6) for j in range(6))
 return c
lifts=[lift_lag(bs) for bs in q.BASES]
linears=[q.symp(u) for u in q.LINEARS]
lininv=[q.inverse(c) for c in linears]
bases=[('I',q.mono(0,0,0,0,0,0)),('CS',q.mono(0,1,0,0,0,0)),
 ('Q',q.mono(0,0,0,0,2,0)),('CSQ',q.mono(0,1,0,0,2,0))]
orbits={}
for name,u in bases:
 eu=q.m.group(u);found=set()
 for cl in lifts if name!='I' else [q.I]:
  cli=q.inverse(cl)
  for cg,cgi in zip(linears,lininv) if name!='I' else [(q.I,q.I)]:
   c=q.prod(cl,cg);ci=q.prod(cgi,cli)
   eg=tuple(sorted(q.prod(q.prod(c,e),ci) for e in eu))
   found.add(eg)
 for e in found:orbits.setdefault(e,name)
 print(name,'conjugates',len(found),'total',len(orbits),'seconds',round(time.perf_counter()-t,2),flush=True)
@lru_cache(None)
def lag_action(e):
 lookup={l:i for i,l in enumerate(q.LAGS)}
 return tuple(lookup[frozenset(q.app(e,v) for v in l)] for l in q.LAGS)
def invariant_charts(es):return tuple(i for i in range(135) if all(lag_action(e)[i]==i for e in es))
def fixed(es):return frozenset(v for v in range(1,64) if all(q.app(e,v)==v for e in es))
alist=[(n,u) for n,u in bases]
for name,quad,cubic in (('Tz',0,0),('CSyzTz',6,0),('CSxyTz',3,0),('QTz',0,1),('CSyzQTz',6,1),('CSxyQTz',3,1)):
 alist.append((name,q.m.core(cs=quad,cubic=cubic,mark=4)))
cases=[];excluded=0;defs=[]
for ed,kind in sorted(orbits.items()):
 xd=invariant_charts(ed);rd=fixed(ed);defs.append({'kind':kind,'group':ed,'charts':xd,'fixed':sorted(rd)})
for an,a in alist:
 ea=q.m.group(a);xa=invariant_charts(tuple(ea));ra=fixed(ea)
 assert tuple(i for i,j in q.charts(a))==xa and all(i==j for i,j in q.charts(a))
 # Returned axes of the selected stationary cores are fixed pointwise as rays.
 rr=q.rows(a)
 assert all(q.one_label(rr[r])==r for r in ra)
 for di,d in enumerate(defs):
  ed=d['group'];xd=d['charts'];sx=set(xd)
  if any(set(lag_action(e)[i] for i in xd)!=sx for e in ea):excluded+=1;continue
  tests=sorted({q.prod(e,f) for e in ed for f in ea})
  tests.remove(q.I);tests.insert(0,q.I)
  edges=[(q.LMASK[i],q.BASES[i]) for i in set(xa)&sx]
  returns=[(r,r) for r in ra&set(d['fixed'])]
  cases.append({'name':an+'_'+str(di)+'_'+d['kind'],'tests':tests,
  'station':[(f,xd) for f in ea],'edges':edges,'returns':returns,'balanced':[]})
 print(an,'cases',len(cases),'excluded',excluded,'seconds',round(time.perf_counter()-t,1),flush=True)
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
(out/'third-level.bin').write_bytes(raw)
(out/'third-level.json').write_text(json.dumps({'cases':cases,'conjugate_groups':defs,'excluded':excluded,'elapsed':time.perf_counter()-t},indent=2)+'\n',encoding='utf-8')
print('DONE',len(cases),len(raw),time.perf_counter()-t,flush=True)
