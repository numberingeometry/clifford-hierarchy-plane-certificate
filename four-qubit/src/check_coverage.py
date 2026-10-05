"""Exact checks of the coverage of the finite domains of (C6) and (C7).
No orbit search or generator implementation is imported.
Arithmetic is integer / Z[zeta_8]; numpy is used only for finite arrays.
"""
from pathlib import Path
import hashlib,json,re,struct,time
import numpy as np

PUB=Path(__file__).resolve().parents[1]
DATA=PUB/'data'
def need(ok,msg):
    if not ok: raise ValueError(msg)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def pack(p): return sum(int(y)<<(4*x) for x,y in enumerate(p))
def perm(p): return len(p)==16 and sorted(map(int,p))==list(range(16))
def affine(p):
    return perm(p) and all(p[x^y]==(p[x]^p[y]^p[0]) for x in range(16) for y in range(16))
def inv(p):
    r=[0]*16
    for x,y in enumerate(p):r[y]=x
    return r
def pw(v,x):return (2*((v&15)&(v>>4)).bit_count()+4*((v>>4)&x).bit_count())%8
def conjugate(m,v):
    p,w=m;ip=inv(p);q=[];z=[]
    for x in range(16):
        y=ip[x];q.append(p[y^(v&15)])
        z.append((-w[y]+pw(v,y)+w[y^(v&15)])%8)
    return q,z
def pauli_ray(m):
    p,w=m;a=p[0]
    if p!=[x^a for x in range(16)]:return False
    b=sum((((w[1<<j]-w[0])%8)//4)<<j for j in range(4))
    return b<16 and all((w[x]-w[0])%8==(4*((b&x).bit_count()%2)) for x in range(16))
def clifford(m):return all(pauli_ray(conjugate(m,1<<j)) for j in range(8))
def tier(m):
    if pauli_ray(m):return 1
    if clifford(m):return 2
    return 3 if all(clifford(conjugate(m,1<<j)) for j in range(8)) else 4
def omega(u,v):return (((u&15)&(v>>4)).bit_count()+((u>>4)&(v&15)).bit_count())%2
def coeffs(m):
    p,w=m;r=[]
    for v in range(256):
        z=[0]*4
        for x in range(16):
            if p[x] != x^(v&15):continue
            e=(w[x]+pw(v,p[x]))%8
            z[e%4]+=1 if e<4 else -1
        r.append(tuple(z))
    return r
def invariants(m):
    co=coeffs(m);s=[v for v,z in enumerate(co) if any(z)]
    need(s and 0 not in s,"root trace")
    nr=[v for v in range(1,256) if all(omega(v,u^s[0])==0 for u in s)]
    nf=[v for v in range(1,256) if all(omega(v,u)==0 for u in s)]
    fp=tuple(sorted(min(co[v],tuple(-c for c in co[v])) for v in s))
    return len(nr),len(nf),fp,nr
def decode(pi,k):
    ids=np.arange(65536,dtype=np.uint32)
    bits=((ids[:,None]>>np.arange(16,dtype=np.uint32))&1).astype(np.int16)
    f=np.empty((65536,16),dtype=np.int16)
    for x,y in enumerate(pi):
        if x==y:f[:,x]=2*bits[:,x]
        elif x<y:f[:,x]=bits[:,x]+2*bits[:,y]
        else:f[:,x]=(-k-bits[:,y]-2*bits[:,x])%4
    return f
def check_forest(rec,words,gens,roots):
    pi=rec['pi'];k=rec['lambda'];local=rec['local_generators']
    need(perm(pi) and all(pi[pi[x]]==x for x in range(16)),"forest permutation")
    need(k in (0,1) and (k==0 or all(pi[x]!=x for x in range(16))),"forest sector")
    need(words.shape==(65536,),"forest length")
    f=decode(pi,k)
    need(bool(np.all((f+f[:,pi]+k)%4==0)),"Hermitian phase domain")
    fixed=[x for x in range(16) if pi[x]==x]
    balanced=(np.sum(np.where(f[:,fixed]==0,1,-1),axis=1)==0)
    dep=words>>24;edge=(words>>16)&255;parent=words&65535
    need(bool(np.array_equal(dep<255,balanced)),"balanced domain coverage")
    root_rows=np.flatnonzero(balanced & (dep==0))
    need(len(root_rows)>0 and bool(np.all(edge[root_rows]==255)),"forest roots")
    seen=[]
    for i in root_rows:
        rid=int(parent[i]);need(rid in roots,"root ID")
        rp,rw=roots[rid]
        need(pi==rp and np.array_equal((k+2*f[i])%8,rw),"root operator equality")
        seen.append(rid)
    non=np.flatnonzero(balanced & (dep>0))
    need(bool(np.all(edge[non]<len(local))),"forest edge bounds")
    need(bool(np.all(dep[parent[non]]<dep[non])),"strict depth decrease")
    for j,gid in enumerate(local):
        need(0<=gid<len(gens),"generator ID")
        rows=non[edge[non]==j]
        if not len(rows):continue
        q,z=gens[gid];iq=inv(q)
        need([q[pi[iq[x]]] for x in range(16)]==pi,"edge permutation semantics")
        # g M(parent) g^-1 at input x, including exact eighth-root phase.
        expected=(k+2*f[parent[rows]][:,iq]+np.array([z[pi[y]]-z[y] for y in iq]))%8
        need(bool(np.array_equal(expected,(k+2*f[rows])%8)),"edge phase semantics")
    return len(non)+len(root_rows),seen,int(max(dep[balanced]))
def main():
    start=time.perf_counter()
    inputs=[DATA/'roots-n4-all.txt',DATA/'lock-n4.txt',DATA/'fingerprints2-n4.txt',
            DATA/'pointdomain-n4.bin',DATA/'atlas-witnesses.json',DATA/'phase-forests/manifest.json']
    roots={};records={}
    for line in inputs[0].read_text().splitlines():
        a=list(map(int,line.split()));need(len(a)==36 and a[0] not in roots,"root format")
        rid,_,sector,k=a[:4];p=a[4:20];f=a[20:36]
        need(perm(p) and k in (0,1) and sector==k and all(0<=x<4 for x in f),"root fields")
        w=[(k+2*x)%8 for x in f]
        need(all(p[p[x]]==x and (w[x]+w[p[x]])%8==0 for x in range(16)),"root involution")
        roots[rid]=(p,w)
    need(set(roots)==set(range(1020)),"root indices")
    a=json.loads(inputs[4].read_text());classes=a['classes'];wits=a['witnesses']
    need(len(classes)==36 and len({tuple(p) for p in classes})==36,"permutation representatives")
    wd={tuple(x['pi']):x for x in wits};need(len(wd)==len(wits),"duplicate witness")
    tri=[]
    for key in range(32768):
        p=tuple(sum((((x>>j)&1)^((key>>((1<<j)-1+(x&((1<<j)-1))))&1))<<j for j in range(4)) for x in range(16))
        if all(p[p[x]]==x for x in range(16)):tri.append(p)
    need(len(tri)==2064 and set(tri)==set(wd),"all triangular involutions covered")
    need({tuple(w["canon"]) for w in wits}=={tuple(c) for c in classes},"every representative has triangular provenance")
    for p in tri:
        w=wd[p];q=w['alpha'];c=w['canon'];need(affine(q) and c in classes,"affine witness")
        iq=inv(q);need([q[p[iq[x]]] for x in range(16)]==c,"conjugacy equation")
    print("PERMUTATION_WITNESSES_PASS seeds=32768 involutions=2064 classes=36",flush=True)

    raw=inputs[3].read_bytes();n=struct.unpack_from('<I',raw)[0]
    need(n==163216 and len(raw)==4+20*n,"point file format")
    words=np.frombuffer(raw,dtype='<u8',count=n,offset=4)
    need(bool(np.all(words[1:]>words[:-1])),"point list sorted unique")
    pp=((words[:,None]>>np.arange(0,64,4,dtype=np.uint64))&15).astype(np.int16)
    need(bool(np.all(np.sort(pp,axis=1)==np.arange(16))),"point list permutations")
    need(bool(np.all(np.take_along_axis(pp,pp,axis=1)==np.arange(16))),"point list involutions")
    def contains(keys):
        pos=np.searchsorted(words,keys)
        return bool(np.all(pos<n)) and bool(np.all(words[pos]==keys))
    need(contains(np.array([pack(p) for p in tri],dtype=np.uint64)),"triangular seeds in point list")
    elementary=[[x^(1<<j) for x in range(16)] for j in range(4)]
    elementary += [[x^(((x>>j)&1)<<k) for x in range(16)] for j in range(4) for k in range(4) if j!=k]
    for q in elementary:
        moved=np.array(q)[pp[:,q]] # all generators here are involutions
        keys=np.sum(moved.astype(np.uint64)<<np.arange(0,64,4,dtype=np.uint64),axis=1,dtype=np.uint64)
        need(contains(keys),"affine point closure")
    print("POINT_COVERAGE_PASS list=163216 elementary_generators=16",flush=True)

    manifest=json.loads(inputs[5].read_text());gens=[]
    for g in manifest['generators']:
        q,z=g['perm'],g['weights'];need(perm(q) and len(z)==16 and all(0<=x<8 for x in z),"gate format")
        need(clifford((q,z)),"edge gate must be Clifford");gens.append((q,z))
    expected={(ci,k) for ci,p in enumerate(classes) for k in range(2) if k==0 or all(p[x]!=x for x in range(16))}
    actual=[(r['class'],r['lambda']) for r in manifest['cases']]
    need(set(actual)==expected and len(actual)==len(expected)==43,"all permutation sectors")
    total=0;ids=[];maxdepth=0
    for r in manifest['cases']:
        need(r['pi']==classes[r['class']],"forest class")
        fpath=DATA/'phase-forests'/Path(r['file']).name;inputs.append(fpath)
        need(sha(fpath)==r['sha256'],"forest byte integrity")
        w=np.fromfile(fpath,dtype='<u4');balanced,rr,d=check_forest(r,w,gens,roots)
        total+=balanced;ids+=rr;maxdepth=max(maxdepth,d)
    need(len(ids)==1020 and set(ids)==set(roots),"atlas roots exactly covered")
    print(f"PHASE_COVERAGE_PASS cases=43 raw=2818048 balanced={total} roots=1020 maxdepth={maxdepth}",flush=True)

    # Reconstruct the invariant table; neither recorded levels nor atlas statistics are trusted.
    fp_expected={};large=[];tiers={}
    for rid,m in roots.items():
        nn,nf,fp,rays=invariants(m);t=tier(m);tiers[t]=tiers.get(t,0)+1
        fp_expected.setdefault(fp,set()).add((nf,nn))
        records[rid]=(nn,nf,t,rays)
        if nn>=63:
            need(t<=2,"large normalizer must be Clifford");large.append(rid)
    tok=iter(map(int,inputs[2].read_text().split()));nfps=next(tok);fps={}
    for _ in range(nfps):
        nco=next(tok);fp=tuple(tuple(next(tok) for _ in range(4)) for _ in range(nco))
        pairs={(next(tok),next(tok)) for _ in range(next(tok))}
        need(fp not in fps,"duplicate fingerprint");fps[fp]=pairs
    need(next(tok,None) is None and fps==fp_expected,"fingerprint table exact agreement")
    modes={};count=0
    for line in inputs[1].read_text().splitlines():
        f=line.split();rid=int(f[0]);count+=1
        need(rid==count-1,"lock ID ordering")
        nn,nf,t,rays=records[rid];mode=f[f.index('mode')+1]
        need(nn==int(f[f.index('N')+1]),"lock root invariant")
        modes[mode]=modes.get(mode,0)+1
        if mode=='largeN':need(rid in large,"largeN disposition")
        else:
            need(mode in ('pair','exact') and rid not in large,"lock mode")
            match=re.search(r'\((\d+),(\d+),(\d+)\)',line);need(match is not None,"lock plane")
            u,v,w=map(int,match.groups());need(0<u<256 and 0<v<256 and u!=v and u^v==w,"plane relation")
            need(any(x in rays for x in (u,v,w)),"plane locks its root")
    need(count==1020 and len(large)==10,"disposition coverage")
    print(f"ROOT_INVARIANTS_PASS roots=1020 fingerprints={len(fps)} largeN=10 tiers={tiers}",flush=True)
    # Semantic negative controls: corrupt an actual forest edge while preserving dimensions.
    rec=next(r for r in manifest['cases'] if r['balanced_rows']>r['atlas_roots'])
    words=np.fromfile(DATA/'phase-forests'/Path(rec['file']).name,dtype='<u4')
    i=int(np.flatnonzero(((words>>24)>0)&((words>>24)<255))[0])
    bad=words.copy();bad[i]=(bad[i]&np.uint32(0xffff0000))|np.uint32(i)
    try:check_forest(rec,bad,gens,roots)
    except ValueError as e:need('depth' in str(e),"negative control cause")
    else:raise ValueError("self-parent corruption accepted")
    # A root edge must denote its precise operator, not merely an existing atlas ID.
    i=int(np.flatnonzero((words>>24)==0)[0]);bad=words.copy()
    bad[i]=(bad[i]&np.uint32(0xffff0000))|np.uint32((int(bad[i]&65535)+1)%1020)
    try:check_forest(rec,bad,gens,roots)
    except ValueError as e:need('operator equality' in str(e),"negative control root cause")
    else:raise ValueError("wrong-root corruption accepted")
    report={'status':'PASS','scope':'finite-domain coverage and root data; the per-gate census is separate',
       'seconds':time.perf_counter()-start,'triangular_seeds':32768,'triangular_involutions':2064,
       'permutation_classes':36,'point_list':163216,'phase_cases':43,'raw_phase_rows':2818048,
       'balanced_phase_rows':total,'roots':1020,'forest_generators':len(gens),'max_forest_depth':maxdepth,
       'largeN':large,'tiers':tiers,'dispositions':modes,'semantic_negative_controls':2,
       'checker_sha256':sha(Path(__file__)),
       'inputs':{p.relative_to(PUB).as_posix():sha(p) for p in inputs}}
    (PUB/'evidence/coverage.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inputs'}),flush=True)
if __name__=='__main__':main()

