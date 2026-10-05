"""Run the proposed-plane checker, with strict scope and result validation."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,hashlib,json,re,subprocess,time
P=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(v,m):
    if not v:raise ValueError(m)
def validate(text,rid):
    lines=[x for x in text.splitlines() if x.startswith('ROOT ')]
    need(len(lines)==1,'exactly one final root record required')
    s=lines[0];f=s.split()
    need(f[1]==str(rid) and f[6]=='survivors' and f[7]=='1','root/acceptance mismatch')
    for k in ('incomplete_families','transport_failures','ray_mismatch'):
        need(k in f and f[f.index(k)+1]=='0','nonzero or missing '+k)
    for k,v in [('kernel','publication-lowtier-v1'),('scope','full'),('guided','1'),('charts','0:2294')]:
        need(k in f and f[f.index(k)+1]==v,'scope mismatch '+k)
    need('FIBER_TOO_LARGE' not in text and 'INCOMPLETE' not in text,'unfinished family')
    need('skel_checked' in f and f[f.index('skel_checked')+1]=='374580720','incomplete chart/permutation visits')
    return s
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--mode',choices=['all','exact','pair'],default='all')
    args=ap.parse_args();need(1<=args.workers<=64,'worker count')
    subprocess.run([__import__('sys').executable,str(P/'src/check_coverage.py')],check=True)
    exe=P/'src'/('ge4_check.exe' if __import__('os').name=='nt' else 'ge4_check')
    files=[P/'src/ge4_check.c',*sorted((P/'src').glob('*.h')),
           P/'data/roots-n4-all.txt',P/'data/lock-n4.txt',P/'data/pointdomain-n4.bin',P/'data/fingerprints2-n4.txt']
    signature={x.relative_to(P).as_posix():sha(x) for x in files};binary=sha(exe)
    jobs=[]
    for line in (P/'data/lock-n4.txt').read_text().splitlines():
        f=line.split();mode=f[f.index('mode')+1]
        if mode=='largeN' or (args.mode!='all' and mode!=args.mode):continue
        jobs.append(int(f[0]))
    out=P/'evidence/recheck';out.mkdir(exist_ok=True);start=time.perf_counter()
    def run(rid):
        base=out/f'root-{rid:04d}';jp=base.with_suffix('.json');rp=base.with_suffix('.txt')
        if jp.exists():
            old=json.loads(jp.read_text())
            if old.get('inputs')==signature and old.get('binary_sha256')==binary and old.get('status')=='PASS':
                need(sha(rp)==old['result_sha256'],'cached result hash')
                validate(rp.read_text(),rid);return dict(old,reused=True)
        if rp.exists():rp.unlink()
        command=[str(exe),'lock','roots-n4-all.txt','pointdomain-n4.bin',str(rp),str(rid),str(rid),
                 'fingerprints2-n4.txt','rigid3','inside','guided']
        t=time.perf_counter()
        with base.with_suffix('.stdout').open('w') as stdout,base.with_suffix('.stderr').open('w') as stderr:
            result=subprocess.run(command,cwd=P/'data',stdout=stdout,stderr=stderr)
        need(result.returncode==0,f'root {rid} process failed')
        line=validate(rp.read_text(),rid)
        expected=re.search(r'\(\d+,\d+,\d+\)',(P/'data/lock-n4.txt').read_text().splitlines()[rid]).group()
        need(expected in line,'proposed plane mismatch')
        r={'root':rid,'status':'PASS','wall_seconds':time.perf_counter()-t,'inputs':signature,
           'binary_sha256':binary,'result_sha256':sha(rp),'result':rp.relative_to(P).as_posix()}
        jp.write_text(json.dumps(r,indent=2)+'\n');return dict(r,reused=False)
    records=[];failures=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run,r):r for r in jobs}
        for f in as_completed(futures):
            try:
                r=f.result();records.append(r)
                print(json.dumps({'done':len(records)+len(failures),'total':len(jobs),'root':r['root'],'seconds':r['wall_seconds'],'reused':r['reused']}),flush=True)
            except Exception as e:
                failures.append({'root':futures[f],'error':str(e)})
                print(json.dumps({'failed':failures[-1]}),flush=True)
    report={'status':'PASS' if not failures else 'FAIL','scope':args.mode,'failures':failures,
        'new_roots_completed':sum(not r['reused'] for r in records),'cached_roots_reused':sum(r['reused'] for r in records),'roots':sorted(r['root'] for r in records),
        'workers':args.workers,'elapsed_seconds':time.perf_counter()-start,
        'sum_worker_elapsed_seconds':sum(r['wall_seconds'] for r in records),'inputs':signature,
        'binary_sha256':binary}
    (P/f'evidence/recheck-{args.mode}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('inputs','roots')}),flush=True)
    if failures:raise SystemExit(2)
if __name__=='__main__':main()


