#!/usr/bin/env python3
"""Command-line driver: statements against the paper, file integrity, and exact recomputation of the finite facts."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT/'src'))


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def pairs(values):
    result={}
    for k,v in values:
        need(k not in result, 'duplicate JSON key: '+k)
        result[k]=v
    return result


def read(p):
    return json.loads(p.read_text(encoding='utf-8'), object_pairs_hook=pairs)


def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')


def check_dependency_dag(graph):
    state={}
    def visit(label):
        need(state.get(label)!=1,'cyclic proof dependencies: '+label)
        if state.get(label)==2:
            return
        state[label]=1
        for dep in graph[label]:
            need(dep in graph,'unknown proof dependency')
            visit(dep)
        state[label]=2
    for label in graph:
        visit(label)


def validate_mathematical_result(fact,route,report):
    need(isinstance(report,dict),'checker returned no mathematical verdict')
    need(report.get('status') in ('accepted','PASS'),'checker did not accept '+fact['navigation'])
    need(report.get('fact')==fact['label'].split(':',1)[1],'checker returned the wrong fact')
    need(report.get('route')==route,'checker returned the wrong route')
    need(report.get('exact') is True,'checker did not establish exact arithmetic')


BOX=re.compile(r'\\begin\{(finitefact|evidencefact)\}')
COLOUR_GROUP=re.compile(r'^[ \t]*(?:\\begingroup\\color\{(?:violet|cyan)\}|\\endgroup)[ \t]*$',re.M)


def strip_comments(text):
    return '\n'.join(re.sub(r'(?<!\\)%.*','',line) for line in text.split('\n'))


def braced(text,i):
    """Return the content of the balanced group starting at text[i]=='{' and the index after it."""
    need(i<len(text) and text[i]=='{','expected a braced argument')
    depth=0
    for j in range(i,len(text)):
        if text[j]=='{' and text[j-1]!='\\':
            depth+=1
        elif text[j]=='}' and text[j-1]!='\\':
            depth-=1
            if depth==0:
                return text[i+1:j],j+1
    raise ValueError('unbalanced braces')


def uncolour(text):
    text=COLOUR_GROUP.sub('',text)
    while (m:=re.search(r'\\textcolor\{(?:violet|cyan)\}\{',text)):
        inner,end=braced(text,m.end()-1)
        text=text[:m.start()]+inner+text[end:]
    return text


def extract_boxes(text):
    """Every fact box of the paper: label -> environment, title, consumers and the statement above the dashed line."""
    text=strip_comments(text.replace('\r\n','\n'))
    out={}
    for m in BOX.finditer(text):
        label,i=braced(text,m.end())
        title,i=braced(text,i)
        consumers,i=braced(text,i)
        need(re.fullmatch(r'(?:cert|ev):[a-z0-9-]+',label) is not None,'malformed fact label: '+label)
        need(label not in out,'duplicate paper label: '+label)
        end=text.index('\\end{'+m[1]+'}',i)
        body=text[i:end]
        lower=body.find('\\tcblower')
        if lower>=0:
            body=body[:lower]
        lines=[line.rstrip() for line in uncolour(body).split('\n')]
        out[label]={'environment':m[1],'title':' '.join(title.split()),'consumers':' '.join(consumers.split()),
                    'statement':'\n'.join(lines).strip('\n')+'\n'}
    return out


def extract_statements(text):
    return {label:box['statement'] for label,box in extract_boxes(text).items()}


def statement_digest(label,statement):
    """Statements are compared up to colour groups, comments and line breaks: the label and the words are hashed."""
    return hashlib.sha256((label+'\n'+' '.join(statement.split())).encode('utf-8')).hexdigest()


def registry(root=ROOT):
    data=read(root/'paper-map.json')
    facts=data['facts']
    need(len(facts)>0,'empty fact registry')
    labels=[f['label'] for f in facts]
    folders=[f['folder'] for f in facts]
    codes=[f['navigation'] for f in facts]
    need(len(set(labels))==len(set(folders))==len(set(codes))==len(facts),'nonbijective fact registry')
    in_paper=[f for f in facts if f.get('in_paper',True)]
    proofs=sum(f['classification']=='proof' for f in in_paper)
    extra=len(facts)-len(in_paper)
    need(set(codes)=={f'C{i}' for i in range(1,proofs+1)}|{f'E{i}' for i in range(1,len(in_paper)-proofs+1)}
         |{f'A{i}' for i in range(1,extra+1)},'changed navigation/classification')
    need(all(f.get('in_paper',True) or f['navigation'].startswith('A') for f in facts),'additional checks are numbered A1, A2, ...')
    expected_dirs={p.relative_to(root).as_posix() for p in (root/'certificates').iterdir() if p.is_dir()}
    need(expected_dirs==set(folders),'certificate folders and registry disagree')
    for f in facts:
        need(f['folder']=='certificates/'+f['label'].split(':',1)[1],'folder does not match permanent label')
        need(f['classification']==('proof' if f['label'].startswith('cert:') else 'evidence'),'classification drift')
        spec=read(root/f['folder']/'spec.json')
        statement=root/f['folder']/'statement.tex'
        need(spec['paper_label']==f['label'] and spec['fact']==f['label'].split(':',1)[1],'spec identity drift')
        text=statement.read_text(encoding='utf-8')
        need(statement_digest(f['label'],text)==f['statement_sha256']==spec['statement']['sha256'],'statement hash drift: '+f['navigation'])
        need(text==spec['statement']['full_tex'],'spec full statement drift')
        need(set(f['dependencies'])<=set(labels),'unresolved fact dependency')
        need(spec['dependent_certificates']==f['dependencies'],'proof dependency drift')
        need(set(spec['routes'])=={'recompute','verify_evidence'},'both routes required')
        for route in spec['routes'].values():
            need(route['availability'] in ('available','partial','not-yet-available'),'unknown route availability')
    check_dependency_dag({f['label']:f['dependencies'] for f in facts})
    frozen=extract_statements((root/'paper/finite-certificate.tex').read_text(encoding='utf-8'))
    need(set(frozen)=={f['label'] for f in in_paper},'paper, registry and folders disagree')
    for f in in_paper:
        need(statement_digest(f['label'],frozen[f['label']])==f['statement_sha256'],'frozen statements changed')
    return data


def check_paper(manuscript=None,root=ROOT):
    data=registry(root)
    if manuscript:
        supplied=extract_statements(Path(manuscript).read_text(encoding='utf-8'))
        expected={f['label']:f['statement_sha256'] for f in data['facts'] if f.get('in_paper',True)}
        need(set(supplied)==set(expected),'supplied manuscript has missing or extra facts')
        for label,value in supplied.items():
            need(statement_digest(label,value)==expected[label],'supplied statement changed: '+label)
    in_paper=[f for f in data['facts'] if f.get('in_paper',True)]
    proofs=sum(f['classification']=='proof' for f in in_paper)
    return {'status':'matched','kind':'paper-correspondence','proof_facts':proofs,'evidence_items':len(in_paper)-proofs,
            'additional_checks':len(data['facts'])-len(in_paper),
            'paper_revision':data['paper']['revision'],'mathematical_acceptance':False}


def integrity(root=ROOT):
    manifest=read(root/'provenance/integrity.json')
    expected=manifest['files']
    protected={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()
               and p.relative_to(root).parts[0] not in ('.git','runs','.venv','.tools')
               and '__pycache__' not in p.relative_to(root).parts
               and p.suffix not in ('.pyc','.pyo') and p.relative_to(root).as_posix()!='provenance/integrity.json'}
    need(protected==set(expected),'missing or extra release file')
    for rel,digest in expected.items():
        p=(root/rel).resolve()
        need(p.is_relative_to(root.resolve()),'nonportable integrity path')
        need(sha(p)==digest,'integrity mismatch: '+rel)
    return {'status':'matched','kind':'artifact-integrity','checked_files':len(expected),'mathematical_acceptance':False}


def select(code,data):
    matches=[f for f in data['facts'] if code in (f['navigation'],f['label'],f['label'].split(':',1)[1])]
    need(len(matches)==1,'unknown fact: '+code)
    return matches[0]


def c1(route):
    owner=ROOT/'four-qubit'
    with tempfile.TemporaryDirectory(prefix='paper-i-c1-') as temp:
        work=Path(temp)
        shutil.copytree(owner/'data',work/'data')
        (work/'src').mkdir(); (work/'evidence').mkdir()
        shutil.copyfile(owner/'src/check_coverage.py',work/'src/check_coverage.py')
        result=subprocess.run([sys.executable,'-B',str(work/'src/check_coverage.py')],cwd=work,
                              capture_output=True,text=True,encoding='utf-8')
        need(result.returncode==0,'C1 coverage mathematics rejected: '+result.stderr[-2000:])
        report=read(work/'evidence/coverage.json')
        expected={'triangular_seeds':32768,'triangular_involutions':2064,'permutation_classes':36,
                  'point_list':163216,'phase_cases':43,'raw_phase_rows':2818048,'roots':1020,'semantic_negative_controls':2}
        need(report['status']=='PASS' and all(report[k]==v for k,v in expected.items()),'C1 required domain changed')
        need(len(set(report['largeN']))==10,'C1 ten-root threshold failed')
        return {'status':'accepted','route':route,'exact':True,'result':report,
                'fact':'r4-representatives-up-to-clifford-conjugacy',
                'scope':'complete finite row coverage and all 1020 root invariants, with exact supplied affine/phase-forest semantics; analytic entrance retained',
                'independent_implementation':False}


# Each fact is checked by one backend, named in its spec.json as "backend": [module, key]. The keys are
# internal; the numbering of the paper is in paper-map.json.
BACKEND_MODULES=('c7','c8','c14','c15','c16','c17','conventional','conventional_large','e6_level',
                 'exact_audit','n3_descent','root_words','row_coverage')


def backend_of(fact):
    entry=read(ROOT/fact['folder']/'spec.json').get('backend')
    need(isinstance(entry,list) and len(entry)==2 and entry[0] in BACKEND_MODULES,
         'no maintained checker for '+fact['navigation'])
    return entry[0],entry[1]


def execute(fact,route,args=None,workspace=None):
    slug=fact['label'].split(':',1)[1]
    backend,key=backend_of(fact)
    if backend=='row_coverage':
        return {**c1(route),'fact':slug}
    if route=='verify_evidence':
        raise ValueError('no complete supplied-evidence verifier for '+fact['navigation'])
    if backend in ('exact_audit','n3_descent'):
        return importlib.import_module('checking.'+backend+'.runner').run(key,output_dir=None)
    if backend=='root_words':
        report=importlib.import_module('checking.root_words.runner').run(key,output_dir=None)
        need(isinstance(report,dict) and report.get('fact')==key,'root-word checker returned the wrong fact')
        return {**report,'fact':slug}
    if backend=='e6_level':
        return importlib.import_module('checking.e6_level').run(workspace=workspace)
    if backend=='conventional':
        report=importlib.import_module('checking.conventional').run(key)
        need(isinstance(report,dict),'exact conventional checker returned no result')
        return {'exact':True,**report}
    if backend=='conventional_large':
        report=importlib.import_module('checking.conventional_large').run(key,mode='full',allow_expensive=bool(args and args.allow_expensive),
             compiler=args.compiler if args else None,workers=args.workers if args else 1,work_directory=workspace)
        need(isinstance(report,dict) and report.get('fact')==key,'large checker returned the wrong fact')
        return {'exact':True,**report,'fact':slug}
    report=importlib.import_module('checking.'+backend).run()
    if backend in ('c14','c15','c16'):
        # These checkers do not report the 'exact' field themselves; their arithmetic is exact,
        # and their verdicts are passed on unchanged.
        need(isinstance(report,dict),'exact checker returned no result')
        report={'exact':True,**report}
    return report


def mathematical(args):
    need(sys.flags.optimize==0,'optimized Python (-O) turns off the assertions the checks rely on; run Python without -O')
    need(sum(bool(x) for x in (args.all,args.fast,args.fact))==1,'choose exactly one of --fact, --fast, --all')
    data=registry()
    chosen=data['facts'] if args.all else [select(c,data) for c in args.fact]
    if args.fast:
        chosen=[f for f in data['facts'] if read(ROOT/f['folder']/'spec.json')['routes'][args.route].get('integrated')
                and not read(ROOT/f['folder']/'spec.json')['routes'][args.route].get('requires_explicit_expensive')]
    need(chosen,'select --fact, --fast or --all')
    # Preflight the entire requested set. No full acceptance through an implicit skip.
    pending=[]
    for f in chosen:
        route=read(ROOT/f['folder']/'spec.json')['routes'][args.route]
        if route['availability']!='available' or not route.get('integrated'):
            pending.append({'fact':f['navigation'],'route':args.route,'availability':route['availability'],'remaining_work':route['remaining_work']})
        elif route.get('requires_explicit_expensive') and not args.allow_expensive:
            pending.append({'fact':f['navigation'],'route':args.route,'availability':'available-expensive-command',
                            'remaining_work':['Use --allow-expensive for the complete census; see runtime_seconds in spec.json.']})
    if pending:
        print(json.dumps({'status':'incomplete','kind':'mathematical-check','route':args.route,'pending':pending,'mathematical_acceptance':False},indent=2))
        return 2
    integrity()
    for f in chosen:
        identity=read(ROOT/f['folder']/'spec.json').get('checker_identity')
        need(isinstance(identity,dict) and identity,'missing checker identity: '+f['navigation'])
        need(all(sha(ROOT/rel)==digest for rel,digest in identity.items()),'stale maintained checker: '+f['navigation'])
    cost=[{'fact':f['navigation'],'finite_domain':read(ROOT/f['folder']/'spec.json')['finite_domain'],
           'runtime_seconds':read(ROOT/f['folder']/'spec.json')['routes'][args.route]['runtime_seconds']} for f in chosen]
    print(json.dumps({'checking':cost,'route':args.route},indent=2),flush=True)
    accepted=[]
    for f in chosen:
        spec=read(ROOT/f['folder']/'spec.json')
        start=time.perf_counter()
        runid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+'-'+f['navigation']
        run=ROOT/'runs'/runid
        run.mkdir(parents=True,exist_ok=False)
        report=execute(f,args.route,args,workspace=run/'workspace')
        # Each implementation must itself reject failures. Driver binds its scoped verdict.
        validate_mathematical_result(f,args.route,report)
        bound={'schema_version':1,'status':'accepted','kind':'mathematical-check','exact':True,
               'fact':f['label'],'navigation':f['navigation'],'route':args.route,'statement_sha256':f['statement_sha256'],
               'paper_revision':data['paper']['revision'],'elapsed_seconds':time.perf_counter()-start,
               'environment':{'python':sys.version,'platform':sys.platform},
               'analytic_premises':spec['analytic_entrance_exit_labels'],'formal_verification':False,
               'conventions':spec.get('checker_conventions','paper/conventions.json'),
               'checker_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in (ROOT/'src/checking').rglob('*.py')},'result':report}
        write(run/'result.json',bound)
        accepted.append({'fact':f['navigation'],'run':run.relative_to(ROOT).as_posix(), 'elapsed_seconds':bound['elapsed_seconds']})
        print(json.dumps({'accepted':accepted[-1],'route':args.route}),flush=True)
    unchecked=[f['navigation'] for f in data['facts'] if f not in chosen]
    print(json.dumps({'status':'accepted-requested-scope' if not args.fast else 'partial-repository-acceptance',
                      'route':args.route,'accepted':accepted,'pending_facts':unchecked,
                      'all_proof_facts_accepted':all(f['classification']!='proof' or f in chosen for f in data['facts'])},indent=2))
    return 3 if args.fast and unchecked else 0


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest='command',required=True)
    p=subs.add_parser('paper');p.add_argument('--manuscript',type=Path)
    subs.add_parser('integrity')
    subs.add_parser('status')
    p=subs.add_parser('math')
    p.add_argument('--fact',nargs='+',default=[])
    p.add_argument('--route',choices=('recompute','verify_evidence'),default='recompute')
    p.add_argument('--fast',action='store_true')
    p.add_argument('--all',action='store_true')
    p.add_argument('--allow-expensive',action='store_true')
    p.add_argument('--compiler',help='installed C compiler path or supported command')
    p.add_argument('--workers',type=int,default=1)
    args=parser.parse_args(argv)
    try:
        if args.command=='math':
            rc=mathematical(args)
        elif args.command=='paper':
            print(json.dumps(check_paper(args.manuscript),indent=2));rc=0
        elif args.command=='integrity':
            print(json.dumps(integrity(),indent=2));rc=0
        else:
            data=registry()
            rows=[]
            for f in data['facts']:
                spec=read(ROOT/f['folder']/'spec.json')
                rows.append({'fact':f['navigation'],'name':spec['fact'],'classification':spec['classification'],
                             'recompute':spec['routes']['recompute'],'verify_evidence':spec['routes']['verify_evidence']})
            recorded=all(r['recompute'].get('last_accepted_run') for r in rows if r['classification']=='proof')
            print(json.dumps({'facts':rows,'all_proof_facts_have_recorded_runs':recorded},indent=2));rc=0
    except (ValueError,FileNotFoundError,KeyError,ImportError) as exc:
        print(json.dumps({'status':'rejected','error':str(exc),'mathematical_acceptance':False},indent=2),file=sys.stderr)
        rc=1
    raise SystemExit(rc)


if __name__=='__main__':
    main()
