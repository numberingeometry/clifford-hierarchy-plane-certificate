"""Fresh census of the 1,010 four-qubit gates of (C7) in a new workspace.

Copies src/, data/ and tools/ into runs/<directory>/, compiles the exact C kernel, runs its regression tests,
recomputes the coverage of the finite domain, and recomputes every gate with tools/recheck_locks.py. Each
gate's result must reproduce the plane, level and sector listed for it in data/lock-n4.txt.
"""
from pathlib import Path
import argparse, datetime, hashlib, importlib.util, json, os, shutil, subprocess, sys, time

if not __debug__:
    raise SystemExit('Run without -O: the exact checkers use assertions.')
ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
REGRESSIONS = {
    'testfamily': 'FAMILY TEST PASSED: repaired invariant correct on all 8 cases, legacy accumulator wrong where expected',
    'fasttests': 'FAST_TESTS_PASS modular_systems=12000 affine_families=2160',
    'testtiers': 'LOW_TIER_TEST_PASS roots=1020 generator_vs_all_Paulis=equal',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compiler():
    cc = os.environ.get('CC')
    if cc:
        return [cc]
    for name in ('clang', 'gcc'):
        if shutil.which(name):
            return [name]
    if importlib.util.find_spec('ziglang'):
        return [sys.executable, '-m', 'ziglang', 'cc']
    raise RuntimeError('Install GCC, Clang, or ziglang; alternatively set CC to a C compiler executable.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['locks'])
    p.add_argument('--workers', type=int, default=6)
    p.add_argument('--directory', help='New directory name under runs/, never an existing directory')
    a = p.parse_args()
    if not 1 <= a.workers <= 64:
        p.error('--workers must be between 1 and 64')
    if sys.byteorder != 'little':
        raise SystemExit('The binary input format requires a little-endian host.')
    spec = importlib.util.spec_from_file_location('verify', ROOT/'scripts/verify.py')
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    verifier.validate_integrity(ROOT)
    frozen = json.loads((ROOT/'manifest.json').read_text())['files']
    name = a.directory or datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+a.mode
    if Path(name).name != name or name in ('.', '..'):
        raise SystemExit('Directory must be a single new name.')
    work = ROOT/'runs'/name
    if work.exists():
        raise SystemExit('Refuse to reuse a run directory.')
    work.mkdir(parents=True)
    for d in ('src', 'data', 'tools'):
        shutil.copytree(ROOT/d, work/d, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    (work/'evidence').mkdir()
    (work/'logs').mkdir()
    started = time.perf_counter()
    steps = []
    report = {'status': 'RUNNING', 'mode': a.mode, 'fresh_all_roots': False, 'reused_canonical_records': False,
              'steps': steps, 'runner_sha256': sha(Path(__file__).resolve()),
              'source_sha256': {f.relative_to(ROOT).as_posix(): sha(f) for f in (ROOT/'src').rglob('*')
                                if f.is_file() and '__pycache__' not in f.parts},
              'manifest_sha256': sha(ROOT/'manifest.json'),
              'package_inputs': {name: sha(ROOT/name) for name in frozen},
              'workspace_inputs': {f.relative_to(work).as_posix(): sha(f) for f in work.rglob('*') if f.is_file()},
              'python': sys.version.split()[0]}

    def save():
        report['elapsed_seconds'] = time.perf_counter()-started
        (work/'run.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')

    def run(label, args, cwd=None):
        t = time.perf_counter()
        env = dict(os.environ, PYTHONHASHSEED='0', PYTHONUTF8='1')
        cp = subprocess.run(args, cwd=cwd or work, capture_output=True, text=True, encoding='utf-8', env=env)
        (work/'logs'/(label+'.txt')).write_text(cp.stdout+'\n'+cp.stderr, encoding='utf-8')
        steps.append({'check': label, 'exit': cp.returncode, 'seconds': time.perf_counter()-t})
        save()
        if cp.returncode:
            raise RuntimeError(label+' failed: '+cp.stderr[-1200:])
        print(label+' PASS', flush=True)
        return cp.stdout

    save()
    try:
        cc = compiler()
        exe = work/'src'/('ge4_check.exe' if os.name == 'nt' else 'ge4_check')
        report['compiler'] = Path(cc[0]).name
        run('compiler-version', cc+['--version'])
        abi = work/'abi.c'
        abi.write_text('#include <stdint.h>\nint main(void){unsigned x=1;return sizeof(unsigned)!=4 || '
                       'sizeof(uint64_t)!=8 || *(unsigned char*)&x!=1;}\n')
        abiexe = work/('abi.exe' if os.name == 'nt' else 'abi')
        run('compile-abi', cc+[str(abi), '-o', str(abiexe)])
        run('binary-abi', [str(abiexe)])
        run('compile-lock', cc+['-O2', '-std=c99', '-DFAST_TESTS', '-o', str(exe), str(work/'src/ge4_check.c')])
        report['lock_binary_sha256'] = sha(exe)
        run('coverage', [sys.executable, '-X', 'utf8', str(work/'src/check_coverage.py')])
        for test in ('selftest', 'testfamily', 'fasttests', 'testtiers'):
            args = [str(exe), test]
            if test == 'testtiers':
                args += [str(work/'data/roots-n4-all.txt')]
            output = run(test, args, cwd=work/'data')
            if test in REGRESSIONS and REGRESSIONS[test] not in output:
                raise RuntimeError('Missing exact regression result: '+test)
        run('all-root-locks', [sys.executable, str(work/'tools/recheck_locks.py'), '--mode', 'all',
                               '--workers', str(a.workers)])
        result = json.loads((work/'evidence/recheck-all.json').read_text())
        if (result.get('status') != 'PASS' or result.get('scope') != 'all' or
                result.get('roots') != sorted(verifier.ROOTS) or result['new_roots_completed'] != 1010 or
                result['cached_roots_reused'] != 0 or result['failures']):
            raise RuntimeError('The fresh census did not cover all 1010 gates.')
        signature = {name: sha(work/name) for name in verifier.LOCK_INPUTS}
        if result['inputs'] != signature or result['binary_sha256'] != sha(exe):
            raise RuntimeError('Census source/binary binding mismatch.')
        lines = (work/'data/lock-n4.txt').read_text().splitlines()
        for rid in sorted(verifier.ROOTS):
            base = work/'evidence/recheck'/f'root-{rid:04d}'
            record = json.loads(base.with_suffix('.json').read_text())
            text = base.with_suffix('.txt').read_text()
            if (record['root'] != rid or record['status'] != 'PASS' or record['inputs'] != signature or
                    record['binary_sha256'] != sha(exe) or record['result_sha256'] != sha(base.with_suffix('.txt'))):
                raise RuntimeError('Per-gate binding mismatch: '+str(rid))
            fields = lines[rid].split()
            verifier._root_log(text, rid, fields[-1], int(fields[2]), int(fields[4]))
        report['fresh_all_roots'] = True
        for name, digest in report['package_inputs'].items():
            if sha(ROOT/name) != digest:
                raise RuntimeError('Package changed during execution: '+name)
        report['status'] = 'PASS'
        save()
        print(json.dumps({'status': 'PASS', 'mode': a.mode, 'fresh_all_roots': report['fresh_all_roots'],
                          'run': work.relative_to(ROOT).as_posix(), 'seconds': report['elapsed_seconds']}, indent=2))
    except Exception as e:
        report['status'] = 'FAIL'
        report['error'] = str(e)
        save()
        raise


if __name__ == '__main__':
    main()
