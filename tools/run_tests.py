#!/usr/bin/env python3
"""Run the tests. The fixture controls of (C15) need the workspace of a fresh run (--c6-workspace)."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT/'src'))
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--c6-workspace', type=Path)
ap.add_argument('--require-c6-controls', action='store_true')
a = ap.parse_args()
if a.c6_workspace:
    work = a.c6_workspace.resolve()
    if not (work/'physical-bindings.json').is_file():
        ap.error('C6 workspace lacks physical bindings')
    os.environ['C6_TEST_WORKSPACE'] = str(work)
elif a.require_c6_controls and not os.environ.get('C6_TEST_WORKSPACE'):
    ap.error('--require-c6-controls needs --c6-workspace or C6_TEST_WORKSPACE')
stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
output = ROOT/'runs'/('tests-'+stamp)
output.mkdir(parents=True, exist_ok=False)
start = time.perf_counter()
with (output/'tests.txt').open('w', encoding='utf-8', newline='\n') as stream:
    suite = unittest.defaultTestLoader.discover(str(ROOT/'tests'))
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
summary = {
    'kind': 'implementation-and-rejection-controls',
    'status': 'PASS' if result.wasSuccessful() else 'FAIL',
    'tests_run': result.testsRun, 'failures': len(result.failures),
    'errors': len(result.errors), 'skipped': len(result.skipped),
    'complete_C6_fixture_controls': bool(os.environ.get('C6_TEST_WORKSPACE')),
    'elapsed_seconds': time.perf_counter()-start,
    'environment': {'python': sys.version, 'platform': sys.platform},
    'mathematical_acceptance': False,
    'output': output.relative_to(ROOT).as_posix(),
}
(output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps(summary, indent=2))
raise SystemExit(0 if result.wasSuccessful() else 1)
