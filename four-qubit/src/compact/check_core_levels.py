"""Exact hierarchy level of the thirty canonical block cores of the block cover, and of all 512 rooted block forms.

The block cover lists twenty non-third-level cores (fourteen marked, six unmarked) and ten
third-level cores.  A mark can change membership in C_3 (pi is third-level, pi T_z is not), so membership is recorded
for marked and unmarked forms alike.  This script records, with the same exact primitives as the other compact checks (`compact_core.is_c3`:
every Pauli correction is Clifford):
  * the twenty listed non-third-level cores are outside C_3 and the ten listed third-level cores are inside;
  * `is_c3` for every rooted block form mono(t,a,b,c,d,p) of the displayed normal form: 2 permutations x 2^3 quadratic
    coefficients x 4 cubic values x 8 marks = 512 forms, with the count of (form, mark) pairs whose membership differs
    from the unmarked form (mark sensitivity).
Run without -O from the package root:  python src/compact/check_core_levels.py
Writes evidence/compact/core-levels.json (status PASS only if every assertion holds).
"""
import sys
if not __debug__: raise RuntimeError("Run without Python -O: assertions are certificate checks")
from pathlib import Path
from itertools import product
import json, hashlib, time
sys.path.insert(0, str(Path(__file__).resolve().parent))
import compact_core as q

t0 = time.perf_counter()
here = Path(__file__).resolve().parents[2] / 'evidence/compact'
# the thirty cores of the cover table, in the parametrization of compact_core.mono
# mono(t, a, b, c, d, p): t permutation (0 affine, 1 Toffoli z -> z + xy); a,b,c the xy, xz, yz coefficients of i^(...);
# d in Z/4 the cubic coefficient of i^(d xyz) (1,3: S_pm; 2: Q = CCZ); p the mark (1 = x, 2 = y, 4 = z; 0 unmarked)
non_third = []
for sign in (1, 3):
    non_third += [(f'S{sign}_Tz', (0, 0, 0, 0, sign, 4)), (f'piS{sign}_Tx', (1, 0, 0, 0, sign, 1)), (f'piS{sign}_Tz', (1, 0, 0, 0, sign, 4))]
for eta in (0, 1):
    non_third += [(f'piCSyz_Q{eta}_Tx', (1, 0, 0, 1, 2 * eta, 1)), (f'piCSyz_Q{eta}_Ty', (1, 0, 0, 1, 2 * eta, 2)),
                  (f'pi_Q{eta}_Tz', (1, 0, 0, 0, 2 * eta, 4)), (f'piCSyz_Q{eta}_Tz', (1, 0, 0, 1, 2 * eta, 4))]
for sign in (1, 3):
    non_third += [(f'S{sign}', (0, 0, 0, 0, sign, 0)), (f'piS{sign}', (1, 0, 0, 0, sign, 0))]
non_third += [(f'piCSyz_Q{eta}', (1, 0, 0, 1, 2 * eta, 0)) for eta in (0, 1)]
third = []
for eta in (0, 1):
    third += [(f'Q{eta}_Tz', (0, 0, 0, 0, 2 * eta, 4)), (f'CSyz_Q{eta}_Tz', (0, 0, 0, 1, 2 * eta, 4)), (f'CSxy_Q{eta}_Tz', (0, 1, 0, 0, 2 * eta, 4))]
third += [(f'CSyz{c}_Q{eta}', (0, 0, 0, c, 2 * eta, 0)) for c in (0, 1) for eta in (0, 1)]
assert len(non_third) == 20 and sum(1 for n, _ in non_third if not n.startswith(('S', 'piS', 'piCSyz_Q0', 'piCSyz_Q1')) or '_T' in n) == 14
assert len(third) == 10
marked = [n for n, prm in non_third if prm[5]]
unmarked = [n for n, prm in non_third if not prm[5]]
assert len(marked) == 14 and len(unmarked) == 6
records = []
for name, prm in non_third:
    u = q.mono(*prm); c3 = q.is_c3(u)
    assert not c3, name
    records.append({'core': name, 'parameters': list(prm), 'marked': bool(prm[5]), 'third_level': c3})
for name, prm in third:
    u = q.mono(*prm); c3 = q.is_c3(u)
    assert c3, name
    records.append({'core': name, 'parameters': list(prm), 'marked': bool(prm[5]), 'third_level': c3})
# all 512 rooted block forms and the mark sensitivity of C_3 membership
forms, sensitive = [], 0
for t, a, b, c, d in product(range(2), range(2), range(2), range(2), range(4)):
    base = q.is_c3(q.mono(t, a, b, c, d, 0))
    for p in range(8):
        c3 = q.is_c3(q.mono(t, a, b, c, d, p))
        forms.append({'parameters': [t, a, b, c, d, p], 'third_level': c3})
        if p and c3 != base:
            sensitive += 1
assert len(forms) == 512
assert q.is_c3(q.mono(1, 0, 0, 0, 0, 0)) and not q.is_c3(q.mono(1, 0, 0, 0, 0, 4))      # pi versus pi T_z
core = Path(__file__).resolve().parent / 'compact_core.py'
r = {'status': 'PASS', 'listed_cores': records, 'non_third_level': 20, 'marked_non_third_level': 14, 'unmarked_non_third_level': 6,
     'third_level_cores': 10, 'rooted_forms': 512, 'rooted_forms_third_level': sum(f['third_level'] for f in forms),
     'mark_sensitive_pairs': sensitive, 'forms': forms,
     'inputs': {'src/compact/compact_core.py': hashlib.sha256(core.read_bytes()).hexdigest()}, 'elapsed_seconds': round(time.perf_counter() - t0, 2)}
(here / 'core-levels.json').write_text(json.dumps(r, indent=1) + '\n', encoding='utf-8')
print(json.dumps({k: v for k, v in r.items() if k not in ('forms', 'listed_cores')}, indent=1), flush=True)
