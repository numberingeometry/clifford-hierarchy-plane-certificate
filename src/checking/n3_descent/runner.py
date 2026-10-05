"""Runs the exact checks of (C4) and (C5) and validates their verdicts (keys C11 and C12).

Paulis use v = a + 8b, i^(a.b) X^a Z^b, with qubit 0 at the low bit. (C4) uses signed Pauli actions with
their actual relative signs; (C5) uses integer row lattices.
"""
from pathlib import Path
import argparse
import contextlib
import hashlib
import io
import json
import time

if not __debug__:
    raise RuntimeError("Run without -O: exact descent assertions must be enabled.")

FACTS = {
    "C11": "clifford-involution-groups-in-no-pauli-masa",
    "C12": "second-difference-equations-for-phases",
}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def validate_c11(result):
    counts = {
        'lagrangians': 135, 'non_Pauli_rays': 3780, 'Pauli_rays': 63,
        'total_candidate_rays': 3843, 'Clifford_generators': 12,
        'first_orbit': 3780, 'first_stabilizer_generators': 4150,
        'second_candidates': 48, 'second_orbit': 48,
        'third_candidates': 12, 'third_Pauli': 1, 'third_non_Pauli': 11,
        'third_subgroups': 3, 'pauli_subgroup_orbit': 1,
        'nonpauli_subgroup_orbit': 2, 'cubic_difference_matrices': 1024,
        'cubic_difference_rank': 7, 'nonaffine_relative_action_order': 3,
    }
    for field, expected in counts.items():
        need(result[field] == expected, f"C11 {field} changed or domain incomplete")
    need(result['third_by_Pauli_kernel'] == {0: 8, 1: 4}, "C11 third-group partition changed")
    need(result['signed_action_orbit_equalities_checked'] is True and
         result['CCZ_normal_form_label_sets_checked'] is True,
         "C11 signed action or normal-form equality failed")


def validate_c12(result):
    need(result['families'] == [
        {'model': 'N', 'systems': 28, 'directions': [0, 1, 2]},
        {'model': 'P', 'systems': 4, 'directions': [0, 1]},
        {'model': 'P', 'systems': 28, 'directions': [0, 1]},
    ], "C12 phase-system family incomplete")
    need(result['phase_systems'] == 60 and result['target_rows_checked'] == 1200 and
         result['equation_rows_generated'] == 10112, "C12 equation or target domain incomplete")
    need(result['all_integer_lattice_inclusions'] is True, "C12 lattice implication failed")


def run(fact, output_dir=None):
    """Recompute the selected finite fact on its complete domain; optionally save a new report."""
    if not __debug__:
        raise RuntimeError("Run without -O: exact descent assertions must be enabled.")
    aliases = {**{code: code for code in FACTS}, **{name: code for code, name in FACTS.items()}}
    need(fact in aliases, "Unsupported descent fact: " + str(fact))
    code = aliases[fact]
    from . import n3_descent_check as owner
    started = time.perf_counter()
    trace = io.StringIO()
    with contextlib.redirect_stdout(trace):
        if code == "C11":
            result = owner.check_c11()
            validate_c11(result)
        else:
            result = owner.phase_check()
            validate_c12(result)
    here = Path(__file__).parent
    report = {'status': 'accepted', 'route': 'recompute', 'fact': FACTS[code],
              'navigation': code, 'exact': True, 'result': result,
              'elapsed_seconds': time.perf_counter() - started,
              'trace': trace.getvalue(),
              'source_sha256': {name: hashlib.sha256((here / name).read_bytes()).hexdigest()
                                for name in ('runner.py', 'n3_descent_check.py')}}
    if output_dir is not None:
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        with (output / 'result.json').open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fact', choices=list(FACTS) + list(FACTS.values()))
    parser.add_argument('--output-dir')
    args = parser.parse_args()
    print(json.dumps(run(args.fact, args.output_dir), indent=2))


if __name__ == '__main__':
    main()
