"""Fact (C18): the gate U_{E6} lies in C_5(5). Exact recomputation.

U_{E6} = sum_{c,t} |ct><ct| (x) A^c B^t, with A = R_{p1} ... R_{p5} and B = R_{p0} R_{p2} the products of root
lifts of Appendix G. The hierarchy recursion of e2_hierarchy.py builds this gate from its root words over
Z[i,1/2] and decides its membership in C_5 exactly: all 1,024 Pauli conjugations at levels 5 and 4, and the
ten generators at levels 3 and 2. It also records the failure of membership in C_4 with an independently
verified chain. Every Clifford block used by the recursion is then audited against direct matrix arithmetic.
"""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

from . import e2_hierarchy as hierarchy

FACT = 'e6-gate-level-five'
WORDS = {'E6': [[1, 24, 2, 48, 6], [8, 24]]}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def count(value, expected, message):
    need(type(value) is int and value == expected, message)


def write(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


def validate_hierarchy(report):
    need(type(report) is dict and report.get('status') == 'HIERARCHY_PASS'
         and report.get('exact') is True, 'hierarchy verdict missing or inexact')
    results = report.get('results')
    need(type(results) is dict and set(results) == {'E6'}, 'hierarchy domain is not the gate U_E6')
    item = results['E6']
    need(type(item) is dict and item.get('in_C5') is True and item.get('in_C4') is False
         and item.get('exact') is True, 'U_E6 exact level-five conclusion missing')
    need(item.get('lower_failure_chain_independently_verified') is True, 'lower-level failure chain not verified')
    chain = item.get('lower_failure_chain')
    need(type(chain) is dict and chain.get('level') == 4, 'level-four failure-chain record missing or misbound')
    count(item.get('complete_root_Pauli_count'), 1024, 'root Pauli domain incomplete')
    count(item.get('full_Pauli_domain_at_levels_at_least_four'), 1024, 'recursive Pauli domain incomplete')
    count(item.get('generator_count_at_levels_two_and_three'), 10, 'lower group-generator domain changed')
    need(item.get('root_words') == WORDS['E6'] and item.get('control_block_order') == ['I', 'A', 'B', 'AB'],
         'root words or coherent block order changed')
    need(item.get('unitarity_checked_exactly') is True and item.get('root_lifts_square_to_declared_Paulis') is True,
         'exact construction binding missing')
    coverage = item.get('upper_coverage')
    need(type(coverage) is dict and coverage.get('every_distinct_upper_C4_node_checked_all_1024_Paulis') is True
         and coverage.get('every_upper_C3_and_C2_node_checked_all_ten_generators') is True,
         'complete upper hierarchy coverage missing')

    def levels(field):
        values = coverage.get(field)
        need(type(values) is dict, 'hierarchy coverage counters missing')
        out = {}
        for key, value in values.items():
            need(type(key) in (int, str) and key in (1, 2, 3, 4, 5, '1', '2', '3', '4', '5'),
                 'invalid hierarchy counter level')
            level = int(key)
            need(level not in out and type(value) is int and value > 0, 'duplicate or malformed hierarchy counter')
            out[level] = value
        return out
    visits = levels('distinct_memoized_operators_by_level')
    checks = levels('computed_conjugations_by_level')
    need(set(visits) == {1, 2, 3, 4, 5} and set(checks) == {2, 3, 4, 5},
         'upper hierarchy counter-level domain incomplete or changed')
    count(visits.get(5), 1, 'upper root domain changed')
    count(checks.get(5), 1024, 'upper root coverage incomplete')
    count(visits.get(4), 1024, 'level-four descendant domain incomplete')
    count(checks.get(4), 1048576, 'recursive level-four coverage incomplete')
    for level in (2, 3):
        need(level in visits, 'lower group domain omitted')
        count(checks.get(level), 10*visits[level], 'lower group domain incomplete')


def validate_block_audit(report):
    need(type(report) is dict and report.get('status') == 'EXACT_BLOCK_AUDIT_PASS'
         and report.get('exact') is True, 'compressed block arithmetic lacks exact dense audit')
    blocks = report.get('blocks_checked')
    need(type(blocks) is int and blocks > 0 and report.get('all_registered_blocks_covered') is True,
         'compressed block audit incomplete')
    count(report.get('Pauli_action_domain_per_block'), 64, 'block Pauli-action domain incomplete')
    count(report.get('signed_Pauli_images_checked'), 64*blocks, 'signed Pauli cache audit incomplete')
    count(report.get('unitarity_blocks_checked'), blocks, 'block unitarity audit incomplete')


def run(workspace=None):
    need(not sys.flags.optimize, 'the check must run with assertions enabled')
    if workspace is None:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        workspace = Path(__file__).resolve().parents[2]/'runs'/('direct-'+stamp+'-C18')/'workspace'
    work = Path(workspace)
    work.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    report = hierarchy.run_hierarchy(output_dir=work/'hierarchy')
    validate_hierarchy(report)
    gate = hierarchy.gates()['E6']
    need(hierarchy.verify_failure(gate, report['results']['E6']['lower_failure_chain'], 4) is True,
         'lower failure-chain semantic verifier did not accept')
    arithmetic = hierarchy.audit_exact_blocks(output_dir=work/'arithmetic')
    validate_block_audit(arithmetic)
    item = report['results']['E6']
    result = {'status': 'accepted', 'fact': FACT, 'route': 'recompute', 'exact': True,
              'gate': 'U_E6', 'root_words': item['root_words'], 'control_block_order': item['control_block_order'],
              'in_C5': True, 'in_C4': False, 'qubits': 5, 'field': 'Q(i)', 'entry_ring': 'Z[i,1/2]',
              'conventions': 'paper/conventions.json; labels x+(z<<5); c,t,d0,d1,d2 from least to most significant; '
                             'actual UP_vUdagger; only a common scalar is quotiented',
              'hierarchy': report, 'exact_block_audit': arithmetic,
              'seconds': time.perf_counter()-started, 'formal_verification': False}
    write(work/'C18-check.json', result)
    return result
