"""Runs the exact checks of (C8), (C14) and (C2) and validates their verdicts (keys C9, C10, C13).

All gates use x = x0 + 2 x1 + 4 x2 with qubit 0 at the low bit, phases modulo 8, Hermitian
P_(a+8b) = i^(a.b) X^a Z^b, and products with the right factor applied first.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

if not __debug__:
    raise RuntimeError("Run without -O: exact audit assertions must be enabled.")

FACTS = {
    "C9": "pauli-terms-and-cubic-coefficient",
    "C10": "reduction-to-common-input-pair-domain",
    "C13": "three-qubit-flag-preserving-involution-rays",
}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def validate_c9(result):
    need(result.get("holds") is True, "C9 mathematical conclusion failed")
    a, b = result["a"], result["b"]
    need((a["gates"], a["paulis_per_gate"], a["cases"]) == (32, 56, 1792), "C9 odd domain incomplete")
    need(a["holds"] is True and a["minimum_terms"] == 8 and a["minimum_exactly_eight"] is True,
         "C9 odd Pauli-support bound failed")
    need(a["histogram"] == {"8": 1024, "10": 768}, "C9 odd support census changed")
    need((b["forms"], b["labels_per_form"]) == (256, 64), "C9 even domain incomplete")
    need(b["holds"] is True and not b["non_subspace_forms"] and b["minimum_dimension"] == 5,
         "C9 even-label subspace conclusion failed")
    need(b["dimension_histogram"] == {"5": 112, "6": 144}, "C9 even dimension census changed")


def validate_c10(result):
    need(result.get("holds") is True and result["cores_all_non_third_level"] is True,
         "C10 mathematical conclusion failed")
    one, two, three, four = (result[f"claim{i}"] for i in range(1, 5))
    need(all(item["holds"] is True for item in (one, two, three, four)), "C10 subclaim failed")
    need((one["affine"], one["toffoli_type"], one["total"]) == (1344, 9408, 10752),
         "C10 permutation domain incomplete")
    need(one["every_one_in_a_left_affine_coset"] is True and one["cosets_disjoint"] is True,
         "C10 affine-coset classification failed")
    by_direction = two["by_radical_direction(0=identity; k=x+2y+4z)"]
    need(set(map(int, by_direction)) == set(range(8)), "C10 radical directions incomplete")
    for raw_key, value in by_direction.items():
        key = int(raw_key)
        need(value == {"classes": 32, "listed_classes_hit": 16 if key == 7 else 32},
             "C10 phase quotient census changed")
    need(two["complete_directions"] == list(range(7)) and
         {int(k): v for k, v in two["incomplete_directions"].items()} == {7: 16},
         "C10 exceptional direction lost")
    need((three["direction"], three["classes"], three["cores"], three["counterexamples"]) ==
         ([1, 1, 1], 32, 20, []), "C10 exceptional exclusion incomplete")
    need((four["domain"], four["qualifying_pairs"]) == (5120, 360), "C10 pair domain changed")


def ray_key(permutation, phase):
    phase = tuple(phase)
    return tuple(permutation), min(phase, tuple((v + 4) % 8 for v in phase))


def validate_c13_conjugators(supplied, generated, H):
    """Validate supplied conjugators against the regenerated exceptional domain.

    This is part of the recomputation: the complete ray domain is regenerated, and the
    supplied conjugators are checked against it. It is not a separate verifier of supplied evidence.
    """
    from . import mono3 as q
    need(isinstance(supplied, list) and len(supplied) == 256, "C13 supplied conjugator count")
    need(len(generated) == 256, "C13 regenerated exceptional-ray count")
    generated_keys = {ray_key(item["pi"], item["phase_zeta"]) for item in generated}
    need(len(generated_keys) == 256, "C13 regenerated exceptional rays duplicated")
    supplied_keys = set()
    targets = {H, (H[0], tuple((v + 4) % 8 for v in H[1]))}
    for item in supplied:
        need(set(item) == {"pi", "phase_zeta", "alpha", "E_phase_zeta"}, "C13 conjugator record fields")
        p, phase = item["pi"], item["phase_zeta"]
        alpha, E = item["alpha"], item["E_phase_zeta"]
        need(all(isinstance(row, list) and len(row) == 8 for row in (p, phase, alpha, E)),
             "C13 conjugator dimensions")
        need(all(type(v) is int for row in (p, phase, alpha, E) for v in row),
             "C13 conjugator integer fields")
        need(sorted(p) == list(range(8)) and tuple(alpha) in q.AFFINE, "C13 conjugator permutations")
        need(all(0 <= v < 8 for v in phase + E), "C13 conjugator phase range")
        key = ray_key(p, phase)
        need(key in generated_keys and key not in supplied_keys, "C13 omitted or duplicate exceptional ray")
        supplied_keys.add(key)
        M = tuple(p), tuple(phase)
        C = q.mul(q.diag(E), q.perm(alpha))
        need(q.is_clifford(C), "C13 supplied conjugator is not Clifford")
        need(q.conj(C, M) in targets, "C13 supplied conjugation equality failed")
    need(supplied_keys == generated_keys, "C13 supplied exceptional domain incomplete")
    return {"supplied_conjugators_checked": 256, "unique_exceptional_rays": 256,
            "Cliffordness_and_exact_equality": True,
            "coverage_basis": "the complete exceptional-ray domain is regenerated by the recomputation"}


def validate_c13(result, generated, H, supplied):
    from . import mono3 as q
    expected = {
        "8": {"pi": 1, "rays": 35, "in_C3": 35, "conj_to_Hstar": 0},
        "6": {"pi": 4, "rays": 160, "in_C3": 160, "conj_to_Hstar": 0},
        "4": {"pi": 10, "rays": 480, "in_C3": 480, "conj_to_Hstar": 0},
        "2": {"pi": 12, "rays": 768, "in_C3": 512, "conj_to_Hstar": 256},
        "0": {"pi": 17, "rays": 4352, "in_C3": 4352, "conj_to_Hstar": 0},
    }
    need(result.get("holds") is True and result["table_matches_manuscript"] is True,
         "C13 ray classification failed")
    need(result["flag_involutions"] == 44 and result["table"] == expected,
         "C13 full ray census changed")
    need(result["non_C3_rays"] == result["conjugators_found"] == 256,
         "C13 exceptional partition incomplete")
    need(result["Hstar"] == {"in_C3": False, "in_C4": True}, "C13 Hstar level control failed")
    need(not q.is_c3(H) and all(q.is_c3(q.conj(H, P)) for P in q.PAULIS),
         "C13 Hstar C4 membership failed on the full 64-Pauli domain")
    result["Hstar"]["full_Pauli_domain_checked"] = 64
    result["conjugator_validation"] = validate_c13_conjugators(supplied, generated, H)


def run(fact, output_dir=None):
    """Recompute a complete fact and optionally write a new result.json."""
    if not __debug__:
        raise RuntimeError("Run without -O: exact audit assertions must be enabled.")
    aliases = {**{code: code for code in FACTS}, **{name: code for code, name in FACTS.items()}}
    need(fact in aliases, "Unsupported exact-audit fact: " + str(fact))
    code = aliases[fact]
    started = time.perf_counter()
    if code == "C9":
        from . import check_c9
        result = check_c9.recompute()
        validate_c9(result)
        checker = "check_c9.py"
    elif code == "C10":
        from . import check_c10
        result = check_c10.recompute()
        validate_c10(result)
        checker = "check_c10.py"
    else:
        from . import check_n3_exceptional
        result, generated, H = check_n3_exceptional.recompute()
        supplied = json.loads((Path(__file__).parent / "data/n3-exceptional-conjugators.json").read_text(encoding="utf-8"))
        validate_c13(result, generated, H, supplied)
        checker = "check_n3_exceptional.py"
    # The statements are bound to the paper through paper-map.json, not through the checkers.
    result.pop("manuscript")
    owner = Path(__file__).parent
    source_names = ["runner.py", "mono3.py", checker]
    if code == "C13":
        source_names.append("data/n3-exceptional-conjugators.json")
    report = {"status": "accepted", "route": "recompute", "fact": FACTS[code],
              "navigation": code, "exact": True, "result": result,
              "elapsed_seconds": time.perf_counter() - started,
              "source_sha256": {name: hashlib.sha256((owner / name).read_bytes()).hexdigest()
                                 for name in source_names}}
    if output_dir is not None:
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        with (output / "result.json").open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fact", choices=list(FACTS) + list(FACTS.values()))
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    print(json.dumps(run(args.fact, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
