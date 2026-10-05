"""The claims of (C11), (C12) and (C13) that the package generators do not state, checked on their output.

Only records generated in the temporary workspace are read; no saved count is accepted as a result. The
checks reuse the exact monomial and F_2 arithmetic of the package, so they are not an independent
implementation. Keys as in conventional.py: C3 = (C12), C4 = (C11), C5 = (C13).
"""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sys


class CheckFailure(ValueError):
    pass


def require(condition,message):
    if not condition:
        raise CheckFailure(message)


def read_record(path):
    return json.loads(path.read_text(encoding="utf-8"))


def audit(code,workspace):
    require(__debug__,"run without -O: archived exact primitives use assertions")
    sys.path.insert(0,str(workspace/"src/compact"))
    import compact_core as q
    out=workspace/"evidence/compact"
    levels=read_record(out/"core-levels.json")
    require(levels.get("status")=="PASS","fresh core-level computation failed")
    forms=levels["forms"]
    require(len(forms)==512 and {tuple(r["parameters"]) for r in forms}
            ==set(product(range(2),range(2),range(2),range(2),range(4),range(8))),"rooted form coverage is incomplete")
    cores=levels["listed_cores"]
    require(len(cores)==30 and len({tuple(r["parameters"]) for r in cores})==30,"canonical core domain is incomplete")
    non_third=[r for r in cores if not r["third_level"]]
    third=[r for r in cores if r["third_level"]]
    require(len(non_third)==20 and len(third)==10,"canonical core level partition changed")
    for record in cores:
        require(q.is_c3(q.mono(*record["parameters"])) is record["third_level"],"fresh exact core classification mismatch")
    require(len(q.LAGS)==135 and len(set(q.LAGS))==135,"residual Lagrangian domain is incomplete")
    require(all(len(L)==8 and len(q.BASES[i])==3 and q.span(q.BASES[i])==L
                and not any(q.om(a,b) for a in L for b in L) for i,L in enumerate(q.LAGS)),"invalid residual Lagrangian")

    if code=="C3":
        cover=read_record(out/"block-cover.json")
        edges=read_record(out/"block-edge-audit.json")
        require(cover.get("status")==edges.get("status")=="PASS","fresh cover or exact edge audit failed")
        rows=cover["forms"]
        require(len(rows)==512 and {tuple(r["parameters"]) for r in rows}=={tuple(r["parameters"]) for r in forms},"cover omitted rooted forms")
        orbits=cover["orbits"]
        require(len(orbits)==30 and sum(r["states"] for r in orbits)==512,"cover orbits do not partition the domain")
        require(all(len(r["cores"])==1 for r in orbits),"a cover orbit has no unique listed core")
        actual=derivatives=0
        for record in non_third:
            A=q.mono(*record["parameters"])
            for P in q.PAULIS:
                correction=q.conj(A,P)
                require(q.is_c3(correction),"a non-third-level core has a correction outside C3")
                actual+=1
                right_derivative=q.mul(correction,q.inv(P))
                require(q.is_c3(right_derivative),"a non-third-level core has a right derivative outside C3")
                derivatives+=1
        # The third-level classes meeting the forms with p'=0 are those of CS_yz^a CCZ^b, that is,
        # mono(0,0,0,a,2b,0): the coefficient c multiplies yz. Orbits of Clifford moves lie in classes.
        orbit_of={tuple(r["parameters"]):r["orbit"] for r in rows}
        gates=[(0,0,0,a,2*b,0) for a in range(2) for b in range(2)]
        require(all(q.is_c3(q.mono(*g)) for g in gates),"a gate CS_yz^a CCZ^b is not third-level")
        gate_classes={orbit_of[g] for g in gates}
        require(len(gate_classes)==4,"the gates CS_yz^a CCZ^b do not lie in four distinct classes")
        zero_mark=[p for p in orbit_of if p[5]==0]
        require(len(zero_mark)==64,"the forms with p'=0 are incomplete")
        zero_mark_third=[p for p in zero_mark if q.is_c3(q.mono(*p))]
        require({orbit_of[p] for p in zero_mark_third}==gate_classes,
                "the third-level classes meeting p'=0 are not those of CS_yz^a CCZ^b")
        return {"status":"accepted","fact":"block-form-clifford-classes","rooted_forms":len(forms),
                "zero_mark_forms_checked":len(zero_mark),"zero_mark_third_level_forms":len(zero_mark_third),
                "zero_mark_third_level_classes":len(gate_classes),
                "cover_move_orbits":len(orbits),"listed_third_level_cores":len(third),"listed_non_third_level_cores":len(non_third),
                "actual_corrections_checked":actual,"right_derivatives_checked":derivatives,
                "normalized_edges_checked":edges["applicable_edges_including_self"],
                "odd_scalar_edges":edges["odd_scalar_edges"],"transported_mark_changes":edges["mark_coordinate_changes"],
                "analytic_scope":"exact orbit cover under the proof's licensed Clifford moves; the paper's normal-form/class reduction remains an analytic premise"}

    data=read_record(out/"third-level.json")
    defs=data["conjugate_groups"]
    require(len(defs)==1396 and len({tuple(tuple(e) for e in d["group"]) for d in defs})==1396,"conjugate-group domain is incomplete or duplicated")
    chartsets={tuple(d["charts"]) for d in defs}
    require(len(chartsets)==451,"distinct input-Lagrangian set count changed")
    lookup={L:i for i,L in enumerate(q.LAGS)}

    @lru_cache(None)
    def lag_action(matrix):
        require(len(matrix)==6,"invalid symplectic matrix dimension")
        require(all(q.om(matrix[i],matrix[j])==q.om(1<<i,1<<j) for i in range(6) for j in range(6)),"matrix does not preserve the symplectic form")
        permutation=tuple(lookup[frozenset(q.app(matrix,v) for v in L)] for L in q.LAGS)
        require(len(set(permutation))==135,"symplectic Lagrangian action is not bijective")
        return permutation

    for d in defs:
        E=tuple(tuple(e) for e in d["group"])
        require(len(E) in (1,4,8) and q.I in E,"covered correction group has invalid order or lacks identity")
        require(all(q.prod(e,f) in E for e in E for f in E),"covered correction group is not closed")
        require(all(q.inverse(e) in E for e in E),"covered correction group omits inverses")
        actual_charts=tuple(i for i in range(135) if all(lag_action(e)[i]==i for e in E))
        require(actual_charts==tuple(d["charts"]),"generated input charts do not match the exact correction group")
        fixed=tuple(v for v in range(1,64) if all(q.app(e,v)==v for e in E))
        require(fixed==tuple(d["fixed"]),"generated returned labels do not match exact fixed vectors")
    # Completeness: the recorded groups are exactly the groups E_D of the Clifford conjugates D of
    # CS_yz^a CCZ^b. Conjugating D by C conjugates E_D by gamma_C, and the transvections
    # t_v(x)=x+omega(x,v)v, which are involutions, generate Sp(6,2).
    recorded={tuple(tuple(e) for e in d["group"]) for d in defs}
    transvections=[tuple((1<<j)^(v if q.om(1<<j,v) else 0) for j in range(6)) for v in range(1,64)]
    require(all(q.prod(t,t)==q.I for t in transvections),"a transvection is not an involution")
    orbit=set()
    frontier=[]
    for a in range(2):
        for b in range(2):
            E=tuple(sorted(tuple(e) for e in q.m.group(q.mono(0,0,0,a,2*b,0))))
            if E not in orbit:
                orbit.add(E)
                frontier.append(E)
    while frontier:
        E=frontier.pop()
        for t in transvections:
            F=tuple(sorted(q.prod(q.prod(t,e),t) for e in E))
            if F not in orbit:
                require(F in recorded,"a conjugate group E_D is missing from the recorded groups")
                orbit.add(F)
                frontier.append(F)
    require(orbit==recorded,"a recorded group is not a group E_D of a Clifford conjugate of CS_yz^a CCZ^b")
    bases=[("I",q.mono(0,0,0,0,0,0)),("CS",q.mono(0,1,0,0,0,0)),
           ("CCZ",q.mono(0,0,0,0,2,0)),("CS_CCZ",q.mono(0,1,0,0,2,0))]
    base_checks=[]
    for name,D in bases:
        require(q.is_c3(D),"a named base gate is not third-level")
        E=tuple(q.m.group(D))
        exact_input=tuple(i for i,j in q.charts(D))
        expected=tuple(i for i in range(135) if all(lag_action(e)[i]==i for e in E))
        require(exact_input==expected and all(i==j for i,j in q.charts(D)),"base input charts fail correction-group classification or stationarity")
        fixed={v for v in range(1,64) if all(q.app(e,v)==v for e in E)}
        pauli_rows=q.rows(D)
        returned={v for v in range(1,64) if q.one_label(pauli_rows[v]) is not None}
        require(fixed==returned and all(q.one_label(pauli_rows[v])==v for v in fixed),"base Pauli-return characterization failed")
        base_checks.append({"base":name,"input_lagrangians":len(exact_input),"returned_nonzero_labels":len(returned),"correction_group_order":len(E)})

    common={"residual_lagrangians":135,"linear_changes":len(q.LINEARS),"base_representatives":len(bases),
            "covered_conjugate_groups":len(defs),"conjugation_closed_groups":len(orbit),"distinct_input_lagrangian_sets":len(chartsets),
            "exact_base_classification":base_checks,
            "analytic_scope":"the recorded groups are checked to be exactly the Sp(6,2)-conjugates of the groups of CS_yz^a CCZ^b; Clifford equivariance, X(C^dagger D C)=gamma_C^{-1}X(D) and E_{C^dagger D C}=gamma_C^{-1}E_D gamma_C, transfers the exact chart classification of the four gates to them"}
    if code=="C4":
        spectrum=[]
        for parameters in product(range(2),range(2),range(2),range(2),range(4)):
            D=q.mono(*parameters,0)
            charts=q.charts(D)
            third_level=q.is_c3(D)
            require(charts,"an unmarked form has no input Lagrangian")
            require(len(charts)<2 or third_level,"multichart unmarked form is not third-level")
            spectrum.append({"parameters":list(parameters),"input_lagrangians":len(charts),"third_level":third_level})
        return {**common,"status":"accepted","fact":"mixed-blocks-with-two-input-lagrangians",
                "unmarked_forms_checked":len(spectrum),"classification_failures":0,"unmarked_chart_spectrum":spectrum}

    require(code=="C5","unsupported conventional fact")
    lz=lookup[frozenset(range(0,64,8))]
    groups=set()
    secondary_records=[]
    pairs=actual=invariant=0
    for record in non_third:
        A=q.mono(*record["parameters"])
        generators=set()
        for P in q.PAULIS:
            correction=q.conj(A,P)
            require(q.is_c3(correction),"secondary-group derivative premise failed")
            actual+=1
            generators.update(q.m.group(correction))
        J={q.I}
        frontier=[q.I]
        while frontier:
            h=frontier.pop()
            for g in generators:
                gh=q.prod(g,h)
                if gh not in J:
                    J.add(gh)
                    require(len(J)<=64,"secondary group exceeds its declared finite bound")
                    frontier.append(gh)
        require(len(J) in (8,64),"secondary group order changed")
        require(all(q.prod(g,h) in J for g in generators for h in J),"secondary-group closure is incomplete")
        local_invariant=0
        for d in defs:
            pairs+=1
            X=tuple(d["charts"])
            sx=set(X)
            if all({lag_action(g)[i] for i in X}==sx for g in J):
                require(lz in sx,"an arising invariant input-Lagrangian set omits L_A")
                invariant+=1
                local_invariant+=1
        groups.add(tuple(sorted(J)))
        secondary_records.append({"core":record["core"],"parameters":record["parameters"],
                                  "secondary_group_order":len(J),"invariant_chart_domains":local_invariant})
    distribution=Counter((r["secondary_group_order"],r["invariant_chart_domains"]) for r in secondary_records)
    require(distribution==Counter({(64,30):10,(8,32):10}),"secondary-group invariant-domain distribution changed")
    return {**common,"status":"accepted","fact":"invariant-sets-of-input-lagrangians",
            "non_third_level_cores":len(non_third),"actual_corrections_checked":actual,
            "core_group_pairs_checked":pairs,"invariant_pairs":invariant,"distinct_secondary_groups":len(groups),
            "implication_failures":0,"secondary_group_records":secondary_records}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code",required=True,choices=("C3","C4","C5"))
    parser.add_argument("--workspace",required=True,type=Path)
    args=parser.parse_args()
    print(json.dumps(audit(args.code,args.workspace),sort_keys=True))
