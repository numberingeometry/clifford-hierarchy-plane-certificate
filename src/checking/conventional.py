"""Facts (C11), (C12) and (C13): exact recomputation with the scripts of the four-qubit package.

The scripts use assertions as checks, so optimized Python (-O) is rejected. Their sources are copied
into a temporary workspace and run there; the package folder is never written. The generator of
the conjugate groups is run afresh for (C11) and (C13), without the symplectic census of (C15) or the
census of (C7). Internal keys: C3 = (C12), C4 = (C11), C5 = (C13).
"""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


class CheckFailure(ValueError):
    pass


def require(condition,message):
    if not condition:
        raise CheckFailure(message)


SOURCE_FILES = (
    "src/completion_core.py",
    "src/compact/compact_core.py",
    "src/compact/check_core_levels.py",
    "src/compact/check_block_cover.py",
    "src/compact/build_third_level.py",
    "scripts/audit_block_edges.py",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def default_package():
    return Path(__file__).resolve().parents[2]/"four-qubit"


def execute(command,workspace):
    completed = subprocess.run([sys.executable,"-B",*map(str,command)],cwd=workspace,
                               capture_output=True,text=True,encoding="utf-8")
    require(completed.returncode==0,
            "exact source command failed: "+" ".join(map(str,command))+"\n"+completed.stdout[-12000:]+completed.stderr[-12000:])
    return completed.stdout


def validate_summary(code,result):
    require(result.get("status")=="accepted","conventional mathematical audit did not accept")
    expected = {
        "C3":{"rooted_forms":512,"cover_move_orbits":30,"listed_third_level_cores":10,
              "listed_non_third_level_cores":20,"actual_corrections_checked":1280,
              "right_derivatives_checked":1280,"normalized_edges_checked":52320,
              "odd_scalar_edges":768,"transported_mark_changes":36096,"zero_mark_forms_checked":64,
              "zero_mark_third_level_forms":20,"zero_mark_third_level_classes":4},
        "C4":{"residual_lagrangians":135,"linear_changes":168,"base_representatives":4,
              "covered_conjugate_groups":1396,"conjugation_closed_groups":1396,"distinct_input_lagrangian_sets":451,
              "unmarked_forms_checked":64,"classification_failures":0},
        "C5":{"residual_lagrangians":135,"covered_conjugate_groups":1396,"conjugation_closed_groups":1396,
              "distinct_input_lagrangian_sets":451,"non_third_level_cores":20,
              "actual_corrections_checked":1280,"core_group_pairs_checked":27920,
              "invariant_pairs":620,"distinct_secondary_groups":2,"implication_failures":0},
    }
    require(code in expected,"unknown conventional certificate")
    for key,value in expected[code].items():
        require(result.get(key)==value,f"{code} {key}: expected {value}, received {result.get(key)!r}")


def run(code,package_root=None):
    require(__debug__,"optimized Python disables the assertions the checks rely on; run without -O")
    require(code in ("C3","C4","C5"),"certificate must be C3, C4 or C5")
    started = time.perf_counter()
    package = Path(package_root) if package_root is not None else default_package()
    require(package.is_dir(),"the four-qubit package is unavailable")
    for name in SOURCE_FILES:
        require((package/name).is_file(),"missing package source "+name)
    bindings = {name:digest(package/name) for name in SOURCE_FILES}
    helper = Path(__file__).with_name("conventional_audit.py")
    require(helper.is_file(),"maintained conventional mathematical audit is missing")
    helper_hash = digest(helper)
    with tempfile.TemporaryDirectory(prefix="paper_i_conventional_") as temporary:
        workspace = Path(temporary)
        for name in SOURCE_FILES:
            target = workspace/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(package/name,target)
            require(digest(target)==bindings[name],"source copy mismatch for "+name)
        (workspace/"evidence/compact").mkdir(parents=True,exist_ok=True)
        execute([workspace/"src/compact/check_core_levels.py"],workspace)
        if code=="C3":
            execute([workspace/"src/compact/check_block_cover.py"],workspace)
            # Run the audit script from the package folder: it refuses to use its own folder as the workspace.
            execute([package/"scripts/audit_block_edges.py","--workspace",workspace],workspace)
        else:
            execute([workspace/"src/compact/build_third_level.py"],workspace)
        stdout = execute([helper,"--workspace",workspace,"--code",code],workspace)
        result = json.loads(stdout)
        validate_summary(code,result)
        require(all(digest(package/name)==value for name,value in bindings.items()),"package sources changed during recomputation")
        require(all(digest(workspace/name)==value for name,value in bindings.items()),"working source bytes changed during recomputation")
        require(digest(helper)==helper_hash,"maintained audit changed during recomputation")
        result.update({"code":code,"route":"recompute","input_source_sha256":bindings,
                       "maintained_audit_sha256":helper_hash,"seconds":time.perf_counter()-started,
                       "packed_conventions":"qubit x is bit0, y bit1, z bit2; label v=x_bits+8*z_bits; Hermitian Pauli phase i^(x dot z); products act right to left",
                       "implementation_independence":"fresh complete recomputation with archived exact primitives; additional maintained obligation checks share those primitives"})
        return result


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("code",choices=("C3","C4","C5"))
    parser.add_argument("--package-root",type=Path)
    args=parser.parse_args()
    print(json.dumps(run(args.code,args.package_root),sort_keys=True))
