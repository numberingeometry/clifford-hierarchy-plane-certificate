# Finite certificates for *When is the Clifford hierarchy generalized semi-Clifford?*

This repository contains the programs that check the finite facts of

> Rongbiao Thomas Wang and Bobby Zixuan Zhang, *When is the Clifford hierarchy generalized semi-Clifford?* (2026).

The paper states 19 finite facts, (C1)–(C19). Each appears in a box where a proof uses it, after the argument that reduces the claim there to a finite computation; that reduction is part of the written proof and is not checked here. For each fact, this repository recomputes the boxed claim from scratch, in exact arithmetic, on its complete finite domain.

## Requirements

- Python 3.11 or later. The release run used CPython 3.13.15.
- NumPy, for the row coverage used by (C6) and (C7), with integer arrays only: `pip install -r requirements.txt`.
- A C compiler, for (C7) and (C15): Clang, GCC, or Zig via `pip install -r requirements-compiler.txt`. Set `CC` or pass `--compiler` to choose one.

Run Python without `-O`, since the checks use assertions.

## Usage

```sh
python check.py paper                          # the frozen statements are consistent
python check.py paper --manuscript main.tex    # and agree with the boxes of a source of the paper
python check.py integrity                      # every file matches provenance/integrity.json
python check.py math --fact C8                 # one check, by number, label or folder name
python check.py math --all --allow-expensive --workers 8
python tools/run_tests.py                      # rejection and control tests
```

`check.py math` writes each accepted result to `runs/<time>-<check>/result.json`. A check is accepted only if its checker returns an exact verdict for the stated claim on the complete domain; any failure stops the run with a nonzero exit code. (C7) is the only expensive check, about 8.2 hours with one worker; `--workers N` distributes its 1,010 gates over N processes.

## The facts of the paper

| Fact | Statement | Folder | Release run |
|---|---|---|---|
| (C1) | Stationary gates $\mathrm{CS}_{yz}^a\mathrm{CCZ}^b$ and complete affine flags | [stationary-cs-ccz-gates-and-complete-affine-flags](certificates/stationary-cs-ccz-gates-and-complete-affine-flags/) | under 1 s |
| (C2) | Three-qubit involutions with a flag-preserving permutation | [three-qubit-flag-preserving-involution-rays](certificates/three-qubit-flag-preserving-involution-rays/) | 4 s |
| (C3) | Cliffords commuting with $H_*$ | [cliffords-commuting-with-hstar](certificates/cliffords-commuting-with-hstar/) | under 1 s |
| (C4) | Clifford involution groups in no Pauli MASA | [clifford-involution-groups-in-no-pauli-masa](certificates/clifford-involution-groups-in-no-pauli-masa/) | 5 s |
| (C5) | Second-difference equations for phases | [second-difference-equations-for-phases](certificates/second-difference-equations-for-phases/) | under 1 s |
| (C6) | Representatives of $\mathscr R_4$ up to Clifford conjugacy | [r4-representatives-up-to-clifford-conjugacy](certificates/r4-representatives-up-to-clifford-conjugacy/) | 5 s |
| (C7) | A plane meeting the fixed Pauli labels | [plane-meeting-fixed-pauli-labels](certificates/plane-meeting-fixed-pauli-labels/) | 84 min (6 workers) |
| (C8) | Pauli terms and the cubic coefficient | [pauli-terms-and-cubic-coefficient](certificates/pauli-terms-and-cubic-coefficient/) | under 1 s |
| (C9) | Normal form of the triples $(\rho,u,D)$ | [normal-form-of-exceptional-balanced-configurations](certificates/normal-form-of-exceptional-balanced-configurations/) | 2 s |
| (C10) | The case $\mathrm{rank}\,E=2$ | [balanced-rank-two-branch](certificates/balanced-rank-two-branch/) | 2 s |
| (C11) | Mixed blocks with two input Lagrangians | [mixed-blocks-with-two-input-lagrangians](certificates/mixed-blocks-with-two-input-lagrangians/) | 12 s |
| (C12) | Clifford classes of the block forms | [block-form-clifford-classes](certificates/block-form-clifford-classes/) | 47 s |
| (C13) | Invariant sets of input Lagrangians | [invariant-sets-of-input-lagrangians](certificates/invariant-sets-of-input-lagrangians/) | 18 s |
| (C14) | Reduction of $B_0$ to $256$ gates | [reduction-to-common-input-pair-domain](certificates/reduction-to-common-input-pair-domain/) | 11 s |
| (C15) | Pair domains for blocks | [block-pair-domains](certificates/block-pair-domains/) | 2 min |
| (C16) | A $2$-group criterion in $S_6$ | [s6-two-group-criterion](certificates/s6-two-group-criterion/) | 2 s |
| (C17) | Sylow cosets in $S_6$ and three-element subsets | [s6-sylow-cosets-and-three-element-subsets](certificates/s6-sylow-cosets-and-three-element-subsets/) | under 1 s |
| (C18) | The gate $U_{E_6}$ lies in $\mathcal{C}_5(5)$ | [e6-gate-level-five](certificates/e6-gate-level-five/) | 31 s |
| (C19) | Root words for $E_6$ | [e6-root-words](certificates/e6-root-words/) | under 1 s |

Each folder contains `statement.tex`, `spec.json` (finite domain, dependencies, checker, command and hashes of the checker sources) and a README describing the check.

## Release run

All 19 checks were accepted in one run of `python check.py math --all --allow-expensive --workers 6` on a fresh copy of this repository, at paper revision `fbfe343`, in 89 min, with CPython 3.13.15 on Windows 11 (8 cores). (C2), (C4), (C6), (C10), (C18) were run again at `830fed2`, the revision of the frozen statements, after changes of wording or typography in their statements; the other 14 statements are identical at the two revisions. The record of each check is in [results/](results/): it contains the statement hash, the paper revision, the hashes of the checker sources and the verdict. For (C7), the census summary in [results/plane-meeting-census/](results/plane-meeting-census/) lists the 1,010 gates, all computed afresh, with no failure. The 85 tests, including the fixture controls of (C15) on the workspace of that run, passed; see [results/tests/](results/tests/).

## Scope

- The statements in [paper/finite-certificate.tex](paper/finite-certificate.tex) are those of revision `830fed2` of the paper. `check.py paper --manuscript` compares them with a source of the paper after removing comments and colour and normalizing whitespace, so any change of a statement is detected.
- Verdicts use exact arithmetic only: integers, binary matrices and exact monomials over cyclotomic integers.
- The programs are not formally verified. The tests in `tests/` include rejection controls: altered inputs that the checkers must reject.
- Every verdict is recomputed; no earlier verdict or count is reused. Where a check uses supplied witnesses (the conjugators of (C2), the affine conjugators and phase forests of (C6)), each witness is verified exactly. The route `--route verify_evidence` exists only for (C6).

## Layout

- `check.py`: command-line driver, statement comparison and integrity check.
- `certificates/<fact>/`: one folder per fact.
- `src/checking/`: the checkers. Each `spec.json` names its checker (`backend`) and an internal key.
- `four-qubit/`: the sources and input data of the four-qubit computations used by (C6), (C7) and (C10)–(C15), including the C kernels of (C7) and (C15). Its files are listed with their SHA-256 hashes in `four-qubit/manifest.json` and checked before every run.
- `paper/`: the frozen statements, the paper revision (`pin.json`) and the conventions (`conventions.json`).
- `provenance/`: the history of statement changes and the integrity manifest.
- `tests/`, `tools/`: tests and maintenance scripts (`tools/sync_paper.py` updates the statements from a source of the paper; `tools/freeze_manifest.py` writes the integrity manifest).

## Conventions

Pauli labels are v = (x, z) in F_2^{2n} with P_v = i^{x·z} X^x Z^z, qubit 0 is the least significant bit, and products act from right to left. Details are in `paper/conventions.json` and in the docstrings of the checkers.

## Citation and license

Please cite the paper. `CITATION.cff` describes this repository. The code is released under the MIT License.
