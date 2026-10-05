# (C2) Three-qubit involutions with a flag-preserving permutation

Fact (C2) of the paper, label `cert:three-qubit-flag-preserving-involution-rays`, stated in a box at its first use (used in `gauss-lem:exceptional`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $M$ be a gate of the form~\eqref{core:eq:row} on three qubits whose permutation preserves the flag $0<\Span\{e_z\}<\Span\{e_y,e_z\}<\F_2^3$, and let $H_*$ be as in the proof of Proposition~\ref{gauss-lem:exceptional}. Then $[M]\in\C_3(3)$, or $[M]$ is conjugate to $[H_*]$ by an affine--diagonal Clifford gate.
```

## Check

```sh
python check.py math --fact three-qubit-flag-preserving-involution-rays
```

- Finite domain: The 44 involutive permutations preserving the flag and the 5,795 classes [M] allowed by the involution and trace conditions.
- Checked: Clifford membership of each class on six generating Paulis; for each of the 256 classes outside C_3(3), the conjugator to [H_*]; the counts of the table in the paper.
- Runtime: 4 s.
- Implementation: `src/checking/exact_audit/`.
- Release run: accepted in 4 s; record in [results/three-qubit-flag-preserving-involution-rays.json](../../results/three-qubit-flag-preserving-involution-rays.json).
