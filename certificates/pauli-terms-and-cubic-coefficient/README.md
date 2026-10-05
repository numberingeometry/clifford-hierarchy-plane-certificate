# (C8) Pauli terms and the cubic coefficient

Fact (C8) of the paper, label `cert:pauli-terms-and-cubic-coefficient`, stated in a box at its first use (used in `app:return-single`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $D=\pi_0\,i^{axy+bxz+cyz+dxyz}$ with $\pi_0\in\{I,\pi\}$, $a,b,c\in\F_2$ and $d\in\Z/4$. (a) If $d$ is odd, then $D^\dagger PD$ has at least eight Pauli terms for every non-diagonal three-qubit Pauli $P$. (b) If $d$ is even, then for every $p\in\F_2^3$ the labels $v$ with $DT_pP_vT_p^\dagger D^\dagger\in\C_2(3)$ form a subspace of dimension at least five.
```

## Check

```sh
python check.py math --fact pauli-terms-and-cubic-coefficient
```

- Finite domain: The 32 forms with d odd times the 56 non-diagonal Paulis; the 256 forms DT_p with d even times the 64 Pauli labels.
- Checked: (a) At least eight Pauli terms in every case; (b) for every p, the labels v with DT_pP_vT_p^dagger D^dagger Clifford form a subspace of dimension at least five.
- Runtime: under 1 s.
- Implementation: `src/checking/exact_audit/`.
- Release run: accepted in under 1 s; record in [results/pauli-terms-and-cubic-coefficient.json](../../results/pauli-terms-and-cubic-coefficient.json).
