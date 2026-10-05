# (C3) Cliffords commuting with $H_*$

Fact (C3) of the paper, label `cert:cliffords-commuting-with-hstar`, stated in a box at its first use (used in `gauss-lem:three-correction-compatibility`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
With $H_*$ as above, a Hermitian gate $K\in\C_2(3)$ commutes with $H_*$ if and only if $K=\pm(-1)^{ax+bxy+c(1+x)z}$ with $a,b,c\in\F_2$. The traceless ones are $\pm Z_x$ and $\pm$ two non-Pauli gates whose supports span $L_Z$.
```

## Check

```sh
python check.py math --fact cliffords-commuting-with-hstar
```

- Finite domain: The 1,344 affine permutations times the 512 diagonal Clifford phases, modulo scalars.
- Checked: Exactly 64 classes commute with H_*, and exactly eight of them have a Hermitian representative, namely plus or minus the listed gates.
- Runtime: under 1 s.
- Implementation: `src/checking/c14.py`.
- Release run: accepted in under 1 s; record in [results/cliffords-commuting-with-hstar.json](../../results/cliffords-commuting-with-hstar.json).
