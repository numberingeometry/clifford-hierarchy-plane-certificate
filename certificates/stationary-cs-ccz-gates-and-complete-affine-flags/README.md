# (C1) Stationary gates $\mathrm{CS}_{yz}^a\mathrm{CCZ}^b$ and complete affine flags

Fact (C1) of the paper, label `cert:stationary-cs-ccz-gates-and-complete-affine-flags`, stated in a box at its first use (used in `lem:tool-bridge`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
The gates $G=\mathrm{CS}_{yz}^a\CCZ^b$ with $a,b\in\F_2$, that is, $I$, $\mathrm{CS}_{yz}$, $\CCZ$ and $\mathrm{CS}_{yz}\CCZ$, are stationary, with $135$, $3$, $15$ and $3$ GSC pairs. Let $\pi$ be the permutation of $CGC^\dagger$, where $C\in\C_2(3)$ maps the Pauli MASA of an input Lagrangian of $G$ to the diagonal matrices. Then for every affine permutation $\alpha$ of $\F_2^3$, if $\alpha\pi$ and the translations preserve a common complete affine flag, so do $\alpha$, $\pi$ and the translations.
```

## Check

```sh
python check.py math --fact stationary-cs-ccz-gates-and-complete-affine-flags
```

- Finite domain: The four gates I, CS_yz, CCZ and CS_yz CCZ; the 135 Lagrangians of V_3; the 21 complete affine flags and the 1,344 affine permutations of F_2^3.
- Checked: The witnesses of each gate (135, 3, 15 and 3, all of the form (L,L)); for each of the 156 input Lagrangians, the flag implication for all 1,344 affine permutations against the 21 flags; the counterexample D of the paper.
- Runtime: under 1 s.
- Implementation: `src/checking/c16.py`.
- Release run: accepted in under 1 s; record in [results/stationary-cs-ccz-gates-and-complete-affine-flags.json](../../results/stationary-cs-ccz-gates-and-complete-affine-flags.json).
