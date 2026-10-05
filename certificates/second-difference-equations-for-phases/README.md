# (C5) Second-difference equations for phases

Fact (C5) of the paper, label `cert:second-difference-equations-for-phases`, stated in a box at its first use (used in `lem:m-M-1`). Facts it depends on: none.

## Statement

As in the paper, above the dashed line of the box ([statement.tex](statement.tex)):

```latex
Let $f:\F_2^3\to\mathbb R/\mathbb Z$, let $T$ be a symmetric binary $3\times3$ matrix, and let $\Theta$ be as in the proof of Lemma~\ref{lem:m-M-1}. Then:
\[
\begin{array}{c|c|c}
T&\text{equations}&\text{consequence}\\\hline
\operatorname{rank}T=3&\Delta_k\Delta_af=0\ (k\cdot Ta=0)&4\Delta_jf=0\ (j=1,2,3)\\
\operatorname{rank}T=3\text{ or }\operatorname{im}T=\Span\{e_1,e_2\},\ T\ne\Theta&\Delta_{\Theta Ta}\Delta_af=0\ (\text{all }a)&4\Delta_jf=0\ (j=1,2)\\
T=\Theta&\Delta_{\Theta Ta}\Delta_af=0\ (\text{all }a)&2\Delta_jf=0\ (j=1,2)
\end{array}
\]
```

## Check

```sh
python check.py math --fact second-difference-equations-for-phases
```

- Finite domain: The 28 nonsingular symmetric binary 3x3 matrices of the first system; the 28 nonsingular and 4 rank-two matrices of the others (60 systems).
- Checked: For each system and each of its target rows, membership in the integer lattice spanned by the rows of the equations.
- Runtime: under 1 s.
- Implementation: `src/checking/n3_descent/`.
- Release run: accepted in under 1 s; record in [results/second-difference-equations-for-phases.json](../../results/second-difference-equations-for-phases.json).
