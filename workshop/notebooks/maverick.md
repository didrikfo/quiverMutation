# Maverick's notebook

## What I now believe (after round 109)
- E-269 (iso140 + E-276 bound) says YES on within-quipu pairs in the 4 testable n = 12 low-lam groups (lam 1.2612, 1.3836, 1.4250, 1.4886). Low-power YES: 2-11 short vectors. The other 2 of the 6 groups have singular S = C+C^T (the "1.000136" is a Jordan block at 1; the 47-group is Salem 1.4013 with null S): outside E-269.
- 15 low-lam groups (lam<1.5, S nonsingular) at n = 8..12 (1,2,4,4,4). |det S| is fixed by the key (det C * charpoly_Phi(1)) -- NOT evidence. Real content: SNF(S) and q_S constant in samples (<=300, q_S on 40). One quipu observed per n = 12 group.
- "No certified pair at lam < 1.5" is E-276 restated.
- Planted-error control (group 5, lam 1.4886): 8 single- and 149 two-entry +-1 perturbations in the domain (same Phi charpoly, det S) are ALL Z-congruent to C -> YES. No NO obtained; Power stays "none". Either low-lam classes are rigid or the test is permissive; unseparated.
- Earlier: r105 six mixed pairs YES (lam 2.24..2.63), 4/4 cross NO at n = 12 group of 140 (the only NO controls, lam >= 1.66); necessity step Ladkani math/0610685 Cor 3.13/Lemma 3.14/Cor 3.15 (by count).
- Convention Phi = -C^-T C, n = len + 2; digit strings are relLengths.

## What I tried
- r109: rounds/109/maverick_planted.py, maverick_planted2.py. r107: maverick_n12low.py, maverick_nopair.py. r105: maverick_mixed, maverick_lowlam. r099: iso140. r097/094/091/089: cartanZ, genusscan, quipupairs.
- Always run a within-class control with nontrivial P before trusting a NO; a within-quipu YES is weak.

## Watch for
- Venv loses cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX: theorem numbers by counting. Print lam and check S nonsingular per group.
- Piping a loop into grep hides output until the loop ends; group 4 planted run is slow (>5 min).

## Next
- Find an in-domain NON-congruent C' (change q_S/SNF while keeping det S) and see whether E-269 says NO; or apply the same perturbation recipe to a lam >= 1.66 cross-quipu group to see if the recipe ever yields NO.
- Planted controls on groups 2 and 4; extend E-269 to singular S (47-group); two-Salem-pair cases N(12,5)~N(12,4).
