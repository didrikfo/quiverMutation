# Maverick's notebook

## What I now believe (after round 107)
- E-269 Z-congruence search (iso140 + E-276 closed-form bound) says YES on within-quipu pairs in the 4 testable n = 12 low-lam groups (lam 1.2612, 1.3836, 1.4250, 1.4886; slack 0.5 and 30). The other 2 of the "6" are NOT testable: S = C+C^T singular (the "lam 1.000136" group is a Jordan block at 1, not Salem; the 47-group has Salem 1.4013 but null S).
- No certified-inequivalent pair at lam < 1.5 in n = 8..12 samples: in every group (lam<1.5, S nonsingular) all sampled members share |det S|, SNF and q_S; each group has one quipu. NO power at lam < 1.5 is unavailable by Cartan invariants. Planted-error control (perturb C, expect NO) not yet run -- cheapest way to give the search a NO at low lam.
- Earlier: r105 six mixed pairs YES (lam 2.24..2.63), 4/4 cross NO at n = 12 group of 140; necessity step Ladkani math/0610685 Cor 3.13/Lemma 3.14/Cor 3.15 (by count). No independent second invariant for LNAs (HH* = k, q-Cartan det 1).
- Convention Phi = -C^-T C, n = len + 2; digit strings are relLengths.

## What I tried
- r107: rounds/107/maverick_n12low.py, maverick_nopair.py.
- r105: maverick_mixed.py, maverick_lowlam.py; r099: iso140 (+_ctl, _transp, _setup); r097/094/091/089: cartanZ, genusscan, quipupairs.
- ALWAYS run a within-class control with nontrivial P before trusting a NO. A within-quipu YES at low lam is weak (2-11 short vectors).

## Watch for
- Venv loses cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX engine: theorem numbers by counting. Print lam and check S nonsingular before claiming anything per group.
- Long combined runs hit the 2-min tool limit: run one group at a time.

## Next
- Planted-error control at lam < 1.5: perturb C_b by an entry that changes Phi-commutant structure and check NO (honest power test, not a derived-equivalence NO).
- Extend E-269 to singular S (quotient by radical) for the 47-group (Salem 1.4013).
- Two-Salem-pair cases N(12,5)~N(12,4) need a different method.
- Sweep the other mixed algebras through the test.
