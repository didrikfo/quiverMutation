# Maverick's notebook

## What I now believe (after round 101)
- r101 (literature): the necessity step of H-024 has a source: Ladkani math/0610685 Cor. 3.13 + Lemma 3.14 + Cor. 3.15 (numbers by counting LaTeX environments). Stated for posets; the proof is general for finite global dimension, so the extension to LNAs is ours (one sentence). P Phi_a = Phi_b P is the proof of Lemma 3.14, not its statement. Bessenrodt-Holm Thm 2.2 gives only unimodular equivalence of C (weaker). Nikulin 1.13.3 still unsourced and not needed for NOs.
- No independent second invariant for LNAs: HH* = k (F-019 derivation, inputs rechecked in 2312.14699), AAG gentle only, q-Cartan det = 1. The E-259 pairs already have an independent NO: the quipu theorem. So the Cartan NOs are cross-checked there.
- 10 literature pairs (Ueda Cor 1.3, 1.6; Ladkani 0911.5137 Cor 1.2), n = 9..12, all pass char Phi, Z[Phi]-module, q_S (3 non-vacuous); E-269 search says YES on the 3 with signature (n-1,1) (n = 10, 11, 12). All radical powers, so no mixed-digit control yet.
- r099: n = 12 key group of 140 pair is not Z-congruent (Phi-equivariance, E-269); r097: q_S separates 7 of 8 pairs. Convention Phi = -C^-T C, n = len + 2.

## What I tried
- r101: rounds/101/maverick_litpairs.py (invariant table), maverick_litiso.py (E-269 `test` on lit pairs; needs pq threshold 1.05 and M.N = n: the 099 code assumes Salem number > 1.5 and N = 12).
- r099: maverick_iso140.py, _ctl, _transp, _setup. ALWAYS run a within-class control with nontrivial P before trusting a NO. qfminim returns 2 x vectors; pari.allocatemem(3e9).
- r097/r094/r091/r089: cartanZ, genusscan, quipupairs, earlier T4/H1.

## Watch for
- Venv lacks cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX engine: theorem numbers only by counting `\newtheorem` environments (shared counters); say "by count".
- E-245/E-246 string labels: n = len + 2. N(n,l) = "l"*(n-l)+"0"*(l-2).
- 099 `test` bound formula depends on lam; skeptic should check it for lam < 1.5.

## Next
- Extend the exhaustive test to two Salem pairs (signature (n-2,2)): N(12,5)~N(12,4), N(12,7)~N(12,6) as controls; also N(12,3)~N(12,6) by transitivity.
- Get a mixed-digit positive control (Happel-Seidel / Lenzing-Meltzer-Ruan, or the 1 002 pairs of E-238 if the experimentalist ran them).
- Z[Phi]-class invariant for the 140 pair vs q_S (cheap second invariant); a split of an E-117/F-010 class would be a bug.
