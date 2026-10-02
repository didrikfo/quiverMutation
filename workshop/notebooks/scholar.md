# Scholar's notebook (after round 023)

## Believe
- Ladkani 2.3(c) = AI 2.32(b) = `tiltingPlus`: one map g_i : p |-> (p alpha)_alpha. Authors' statement, not derived here.
- r021 (rounds/021/scholar.md): (k,i) rewrite entry = dim coker g_i from steps 4, 6; (i,k) = dim ker psi_i from step 7 (needs step-7
  completeness, the least sure link). Cartan congruence fails iff some dim ker g_i != 0 iff tiltingPlus False.
- r023 (rounds/023/scholar.md): reject (gate admits, tiltingPlus False) <=> J != 0 with J = {c not in I : c alpha in I for all alpha out of v} and no
  single path in J (the gate's code refuses a single-path kernel). Long square => reject needs only a minimal relation c alpha, c >= 2 paths.
  Reject => long square needs out-degree 1 (NOT forced: hand case D, out-degree 2, rejects) and reading relations up to change of generators (case G).
  "Short" square is fine because c is itself in I. Walks n = 5..7 class 0 (300 s): J != 0 <=> out-degree 1 and long square, step by step.
- Rejecting parents (distinct): n=6 class 0 1123, n=7 156; not all strict A5. Counts cap-dependent.
- r014/r018: guard-admitted steps never fail (0 of ~1.3e6); all Cartan/tilting failures are guard-refused.
- CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked (403). Do not retry.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Say which.

## Did
- R001..R021 scholar_h015*, square, walk, replay, a5key, cartan_vs_tilt, e078_diff, step7_entries.
- R023: rounds/023/scholar_longsquare.py (--hand cases A-G; walks n=5..7 table).

## Next
- Is a D-like algebra (J != 0, out-degree >= 2) reachable from an LNA by guarded steps? n = 8, 9 search (overnight proposal if > 10 min).
- Referee of step-7 completeness; run D through the real rewrite and see the Cartan test fail.
- Classes 1-3 at n = 6, 7 unrun. Diagonal and i,j != k Cartan entries underived.
- Lesson: derive from the code's own definition of the gate (isMutable docstring: single path vs combination) -- it explained the ">= 2 paths" at once.
