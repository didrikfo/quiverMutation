# Scholar's notebook (after round 025)

## Believe
- Ladkani 2.3(c) = AI 2.32(b) = `tiltingPlus`: one map g_i : p |-> (p alpha)_alpha. Authors' statement, not derived here.
- r021: (k,i) rewrite entry = dim coker g_i; (i,k) = dim ker psi_i (needs step-7 completeness, least sure link). Cartan congruence fails iff some dim ker g_i != 0 iff tiltingPlus False.
- r023: reject (gate admits, tiltingPlus False) <=> J != 0 and no single path in J. Long square with minimal relation => reject. Converse needs out-degree 1 and generators fixed.
- r025: the "J != 0 <=> out-degree 1 and long square" iff holds on class-0 walks n = 5..7 (300 s caps) and FAILS at n = 8: 55/59 J != 0 steps have out-degree 2;
  61 distinct rejects all commute into one arrow only ("D-part 1/2", killed into the other by zero relations), all fail the Cartan test through the real rewrite.
  Only one read by hand. Counts cap-dependent (referee 42, me 55).
- Rejecting parents (distinct): n=6 class 0 1123, n=7 156.
- r014/r018: guard-admitted steps never fail; Cartan/tilting failures are guard-refused (verify: rejects here are gate-admitted yet the guard... not re-checked at n=8).
- CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked (403). Do not retry.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Say which.

## Did
- R001..R023 (see rounds/*/scholar_*.py). R025: scholar_n8rejects.py (classify + Cartan test of out-degree >= 2 rejects), n = 8 row via r023 script.

## Next
- n = 9 class 0 walk with the classifier (overnight proposal; needs checkpointing).
- Hand-verify a second D' example (esp. a dim J = 2 one); extend classifier to test c.beta in I via zero relations explicitly.
- Classes 1-3 at n = 6..8 unrun. Diagonal and i,j != k Cartan entries underived.
- Lesson: a regularity checked to n = 7 is not a theorem; always run one size further before writing "iff". Referee found it in 200 s.
