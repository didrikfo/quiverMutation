# Scholar's notebook (after round 018)

## Believe
- Ladkani 2.3(c) = AI 2.32(b) = `tiltingPlus`: one map g_i : p |-> (p beta)_beta. Authors' statement, not derived here.
- r018 (T5): Cartan congruence (Ladkani Prop 3.6 = Lemma 3.5 specialised) vs tiltingPlus: 0 xor in 98k gate-admitted steps
  (E-078 family, n=5 closed, n=6 and n=7 class 0 partial, fixed library). In all 807 failures the diff X-Y lives in row k
  off-diagonal and equals -dim ker g_i. Reason: X_{k,i} = dim coker - dim ker, rewrite gives coker. So congruence is the
  same condition read through the rewrite, plus a check on the rewrite itself (caught E-085's defect). Conditional step:
  rewrite's (k,i) entry = dim coker g_i (data only).
- r014 (T5): guarded BFS from LNAs reaches gate-admitted tiltingPlus-False parents (n=6 dist 8; n=7..9 dist 5..7); all
  guard-refused; 0 of ~1.3e6 guard-admitted steps fail. A5 square shape (E-078) is the mechanism.
- n=8 class 2 "10 key-moved steps" were a rewrite defect (E-085, fixed E-089); walk to depth 8 not re-run (E-090).
- CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked (403). Do not retry.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Say which.

## Did
- R001..R014: scholar_h015*, square, walk (rounds/014/scholar_walk.py, --plan, --stop-on-reject), replay, a5key.
- R018: rounds/018/scholar_cartan_vs_tilt.py (crosstab + diff = -ker check), scholar_e078_diff.py.

## Next
- Theorist: read step 7 to prove rewrite entry = coker; diagonal/column-k entries can they differ?
- Toolsmith: Cartan assertion in mutateAtVertex (cheap rewrite check).
- Finish n=6 classes 1-3 and n=7 classes 1-2 of the crosstab (7 min each) if wanted; n=8 class 2 under the fixed rewrite
  still needs the checkpoint in scholar_walk.py.
- Lessons: a "second invariant" that is computed from the code's own output tests the code, not the theory; ask what it
  actually compares before counting it as independent evidence.
