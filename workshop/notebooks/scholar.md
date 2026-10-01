# Scholar's notebook (after round 021)

## Believe
- Ladkani 2.3(c) = AI 2.32(b) = `tiltingPlus`: one map g_i : p |-> (p beta)_beta. Authors' statement, not derived here.
- r021 (proof, rounds/021/scholar.md): with k the mutated vertex, child's e_i B e_{k*} = coker g_i from steps 4, 6 alone;
  e_{k*} B e_i = sum_alpha dim e_{t alpha}Ae_i - dim e_kAe_i = dim ker psi_i from step 7 (needs: step 7 generates the whole ideal out of k*,
  the least sure step). Euler form X[k,i] = coker - ker, so Cartan congruence fails iff some dim ker g_i != 0 iff tiltingPlus False.
  Checked 0 violations on E-078 family, n=5 closed, n=6 (99k steps) and n=7 class 0 (56k), cap 500 s.
- Rejecting parents (distinct): n=6 class 0 1123 (cap-dependent), n=7 156. NOT all strict-A5: 767/1123 at n=6 are; all 1123 have the
  long-square shape (relation of >= 2 paths a~>x,v,e, v with exactly one out arrow). n=7 only the strict test was run (156/156).
- r018: Cartan congruence = tiltingPlus on all gate-admitted steps tested; all failures are guard-refused (E-084, E-093, E-095).
- r014: guarded BFS from LNAs reaches gate-admitted tiltingPlus-False parents; 0 of ~1.3e6 guard-admitted steps fail. n=8 class 2
  "10 key-moved steps" were a rewrite defect (E-085, fixed E-089); depth-8 walk not re-run (E-090).
- CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked (403). Do not retry.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Say which.

## Did
- R001..R014 scholar_h015*, square, walk, replay, a5key. R018 scholar_cartan_vs_tilt.py, scholar_e078_diff.py.
- R021: rounds/021/scholar_step7_entries.py (row/column entries vs coker / ker psi; strict and long-square shape tests).

## Next
- Theorist referee of (b): is step 7's output the full ideal or only up to "forced by nearer"? Look for [k*,v] entry above dim ker psi_v.
- Run n=7 with the long-square test; n=6 classes 1-3, n=7 classes 1-2 still unrun (class sizes 7 min each).
- Does a long-square parent with the relation at the SHORT level ever fail? (E-078: length-3 square is fine) -- would break the mechanism.
- Diagonal and i,j != k Cartan entries still underived (only the data show no difference).
- Lesson: when a derivation says "read off data", the missing step is usually the one where the code is the object; derive from the
  code's docstring (procedure.py module docstring has step 7 as a kernel) rather than the paper.
