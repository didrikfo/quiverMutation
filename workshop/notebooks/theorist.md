# Theorist notebook (rewritten round 039)

## What I believe now
- J_i = Hom(S_v, e_iA) (E-122); gate = "no single nonzero path in J" (E-114); L1: dim J_i <= d_i - 1 (E-126). The gate does not force d = 2 (hand T1, E-132) nor out(i) = 3 (round 039: hand 7-vertex algebra, d_1 = 3, J = 1, out(1) = 2, gate-admitted, key off).
- Round 039 key finding: in the n = 8 c0 walk every row with J_i != 0 (229 rows, 192 steps, incl. all 9 with d >= 3) has a child that FAILS the key guard; all 121 d >= 3, J = 0 rows pass. n = 6, 7: 0 of 123 H != 0 steps pass. So E-129's d = 2 and E-135's out(i) = 3 describe parents of refused steps, not the class. Observation only (one class per n, capped), no proof.
- Exact guard criterion: C_B = C' + e_v j^T, pass iff R(x) = (1+f(x))(1+f(1/x)) - x b c = 1 + t (determinant lemma, checked numerically at x = 2, 3, 5; t = j^T C'^{-1} e_v, t = 0 in every observed step). Pass/fail is a polynomial identity, det alone never decides.
- "Born vs inherited" J_i: born needs a relation across >= 2 out-branches; d >= 3 born needs out 3 (three one-path branches) or out 2 with a 2-path branch. Heuristic, unproved.
- Earlier: L2 conditional on AI 2.31; 16 key-coincidence fans (E-123) all meet an LNA (E-134), so they are in the LNA class if the guard holds.

## What I tried
- 039: `rounds/039/theorist_d3walk.py`, `theorist_d3rows.py`, `theorist_out2.py`, `theorist_guard.py`. 037: `theorist_d3.py`. 034: `theorist_dimji.py`. 031: `theorist_t3.py`, `theorist_layers.py`. 029: `theorist_circuit.py`, `theorist_keys.py`.

## Next
- Prove/refute: J != 0 (j = e_i) implies R(x) != 1 + t; first obstructions: x = -1 value, x^{n-1} coefficient. Search for any passing J != 0 child (skeptic request).
- Verify row 7798 kernel (5-1-7 + 5-4-1-7 + 5-4-3-7) with skeptic_kernel on the saved pkl.
- If J != 0 never passes, restate T5: the d = 2 / out(i) questions dissolve; the real question is which J = 0 steps pass.
- Blind spots: one class, capped tail walk; a formula checked on the repo's own perI; "never observed" is not "never".
