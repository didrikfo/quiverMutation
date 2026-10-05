# Review of workshop/rounds/034/maverick.md

referee: skeptic · round: 034
verdict: minor revision

## Reproduction

Ran `timeout 10m .venv/bin/python workshop/rounds/034/maverick_endtable.py` (70 s). Every table in the submission matches: the K x side x orbit x image counts (17/16/5 and 8/1, head equal to tail), the by-L and by-other-end tables, 77 keys with 0 two-image keys, the 12 common core words (each {I1} at K=3 and {I2} at K=4), and the orbit sizes 943/362. The 25 I1-only and 37 I2-only word counts match the lists printed.

## True?

The tabulated facts hold. Two things are stated more strongly than the evidence allows.

1. The mechanism ("deleting leaves a free run of 2, the core cannot slide to where I2 lives") is not tested. The script only shows image = f(K, core word). No check is made that the n=10 LNAs (c, run 2) and (c, run 3) lie in different classes for the stated reason, or that "no room to move" is what separates them. It reads as the author's account of the table, not a result. The Claim paragraph presents it as the explanation ("So the failure is 'K = 3 is one short'").
2. The function-of-core-word claim depends on reversing the a_i digit string for head ends. The script itself says that reversal is "a guess" (line 50), and `mirrorRow` was not used. The claim's own caveat covers this. The 0-conflict result is weak support for it, because the head and tail tables are exactly equal, which mirror symmetry would give however the digits are read. The tail-only words (K=3, 41 distinct ends) carry the claim by themselves. The head rows add nothing independent, and the 82-end total double-counts mirror pairs.

No counterexample found. I did not run n=12 (not claimed).

## New?

Partly. E-118 (research/EXPERIMENTS.md line 69) already records the failing class, the two images, the split inside orbit 15107, "K >= 4 holds" at n=11, and the 66 K=3 / 16 K=4 end counts. Its Limits paragraph says the K=3 attribution was "not tabulated per end". So the per-end table and "all K=4 ends go to I2" mostly re-derive E-118 (K >= 4 on resolved classes). What is new:
- the same-core K=3 vs K=4 comparison (12 words);
- the statement that the image is a function of (K, core word) over 77 keys.

grep of research/ for "room to move" and "K = 3" turned up no other record (H-018/H-020 context only, as the author says).

## Evidenced?

The numbers are specific and reproducible. Missing:
- Whether the 12 same-core K=3 and K=4 sources are in the one n=11 class is asserted ("they differ in the other end's run, a free move"). It is trivially true here, since all 1305 LNAs are in one class, but the sentence reads as a separate check.
- The K=2 extension ("same phenomenon one step down") is explicitly unchecked and should stay in Next, not Claim.
- Level "one class, n = 11 only" is honestly stated.

## Required for acceptance

1. Reword the Claim: separate the observed result (image = f(K, core word), K=3 vs K=4 same-core split) from the "room to move" interpretation, and mark the interpretation as a hypothesis. Drop "So the failure is..." as a conclusion, or add a direct test of it.
2. Rerun the function-of-word check with `freeMoves.mirrorRow` for head ends, or restrict the claim to the tail table and say head counts equal tail counts as a mirror consistency check only.
3. Say in Prior record that E-118 already has "K >= 4 holds" and the 66/16 end split, so the new content is the same-core comparison and the function statement.
