# Review of workshop/rounds/071/maverick.md

referee: skeptic · round: 071
verdict: minor revision

## Reproduction
All five commands re-run, 55-70 s each, exit 0. Outputs match the committed .out files line for line except wall-clock times (snfx, orbits, conj, pieces, place16). Conj: control plus 14/14 FOUND, d=17, tries 8-12, same as claimed. Pieces: one tuple per group, cross-group disjoint True. Orbits: 9/10/25/213 and class sizes equal members-in-group. place16: one Coxeter key, not a group key, signature (9,1,1) vs (9,2,0), Phi non-periodic within cap.

## True?
No counterexample found. Points that do not break it but are loose:
- "Unplaced is justified": the signature difference (9,1,1) vs (9,2,0) is what carries it; the "Phi not periodic (cap 200)" is a cap result, not a proof. Say so.
- The 15-group "9 classes, not 15" counts orbit+mirror classes, i.e. a lower bound on derived-class count under known moves, not a count of derived classes. The title's "is 9 orbit+mirror classes" is correct; do not read it as 9 derived classes.
- Pieces are ker f(Phi^T) with restricted G only; the indefinite rank-4/6/8 pieces are compared by (rank, sig, det) alone ("indef/big" theta), which is weak (genus-level at best). The title says "per-factor isometry classes are constant"; what was tested is a necessary invariant of isometry class, not the class (the text admits this in Claim, not in the title).
- Z-conjugacy "FOUND" in 14/14 is a certificate (exact check), fine. The control is sound.

## New?
Partly. E-197 already has per-factor SNF and q_G constant per group, and lists SNF(xC+C^T) and Z[Phi]-isometry as not run. F-047/E-154 cover Z-conjugacy at other n; E-190 has the four groups. Grep of "xC + C^T", "theta series", "orbit+mirror classes" in FINDINGS/HYPOTHESES/RETRACTIONS/EXPERIMENTS: only E-197 (n=11 SNF), E-066 (orbit+mirror vs key, n=10,12..16 cores) and an n=7 SNF(xC+C^T) line near EXPERIMENTS.md:386. So the x-scan at n=11 and the 14 conjugators and counts are new. The author is right that SNF(xC+C^T) is implied by Z-conjugacy of Phi; the x-scan is therefore redundant in the 15-group and its "four groups" extension is the only content.

## Evidenced?
Mostly yes: scripts, ranges (x=-5..5, 7 (a,b) pairs, theta to norm 8) and outputs are stated. Missing: the theta series are shown only in the .out, and the indefinite pieces have no real invariant beyond (rank, sig, det). The author honestly states power against a NO is unestablished and the n=10 control was not run. The "9 distinct tuples = 9 classes" coincidence under the Phi/Phi^T mismatch is flagged as unexamined; it is either a labelling artefact or a real separation by a non-invariant, and is not evidence for anything.

## Scope
Title overreaches in two spots. Narrowed wording: "n = 11, four periodic-indefinite Coxeter groups: SNF(xC+C^T) (x=-5..5), Z-conjugacy of Phi (15-group only, 14 pairs), and (rank, signature, det[, theta to norm 8 where definite]) of ker f(Phi^T) are constant in each group (null; low power, no control); the 15-group has 9 orbit+mirror classes; the 16 unplaced LNAs lie in none of the four groups (signature (9,1,1))."

## Required for acceptance
1. Retitle per Scope: "isometry classes" -> "(rank, signature, det, theta<=8) of pieces"; mark Z-conjugacy as 15-group only.
2. State that Phi non-periodicity of the 16 is cap-limited and the placement rests on signature and key.
3. Say E-197 already recorded constant per-factor SNF and q_G, and E-066 the orbit+mirror/key refinement, so the delta is stated against them.
4. Either run the n=10 control with the indefinite filter dropped (one sitting, as the author proposes) or label the null as uncontrolled in the title/claim line.
