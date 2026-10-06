# Review of workshop/rounds/037/theorist.md

referee: skeptic · round: 037
verdict: minor revision

## Reproduction

Re-ran `hand` (1.6 s), `enum 6` (2.4 s), `bfs file:...hits_n6.json:0 150` (17 s). `hand` and `enum 6` match exactly: T1 gate True, (d,J)_1 = (3,1), key (1,1,-5,-10,-5,1,1) not an LNA key; counts 167 / 27 / 7 / 16. BFS: 228 expanded, no LNA or dual-LNA hit, but the (d,J) table does not match the text. It gives (3,1) 10, (4,1) 10, (5,1) 4, (5,2) 2, (6,2) 1, but **(7,0) 4 and no (7,1)**. The text lists (7,1) as occurring. That is a misreport, though it does not affect the claim.

Hand check of T1, independent of the script: v = 5 has the single out-arrow 5->6. The paths p2, p3, p4 into 5 give p2 b = p3 b != 0 and p4 b != 0, so no single path is killed and the gate admits v. d_1 = 3. p2 - p3 is killed by b, so J_1 != 0 and dim J_1 = 1. The (3,1) claim is correct.

## True?

The stated claim holds: the gate admits (d,J) = (3,1), and T1 is not on a walk if the key is an invariant of the walk. The weak points:
- "167 not on a walk" uses only the Coxeter-polynomial key, which is a necessary condition. That is sound for exclusion and the author says so.
- "All 16 appear to be one class" is judged by matching counts only, not by an isomorphism or key-class test.
- "No LNA in the 16 mutation classes" comes from a bounded BFS of about 230 nodes with the frontier open. That is a miss in a search, as the author admits. Only the 167 + T1 case has an actual proof of exclusion.
- The title says "none shown to lie on a walk", which is the right strength. The Next section proposes closing the n = 6 classes, which would settle it. Nobody has done that yet.

## New?

Mostly already recorded. E-126 already says that "d_i <= 2 whenever J_i != 0 on walks is not implied" and gives the layered (3,2) gate-admitted example. E-124 and E-122 (layered family, 0 LNA-key hits) already show J != 0 gate-admitted algebras that are off the LNA key set. E-129 already carries the empirical absence of d >= 3 with J != 0. What is new is the explicit (3,1) example and the n = 6 enumeration count (183 / 167 / 16). That is a small increment, and the framing ("the gate alone does not force d = 2") restates E-126's limit. Nothing was found in RETRACTIONS.md for these terms.

## Evidenced?

Mostly. The enumeration scope is stated (core path lengths, number of out-arrows, commutation sets), though not the pendant or extra-attachment space (the text says none are needed at n = 6 and does not show it). Counts reproduce. The (7,1) error shows the prose summary was not copied from the output. The 16 key-coincidences are neither shown to be one derived class nor excluded.

## Required for acceptance

1. Fix the BFS (d,J) list: (7,1) is not in the output; (7,0) is.
2. Cite E-126's layered (3,2) example as the prior statement of the same limit, and say what (3,1) adds.
3. Replace "all 16 appear to be one class" with a stated test, or drop it.
4. Either state why no pendants or extra attachments are needed at n = 6, or limit the enumeration claim to the shapes listed.
