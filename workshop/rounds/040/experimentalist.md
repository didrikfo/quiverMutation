# Experimentalist, round 040 (conference, position statement)

**Most promising question.** Does any gate-admitted J_i != 0 step keep the LNA key when it is taken off the walks? Concretely: for every J_i != 0 row of the capped walks at n = 6, 7, 8 (classes 0, 1, 2), compute the child's key with the key guard off and tabulate by (d_i, out(i), key kept or not). If the count is zero everywhere, "silting-not-tilting steps leave the key class" has a table behind it and T5's meeting question (E-136, E-139) can be closed in one direction. If one row keeps the key, the law fails and that row is the cheapest counterexample the workshop could get.

**Why this one.** It is the only question where a single positive row changes the state of T5, and it needs no new theory to run. E-140 covers class 0 only, so the column for classes 1 and 2 is still empty.

**Weakest claim the workshop relies on.** "Every J_i != 0 step is refused by the key guard" (E-140), and with it "d_i = 2 whenever J_i != 0" (E-131). Both are statements about capped BFS: levels 7 to 9 of 10, with d >= 3 living only in the last ~8 % of the walk (E-137). The key guard is also the BFS filter, so on walk rows its refusal is close to built in. The number tests the walk, not the step. The same pattern was called vacuous once already (E-107, "the 42 reproduce and are BFS-reached"). Also: the 9 d >= 3 rows are 5 algebras from one class at one depth.

**What I need.**
- skeptic: confirm that the child key in the replay is computed on the child (not inherited from the parent), and which function `coxkey-pres` uses, so the zero count can be read as a test.
- toolsmith: `--plan` for the n = 8 classes 1 and 2 J_i != 0 column with the key guard off, sized under the 560 s cap; a positive control (a hand-built J_i != 0 child that keeps the key) so that a zero means something.
- theorist: when can a J_i != 0 child keep the key? E-138's `C_B = r C_A r^T + H` gives a candidate condition. Without a stated prediction the table has nothing to fail.
- scholar: nothing needed this round.

**Watch.** A time cap is not a verdict, and counts move with machine load (E-127, E-107). Algebra ids are keys, not isomorphism classes. The 5 parallel rows of E-129 are 3 algebras, not 5. Any "zero" should be reported with the level reached.

**This round.** No runs (conference). I have not edited other files.
