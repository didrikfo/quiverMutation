# Review of workshop/rounds/057/toolsmith.md

referee: skeptic · round: 057
verdict: minor revision

## Reproduction
Used the existing /tmp/tsm pickles (did not rebuild). Re-ran, all with `.venv/bin/python`:
- `fail` for CLS=1 (16 steps) and CLS=2 (9 steps): diff against stored `*_fail_c*_out.txt` differs only in the per-step timing column. All iso, 8 c1 steps with parallel pairs (maxpar 2, step 12 has 2 pairs).
- `paths 0:40` CLS=1: TOTAL iso 214. CLS=2: iso 156. Matches.
- `n10 run 0:5`: iso 35, 16 s. Matches.
- `ctrl` CLS=1: A_no 0, A_iso 9, B_no 9, B_iso 0, n 8. Matches the claimed 9/9 and 9/9.
- Counted in the stored c1 paths output: 214 edges, 56 with a parallel pair. Matches.

## True?
I found no error. I read `symcheck2`: per-block unknown m x m matrix M plus rad^2 corrections, all relations of c required to vanish, Rabinowitsch variable per block for det M != 0. The logic is sound: Groebner basis != [1] gives a homomorphism KQ_c/I_c -> End(T) over the closure that is surjective (arrows hit rad/rad^2 with invertible M, Nakayama) and, with equal dims, an isomorphism. Two points rest on assumptions, not tests:
- The dimension of KQ_c/I_c is taken from `child_info` `d`, and `crels` is taken to generate the full ideal of c. If `reduction` ever dropped a relation, a wrong algebra of equal dims could pass; the test cannot detect that.
- The author concedes there is no wrong-algebra control with equal dims, equal arrows and a parallel pair. The matrix-valued path is therefore tested only against the binomial-to-monomial perturbation (9/9), which cannot show that the matrix machinery rejects a wrong GL_m-twisted algebra. The "absorbed" doubling result is expected. I did not build such a control (not doable here); the power claim for parallel blocks is weaker than the 16/16 headline suggests.
- n = 10: all 35 edges have maxdim 1, no parallel arrows, so the verdict reduces to dims and arrows plus a few scalar equations (6 edges' worth of "no equations"). Agreement is nearly forced there. The author says this.

## New?
Extends E-167 (13 edges, 8/16 decided; the 8 undecided are exactly the parallel-arrow ones). E-168 scope line (EXPERIMENTS.md:35) already records the 44 paths, 26 c1 / 18 c2, 214 / 156 edges, and the 16 failing steps. Nothing in RETRACTIONS.md or HYPOTHESES.md on "End(T)". Genuinely new: the 8 parallel decisions, c2 (9 failing steps, 156 edges), and n = 10 (35 edges, E-169 had no End(T) test).
Discrepancy: E-168 states Hom(T,T[-1]) != 0 for the 16 c1 failing steps only. The submission says "25 of 25" failing steps are blind to the premise, but nothing in the record shows Hom(T,T[-1]) != 0 for the 9 c2 failing steps. "Failing J != 0" for c2 is asserted, not cited.

## Evidenced?
Mostly yes: counts per set, outputs named, determinism check on logged keys, stated commands run in about 1 minute. Missing: the 9 c2 steps' J != 0 evidence (above), and the per-edge ideal sizes are only in the output files.

## Scope
Title says "25 failing n = 7 steps ... 405 path edges" and "the comparison cannot test the J = 0 premise". The second half is correct for the 25 steps tested (it accepts a J != 0 algebra); it is an inference, not a theorem, so "cannot test" should read "did not discriminate at any of the 25 steps". The theorist lemma suggestion (End(T) = mutation algebra always) is a conjecture; the sample is n = 7 and n = 10 with one start.

## Required for acceptance
1. Cite or compute Hom(T,T[-1]) != 0 for the 9 c2 failing steps (E-168's check, one command), or drop "25 of 25 J != 0" to "16 c1 + 9 c2 failing-by-key".
2. Reword "cannot test" in the title to the observed "did not discriminate at 25 of 25 failing steps".
3. State in Claim that `crels` completeness (generating set of I_c) and `child_info` dims are inputs, not outputs, of the test.
4. [next round] A wrong-algebra control with equal dims and arrows and a parallel pair (for example a child with one relation altered, same Cartan), to show the matrix-valued path can reject.
