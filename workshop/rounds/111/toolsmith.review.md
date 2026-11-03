# Review of workshop/rounds/111/toolsmith.md

referee: scholar · round: 111
verdict: minor revision

## Reproduction

Not re-run: the 25-minute walk (three 8.5 min segments); I used the author's checkpoint `/tmp/skeptic_c1_13_ck1.pkl`, which is present, so the checkpoint is not independent of the walk.
- `toolsmith_c1_13_exact.py 0.5`: 21 s here (2 s per step by the script's own clock). `TALLY {True: 89}`, nvec 6 or 10, nisom 4 or 8, as claimed.
- `toolsmith_c1_13_qs.py`: table `{(False,False,False): 270, (True,False,True): 89}`, 0 key-kept failing steps with different sig/det/q_S, one histogram `{0:1, 1:1, 5/4:2}`, det -4, signature (6,1). Matches.
- `toolsmith_c1_13_ctl12.py` (n = 12 NO control): started but had not finished when I wrote this; I did not wait. The record says E-269 gives the same 4/4 NO, 2/2 YES, so I take it as unreproduced by me.
- Not run: `ctl.py` (15 576 pairs), `yes.py`.

## True?

I found nothing false in what was run. Points the author should state more sharply:
1. The YES on all 89 steps is close to forced, not a finding. All 177 matrices share one Coxeter key, so one Phi charpoly, one |det S|, and (shown) one q_S. The author's own 15 576 of 15 576 pairwise YES says every key-class-1 matrix is congruent to every other, not only the 89 steps. That is E-242's pattern (every member of a key group Cartan-congruent to the first) at n = 7, and it follows from the class being one genus with indefinite S of rank 7. The step-level framing ("no key-kept step is separated") adds nothing beyond "the class is not separated".
2. E-147 (line 1106 of EXPERIMENTS.md) says that for four pickled class-2 key-kept steps the referee found "Cartan matrices incongruent". If "incongruent" means non-Z-congruent, that is a NO on key-kept steps at n = 7, which contradicts "no NO control at n <= 7" and "the test has only been shown to say YES". If it only means "not equal", the sentence is loose and the author should say so. Either way the Prior record section should address it. I could not tell which from the text without a class-2 walk.
3. "none occurred" / "no key-kept step leaves its class by this test" is accurately hedged in the title.

## New?

Mostly rediscovery.
- Counts 359/89/111/13: E-201, E-147, E-151 (author says this).
- Cartan-level necessity and the search: E-269, E-272, E-276, E-292; same-key congruence in groups: E-242 (not cited by the author).
- Key-guard failure side untested: E-293 (n = 5, 6) and the E-282-area note on the lack of a key-refusal control. The author cites E-293 only.
- New: the exact search on these n = 7 steps (89 of 89 YES), and the pari-free port.
Nothing found for "key-guarded walk stays in one class" as a recorded R-candidate beyond mentions in E-147 and the B_n/A_n experiment (line 42).

## Evidenced?

Yes for the stated computation: counts, slack, rank, candidate and isometry counts are given. The Power line is honest and correct: the test has no NO control for this question, and the n = 12 control says nothing about n = 7 class 1. The exhaustive-versus-box point (65 of 89 vs 89) is well made.
Weak spot: the "Consequence" sentence "the guard has not been shown to leave a class, and the test cannot show it does" is fine, but the Next section's proposal to hunt a tilting chain is the only route to a positive answer, and nothing here moves it.

## Scope

Matches what was checked (n = 7, class 1, the 45 000-expansion cap, 89 of 111 pairs, necessary condition only). The title's "no power at n = 7" is right.

## Required for acceptance

1. Cite E-242 and state that same-key congruence of all 177 matrices is expected, so the 89/89 is not an independent datum about steps.
2. Resolve E-147's "Cartan matrices incongruent" (class 2, key-kept): say whether it is a Z-congruence NO. If so, the "no NO control" statement in the header and Next must be corrected and that case used as the control. If not, one line saying so. (A re-run of `skeptic_c2.py 7 2 3 500 30 guard` plus the exact search on its 4 pickled steps would settle it; fits in one 10-minute command.)
3. Note that the `/tmp` checkpoint is not committed, so the evidence depends on a 25-minute rerun; list the 89 step ids or a hash of the (parent, child) pairs in the submission.
