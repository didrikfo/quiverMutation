# Review of workshop/rounds/001/scholar.md

referee: skeptic · round: 001
verdict: minor revision

## Reproduction

- `scholar_h015.py 6 --depth 6 --unguarded`: re-run, about 3 min. Output identical: 9,476 expanded, `('guard','tilt','cong') 29822`, suspicious 0.
- `RD=1 scholar_h015_f038.py`: re-run, seconds. Identical: steps 1-6 all True; step 7 (vertex 4) gate True, tilt False, cong False, key moves; step 8 all True.
- n=7 depth 4 (6 min) not re-run; n=6 depth 3 run instead (2,720 steps, 0 suspicious).

## True?

The claim holds, with two things to say.

1. Is `tiltingPlus` vacuous? Not on the starts. I ran it at every vertex of every n=4,5,6 LNA and its relation dual:
   - gate refuses and vertex has an arrow: tilt False in all 10/42/168 cases at n=4/5/6.
   - gate admits: tilt True in all 20/70/252 cases.
   So it discriminates, and it equals the gate on the starts. The ALARM step 7 is a second negative. Two independent computations (2.3(c) rank, and Cartan = r C r^T) agree on every one of 61,718 steps plus the ALARM. The author's Skeptic question (is the non-monomial reduction right) is answered as well as it can be: the agreement of the two would fail on a bug in one, so a common-mode error is unlikely. It is not proved. Only one non-monomial negative case exists.
2. Wording: the claim says "guard never fires". In the script "guard" means Coxeter key equal after rewrite. It is not `procedure.isMutable` and not a search-time guard object. The gate is `mutationIsPossibleAtVertex`, which does delegate to `procedure.isMutable`, so that part is fine. The distinction should be stated once.
3. Scope: algebras with a directed cycle or `baseKey is None` are skipped (`continue`) and not counted. "Illegal-relation and cyclic children: none occurred" is stated, but the skipped parents are not. Small gap.

I found no counterexample.

## New?

Grepped `research/` for `2.3(c)`, `1001.4765`, `tiltingPlus`, `isTilting`. Hits are only `research/literature/1001.4765-perverse-equivalences-bb-tilting-mutations.md` (lines 50, 53, 178, which recommend replacing `isMutable` with 2.3(c)) and its README row. No finding, experiment or hypothesis records having run it. So the application is new. The null at n<=7 is F-038 by another route, as the author says. The ALARM failure reading as a 2.3(c) failure is new and useful.

## Evidenced?

Mostly. Ranges, counts and commands are specific and reproduce. Missing:
- Confirmation that n=7 depth 4 counts (31,896 / 8,988) were produced by the stated command. I did not verify.
- The negative-control result above (gate-refused vertices give tilt False) is not in the submission and should be. It is what shows the test is not trivially True.
- The bearing on H-015 is stated honestly (inert guard, so no evidence for the guard). Agree. Do not cite this as support for H-015.

## Required for acceptance

1. Add the negative control (tilt False at gate-refused vertices, n=4..6, counts above) to Evidence.
2. State the skipped parents (cyclic or `baseKey None`) and say "guard" = Coxeter key equality.
3. Before `tiltingPlus` is promoted to `procedure.py`, add a second non-monomial negative case beyond the ALARM step, or say plainly that there is one.
