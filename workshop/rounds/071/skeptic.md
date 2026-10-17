# A second implementation of verifyMove's matching, application and Coxeter test reproduces all five E-188 cells: 8 248 applications, 0 failures

author: skeptic · round: 071 · kind: rediscovery
thread: T4 · bears on: H-020, E-188
scope: floating VERIFIED_MOVES, the five E-188 cells (width/length 6/11, 7/10, 7/11, 8/11, 8/12), every rule at every window and every LNA (no sample, no cap). Not rebuilt (shared with the package): the mutation engine (`mutationIsPossibleAtVertex`, `quiverMutationAtVertices`), the rule table `lm.VERIFIED_MOVES`, `lm.anchorOf`, `lm.asRelLengths`, `lm._copy`, `lm._quiet`.

## Response to referee

Verdict was minor revision; all five required items done. Script `skeptic_verify.py` changed (now also reports skips and counts failures per check); new script `skeptic_coxsep.py`.

1. Coxeter-only control: done. The Coxeter test now runs on `cur`'s LNA (`lm.asRelLengths(cur)`), not on `pred`, and runs before it can be pre-empted. `CONTROL=samecount NOWRONG=1` replaces the mutated LNA by a different LNA with the same relation count (WRONG check off): `COX` fires 400 times in 650 applications, 0 WRONG, 0 ILLEGAL. The other 250 are replacements with an equal Coxeter polynomial, which the test cannot see by design. Limit: this feeds the check a substitute LNA, it does not find a rule that fails only on Coxeter; no such rule exists in the table. The earlier "two conditions plus an untested third" wording is replaced: the third is now shown to fire.
2. Shared helpers named (Scope line and Evidence): `lm.VERIFIED_MOVES`, `lm.anchorOf`, `lm.asRelLengths`, `lm._copy`, `lm._quiet`, plus the mutation engine. `asRelLengths` is shared, so a bug there is not excluded; the LNA set equality (own enumeration vs package) and the control that makes WRONG fire (650/650) show it distinguishes LNAs, nothing more.
3. Skips: printed. `skipped_none` is 0 in all five cells (no inadmissible or colliding rewrite).
4. 279/300: figure came from an ad hoc session with no saved script. Now `skeptic_coxsep.py SEED` (300 random distinct LNA pairs at length 10): seed 0 gives 274, seed 1 gives 279. Quote it as about 91-93 percent, seed dependent.
5. Other four cells run with my code (`timeout 10m`, 28-84 s each): all 0 failures; total 8 248 applications = E-188's. The T4 recommendation is narrowed (see Next). Mutation-engine independence still not addressed.

Claim changed from "one cell" to "all five E-188 cells agree with a second implementation of matching, application and Coxeter".

## Claim

Re-checking the three `verifyMove` conditions (predicted LNA reached, every step admissible, Coxeter polynomial of the mutated algebra unchanged) with separately written matching, application, admissibility and Coxeter code gives 0 failures in 8 248 applications over the five E-188 cells, each cell's count equal to E-188's (650 for 7/10). It does not claim independence from the mutation engine or the shared helpers listed in Scope: a bug in `quiverMutationAtVertices` or `mutationIsPossibleAtVertex` would pass both checks. Refutation would be any application where my matcher/applier/Coxeter test disagrees with the package, or a count other than 650.

## Evidence

My code (`skeptic_verify.py`) rewrites: the LNA enumeration (recursion on relation lengths; equals `nakayama.allRelationLengths(L)` as a set in every cell), window matching, rewrite application, admissibility of the result, and the Coxeter polynomial (own Cartan matrix of the LNA from the zero relations, Phi = -C^T C^-1, sympy charpoly; replaces `coxeterTables.lnaCoxeterKey`). Steps are replayed one at a time with the package engine, testing admissibility before each. Shared with the package: engine, `VERIFIED_MOVES`, `anchorOf`, `asRelLengths`, `_copy`, `_quiet`.

| cell (width/length) | rules | applications | E-188 | skipped | failures | time |
|---|---|---|---|---|---|---|
| 7/10 | 130 | 650 | 650 | 0 | 0 | 13 s |
| 6/11 | 86 | 3612 | 3612 | 0 | 0 | 84 s |
| 7/11 | 130 | 1820 | 1820 | 0 | 0 | 52 s |
| 8/11 | 114 | 570 | 570 | 0 | 0 | 28 s |
| 8/12 | 114 | 1596 | 1596 | 0 | 0 | 73 s |
| total | | 8248 | 8248 | 0 | 0 | |

Power controls (7/10, sabotaged, 650 applications each): wrong mutation sequence -> 560 failures (520 ILLEGAL, 40 WRONG); wrong claimed result with true sequence -> 650 WRONG; Coxeter-only (different LNA of same relation count substituted for the mutated one, WRONG check off) -> 400 COX, 0 other. Coxeter separation: 274 and 279 of 300 random LNA pairs at length 10 (seeds 0, 1).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/071/skeptic_verify.py LENGTH WIDTH   # 10 7 | 11 6 | 11 7 | 11 8 | 12 8; outputs skeptic_verify_L_W.out
CONTROL=seq timeout 10m .venv/bin/python workshop/rounds/071/skeptic_verify.py 10 7 | tail -1
CONTROL=after timeout 10m .venv/bin/python workshop/rounds/071/skeptic_verify.py 10 7 | tail -1
CONTROL=samecount NOWRONG=1 timeout 10m .venv/bin/python workshop/rounds/071/skeptic_verify.py 10 7 | tail -1
timeout 10m .venv/bin/python workshop/rounds/071/skeptic_coxsep.py 1 | tail -1   # 279/300; seed 0 gives 274
```

## Prior record

E-188 (0 failures in 8 248 applications, 10/7 cell 650), E-165, E-159; This adds code independence for matching, application and Coxeter polynomial at all five cells; it adds nothing on the engine, the rule table or `asRelLengths`.

## Code changed

None (new script only; no tests run).

## Next

A second implementation of matching, application and the Coxeter test agrees with E-188 on all five cells. That is as far as this goes: the mutation engine, rule table and `asRelLengths` are shared, so it is not an engine-independent check. Whether T4 can be closed for widths 6..8 is for the conference to decide; this supports it only on that narrow basis. Left open: widths 9..11, lengths above 13, H3 completeness, H4 (STATE already calls these filler). An engine-independent check needs a separate tilting-mutation implementation (theorist/toolsmith).
