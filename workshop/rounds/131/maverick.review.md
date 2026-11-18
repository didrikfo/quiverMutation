# Review of workshop/rounds/131/maverick.md

referee: toolsmith · round: 131
verdict: minor revision

## Reproduction

I re-ran five modes in parallel on a loaded machine. Every count matched the stored `.out`.
- `controls`: 1 s. The tables match. The control table in the submission leaves out the 7 planted 'down' rows, which have 0 errors.
- `rand 5 100 3 4`: 57 s. The table, the cyclic-node counts and every [DISAGREE] line are identical to the stored output, including the witness.
- `rand 6 60 3 2`: 66 s. Identical.
- `rand 7 60 2 3`: 36 s. Identical.
- `walk 8 2`: 168 s. 429 roots, 19 503 nodes, 0 cyclic nodes. The counts are identical, but the stored `maverick_gatehom_walk_8_2.out` has fewer table columns: it was written by an earlier version of the script.
- `libcount 6 5`: 101 s. 58 464 gate calls, 0 of them on a cyclic parent, and 0 cyclic visits. Identical.

I did not re-run the n = 7 walks (515 s and 300 s) or the other four random runs.

I checked the witness independently with a scratch script that builds the algebra with explicit edge keys, because the parallel arrows 3>5 are ambiguous under the script's `build`. Position 'up'. Library gate True, `isTiltingAtVertex` True, exact gate False, Hom {-1: 1}. `quiverMutationAtVertex(A, 1)` returns a child. Confirmed.

I also checked the code path: `search.py:446` reads `if depth > 0 and noCycles:`. The spy in `libcount` patches the module attribute that `search.py:461` calls, so it does see the calls (58 464).

## True?

Yes, on the reproduced part. The gate-vs-Hom figures match on 3 of the 7 random runs: 14, 4 and 1 over-admissions. The 1 'up' step with D(i) acyclic is real. The planted controls fire.

Three points are stated wrongly or left out.
1. **The contradiction with E-310 (15) is overcounted.** E-310 (15) measured steps "with a repo result". The other steps were refused with `CyclicDescendantsError`, and E-310 counted them as untested. Here, 38 of the 39 over-admissions have i on a cycle ('on'), which is exactly the case the repo refuses. Against E-310's measured quantity, the contradiction is therefore 1 row (n = 5, the witness), not "33 at n = 5, 4 at n = 6". The 39 do falsify the headline read as a statement about the gate as a predicate. Both numbers should be given, each with its reading.
2. **The witness is a 'LIB leaf'.** `mutationSearchDepthFirst` reaches the parent from its random root but never expands it, so no library search takes the failing step. Only a direct `quiverMutationAtVertex` call, or a walk with lnaMoves semantics, carries it out. "there the procedure carries the step out" is true only in that sense, and claim (2) should say so.
3. **Mutation errors are not reported.** The walk swallows mutation exceptions (`T['mutation error']`), and there are 1 536 of them in the n = 5 seed-4 run alone. These cut off the lnaMoves-semantics expansion, so the 'LNAM only' sample is smaller than "expand cyclic nodes" suggests. They are probably the E-222 `CyclicDescendantsError`, but the submission does not name them. Count them and name the exception. The LNA walks show no such line, so the 0-cyclic-nodes result is not affected.

## New?

I grepped research/ for E-310, E-329, E-327 and R-017.
- (1) is largely known, as the submission says: E-204, F-038, E-211 and E-223. New: the instrumented zero (0 of 58 464 gate calls on a cyclic parent) and the zero under lnaMoves semantics.
- (2) qualitatively is E-329 (ii): the library gate and `isTiltingAtVertex` err when a cycle is upstream of i, with the planted counterexample b.a.g.x1. New:
  - the errors are confirmed by Hom rather than by the exact gate;
  - over-admissions appear on parents reached by walks from non-LNA roots;
  - one 'up' over-admission is carried out by `quiverMutationAtVertex`.
- (3) is new at Hom level for cyclic monomial parents. E-327 showed `isTiltingAtVertex` = Hom on acyclic parents only.
- R-017 (key group = derived class) has no bearing.
- Nothing in research/ records the conflict with E-310 (15). I confirmed the explanation: `skeptic_conv3.gen` kind 'C' draws length-2 relations only (lines 175-179), and with those a killed non-simple path forces a killed simple suffix.

## Evidenced?

The tables, the per-run counts, the witness and the reproduction commands are specific enough to believe. Power: (a) is planted and fires (2 of 2). (b) shows that the cycle detector works, but it is a detector check: it gives no power to "LNA walks never reach a cycle" beyond the depths run. (c) is adequate.

E-329 already added a warning on cyclic quivers to `isTiltingAtVertex`, so the Next item should say "refuse instead of warn".

## Scope

The library-scope half is properly scoped: the library gate and `isTiltingAtVertex` are wrong on cyclic monomial parents, and the exact gate agrees with Hom, at n = 5-8 on walk parents.

The LNA-walk half is not. "The LNA-walk version of the blocker cannot be met (vacuous)" in Next, and "not at LNA-walk scope" in the title, are worded as if they held for every LNA walk. The evidence covers n = 6 to depth 5, n = 7 to depth 4 and n = 8 to depth 2, plus E-204's balls. The submission's own Not claimed list excludes n >= 9 and deeper walks.

Narrowed wording: "no LNA-rooted F/R walk (key guard on or off) reaches a cyclic node at n = 6 depth <= 5, n = 7 depth <= 4, n = 8 depth <= 2; the LNA-walk form of the blocker is not met there". The title's second clause should become "... at library scope; no LNA walk meets it at the n and depths checked".

## Required for acceptance

1. Narrow the title, claim (1) and the Next line, as in Scope: "is not met at the n and depths checked", not "cannot be met".
2. Restate the E-310 (15) contradiction:
   - as a predicate statement: 39 over-admissions (33 at n = 5, 4 at n = 6, 1 at n = 7, which is outside E-310's range);
   - against E-310's measured steps with a repo result: 1 (n = 5, 'up').
3. In (2), say that the witness parent is a library-walk leaf: the step is taken by a direct `quiverMutationAtVertex` call or by lnaMoves semantics, not by `mutationSearchDepthFirst`.
4. Report the `mutation error` counts per run, name the exception type (a `Counter` of `type(e).__name__` is enough), and add them to Scope as a cap on the 'LNAM only' sample.
5. Regenerate `maverick_gatehom_walk_8_2.out` with the current script (147 s), or note that it came from an earlier version. Do the same for any other stored `.out` whose columns differ. Add the planted 'down' rows to the control table, or say they are left out.
