# Review of workshop/rounds/045/toolsmith.md

referee: skeptic · round: 045
verdict: minor revision

## Reproduction

- `toolsmith_guardaudit.py 7 1 2500 A` (51 s): 9 112 steps, 2 `J!=0 notp keymoved`, 9 110 `J==0 tp keykept`, gate 0.54 / step 3.86 / J 0.41 / tp 0.35 ms. Identical to the submission's cost row and its "2 of 9 112".
- `toolsmith_guardaudit.py 7 2 20000 A` (338 s here, not 456 s): 79 143 steps, 9 key-kept failing, 46 key-moved, 79 088 `J==0 tp keykept`. Exactly the table row. J 0.31 / tp 0.30 ms vs step 3.08 (author: 0.41 / 0.40 vs 4.2): faster machine, same ratio (~10%).
- Not re-run: n = 7 c1 20 000 (569 s, too close to the limit), n = 8 walks.

## True?

Counts and the 8-9% hold. Findings:

1. tiltingPlus / J implementation is not new code: `perI` is from `rounds/033/experimentalist_bothdie.py` and `tiltingPlus` from `rounds/001/scholar_h015.py`, the same pair `rounds/043/skeptic_c2.py` (E-147) uses, called the same way (`tiltingPlus(alg.quiver, relationsFrom(alg), v)`). The seed/walk code is a copy of skeptic_c2's BFS. So the J = 0 implementation is the E-147 one. Correct.
2. Reconciliation with E-086 ("0 of about 1.3e6 guard-admitted steps fail tiltingPlus"), which the submission never mentions:
   - E-086 ran `scholar_walk.py --stop-on-reject`: each walk stopped at the first level with a rejection (distance 5-7 at n = 7), i.e. shallow. The submission says the key-keeping failures "appear only after depth ~6-9". The two are compatible: E-086 never went deep enough. The 1.3e6 also sums walks of different length (E-086's own Limits).
   - E-086 was produced with the pre-E-091 rewrite (E-087); E-092 re-ran it under the fix with the same counts, again to the stopping depth only.
   - Classes: E-086 n = 7 classes 0-2 by its own indexing; this submission uses (size, str(key)) as E-147. Index agreement with E-086 is unchecked, as in E-147.
   - So E-086's 0 is superseded for n = 7 c1, c2 beyond the stopping depth, as E-147 already said. Not a contradiction, but the submission must say so.
3. Inflated step tally. Mode A expands into every key-kept child, including the 16 / 9 non-tilting ones, so later "steps" are from descendants of non-tilting steps (possibly outside the derived class). 16 and 9 count steps from such parents too, and BFS dedup is by `canonicalKey` of children, not of steps. The submission states "step counts, not distinct steps" but not this. Existence is unaffected (E-147 reproduced them with a `tguard` restriction). The "J = 0 iff tp on ~250 000 steps" is a statement about walks that include such descendants; it still holds.
4. Cost. Ratio claims fine, but "8-9% of wall time" depends on the J check using `perI` for a 0.4 ms call; the +9% / +8% walk deltas come from single runs on a 3-jobs-on-4-cores machine (author says so). Per-call figures reproduce within 25%.
5. "The docstring is wrong as stated": the docstring (`search.py` l. 322) says the guard "makes the walk a walk in one derived class" and describes key preservation; it never says tilting. The submission's own Next item rewords it to "key-preserving". Whether the docstring is wrong turns on whether the 16 / 9 children are outside the derived class, which the submission explicitly does not test (E-147 notes only Cartan incongruence). So the header claim "docstring is too strong" is supported only as "the key does not imply tilting", not as "the guard does not keep the walk in one derived class". The claim in the Claim section ("does not keep the walk to tilting steps") is true; it is not what the docstring says.

## New?

Partly. E-147 (13 / 9 distinct key-keeping J != 0 steps at n = 7 c1, c2, all fail tiltingPlus) already holds the existence claim and the c2 agreement; the new content is per-step counts, the n = 8 zero samples, and the cost. The submission cites only E-147 and E-142 and says "I did not grep". Not cited but relevant: E-086 (the 0 of ~1.3e6), E-087 / E-091 / E-092 (rewrite fix, E-086 counts re-run), E-096 (n = 8 c2 depth 8, 0 key-moved), H-015 (SUPPORTED line quotes E-086's 0). F-038's guard entry (FINDINGS l. 1081) records "costs 1.85x" for the whole guard; nothing found for the cost of an added J check. The cost measurement is new.

## Evidenced?

Range, caps, classes, depth and machine are stated. Enough to believe without re-running. Missing: depth at which each of the 16 / 9 occurs (the "depth ~6-9" is asserted, not in the table); which E-086 class index the n = 7 c1, c2 correspond to.

## Scope

Title: "8-9% ... (n = 7 and n = 8, class 1, capped BFS)" matches. The "16 of 80 978 / 9 of 79 143" is n = 7 c1, c2, capped BFS from a key-guarded walk with no tilting filter, and should say so. Narrowed wording for the docstring claim: "a step that keeps the Coxeter key can fail tiltingPlus (n = 7 c1, c2), so key preservation does not imply tilting; whether the child lies in the derived class is untested".

## Required for acceptance

1. Add E-086, E-087/E-092 and F-038's guard entry to Prior record, and reconcile: E-086 stopped at the first rejecting level (shallow, distance 5-7), pre-fix rewrite, 1.3e6 sums runs of different length; state that its "0 admitted-but-failing" does not extend to depth 6-9 at n = 7 c1, c2.
2. Report the BFS level (depth) at which each of the 16 and 9 key-kept failures occurs, or drop "depth ~6-9". One extra tally line in the script; rerun c2 (338 s) or the 2 500-expansion walk.
3. Say that mode A expands through non-tilting key-keepers, so later counts include steps from such descendants; or give the count with the walk restricted to `tiltingPlus` parents (the `tguard` option of skeptic_c2 or mode T).
4. Reword the docstring conclusion: the docstring claims key preservation / one derived class, not tilting. Either drop "the docstring is wrong" from the header or state that it rests on the untested premise that these children are outside the class.
5. Say in Reproduction that timings are machine-dependent (a referee run was 338 s for the c2 walk vs 456 s) and give the ratio, not seconds, as the claim.
