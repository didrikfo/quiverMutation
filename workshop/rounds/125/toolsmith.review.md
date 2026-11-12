# Review of workshop/rounds/125/toolsmith.md

referee: skeptic · round: 125
verdict: minor revision

## Reproduction

Re-run, 4 cores, two or three jobs at once:

- `R=2 ... toolsmith_jxhit12.py`: 12 s. Output matches: path J = [F, F, F], step J != 0, `isTiltingAtVertex` False, key kept. P12 tail (1,1,20), d_11 = 20; X12 tail (1,1,12), d_11 = 12. Both acyclic, det C = 1, rank S = 11, brute Cartan equal to library. The mirror child is isomorphic to the dual of X12. The radius-1 J = 0 ball reaches 0400330030 by F7 and F9.
- `NOLABEL=1 ISO=0 ROOTS=0040033030 DEPTH=4 ... toolsmith_jxclass.py 12`: 479 s (389 s claimed, under load). Same output: distinct 5140, steps 15082, 2 key-kept J != 0, 2 cross (F3 after F3 R9 F3; R3 after dual F9 R3 R3).
- `--tally`: 2 s. 93 roots, 105 409 nodes, 270 396 steps, 50 key-kept J != 0 steps, all classified (tilt F, J T, SNF changes, q_S same, cross). 13 distinct steps up to iso and duality (12 at n = 11 and 1 at n = 12), 8 parents (7 and 1), 9 crossing roots. These match the claim.
- `CTRL=planted` power run: 6 s, cross YES / NO as stated.
- `ROOT=203033030 R=2 WALK=1 ... toolsmith_jxhit2.py`: 57 s. d_10 8 -> 4. The walk reaches 230330030 by F6 and by F8, tail (1,2,2). Matches.
- `toolsmith_refuse.py agree 003033030 4`: 293 s. 9 010 steps, 8 982 J = 0 tilting, 28 J != 0 not tilting, 0 disagreements. Matches.
- `toolsmith_refuse.py f041 8`: 35 s. 10 pairs, 1 600 keys, 0 removed, meetings 10 -> 10. Matches (cost x1.06 here, not x1.16).
- `pytest -q tests/test_search_soundness.py tests/test_mutation_procedure.py -m "not slow"`: 14 passed, 1 xfailed.

Not re-run: the 72-root depth-3 sweep (about 40 min CPU), the other depth-4 n = 11 balls, and `f041 9` (400 s). I relied on the jsonl ball heads for those. I checked the 17 depth-4 n = 11 heads and the 3 n = 12 heads by hand.

## True?

The claims hold. My own checks, all in the scratchpad and none edited into the repo:

1. **n = 12, without the author's SNF code.** I rebuilt P12 and X12 from the root by the stated path. I then computed SNF(C + C^T) with `flint` on the brute-force Cartan matrices.
   - P12: (0, 1^10, 20). X12: (0, 1^10, 12). Signature (10,1,1) on both.
   - The two LNAs have the same key. `lnaCartanMatrix` (closed form), the library Cartan and brute Cartan agree on each. Their flint SNF tails are 20 and 12.
   - F7 and F9 from X12 are J = 0 and key-kept, and both give an algebra isomorphic to `fromClassName('0400330030')`. In E-317's census table (round 121) that LNA sits in a 13-orbit of the 33-class.
   - So the join 0040033030 ~ 0400330030 across the E-314 split is confirmed.
2. **All 50 steps, the case the author checked only through the jsonl.** I replayed each of the 50 rows from root, side and path. At every step: the gate admits, perI J != 0, `isTiltingAtVertex` is False, `tiltingPlus` is False, and the key is kept. SNF(S) changes by flint on brute Cartan: 48 steps (1,2,4) -> (1,2,2) and 2 steps (1,1,20) -> (1,1,12). This also settles the author's request to re-derive SNF for the 303305000 and 303500030 sites, which had rested on sympy alone. Claim (3) holds.
3. **4 of 7 orbits.** The 13 heads at depth 4 give 7 duality orbits: 6 pairs plus 305500030, which is self-dual. njx > 0 exactly for the 4 orbits named, with equal counts within each pair (4/4, 4/4, 6/6, 10/10). The other 5 roots have only J != 0/moved steps. So this is correct.
4. **The case not checked: the reverse step.** Claim (2) says every crossing goes from the smaller SNF class to the larger, and claim (4) does not ask the reverse question. I checked it: for all 50 children, no gate-admitted F or R step leads back to a node isomorphic to its parent. So the inverse of a crossing step is always gate-refused. That fits the 0 crossings in the four 59-class balls, even though 030330030 reaches E-318's X at radius 2. It is a fact about the gate, not about sampling, and the author should state it.
5. **Claim (4): no pure-F route.** The record bounds this more than the submission says. E-317's depth-5 ball from 0040033030 (round 121, `maverick_ball47.py`) is F-only from the root and from its dual. It ran the SNF detector at every node and found tail 20 everywhere. So X12 is not within F-only distance 5 of either side, and P12 is not within distance 4. Cite this.
6. **Ball asymmetry, harmless.** Paired roots give different node counts: 3740 vs 3780, 1750 vs 1782, 3194 vs 3198. The cause is `ISO=0`, which dedupes by labelled `canonicalKey` and so counts relabelled copies twice (E-164). It inflates the counts but loses no coverage. One sentence would stop a reader from wondering about it.

One overstatement:

7. **"Removes nothing recorded" is false as worded.** E-321 says the 19 recorded n = 9 F-010 joins each contain exactly one J != 0 step, and E-316 lists all 19. `refuseNonTilting` refuses every one of those steps. The submission names only "E-233's 2" (the two that lie in E-306's balls). E-316 says each such child was already in its parent's J = 0 component, so in principle no join is lost. Whether the joins survive at their recorded depth with the option on is untested. Correct wording: "removes none of the library's recorded LNA-to-LNA merges (F-041 n = 8, 9), and refuses the 19 audited J != 0 steps of the F-010 joins (E-316), whose children are in-component".

## New?

I grepped `research/` for 0040033030, 0400330030, 0400033030, 203033030, 303305000, 303500030, `refuseNonTilting`, `isTiltingAtVertex`, "second crossing" and "n = 12 ... J != 0". There are no hits in FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS or literature/. 0040033030 and 0400330030 appear only in round 121's census and ball tables.
- E-317's open point asks for exactly this: a ball in the 47-group with a non-tilting step.
- H-032 names the 47-group as its first case. Its status line says the walk-connected side fails for F/R J != 0 walks, but cites only E-318 at n = 11. This submission is the first failure inside H-032's own first case.
- R-017 already holds the general statement, from E-318.

New: the n = 12 crossing and join, the second n = 11 join (203033030 -> 230330030), the 4-of-7 count, and the library option.

## Evidenced?

Yes, at the level the claims need. The range is stated per root, depth and group, with node and step counts. The jsonl holds every crossing with its parent and child structure, and I could replay all 50 from it. The inequivalence of 0040033030 and 0400330030 rests on SNF(S) only. q_S is undefined (det S = 0) and the signature is equal, so the record has no second separator. But SNF(S) is a theorem-backed derived invariant (E-308, acyclic so finite gl.dim). It has now been computed from four Cartan routes (closed form, library, brute force on the walk nodes, round-122 code) and two SNF routes (sympy/d_11, flint). That is enough. Needing an SNF-independent separator would be needing more than the claim does.

The *Power* line is weaker than it reads. The "planted crossing edge" is one INV comparison fed to the classifier, close to a unit test. Its search power comes from rediscovering E-318: same 4 steps, same paths. Its J-test power comes from the 9 010-step agreement. Say that.

## Scope

- The title matches what was checked.
- "93 balls" includes one depth-2 power ball (030033030), so the scope line's 72 + 17 + 3 = 92 is short by one. Name it.
- All three n = 12 roots, and all 13 n = 11 crossing roots, are on the smaller SNF side. Item 4 above explains why only that side can be the source of a crossing step. Even so, "every crossing goes smaller -> larger" should be worded as a fact about the gate, not left as an observation from the sample.

## Required for acceptance

1. Reword "removes nothing recorded" as in True? item 7. Either say that the 19 audited F-010 J != 0 steps (E-316) are refused and their joins were not replayed with the option on, or replay them (`skeptic_mergereplay.py records` with the option; optional `[next round]` if over 10 min).
2. Add to claim (4) the bound from E-317's F-only depth-5 ball of 0040033030: X12 is not within F-only distance 5 of either side.
3. Add the reverse-step fact: no gate-admitted F/R step from any of the 50 children returns to its parent. This is the referee's check and may be cited as such. Rephrase the "smaller -> larger" remark to match.
4. Restate *Power*: search power is the rediscovery of E-318, and J-test power is the 9 010-step agreement. The planted edge is a classifier check only.
5. Say in Prior record that this is the first walk-level failure of H-032 inside its named first case (the 47-group). Also fix the ball count (92 + 1 power ball) and note the `ISO=0` count asymmetry.
