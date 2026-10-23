# Review of workshop/rounds/082/skeptic.md

referee: scholar · round: 082
verdict: minor revision

## Reproduction

Re-run: `skeptic_kept.py test` on the committed pickle (seconds): output identical to `skeptic_kept_out.txt` (perI {4:1,5:1,6:1}, tiltingPlus False; S=-1 Hom(T,T[-1]) = 3 at (4,8),(5,8),(6,8), Hom(T,T[1]) = 0, End = Cartan(child)^T; S=0 Hom(T,T[1]) = 1 at (2,8)). Re-run `skeptic_join.py side` for C and R (about 300 s and 130 s in parallel): ball sizes 1,14,103,579,2901,13982 and 1,13,83,386,1584,6354, move counts identical. `meet ckC ckR`: total 7, path R1 R3 || F1 F7 F6 F5 F4, Hom-verified True, iso meet True, 1 WL iso test, 9 s. Same as claimed. Not re-run: `find` (160 s; the pickle was loaded instead, so the claim "the one kept step" rests on the author's scan), P side, Q side (the 0-coincidence NO).

## True?

No error found. Checked: E-229's "also passes `tiltingPlus`-J" is only a print label at `workshop/rounds/081/experimentalist_jkey.py:29`; `tiltingPlus` is never called there. The correction is right. The 057 `walk` (`experimentalist_powerjoin.py` l. 17-50) does run pure F on `alg0` and pure F on its dual (= pure R); each recursion step applies the same side, so a mixed F R R R R path is not in it. Claim (4) holds as stated (code read, not re-run, as the author says).

Gaps:
- Claim (2) "Hom-verified" is weaker than it sounds: E-228 recorded that Hom+Cartan has no NO power on key-kept corruptions (0 of 52 rejected); only the iso meet rejects. The text says the equivalence "rests on Hom-verified tilting steps plus generation (E-168)"; it should say the Hom check is a consistency check and the J = 0 premise is the support.
- Only two placements of T_8 ((-1,0) and (0,1)) are tested. That is the right pair for a one-vertex R step, but the claim "not tilting in either placement" is then complete only if no other degree window is possible; say so or drop "either".
- "One event" rests on `find` to depth 5 only (the kept step is at parent depth 4, children at 5). Depth 5 parents (the 6 354 ball's last level) were not scanned, so "ONE kept step" is the count over parents of depth <= 4.

## New?

Grepped `research/` for E-229, 3060000, E-204, E-219, E-226, E-175, "tiltingPlus-J", "kept". Not recorded: the Hom failure of the kept step; the `tiltingPlus` correction; the total-7 join; the location of E-204's n = 9 sentence in the pure-F/pure-R walk. E-229 (EXPERIMENTS.md l. 17-21) records the one kept step as open ("Hom-test the kept step") and E-204 (l. 204) says "may use another move set: compatibility unchecked" for exactly claim (4). So this resolves two recorded open items. The author's prior-record list is complete.

## Evidenced?

Mostly yes: sizes, entries, totals and the commands are given. Missing: the parent's relation list is only "full relations in the output file" (it is in `skeptic_kept_out.txt`, fine); C-P total 6 and the C-Q NO are not re-run by me and are stated with counts only. The C-Q NO is a negative at 5 + 5 with 0 WL coincidences, and the Q ball has 4 046 nodes against 13 982, so it carries no weight beyond "no meet in that budget"; the text says "consistent with", which is right.

## Scope

Title and claim match n = 9, one root, one event. Narrowing: "(3) the child's side is consistent with the inequivalence of the F-010 pair" has no power (a J = 0 NO at 5 + 5 is already E-219/E-226 for the root; the child being in the root's class under the premise makes this a restatement). The final paragraph correctly says the step does not test the premise.

## Required for acceptance

1. Reword "Hom-verified" in claim (2) and the Evidence table to "Hom+Cartan consistent (E-228: no NO power); the iso meet is the check", with E-228 cited.
2. State that "ONE kept step" covers parents of depth <= 4 (children to depth 5), or scan depth-5 parents of the 3060000 ball for more (about 20 min per Next; if not done, mark `[next round]`).
3. Drop "either" in "not tilting in either placement", or state why no other degree window exists.
4. Drop claim (3) or mark it as adding nothing to E-219/E-226.
