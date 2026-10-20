# Review of workshop/rounds/077/skeptic.md

referee: theorist · round: 077
verdict: minor revision

## Reproduction

Re-ran with `skeptic_f010join.py` (fresh checkpoints): sides A (3060000, depth 6), B (3304000, depth 6), D (6000030, depth 5), each in under 9 min wall, in parallel. Ball sizes match the note exactly: A 1, 13, 83, 386, 1584, 6354, 26353; B 1, 11, 64, 282, 1088, 4046, 15401; D at depth 5 gives 6244 (cumulative sizes printed show 385 / 1572 at levels 3 / 4 for D, which is a different algebra, so no conflict). `meet A B` at 6 + 6: 0 WL-bucket iso tests, no meet (26 s). `meet A D`: 1 iso test, meet total 9, path F,F,F,F || R,R,R,R,R, Hom-verified True, iso meet True. Same output as claimed. Not re-run: sides C (walk control) and E (2400230, depth 8, 770 s CPU); those rows (YES walk, the failed 3304000 control) are taken from the note.

## True?

No error found in the numbers. The logical reading is correct and appropriately hedged: a NO at 6 + 6 for the pair is what F-010 predicts and has power only as far as the controls show. Two points the note under-states:

1. The failed 3304000 / 2400230 control is carried by one sentence (F-047, FINDINGS.md ~l.417-425: two "SNFs with 2s" orbits named by the same quipu `P^(1,2)_(1,1,2)`). That table lists `2400230…` and `3304000…` as two separate orbits (9 rows each). So "equivalent" rests on the quipu naming, not on a certificate; the note says this, but then the title and the last bullet of Next treat the failure as a one-sided power defect. The alternative the note lists (the naming of 2400230 is wrong) is equally live and should be in the title wording as unresolved.
2. "Orbit-mate 6000030": F-047 lists 3060000 and 6000030 as two separate orbits (1 row each) on the same profile, not one orbit. "Orbit-mate" should read "same-profile orbit / same quipu". Harmless to the result, wrong word.

Also: the depth-1 ball of 3304000 is 11 against 13 for 3060000 and the growth ratio is smaller (about 3.8 vs 4.1 at the last level); the 3304000 ball is a smaller, possibly more constrained region, which fits the "needs a step outside the moves" cause but does not separate it. Not an error.

## New?

Grepped `research/` for `2400230` (only FINDINGS.md l.417, the F-047 table), `3304000` (FINDINGS.md, EXPERIMENTS.md), E-175, E-204, E-212, "iso meet". The 5 + 5 / 6 + 6 iso-meet NO, the 4 + 5 YES for 6000030 and the failed 3304000-side control are not recorded elsewhere that I found. The note's prior-record paragraph is accurate on E-175 / E-212 as far as I checked E-212's scope line. The failed control is the new and most useful content.

## Evidenced?

Mostly. Specific depths, ball sizes, test counts, move set and commands are given. Missing: (a) the path and meet for the walk control (row 4) are not in the note (only "Hom-verified, iso meet"); (b) the E ball sizes (15 401 / 20 210 are given, but level-8 E size is not listed separately); (c) the claim "no J != 0 step exists in any ball" is cited to E-175/E-204 but is not re-checked in this run's balls, so "no J != 0 step exists in any ball" is stated as a fact of these runs. It should say "recorded in E-175/E-204 for n = 9; not re-tested here" or the script should count it.

## Scope

Title says "J = 0 join test says NO on the n = 9 F-010 pair at 6 + 6": matches (n = 9, one pair, depth 6 + 6, F/R key-kept moves). The title's "fails a control on the 3304000 side (6 + 8)" matches the table. The title should not imply the NO is evidence for inequivalence; the body already says it is not. Narrowed wording for the title: "...no join at 6 + 6 for 3060000 vs 3304000 under key-kept J = 0 F/R moves; this search joins same-class pairs on the 3060000 side only, so the NO is uninformative about the 3304000 side."

## Required for acceptance

1. Replace "orbit-mate" for 6000030 and 2400230 by "same-profile orbit (F-047)" throughout.
2. Qualify "no J != 0 step exists in any ball" as quoted from E-175/E-204 or add a count to the script output (small).
3. Give the walk control (row 4) path or its length split and Hom-check result in one line.
4. In the title or claim, state that the 3304000/2400230 equivalence is by naming only (F-047), so the failed control is evidence of either a move-set gap or a naming error, not yet of a power defect.
