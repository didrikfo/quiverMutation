# Review of workshop/rounds/133/experimentalist.md

referee: maverick · round: 133
verdict: minor revision

## Reproduction

All with `.venv/bin/python -W ignore`, 4 cores shared.
- `experimentalist_band.py 13 2`: 90 s, output **identical** to `experimentalist_band_out.txt` (0 / 5390 (4,3), 0 / 2505 (4,4), 1749 / 4147, 300 / 1547, 1492 / 2134; 8171 / 59917 overlap-2 rows; every table line the same). W = 14 not re-run (364 s, same code).
- `experimentalist_band.py 13 3` (one more free vertex each side, my variant): same 541 450 edges and same counts, so P = 2 padding does not starve the rules.
- `experimentalist_reach.py 13 2`: 19 s, 37 635 rows reached, none at overlap <= 1, 830 in (2, 3, True). Matches.
- `experimentalist_cox.py band 9 2`: 61 s; 28657 rows, 210 polynomials, 6 of 40 clean (4,3) outside S_13, YES control 84 / 0 outside. Matches.
- Engine, `gap2fix.py 4 1 540 "8,4:10,3:11,3;6,2:8,4:10,3" --jobs 2`: 166 s. Control lowers 2 -> 0 by [-13, 10, -13] (71 s, 4 lowered); clean (6,2) cell min 2, 24 reached, 0 lowered. Both match `ctl_out` / `one_b1_out`.

## True?

Reproduces. I tried to break it in two places:

1. **Is the closure exhaustive?** It leaves out `freeMoves`, which is a library derived move ("its value is bridging"). I added it (add or delete a 2-arrow relation, kept in the band) to the W = 13 closure (scratch script, 99 s). Still **0 of 5390** clean (4,3) cells reach overlap <= 1, and 8171 of 59917 overlap-2 rows, as before. So the null holds even under free moves. It still does not hold for routes that leave the band, or for single mutations outside the table. The author's scope says so, but the title's "at any length" and "the closure says none exists" should carry "in the W <= 14 band".
2. **Is "clean" the same as E-336's?** Yes. E-336 `--clean` needs base = pair overlap and no bystander starting at a pair start ± 1, which is {7, 9, 11}. `clean_cells` needs maxOverlap = a - 2 and no start in {s-1, s+1, s+3} = {7, 9, 11}. The band version allows any number of further relations, so it is a superset of E-336's cells, which is fine. It also counts 471 isolated (vacuous) cells, which the text discloses.

A wording error: the claim says the lowering rows have a run of three that "contains a gap-1 (3,3) pair". `band.py` (`pairtypes`) and `reach.py` (`has33`) both test the **whole row** for a (3,3) pair, not the run. The deduction "so it forces (11,3)" is correct for a (3,3) pair inside the run: a predecessor would have to start at 7, which clean forbids, or be a gap-2 (x,4). But that is not the statistic in the table.

There is also a simpler explanation of the characterisation. Of the 414 `VERIFIED_MOVES`, exactly 20 lower a window's own overlap from 2 to <= 1, and **all 20** have a gap-1 (3,3) pair in their before-pattern. Double mutations add none (0 of 23096 overlap-2 rows without a (3,3) pair lower). So "lowering at overlap 2 only goes through a (3,3) pair" is mostly what the rule table already contains. The new content is the reachability statement: a clean (4,3) never reaches such a row. That is not vacuous, because 64 rules act on a (4,3) gap-2 pair.

## New?

The engine cells, the band closure with rules, and the Coxeter sidebar are new. The (3,3)-run characterisation is the band-closure form of **F-022** (a run of three dissolves, a pair is frozen) and of **E-152 / E-230 / E-241** (for (8:3)(9:3), only a bystander sharing >= 2 arrows, i.e. a run of three, lowers). **E-271** is the same "no control found" situation for (8,4)(10,4). None of these is in the Prior record. The (11,3) control is a run-of-three control of exactly the E-152 kind.

## Evidenced?

Yes for the band and Coxeter numbers: band, cells, denominators, isolated counts, time and seed are all stated, and I did not need to re-run anything to read them. The engine table is specific (cells, caps, decided count 8 + 6 + 1 = 15). The *Power* line holds as worded, "harness sense". The engine control shows that the engine sees (3,3)-run dissolution, not gap-2 dissolution, and the author says so. The Coxeter test decides only 6 of 40 clean cells, all with a second overlap-2 pair, so it carries no weight for the claim. It is presented correctly as a sidebar.

## Scope

Narrow two phrases. "none lowers under the library's derived moves at any length" should read "none lowers under doubles + VERIFIED_MOVES (+ free moves, referee) at any length, rows confined to bands W = 13, 14". "no clean control exists (in the closure)" should read "none in those bands".

## Required for acceptance

1. Add F-022, E-152, E-230, E-241, E-271 to the Prior record. Say that the (3,3)-run characterisation restates F-022's run-of-three mechanism on the band, and that the (11,3) control is of the E-152 kind.
2. State that all 20 overlap-2-lowering rules in `VERIFIED_MOVES` contain a gap-1 (3,3) pair, so the characterisation is partly the rule table's content, and the new part is reachability from clean cells.
3. Fix "a run of three ... contains a gap-1 (3,3) pair". Either say "row has a (3,3) pair and longest run 3", or recompute with the pair required inside the run (about 90 s at W = 13).
4. Narrow the title and claim as in Scope. You may cite the referee's free-move closure result (W = 13, 0 / 5390).
