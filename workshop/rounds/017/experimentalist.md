# With the full-reduction patch the E-084 walks give the same counts everywhere they could be compared; the 10 n = 8 class 2 "key moved" steps lie beyond the part of depth 8 that fits in 10 minutes

author: experimentalist · round: 017 · kind: negative (partial re-run)
thread: E-084 / E-085 loose end · bears on: guard-redundancy question, `isTilting` promotion

## Claim

Running the E-084 guarded walk (`scholar_walk.py`, same arguments) with `arrowPaths.reduceAgainstPivots` replaced
by the full reduction of `theorist_fix.py` changes **no count of E-084 that I could re-run**: the gate-admitted
rejections (gate True, `tiltingPlus` False), the guard-admitted-and-failing steps (0) and the walk tables are
identical to the unpatched files for n = 7 classes 0-2, n = 8 classes 0 and 1, and n = 9 classes 0-3. For n = 8
class 2 (the loose end) the patched walk is identical to an unpatched control through 16 500 expansions (same
distinct algebras 42 028, same table, 1 rejection, 0 "M" lines), but the 10 "M" steps of E-084 appear only after
expansion 16 500 of 20 899 at depth 8, and the patched walk did not get that far in 10 minutes (17 058 expanded). So
the claim "the n = 8 class 2 loose end disappears under the patch" is **not** shown by this walk; it rests on
E-085's replay of the 10 saved parents (which I did not repeat). Not claimed: that the patched walk is closed or
complete at any n (all are stopped-on-reject or budget-stopped, as in E-084); nothing about n >= 10.

## Evidence

`experimentalist_walk.py` = exec of `workshop/rounds/014/scholar_walk.py` after setting `AP.reduceAgainstPivots`
to the E-085 full reduction (the library is untouched; `procedure.py` and `arrowPaths.isInIdeal` both call it
through the module, so both see the patch). `--unpatched` is the control; `MAXEXP=N` caps expansions
(deterministic) in place of the time budget. Flags as in E-084: `--stop-on-reject`, 540-585 s budget.

| run | patched | E-084 (unpatched) file | same? |
|---|---|---|---|
| n=7 c0 / c1 / c2 | rej 8 / 4 / 6; guard-tilt 4686 / 18457 / 28375; M 0 | same numbers | yes |
| n=8 c0 | rej 4; guard-tilt 3585; M 0; 2255 algebras | same | yes |
| n=8 c1 | rej 14; guard-tilt 22542; M 0; 13769 algebras | same | yes |
| n=9 c0 / c1 / c2 / c3 | rej 8 / 4 / 12 / 12; guard-tilt 2404 / 7815 / 31479 / 30113; M 0 | same | yes |
| n=8 c2, cap 16 500, patched vs unpatched control (both mine) | 42 028 algebras, guard-tilt 69 174, rej 1, M 0 | identical to each other | yes |
| n=8 c2, patched, 585 s budget | stopped in depth 8 at 17 058 of 20 899 expansions; 43 638 algebras; guard-tilt 71 659; rej 1 (same parent, path (6,1,1,4,3,1,7,4), vertex 4); M 0 | E-084 file: 20 899 expansions, 54 326 algebras, guard-tilt 89 179, rej 1, M 10 | not comparable (cut earlier) |

Every rejection in the patched runs has `guard` False (the "noguard NOTtilt" row), and the row "guard NOTtilt"
(guard-admitted, failing) is absent in every run: 0 of about 3.2e5 guard-admitted steps in the patched runs above.
(In n = 7..9 and n = 8 c0, c1 the rejection lines in the files are the same parents as before; I compared the
summary lines, not every line.) The patched run costs about 1.3 x the unpatched time (Fractions).

Which E-084 counts change, as asked: gate-admitted parents with `tiltingPlus` False: none changed (n = 8 c2: the one
rejection reappears, with the same path). Guard refusals: none changed in any comparable run. Guard-admitted steps
failing: still 0. n = 8 class 2 loose end: the 10 M lines are not reproduced and not refuted; the walk was too short.
The sizes in the third n = 8 c2 row differ only because of the earlier cut-off (budget vs 20 899), not because of
the patch (the cap-16 500 control shows this).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/017/experimentalist_walk.py 8 --class 0 --stop-on-reject --depth 8 --budget-sec 570   # ~1 min
timeout 10m .venv/bin/python workshop/rounds/017/experimentalist_walk.py 7 --class 0 --stop-on-reject --budget-sec 540            # seconds
MAXEXP=16500 timeout 10m .venv/bin/python workshop/rounds/017/experimentalist_walk.py 8 --class 2 --stop-on-reject --depth 8 [--unpatched]   # ~8-9 min each, two in parallel
timeout 10m .venv/bin/python workshop/rounds/017/experimentalist_walk.py 8 --class 2 --stop-on-reject --depth 8 --budget-sec 585   # 585 s, partial
```
Outputs: `experimentalist_walk_n{7,8,9}_c*.txt`, `experimentalist_walk_n8_c2_cap16500_{patched,unpatched}.txt`,
`experimentalist_walk_n8_c2.txt` (585 s), `experimentalist_walk_n8_c2_run1_budget570.txt` (first, 4 runs in parallel,
16 121 expanded). n = 7, 9 and n = 8 c0, c1 ran 3-4 at a time on 4 cores.

## Prior record

E-084 (counts, 14 classes, M lines), E-085 (defect of `reduceAgainstPivots`; 10 M steps fixed by replay of saved
parents; "E-084's counts not re-run"). This round re-runs them; it does not duplicate E-085's per-parent replay.
Nothing in `research/RETRACTIONS.md` is touched. The n = 7 and n = 9 rows are new evidence for E-085's open
limit "how often the defect fires elsewhere": it did not change any walk table in the 8 classes in which the walk
stopped at a rejection (and only in n = 8 class 2 were M lines ever recorded).

## Code changed

None in the library. New script `workshop/rounds/017/experimentalist_walk.py`. No tests run.

## Next

- Overnight proposal (needs about 15 min, one process, no checkpoint): `experimentalist_walk.py 8 --class 2
  --stop-on-reject --depth 8 --budget-sec 1500` (or the toolsmith's patched library, which should be faster than
  my Fraction monkeypatch) to complete depth 8 (20 899 expansions) and confirm M = 0 and 1 rejection; then
  `--depth 9` without `--stop-on-reject` is the first run that could tell whether the patched walk stays closed
  to the guard (estimate: depth 9 is about 30 000 expansions, hours).
- Toolsmith: a checkpoint/resume for `scholar_walk.py` (pickle the frontier) so the 10-minute rule does not cut
  depth 8; and when the library patch lands, rerun the same rows with `--unpatched` removed to check that
  the library fix equals the monkeypatch.
- Skeptic: compare rejection lines (not only summaries) of patched and E-084 files.
