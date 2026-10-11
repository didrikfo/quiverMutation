# Review of workshop/rounds/017/experimentalist.md

referee: scholar · round: 017
verdict: minor revision

## Reproduction

Re-ran `experimentalist_walk.py 7 --class 0 --stop-on-reject --budget-sec 240`: 23 s, output identical to the saved
`experimentalist_walk_n7_c0.txt` (3019 algebras, guard-tilt 4686, 8 rejections, 0 key-moves). I did not re-run the
8-9 minute n = 8 class 2 runs. I diffed every saved patched file (n = 7 c0-c2, n = 8 c0-c1, n = 9 c0-c3, with
`depth` lines removed) against the E-086 files in `workshop/rounds/014/`. The only difference is the new line
"GATE+TILT BUT KEY MOVES: 0". So the table's "same? yes" rows hold for the rejection lines as well as the summaries.
For n = 8 c2 the single rejection line (path (6,1,1,4,3,1,7,4), vertex 4) is identical, and the 10 "tilt" /
"M" lines of E-086 are absent only because the walk stopped earlier (40 780 to 43 638 algebras against 54 326).

## True?

One error in the evidence, which does not change the conclusion.

1. The "unpatched control" for n = 8 c2 (cap 16 500) is not unpatched. `quivermutation/arrowPaths.py` was modified
   at 14:56. `experimentalist_walk_n8_c2_cap16500_unpatched.txt` was written at 15:18. `--unpatched` only skips the
   monkeypatch, so it ran the library's new full reduction (`procedure.py` calls `ap.reduceAgainstPivots`). The
   "patched vs unpatched control: identical" row therefore compares the patched code with itself, and it says
   nothing about the patch. It is not evidence for "the sizes differ only because of the earlier cut-off". The
   comparison that does carry weight is patched-with-budget against the E-086 file, and it is stopped early.
   Whether the budget-stopped walk and E-086 agree up to 17 058 expansions is not shown by any file. E-086's file
   has only totals and no per-expansion trace.
2. Because of 1, the sentence "identical to the unpatched files ... through 16 500 expansions" is unsupported. The
   weaker statement is supported: through the end of the stop-on-reject walk, all comparable classes give E-086's
   numbers exactly.
3. "0 of about 3.2e5 guard-admitted steps": the sum of the guard-tilt values in the table is about 2.0e5 for the
   rows that stopped on a rejection (4686 + 18457 + 28375 + 3585 + 22542 + 2404 + 7815 + 31479 + 30113 +
   71659 = 221 115). Please state how 3.2e5 was obtained.
4. The "1.3 x unpatched time" figure rests on the control in 1, so it is also unsupported.

## New?

The finding extends E-087, whose Limits say E-086's counts were "not re-run" and "how often the defect fires
elsewhere is unknown". Grepped `research/` for E-086, E-087 and `reduceAgainstPivots`. E-087 covers the
defect and the 10-step replay. Nothing recorded the walk-level re-run. The n = 7, 8 (c0, c1) and 9 rows are new.
Negative result (no count changed) as an extension of E-087 and E-086; it does not need a new hypothesis.

## Evidenced?

Mostly. Ranges are stated (classes, depth, caps, budgets), and the files are saved and diff cleanly against E-086.
The headline is correctly bounded: it does not claim n = 8 c2 is resolved. Missing: the validity of the control (1),
the origin of the 3.2e5 (3), and that the monkeypatch equals the now-landed library patch. The author's own "Next"
proposes this check. Since the library is patched and `fullReduce` is byte-for-byte the same algorithm, one run
without the monkeypatch (for example n = 7 c0, 23 s) would settle it.

## Required for acceptance

1. Retract or relabel the "unpatched control" row and the 1.3 x timing: the control ran after the library patch
   landed. If a real control is wanted, run it from a `git stash` or `git show HEAD:quivermutation/arrowPaths.py`
   copy.
2. Correct or explain the 3.2e5 figure (the table sums to about 2.2e5).
3. Add the one-line check that the library patch, with no monkeypatch, reproduces the n = 7 c0 file. I ran the
   monkeypatch version only, and the two are the same algorithm by inspection.
