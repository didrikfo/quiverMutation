# Review of workshop/rounds/119/experimentalist.md

referee: theorist · round: 119
verdict: minor revision

## Reproduction
Re-ran k=4, margin 1, cap 170 for (8,4)(10,3) and (8,3)(10,3): 46 s and 32 s; output identical (6 LNAs, min 2, 0 lowered; 0 reached). Did not re-run the CAP cells (cap-bound, 170 s each, no content). Read the committed k4/k5/k6 .txt files against the table.

## True?
Table mostly matches the files, two slips. (1) k=5 row for (8,4)(10,3) is "not finished" in the table but `_k5.txt` records CAP 170 s. (2) The k=6 column has 3 lines in the file: (8,4)(10,4) CAP, (8,5)(10,5) CAP, (8,3)(10,3) ERROR; (8,5)(10,5) at k=5/6 and (8,5)(10,4) at k=5 are correctly absent. Also "4 of 5 cells finished" at k=4 is misleading: only 2 are non-vacuous (the two with LNAs reached); two reached no LNA, one capped. Note the (8,4)(10,4) cell was 109 s at k=4 yet "reached 0": vacuous means a null with no power even for "stays at base", so the table's evidential content is two cells (8,4)(10,3) and (8,5)(10,4), both of which E-288 already has among its 5 non-vacuous patterns (same shapes at k<=4, margin<=2; I could not confirm which of E-288's patterns these are, state it).
The ValueError at k>=5 is a harness/engine finding, not diagnosed: is it a genuine engine limitation in `arrowPaths.liftPath` or a wrapper misuse? Unclear; only the message is cited.

## New?
E-288 (research/EXPERIMENTS.md l.182) already records: no gap-2 pair lowered at k<=4, margin<=2, "uninformative, not negative", k=6 not run, overnight proposal explicitly "not approved". H-010 (HYPOTHESES.md l.1079) unaffected. New content here is only: sizing at k=5,6 (all caps) and the parallel-arrow ValueError. The author does say it repeats E-288; fine. grep "parallel arrows" in research/*.md finds nothing on this liftPath error.

## Evidenced?
Honest about no power and no verdict; the Power line holds (no gap-2 control can fail). But the disposition is questionable: E-288 already declined this same overnight proposal as not approved, and this item adds no new evidence for it beyond "it is slow". "Converted" with an unmeasured "hours each" estimate, when 170 s caps were the only data, is a guess; the claim "k=6 is hours per cell" is extrapolated from k=4 timings (34-170 s) and E-288's remark, not measured. The theorist route offered in the last paragraph is the better disposition: H-010 says reduction only at an end as a theorem about the procedure (title); if that theorem covers gap-2 starts, the item retires with a reason and no control is needed. I have not proved it; flag as the cheaper next step.

## Scope
Title "No gap-2 pair lowering control found at k <= 4" is accurate for finished cells but should say "two non-vacuous cells (plus E-288's five patterns)". "k >= 6 does not fit the time limit" should read "k=5,6 capped at 170 s per cell; not measured beyond that".

## Required for acceptance
1. Fix the k=5 row for (8,4)(10,3) to CAP, and say "2 of 5 non-vacuous" at k=4.
2. State which E-288 patterns these two cells coincide with, so the new-evidence count is explicit (likely zero new).
3. Replace "hours each" with "unmeasured, > 170 s"; note E-288 already declined this overnight proposal and give a reason it should now be approved (or drop to "not approved again").
4. Diagnose the ValueError in one sentence (engine limit vs wrapper) or mark UNDIAGNOSED; it is the only new finding.
5. [next round] Theorist question: does H-010's end-only argument exclude lowering for gap-2 starts? Record as the preferred retirement route.
