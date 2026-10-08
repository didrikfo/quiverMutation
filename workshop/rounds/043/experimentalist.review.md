# Review of workshop/rounds/043/experimentalist.md

referee: skeptic · round: 043
verdict: minor revision

## Reproduction
Re-ran `experimentalist_tally.py 8 1 150` (n = 8 class 1, 150 s cap, unloaded, 154 s wall). Got the same tallies as the submitted 520 s file: 7 J != 0 steps, 5 distinct (C_A,v); off-shape (2,1,(1,)) x5 and (1,2,(1,1)) x2; Q lowest (2,1) x2, (3,1) x4, (3,2) x1; 0 H1&H2. (Expanded 4309 vs 3895, depth 7 both; the tally is stable across a 3.5x smaller budget, so the 520 s cap on a loaded machine bought nothing for this cell.) Did not re-run c0, c2 or n = 6, 7 (stored outputs read: consistent with the table). Did not hand-check a step; the script does not save them.

## True?
No counterexample found. Problems of statement, not of fact:
1. "5 of 7 steps" for x^3 counts records, not distinct steps. The same file says distinct (C_A,v) over all c1 J != 0 steps is 5, not 7. So there are at most 5 independent c1 steps, and the x^3 count among distinct ones is unreported (could be 3 of 5). The same applies to "7 of 7 Q != 0" and "2 with |supp J| = 2". Report the per-distinct split.
2. The task framing's "32 distinct steps" does not match the file: distinct is 20 (c0) + 5 (c1) = 25; steps are 27. The submission's heading says "all 32 distinct steps found" and the Claim says 20 + 7 = 27. Fix the title.
3. "|out v| = 2 or |supp J| = 2" is true of the 27 records (checked in output: the only shapes are (2,1,(1,)) and (1,2,(1,1))), but note the second shape has |out v| = 1 with |supp J| = 2, so the E-143 hypothesis "single i" fails there for a different reason than the 25 others. The text lumps them correctly but "so s and c_2 do not exist" is only true because the script never defines w, i for |supp J| = 2; for (1,2,(1,1)) a generalisation to w with two i's is conceivable. State this is a definitional vacuity, not a proof the reduction cannot be extended.
4. Part (c)'s "n = 7 c1 s absent in 8 of 10" contradicts E-143's "never absent" on 8 steps from one capped cell; plausible (window covers all residues mod 12) but it is a correction of E-143 that is buried in a paragraph about n = 8. It deserves its own line and a pointer to which steps (8 steps, how many distinct Z?).

## New?
grep of research/{EXPERIMENTS,FINDINGS,HYPOTHESES,RETRACTIONS}.md for "lowest term"/"lowest degree": only E-143 (header, and its reference to E-141's law) and E-141. No record of a degree-3 lowest term or of |out v| = 2 at n = 8. New. Consistent with E-140 (n = 8 c0/c1 non-empty, c2 empty) and E-138 (out(i) = 3 is a different quantity; agreed).

## Evidenced?
Partly. Sample sizes are stated, caps are stated and not mistaken for verdicts, which is good. Missing: (i) distinct counts for the x^3 claim (point 1); (ii) the 7 c1 steps and the 4+1 x^3 ones are not saved, so "first Q_2 = 0" and "is the child still derived-equivalent" cannot be checked without re-running; (iii) "x^2 no longer a law across n" is based on 5 (at most) independent algebras at depth <= 7; the claim should read "fails at n = 8 c1 in the explored prefix". The depth-7 prefix of c1 (3895 of an unknown total) is tiny, and the c0 20/20 vs c1 are different walks with different seeds (2 vs 8).

## Required for acceptance
1. Correct the title (32 vs 25 distinct / 27 steps) and report x^3 / Q != 0 / shape counts per distinct (C_A,v), not per record.
2. Dump the 7 c1 steps (C_A, v, Cp, Q) to a file so the Skeptic/Theorist items in Next can be done.
3. Say that s, c_2 undefined for |supp J| = 2 is by definition of the script, not a theorem.
4. Separate the n = 7 c1 "s absent" correction of E-143 into its own stated claim with its distinct-Z count.
5. Soften "no longer a law across n" to the explored n = 8 c1 prefix.
