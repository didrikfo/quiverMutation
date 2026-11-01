# Review of workshop/rounds/106/scholar.md

referee: experimentalist · round: 106
verdict: minor revision

## Reproduction

`.venv/bin/python -I workshop/rounds/106/scholar_cor46.py`, about 3 s. Output matched the submission: printed variant (2253, 3676 inadmissible, 0 iso failures, R3 in 2253 and R4 in 1423); fixed variant (2253, 0, 0, {}); A_n ~ B_n quipu True for n = 4..12. I re-fetched the 2305.06642 e-print into an empty scratchpad directory. Lines 116-117 (and 549-550) of `Derived_equivalent_Nakayama_algebras.tex` show the right column of rows 1-2 with `n_0+l_0`, and the left column with `n_0+2`. I did not re-derive which rows 3 and 4 the submission means, and I did not re-run the numbering script. I did not fetch the other two sources.

## True?

Nothing found against it. I did not test the hand check (r = 0, n_0 = n-3, l_0 = 3 gives n, but the relation should end at n+1) beyond the script. One gap: "A_n ~ B_n for all n >= 5" is derived from Thm 4.2 plus Lemma 4.4, which cites external theorems that were not re-read. The submission says so. Quipu isomorphism was checked only for n = 4..12, and the claim is "every n >= 5" by argument. Say "n = 5..12 checked, general n by Thm 4.2/4.5".

The submission says E-218's title contains "equivalent only at n = 9, 10". This is accurate: the E-218 header in EXPERIMENTS.md contains that phrase, and so does its Scope ("known only at n = 9 ... and n = 10").

## New?

Mostly already recorded. E-275 has the misprint. E-280 has the theorem numbers and the "superseded" consequence line, which already names E-218 and E-247. The submission acknowledges this. New: the independent Thm 4.2 check of the misprint, and the R-entry proposal. RETRACTIONS.md has R-013 as its latest and no R-014, so the proposal is feasible. Because E-280 already carries the supersession line, the R-entry is a judgment call for the chair. I lean toward adding the one-line pointer on E-218's title.

## Evidenced?

Yes. Counters, script names, sample sizes, seeds and the planted-number power test are all stated. The Power line holds: my run reproduced both the NO (printed variant) and the pass (fixed variant). The 2253-case sample covers r = 0..3 and n = 8..22 only, with l_i = 3..6; the submission says this.

## Scope

The title lists many items; the checked scope is the LaTeX of v1 and the numbers only. The title is accurate on this. Narrow "A_n ~ B_n for all n >= 5" to "follows from Thm 4.2 + Lemma 4.4, quipu-checked for n = 4..12".

## Required for acceptance

1. Reword the A_n ~ B_n sentence as above (n = 4..12 computed, general n by the cited theorems).
2. Chair decision on R-014 versus a pointer on E-218's title, noting that E-280 already carries the supersession line. Not an author item.
