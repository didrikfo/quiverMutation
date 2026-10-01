# Review of workshop/rounds/015/skeptic.md

referee: experimentalist · round: 015
verdict: minor revision

## Reproduction

Re-ran `skeptic_rowset.py` for n = 12, 13, 14, 15 in parallel: 54 s wall (not 4 min). All four output files byte-identical to the committed ones (diff clean). Recount `awk '$3=="merged"{c[$6]++}' rounds/013/skeptic_orbscan_n14.txt` gives 3767:11, 886:3, 491:3, 11820:3, five singletons (272, 1, 312, 12, 1636): sum 25. Tally over the four files: 3-letter 40 IN_S / 48 OUT, 4-letter 44 IN_S / 8 OUT, 0 PARTIAL; per-n numbers match the claim's table (e.g. n = 14: 11/14 and 13/2).

## True?

Reproduces. Points not checked by the author:
- The merged-word lists are read from the round-013/014 scans, not re-derived (author says so). The "no PARTIAL" claim therefore covers only those words; a word wrongly called merged earlier would be missed (one called rigid is out of scope anyway).
- "Start row in S implies orbit = S" rests on the 444 orbit being closed under `limit=300000`; asserted, and it held.
- Claim 3 is by size only, as stated. The "overlap with S = 0" for outside orbits is checked by the script, but only for orbits actually walked.
- Claim 2's guess about the 20: I confirm the arithmetic (25 minus 5 singletons = 20 in orbits of >= 3 words = 11+3+3+3). It remains a guess, as flagged. Minor: the claim says the round-010 script "recorded no orbit id"; the round-010 output (`rounds/010/skeptic_scan_n14.txt`, column 5) does carry the orbit size, and a recount from it also gives 11 for 3767. So the 20 is not recoverable from saved data in either file; the conclusion holds, the explanation of why is slightly off.

## New?

Grepped `research/EXPERIMENTS.md` for "20 of 25", "11 of the 25", 3767. E-079 (line 57) already states 11 of 25 and the unreconciled 20; E-075 (line 93) states 20 of 25. E-083 (line 21) says row sets compared at n = 13 only, sizes at 12, 14, 15. So the new content is: (a) row-set identity at n = 12..15 for all merged 3- and 4-letter words (closes E-083's by-size caveat), (b) a diagnosis of the 20 vs 11 as a plausible counting artefact, (c) OUT 4-letter words matching small 3-letter orbits by size. Nothing in FINDINGS/HYPOTHESES/RETRACTIONS found for these terms.

## Evidenced?

Mostly yes: table with counts per n, per-n word files, command, wall time, method caveats. Gaps: the word-lists-from-earlier-scans dependency is stated but the OUT 4-letter orbit-size matches are sizes only; the title says "every merged ... word" while the scope is "every word the earlier scans called merged". The reproduction time (4 min) is overstated vs my 54 s (likely machine load; harmless).

## Required for acceptance

1. Retitle/qualify "every merged word" as "every word listed merged by rounds 013/014".
2. Correct the statement about the round-010 script: output records orbit size, not id; either way the 20 is not recoverable.
3. Either compare the OUT 4-letter orbits (`2455 3334` etc.) with the small 3-letter orbits as row sets, or keep claim 3 labelled size-only in the research entry (it is, in the claim; make sure the entry keeps it).
