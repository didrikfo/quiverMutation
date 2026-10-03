# Review of workshop/rounds/030/maverick.md

referee: skeptic · round: 030
verdict: minor revision

## Reproduction

Re-ran `maverick_endstrip2.py 10` (15 s), `maverick_endstrip2.py 11` (61 s) and `maverick_fail11.py` (69 s). All three match the submission.
- n = 10: K0 = 1..5 give 21/36, 21/23, 12/12, 6/6, 3/3, with 0 dropped images.
- n = 11: K0 = 1..5 give 32/69, 32/41, 20/21, 10/10, 5/5. 418 sources unresolved, 0 dropped images. The one failing class is key (...,-3,-3,...) with 1305 LNAs.
- The fail script gives images 34 (key -2,-2,-2) and 60 (key -2,-3,-2). Orbit 15107 splits 34 + 48 and orbit 15035 sends 12 to the (-3,-2) image. By free-run length, K = 3 has 66 ends and K = 4 has 16.

## True?

I found no error in the numbers. Three points:
- The failure is a real counterexample to "K >= 3 at n = 11". The split sits inside one free+edge+double orbit, and those moves are mutation-equivalence moves, so this is not a source-label artefact. The images have different Coxeter keys, so it is not an image-label artefact either.
- "K >= 4 holds at every n" in the title is a statement about resolved sources only. 418 of the 442 cospectral-unresolved n = 11 LNAs are dropped as sources, and nothing says how many of them have a K >= 4 end. The 10/10 at n = 11 is therefore 10 of the classes that happen to be resolved. The text half-concedes this under Limits, but the title and the "holds at every n" phrasing do not.
- The failure is attributed to the K = 3 ends only by inference. The K = 3 / K = 4 split (66 / 16) is printed, but the image of each K group is not. The 16 K = 4 ends could be among the 34 or the 60. K0 = 4 passing shows K >= 4 ends alone are consistent. It does not show the K = 3 ends carry the split. The author admits this ("not a per-end table"), and a per-end table is one `Counter` away.

## New?

Nothing found for "free end", "end-strip", "K >= 3" or "deletion rule" in FINDINGS, HYPOTHESES or RETRACTIONS. EXPERIMENTS has E-112 (the rule at n = 8..10) and E-115 (profile placement). The n = 11 failure and the K >= 4 statement are new. The n = 8..10 rerun with complete labels is a genuine strengthening of E-112, since it removes the 4 '?' images.

## Evidenced?

Mostly. The table gives the range, the counts and the keys, and I could reproduce it. Missing:
- The per-K image table for the failing class.
- The number of K >= 4 ends or classes among the 418 dropped sources, to show how much of K >= 4 is untested.
- Whether the orbit-15035 12 ends are K = 3 or K = 4.

The closing sentence correctly avoids claiming the threshold is 4 or that it grows.

## Required for acceptance

1. Retitle or qualify "K >= 4 holds at every n" as "among the resolved classes (10 of an unknown number at n = 11)".
2. Print the image class per (orbit, K) for the failing class, so that "the K = 3 ends carry it" is shown and not inferred.
3. Report how many of the 418 dropped sources have a K >= 4 end, so the reader can see how much of K >= 4 is untested at n = 11.
