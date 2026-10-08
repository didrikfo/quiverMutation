# Review of workshop/rounds/041/experimentalist.md

referee: skeptic · round: 041
verdict: minor revision

## Reproduction

- `experimentalist_keyoff.py 7 1 150 99` (150 s cap, 2.5 min): same as the table. 12 seeds, 3659 expanded, 13637 gate steps, 2 J != 0 steps, 0 vertex-set drops, 0 illegal. Both steps have a parent with the class key and a child that is not the class key. The child key is (1,1,-1,-3,-3,-1,1,1), an LNA key. The expansion count matches exactly, so the time cap is not what decides the walk.
- `experimentalist_control.py 6 1 60` (61 s): 1175 + 53 J != 0 steps with child key == parent key (the writeup has 1264 and 55; the run is time-capped, so the counts differ). The control does say yes, so the test is not constant-no.
- I did not re-run the 480 s n = 8 walks or the off-mode runs.

## True?

I found nothing false. The claim is stated as negative and non-exhaustive.

- The "282 steps" figure does not match the table. 147 + 49 + 2 + 59 + 21 = 278. Fix it, or say what the extra 4 are.
- "In fact every one lands on a different LNA key (1,1,0,0,0,1,1) in the n = 6 examples printed" is a stated example at n = 6, and the claim sentence runs it together with all n. For n = 7 c1 the key is the 8-term one above. The "lands on another LNA key" statement is therefore supported only by the printed examples. The script's `child is an LNA key at n` column is the evidence, so give its tally per cell.
- The control hits are not LNA parents, and the author says so. They are no evidence about LNA parents, in either direction.

## New?

E-138 covers n = 8 c0 and n = 6, 7 c0, so the verdict for c0 is not new, as the author says. I grepped `research/EXPERIMENTS.md` and `HYPOTHESES.md` for the key-guard and child-key terms, and for c1 and c2 at n = 6, 7, 8, and found nothing. The new content is:
- classes 1 and 2;
- n = 7 c1;
- the off-mode expansion;
- the control.

It is a small increment.

## Evidenced?

The method, caps, counts and output files are given. Weak points:
- 4 of 9 cells have 0 J != 0 steps, so they say nothing. The "0 keep" total rests on 5 cells, and on two of those the sample is 2 and 21 steps.
- The off-mode paragraph reports ranges as prose ("depth 6, 26 steps"), not as a table, and it has no expansion counts for n = 7.
- The off-mode sentence "off-class algebras ... never had a gate-admitted J != 0 step" is only as strong as those depths.

## Required for acceptance

1. Reconcile 282 with the table (278).
2. Give the per-cell tally of child key (and whether it is an LNA key) so the "different LNA key" statement covers n = 7 and n = 8 as well as n = 6.
3. Say in the title and claim that the sample is 5 non-vacuous cells, not 9.
4. Put the off-mode runs into a table with expansion counts.
