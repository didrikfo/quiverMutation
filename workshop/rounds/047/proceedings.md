# Round 047 -- proceedings
kind: ordinary. Step 0.5: round-046 questions unanswered; decided (below). Next round (048) is a conference.

## experimentalist (T10 ii) -- accept (response did the required items)
Claim: at n = 8 class 2, E-096's depth-8 guarded walk has 0 of 104 629 key-kept edges failing `tiltingPlus`; the 2 failures are key-refused. Referee (skeptic): minor revision; slice 1 reproduced, slice 2 not. Response: the tree + starts vs distinct off-by-2 explained (2 keyless children), taint bookkeeping demoted, resume stated, claim narrowed. Promoted: **E-155**. Depth 9 and the n = 7 c1, c2 control are not done.

## toolsmith (T10 iii) -- accept, narrowed
Claim: the 10 F-041 n = 8 merges have shortest witness 3 + 3 = 6 inside the key-guarded gate graph; 60 of 60 edges pass; `merges.py --witness` stores paths. Referee (theorist): minor revision; reproduced (19 s, 27 tests). Response: title and scope narrowed; inverse B-half tested by an op-forward definition (30/30), not an independent E-149 reverse step; `merges.py 10 --depths 5 --witness` untested by data. Promoted: **E-156** (adds little beyond E-150 / F-041; the new part is the tool and the per-edge check along stored paths). Code (`merges.py --witness`, new test) committed with the round.

## maverick (T1/T2 breadth slot) -- note
Claim: T1/T2 have no open question with a cheap census; settled and open items listed; reopening test = a rule for k(c) that predicts a held-out core. Referee (scholar): minor revision (omitted E-076..E-100 items, mis-cited fits, over-read E-072); response did all six items. Nothing run, nothing promoted. Decision: T1/T2 stay dormant; **the H-021 header change** (iff REFUTED by E-058; H-021' SUPPORTED as description) is left to the round-048 ledger. T3/T8 is not closed with them.

## Consequences
None new. E-155 and E-156 point the E-150 way (no recorded merge examined uses a J != 0 step); H-015 stays OPEN (E-147/E-151 stand at n = 7). The `mutationSearchDepthFirst` docstring is not reworded; T10 (i) is untouched.

## Questions for the steering committee
1. **Overnight: depth-9 continuation of the n = 8 c2 walk** (experimentalist's proposal, ~4-5 slices, checkpoint >100 MB)? Recommend: no; first the cheap control, the same edge tally on n = 7 classes 1, 2 where E-151's failures lie.
2. **PDFs of 1009.3370 and 2509.12983** (arxiv.org blocked). Recommend: yes if you can supply them.

## Decisions taken for the steering committee
- round 046 q1: `tiltingPlus` not promoted; q2: literature parked; q3: no overnight, deep replay sized first (all in STEERING).
