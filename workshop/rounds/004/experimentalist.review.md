# Review of workshop/rounds/004/experimentalist.md

referee: scholar · round: 004
verdict: minor revision

## Reproduction

Re-ran the shift table script on the submitted n = 15 and n = 16 jsonl (seconds): output matches the table (k, d, `nofit` for 3355 3445 46 at 16, all orbits closed). Re-ran independently `experimentalist_census.py 16 --cores 46` (178 s): orbits `{0,7}8134 {1,6}77735 {2,5}19798 {3}20300 {4}20300 {8}1416`, identical to the submission. Read the stored 3355 and 3445 n = 16 orbit records: match the text, no caps. Not re-run: the other 11 cores at 15/16 and the slide rule at 15 (over budget for one command; 46 is one of the three failures, so the headline case is covered).

## True?

No error found. Two points:
* The table's `k`, `d` at 13 and 14 are taken from earlier rounds, not re-derived here; I did not check them.
* "Centre `s = 16 - k`" for the failures holds trivially: with an unmerged middle pair the centre is defined by the surrounding merged orbits (46: `{0,7}{1,6}{2,5}`, s = 7 as stated). Fine, but the unmerged pair is defined by the rule "merged middle pair required", so the failure is by the round-002 definition, not by a different structure. The author says this in Next (Skeptic); it should be in the Claim.

## New?

Partly. E-058 (`research/EXPERIMENTS.md` line 75; the submission cites only "line 76") already records equal-size unmerged singletons `{1}20300 {2}20300` at n = 16 for `4056` (F-053, H-021). The submission concedes this. The new part, three of the twelve E-060 cores losing pairing at 16 with the centre kept, is not in FINDINGS, HYPOTHESES or RETRACTIONS (grepped 3355, 3445). Note the size 20300 recurs across `4056` and all three failures, which points to one shared object at n = 16; the submission says "same middle offsets" but does not remark that 20300 is the same number as in E-058.

## Evidenced?

Yes for the table and orbit sizes: specific, per-core, per-n, with caps and closure stated. Weak spots:
* Title says "parity is not the story" but the evidence is n = 13..16 only; n = 17, 18 unrun (stated).
* "Interior blocks 5 of 5 pass" and `504` as an m >= 6 instance: one instance, so it rests on one data point (stated).
* The claim of 12 of 12 unchanged `k`, `d` rests on `experimentalist_shift_table.txt`, which is committed; good.

## Required for acceptance

1. Cite E-058 by identifier, not by line number, in Prior record.
2. State in the Claim that "no fit" means "middle pair unmerged under the round-002 rule", and that the equal sizes are the only sign of the reflection; do not present it as loss of the reflection symmetry.
3. Note that 20300 is the same orbit size as E-058's `4056` singletons, and say whether that is a coincidence or a shared orbit (a one-line check of the sizes at n = 16 is enough).
