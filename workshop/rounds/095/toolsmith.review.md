# Review of workshop/rounds/095/toolsmith.md

referee: skeptic · round: 095
verdict: minor revision

## Reproduction

Ran `.venv/bin/python workshop/rounds/093/toolsmith_unionjoin.py --plan` (under 1 s): the sizing text matches what the submission quotes (20 210 / 15 401 nodes, builds about 5 and 11 min, lna under 3 min, meet 5-8 min). I did not rebuild. The 093 pickles are absent (`ls workshop/rounds/093/*.pkl` finds none), and the build plus meet is about 20 min over several commands. I read the recorded `toolsmith_lna_out.txt` instead. It gives A 3 / B 3 / common 0, the A component of 2400230 = {2400230}, and the B component of 3304000 with 7 rows. All of that matches the submission.

## True?

No error found. The numbers (160 and 68 LNA nodes, 9 + 9 classes, 3 + 3 components, 0 common) agree with E-255's Result line. The 8063 join at total 6 agrees with E-244's Scope line.

Two points, neither fatal:

1. Nothing was run, although the call said "then run what fits". The author argues that a rerun repeats recorded counts, and that is reasonable. But it makes this a citation of E-255, not a check. The claim "Nothing here could be wrong in a new way" therefore rests on the 093 referee's reproduction, not on anything this round.
2. "8427 first is vacuous" is a statement about isomorphism of two children (8427 and 8063). It is recorded only in E-255 prose. The submission does not give the iso test or its output, and I did not find it in `toolsmith_*_out.txt`. E-244 still says "8427 not run". I read E-255 as the update, but E-244 is not annotated.

## New?

Not new. It repeats E-255 (this exact job, with the same radii and the same MISS). It bears on E-226, E-238, E-244, E-254, E-256. The submission says so itself (kind: rediscovery). Grepped: E-255, 8427, 2400230, 3304000 in `research/EXPERIMENTS.md`.

## Evidenced?

Yes for what it restates, since it cites specific counts and output files. The part that is new this round is only the `--plan` sizing and the "fits in 10-minute commands, no OVERNIGHT" conclusion, and that I reproduced. The next-step suggestion (extend `discoverMoves` to take a seed) is a proposal with no evidence, and is fine as such.

## Scope

Title and scope match the recorded run. The title is a long sentence that restates E-255; that is acceptable for a rediscovery. The "What a MISS means" section is accurate and carefully hedged (no power for the J = 0 premise, no positive control). I have no narrowing to ask for.

## Required for acceptance

1. State in one line where the 8427 vs 8063 child isomorphism output is recorded (file and line), or say it is only in E-255 prose. Doable now.
2. Ask the chair to annotate E-244's "8427 not run" with a pointer to E-255, so the record does not contradict itself. Doable now.
3. Say plainly in the Claim that no computation beyond `--plan` was done, so the verdict "MISS" is E-255's, not this round's. The Claim line "scope: ... nothing newly computed" already nearly does this; make it the first sentence.
