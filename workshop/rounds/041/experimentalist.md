# With the key guard off, no gate-admitted J_i != 0 step at n = 6, 7, 8 (classes 0-2, time-capped) has a child with the class key; a random-algebra control shows the test can say yes

author: experimentalist · round: 041 · kind: negative
thread: T5 · bears on: E-137, E-138, H-015

## Claim
In guarded walks from the seeds of LNA key classes 0, 1, 2 at n = 6, 7, 8, every gate-admitted step with J_i != 0 that was reached (table below, 282 steps incl. repeats across the 150 s / 480 s caps) has a child whose Coxeter key, computed on the reduced child, differs from the class key (0 keep it). In fact every one of them lands on a different LNA key of the same n: (1,1,0,0,0,1,1) in the n = 6 examples printed. This extends E-138 (n = 8 c0, n = 6/7 c0) to c1, c2 and n = 7 c1, but is NOT a proof: all walks are time-capped and 4 of 9 (n, class) cells reached no J != 0 step at all (vacuous, not a verdict). A positive control (below) shows the "child key == parent key" test returns yes on J != 0 steps for non-LNA parents; I did NOT find a J != 0 child that keeps an LNA key.

## Evidence
Child key = `search._coxeterKeyOrNone(reducePathAlgebra(mutation))`, vertex set checked (0 dropped, 0 illegal children in all runs). Tested, never filtered; expansion follows only key-preserving children (guard on for the walk, off for the test).

| n | class | seeds | expanded | gate steps | J!=0 steps | child==class key |
|---|---|---|---|---|---|---|
| 6 | 0 | 2 | 7361 | 27387 | 147 | 0 |
| 6 | 1 | 24 | 8911 | 29702 | 0 | - |
| 6 | 2 | 26 | 9058 | 29402 | 0 | - |
| 7 | 0 | 8 | 4530 | 17177 | 49 | 0 |
| 7 | 1 | 12 | 3659 | 13637 | 2 | 0 |
| 7 | 2 | 14 | 4164 | 15919 | 0 | - |
| 8 | 0 | 2 | 4818 | 20347 | 59 | 0 |
| 8 | 1 | 8 | 6837 | 29186 | 21 | 0 |
| 8 | 2 | 18 | 7955 | 32352 | 0 | - |

Caps: 150 s (n = 6, 7), 480 s (n = 8) under load 4 to 8; n = 8 c0 expanded 4818 against 7852 in E-138, so fewer J != 0 rows (59 against 192): this is a smaller sample, not a disagreement. c0 is the only class with J != 0 steps early; c2 had none in 8 000 expansions.

Guard OFF expansion (children failing the key are followed too, depth-limited, `off`): n = 6 c0 depth 10, 11253 expansions: 316 J != 0 steps, 0 keep the class key; c1 depth 6 (7026), c3 depth 8 (14976), n = 5 c0, c1 depth 12: 0 J != 0 steps; n = 7 c0 depth 6: 26 steps, c1: 4 steps, 0 keep. In every case the parent of a J != 0 step had the class key: off-class algebras reached within these depths never had a gate-admitted J != 0 step.

Positive control (random acyclic algebras on 6 vertices, random arrows, commutativity and zero relations; `experimentalist_control.py 6 1 60`, 60 s, 28813 algebras): 1264 gate-admitted J != 0 steps, 55 have child key == parent key (4 further children illegal). So the comparison, the J detector and the child-key code can answer yes. Limits of the control: parents are not LNAs and their keys are mostly non-LNA keys, some relation sets are non-homogeneous; I did not check these 55 for being tilting (tiltingPlus is False by J != 0) nor for derived equivalence of the child, so they show the test is not constant-no, not that the child is derived equivalent. The question for an LNA parent stays: none found.

## Reproduction
```
timeout 10m .venv/bin/python -u workshop/rounds/041/experimentalist_keyoff.py 8 1 480 99     # 480 s; args n class budget maxdepth [off]
timeout 10m .venv/bin/python -u workshop/rounds/041/experimentalist_keyoff.py 6 0 400 10 off   # 157 s
timeout 5m .venv/bin/python workshop/rounds/041/experimentalist_control.py 6 1 60             # 60 s
```
Outputs: `workshop/rounds/041/experimentalist_keyoff_n{6,7,8}c{0,1,2}.txt`, `..._keyoff_off_*.txt`.

## Prior record
E-138 (n = 8 c0 all 192, n = 6, 7 c0 123 fail), E-137. New: classes 1, 2, n = 7 c1, the off-mode expansion, the random-algebra control. Not new: the verdict for c0.

## Code changed
None in the library. New scripts `experimentalist_keyoff.py` (reuses rounds/033 `perI`), `experimentalist_control.py`. No tests touched.

## Next
- Overnight proposal: n = 8 c0, c1 to closure or 7 000+ expansions with `experimentalist_keyoff.py` (c2 shows no J != 0 steps in 8 000 expansions; is that a class property?).
- theorist: why every child lands on another LNA key (1,1,0,0,0,1,1) at n = 6; is the key change a shift between classes, a map between classes? Tally child key per class.
- skeptic: check the 55 control hits (relations legal? child derived equivalent?) and try to build an LNA-parent analogue.
