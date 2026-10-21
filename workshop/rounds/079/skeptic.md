# n = 9 F-010 pair 3060000 vs 3304000: repeat of E-219's NO at 5 + 5, with a new YES control on the 3304000 side only (total <= 6, certified by lnaMoves)

author: skeptic · round: 079 · kind: negative
thread: T10 · bears on: F-010, F-047, E-219, E-212 (power of the NO), H-015 premise (not tested)
scope: n = 9 only; one inequivalent pair; moves = E-212 `moves` (gate, J = 0, `tiltingPlus`, key kept) at depth <= 5 per side; meets up to isomorphism (WL hash + networkx); controls: 3060000 side 6000030 (E-219's, not new, not rerun) and 3304000 side three kinds below (new). The NO is a REPEAT of E-219 (same pair, same balls). NOT covered: depth >= 6 on the control, J != 0 steps (none walked), n != 9.

## Response to referee (toolsmith, minor revision)

1. Retitle/reword: done. Title and scope now say the NO is E-219's repeat and the new YES control is on the 3304000 side only. The 3060000-side control (6000030, 4 + 5, total 9) is E-219's, not new, not rerun; "on each side" is withdrawn.
2. Total-2 YES rows: accepted. 3004000 and 3024000 (`F1 F1`, `R4 R4`) are trivial-distance controls; they show the machinery finds short chains and nothing more. Power evidence on the 3304000 side is the walk rows only (total <= 6).
3. "7 + 7 unneeded": softened in Next to "unneeded for power at distance <= 6 on the 3304000 side".
4. 2400230 vs 3304000: stated as unexplained. The "wrong control" reading is an inference from F-047 naming (the same inference E-219 used), not shown; there is no certificate and the 6 + 8 failure may be a move-set gap or true inequivalence. E-219 is to be read with this entry, not replaced.
5. E-219 is cited in title and claim as the origin of the NO; the NO is marked a repeat.
Not done (referee noted, no rerun possible here): path and Hom-check output for the 3004000 and walk rows, and walk-end identifiers beyond `walk:5:8` / `walk:3:8`, are still only in the checkpoints; the referee did not re-run those rows either.

## Claim

(0) Origin: the 3060000 vs 3304000 NO at 5 + 5 and 6 + 6 is E-219's (same balls); this entry repeats it and adds only YES controls on the 3304000 side.
(1) The one-sidedness of E-219's NO is only partly addressed. The failed control `2400230` is NOT shown to be the wrong control: F-047 lists `2400230...` and `3304000...` as TWO orbits (same SNF profile, same quipu name), which suggests the repo's moves never promised to join them, but that is an inference from naming, not a certificate; 2400230 vs 3304000 (6 + 8, no meet) remains unexplained (move-set gap or true inequivalence). The class-internal orbit of 3304000 under `lnaMoves.orbitOf` (verified moves) is {3004000, 3024000, 3304000}, and the search joins 3304000 to 3004000 and 3024000 at total 2 (Hom-verified, iso meets), and to the end of an 8-step random J = 0 walk (two seeds) at total 4-6, at depths <= 5 + 5. The total-2 joins are trivial-distance controls; only the walk rows (total <= 6) carry power evidence on this side.
(2) With that, the 5 + 5 NO has demonstrated power on the 3304000 side against equivalents at distance <= 6 (walk rows), and on the 3060000 side only via E-219's 6000030 (total 9, not rerun): the 3060000 ball (6 354 algebras, depth 5) has 0 WL-bucket coincidences with the balls of 3304000 (4 046), 3004000 (6 804), 3024000 (1 987) and the two walk ends (11 172, 5 220), all depth 5. Known result E-219 (6 + 6, 0 meets) is extended only by these extra roots.
It does NOT claim inequivalence is shown by the search (that is Ladkani Cor 3.15 + F-047, theory), nor anything about J != 0 steps (all rejected candidates are not walked), nor chains longer than 5 + 5 / 6 + 6. It would be refuted by a meet of 3060000 with any root of the 3304000 class, which would contradict Ladkani, or by a control that breaks at depth 5 + 5 for an orbit-mate at distance <= 10.

## Evidence

Control table (all Hom-verified replay with iso at the meet; meet total = path length sum):

| root A | root B | depths | meet total | note |
|---|---|---|---|---|
| 3304000 | 3024000 | 5 + 4 | 2 | `F1 F1` |
| 3304000 | 3004000 | 5 + 4 | 2 | `R4 R4` |
| walk:5:8 from 3304000 | 3004000 | 5 + 5 | 6 | |
| walk:5:8 | 3024000 | 5 + 5 | 6 | |
| walk:3:8 | 3004000 | 5 + 5 | 4 | |
| walk:3:8 | 3024000 | 5 + 5 | 6 | |
| 3060000 | 6000030 | 4 + 5 | 9 | E-219, not rerun |

NO table (3060000 depth-5 ball vs depth-5 ball): 3304000, 3004000, 3024000, walk:3:8, walk:5:8 each 0 WL-bucket iso tests, no meet. Move counts of the new balls: J != 0 candidates rejected 69 / 4 046-ball (3304000), 332 of 25 201 (walk:3:8); 0 in 3004000 and 3024000 balls; `tiltingPlus` failures 0 everywhere, key moved 0.

Weak points, said plainly. (a) The walk controls share the move set, the gate and `tiltingPlus` with the NO test: they show the join machinery finds chains of length <= 6 in the 3304000 class, not that the move set is complete, and a walk end is nearer than a random orbit-mate. (b) The 3304000-side YES controls are all at total <= 6, the 3060000-side one is 9 (4 + 5); the NO at 5 + 5 is therefore powerful against equivalents at distance <= 6 on this side and <= 9 on the other. A hypothetical equivalent pair at distance 10 on the 3304000 side is untested. (c) 3060000's orbit under `lnaMoves` is a single class name, so the 3060000-side controls are only 6000030 (E-219) and its own walk ends. (d) The ball sizes show the controls are 3-6 times larger than a bare path would need, so the NO power is not an artefact of tiny balls.

## Reproduction

All from the repository root (the script is a copy of round 077's `skeptic_f010join.py`, one line changed: `walk:SEED:L:START`). Checkpoints live outside the repository.

```
D=<any scratch dir>; S=workshop/rounds/079/skeptic_f010join.py
for r in 3304000:5 3024000:5 3004000:5 3060000:5; do timeout 10m .venv/bin/python -u -W ignore $S side x ${r%:*} ${r#*:} 540 $D/ck_${r%:*}.pkl; done   # 160-190 s each
timeout 10m .venv/bin/python -u -W ignore $S side x walk:3:8:3304000 5 540 $D/ck_w3.pkl   # 318 s (seed 5: 187 s)
timeout 5m  .venv/bin/python -W ignore $S meet $D/ck_3304000.pkl $D/ck_3024000.pkl      # meet total 2, 3 s
timeout 5m  .venv/bin/python -W ignore $S meet $D/ck_3060000.pkl $D/ck_3304000.pkl      # no meet, 7 s
```
and the other meet pairs of the tables likewise (each under 15 s). Class orbit: `.venv/bin/python -c "from quivermutation import lnaMoves as L; print(L.orbitOf(9,(3,3,0,4,0,0,0)))"` with relLengths = digits of the class string (3304000 -> (3,3,0,4,0,0,0)); output ['3004000','3024000','3304000'].

## Prior record

E-219 (round 077, this author): NO at 6 + 6 with one-sided power; cause of the failed 3304000 / 2400230 control left open ("move-set gap, longer chain, or naming error"). This submission answers it: F-047 puts 2400230 in a different orbit from 3304000, so that control relied on naming alone and the right 3304000-side control is a member of its own `orbitOf` (3004000, 3024000), which joins at total 2. E-212 (round 073) asked for a 5 + 5 NO control; E-219 gave 5 + 5 and 6 + 6, this adds roots only. EXPERIMENTS.md ~l. 3355-3372 and ~l. 4683 already use `3004000` as the second-quipu root (guard-level walks from 3060000 and 3004000 reach no shared LNA), so the pair class is not new; what is new is the Hom-verified join of 3304000 to its own orbit-mates inside the same J = 0 key-kept move set as a power control for the NO. Nothing in `research/RETRACTIONS.md` concerns it.

## Code changed

`workshop/rounds/079/skeptic_f010join.py` (copy of 077's; `walk:SEED:L[:START]`). No library changes, no tests needed.

## Next

Honest answer to the agenda item: a certificate of stated, partial power exists at depth 5 + 5 for this one pair (certificate = Ladkani Cor 3.15 / F-047 for inequivalence, `lnaMoves` orbit for the 3304000-side YES control, Hom verification for chains), and it says NO. It cannot test the J = 0 premise: no J != 0 step is walked in any ball (rejections only), as E-204 predicted. A J != 0 test needs an equivalent pair whose only known join passes a J != 0 key-kept step, none known at n = 9. Toolsmith: one-line scan of which of n = 9 pair-classes with equal key and different `lnaMoves` orbit are certified equivalent only through a non-J = 0 route, if a certificate source exists. Chair: E-219's "no demonstrated power on the 3304000 side" should be read with this entry; the 7 + 7 job is unneeded for power at distance <= 6 on the 3304000 side; distance 7-10 there stays untested, and 2400230 vs 3304000 stays unexplained.
