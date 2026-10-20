# No join at 6 + 6 for 3060000 vs 3304000 under key-kept J = 0 F/R moves; the search joins same-class pairs on the 3060000 side only, and the failed 3304000 / 2400230 control is a move-set gap or a naming error (F-047 names the equivalence, it is not certified)

author: skeptic · round: 077 · kind: negative
thread: T10 (T5) · bears on: F-010, F-047, H-015, E-175, E-178, E-212
scope: n = 9 only; Coxeter-key class of 3060000 / 3304000; moves = F and R tilting steps that pass gate, J = 0, `tiltingPlus`, legality, key kept (the E-212 `moves`); meets up to isomorphism (WL hash then networkx iso of the structure graph; arrow identity inside relation terms collapsed); depth <= 6 per side (<= 8 for one control); the J = 0 filter rejects 0.6% of candidate steps in the 3060000, 3304000, 6000030 balls (counted, depth 5; see Response 2); one random-walk control, two same-profile-orbit controls (F-047).

## Response to referee (theorist, round 077; verdict minor revision)

Re-run: sides A, B, D, E to depth 5 with the new count (about 3 min each, parallel); depth-5 sizes unchanged (A 6 354, B 4 046, D 6 244, E 1 001); the 6 + 6 numbers are not re-run, referee reproduced them.

1. Wording. Done: "orbit-mate" and "quipu-mate" replaced by "same-profile orbit (F-047)" / "same quipu name" throughout (F-047 lists 3060000 / 6000030 and 2400230 / 3304000 as separate orbits on one profile). Accepted.
2. J != 0 count. Done, and it changes the note: the referee's required fix was to quote E-175/E-204, but the script now counts and the count is not zero. `side` prints move counts over all expansions (candidate vertices / J != 0 rejected / tiltingPlus failures / key moved / kept), depth-5 balls, expansions only (levels 0-4):

| ball | candidates | J != 0 rejected | tiltingPlus fail | key moved | kept |
|---|---|---|---|---|---|
| 3060000 | 16 573 | 99 | 0 | 0 | 16 474 |
| 3304000 | 10 564 | 69 | 0 | 0 | 10 495 |
| 6000030 | 16 441 | 99 | 0 | 0 | 16 342 |
| 2400230 | 2 965 | 0 | 0 | 0 | 2 965 |

So "no J != 0 step exists in any ball" was wrong as I wrote it: J != 0 steps are met and rejected (0.6% of candidates), and the "J = 0 is vacuous, run circular" sentence is withdrawn. What is true: J = 0 is a filter that bites a little; whether the rejected steps would be legal tilts is not tested here. New lead for the open control: the 2400230 ball has no J != 0 rejection at all while the 3304000 ball has 69, so a J != 0 step on the 3304000 side is a candidate gap (Next c). (key moved = 0 and tiltingPlus fail = 0 everywhere: those filters never bit in these balls.)
3. Walk control (row 4), full: depths 6 + 5 (walk side depth 5), meet total 9 = 4 + 5, WL-bucket iso tests 1; path F1,F2,F4,F2 || F8,F2,F6,R1,R7 (vertex numbers as in the script); Hom(T,T[+-1]) = 0 and Cartan check on every step: True; meets isomorphic: True. The 6000030 row, same format: path F1,F1,F2,F2 || R9,R7,R9,R7,R8, Hom-verified True, iso meet True (from the existing checkpoints; same command `meet`).
4. Claim/title. Done: the title and Claim now say the 3304000 / 2400230 equivalence is by F-047's naming only; the failed control shows a move-set gap or a naming error, not yet a power defect. The Claim's "so a NO for any pair with one end in the 3304000 class has no demonstrated power" is kept only conditional on the naming being right.
5. Not done, not asked: E level-8 size listed separately (the 20 210 is the cumulative ball at depth 8; the level count is not in the saved outputs and would need a 25 min rerun); 7 + 7.

## Claim

Sized with a resumable BFS (9.7 expansions/s per process; balls 1, 13, 83, 386, 1 584, 6 354, 26 353 for 3060000 and 1, 11, 64, 282, 1 088, 4 046, 15 401 for 3304000, 0 keyless children dropped). The pair 3060000 vs 3304000 (certified derived inequivalent, same Coxeter key, F-010) has no meet at 5 + 5 and none at 6 + 6: 0 WL-hash coincidences between the two balls, so not one isomorphism test was even needed. That NO is what the truth predicts. Its power is shown on one side only: 3060000 joins its own 10-step random-walk end at 4 + 5 (total 9) and joins 6000030 (same-profile orbit, F-047) at 4 + 5 (total 9, E-175 needed 6 + 6 with labelled keys), both chains Hom-verified with isomorphic meets. On the other side, 3304000 does not join 2400230 (same F-047 profile and quipu name P^(1,2)_(1,1,2); equivalence rests on that naming only, no certificate) at 5 + 5, 6 + 7 or 6 + 8 (balls 15 401 and 20 210). So, if F-047's naming is right, a NO for any pair whose one end is in the 3304000 class has no demonstrated power; the 6 + 6 NO for this pair rests on the 3060000-side controls alone and on the other end being unable to meet even its own class.

What a NO shows: no F/R chain of at most 6 + 6 key-kept J = 0 tilting steps links these two algebras, which is consistent with the inequivalence certificate. What it does not show: inequivalence (that comes from Happel's theorem on the two cospectral quipus plus the class naming, F-010/F-047, not from this search); that no longer chain exists (depth bound; the YES control needed 9 of 12); anything about the J = 0 premise (steps with J != 0 are excluded by the move set, never tried; counts in Response 2); or any pair other than this one. What would refute the NO's reading: a meet at any depth (it would be a false positive of the iso test, a bug, or a wrong certificate).

## Evidence

| run | depths | ball sizes | WL-bucket iso tests | verdict |
|---|---|---|---|---|
| NO: 3060000 vs 3304000 | 5 + 5 | 6 354 / 4 046 | 0 | no meet |
| NO: 3060000 vs 3304000 | 6 + 6 | 26 353 / 15 401 | 0 | no meet |
| NO: 3304000 vs 6000030 (inequivalent, other profile) | 5 + 5 | 4 046 / 6 244 | 0 | no meet |
| YES: 3060000 vs end of 10-step walk (seed 7) | 6 + 5 | 26 353 / 19 305 | 1 | meet total 9 (4 + 5), Hom-verified, iso meet |
| YES: 3060000 vs 6000030 | 5 + 5 | 6 354 / 6 244 | 1 | meet total 9 (4 + 5), Hom-verified, iso meet |
| YES: 3304000 vs 2400230 | 5 + 5, 6 + 7, 6 + 8 | 15 401 / 20 210 | 0 | NO MEET (control fails) |

Reading. The two YES controls on the 3060000 side find exactly one WL-equal pair each, an isomorphism, and a chain whose every edge passes the independent Hom(T,T[+-1]) = 0 and Cartan test: the test can say YES, and the same code finds no coincidence for the inequivalent pair. The failed control is not a bug in the pair test; it is the only evidence about the other class: 3304000 and 2400230 are equivalent by the quipu naming of F-047 but not joined by 14 steps of this move set. Possible causes, not separated: the equivalence needs a step outside the J = 0 / key-kept / vertex-preserving moves, a chain longer than 14, or the quipu naming of that orbit (the sentence in F-047 is the only source; grep finds `2400230` nowhere else) is wrong. Sizing: depth 6 per side costs about 15 min of one core (two sides in parallel, 4 cores); depth 7 would be about 20 000 and 12 000 expansions more per side, i.e. about 40 and 20 min each, so 7 + 7 is a 1-hour job, not run (proposal below). The earlier E-212 5 + 5 hit the cap because it hashed and compared in one process; here balls are stored and met separately (meets 8 to 37 s).

## Reproduction

All from the repository root; the script writes a checkpoint (resumable) so a ball can take several < 10 min commands. D = any directory for checkpoints.
```
timeout 10m .venv/bin/python -W ignore workshop/rounds/077/skeptic_f010join.py plan 3060000 120      # 9.7 exp/s, sizes to depth 4
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/077/skeptic_f010join.py side A 3060000 6 540 $D/ck_A.pkl   # run twice (540 s each): depth 6, 1 011 s CPU
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/077/skeptic_f010join.py side B 3304000 6 540 $D/ck_B.pkl   # twice: 822 s CPU
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/077/skeptic_f010join.py side C walk:7:10 5 540 $D/ck_C.pkl # control walk, depth 5 (two runs)
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/077/skeptic_f010join.py side D 6000030 5 480 $D/ck_D.pkl   # 163 s
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/077/skeptic_f010join.py side E 2400230 8 480 $D/ck_E.pkl   # three runs, 770 s CPU
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/077/skeptic_f010join.py meet $D/ck_A.pkl $D/ck_B.pkl      # 29 s
```
(`meet` uses each ball's completed depth; copy a checkpoint first if you will resume it. Take a depth-5 copy of A, B, C by running `side` with depth 5 before 6.)

## Prior record

E-175 (6 + 6, labelled `canonicalKey` meets, 3060000 vs 3304000 and controls up to total 9), E-178 (n = 10), E-212 (iso meets; 4 + 4 no meet, 5 + 5 NO control hit the cap, "NO power thin"), E-204 (7 + 7 labelled, no join). New here: the 5 + 5 and 6 + 6 iso-meet NO result (toolsmith's proposal for OVERNIGHT is now done in under 1 hour of CPU, no overnight needed), the same-profile-orbit YES control 6000030 at 4 + 5 (E-175 had it at 6 + 6 labelled), and the failed 3304000-side control, which is not in `research/` (the E-212 YES control was the 3060000 side only, as here). Not claimed: no J != 0 step is tested (they are excluded by construction).

## Code changed

New file only: `workshop/rounds/077/skeptic_f010join.py` (standalone, resumable; reuses `toolsmith_tiltpath.py` and `skeptic_tilt.py` by exec). No library file touched; no tests apply.

## Next

- toolsmith/skeptic: decide why 3304000 and 2400230 do not join: (a) 7 + 7 on the 3304000 side (about 20 min per side for the missing level 7, E-ball level 8 -> 9 about 25 min), (b) check by an independent route (Cartan/Z-conjugacy profile of both is equal per F-047, Ladkani Cor 3.15 on their quipus) that 2400230 really is in the 3304000 class, (c) test whether a single J != 0 or key-moving step bridges the balls (the J != 0 rejections, 69 in the 3304000 ball vs 0 in the 2400230 ball, are the first candidates).
- OVERNIGHT proposal (not run): `side` to depth 7 for A, B (about 1 h on 2 cores) and meet 7 + 7; only worth it if the 3304000-side control is repaired first.
- Anyone quoting "no join at 6 + 6" must add: one side's equivalents are joined by this test, the other side's are not.
