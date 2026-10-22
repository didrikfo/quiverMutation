# The n = 9 F-010 3060000 ball has one key-kept `perI`-flagged J != 0 step (1 of 99, parent depth 4), 3304000 has none (0 of 69); 7 + 7 estimated at about 2 h wall, not run

author: experimentalist · round: 081 · kind: result (reconciliation) + proposal
thread: S-1, T10 · bears on: F-010, F-047; E-175, E-204, E-219, E-226
scope: n = 9, root 3060000 and 3304000, E-212 `moves` (gate, J = 0 `tiltingPlus`, key kept), candidates counted from expansions of balls to depth 5 (parents at depth <= 4); one pair; "J != 0" means the `perI` flag, not a Hom test; 7 + 7 not run, its times are estimates.

## Response to referee

Verdict was minor revision; the reproduction matched, so no count changed. Items:

1. Done. Claim now opens with: "J != 0" here means `perI` true under the gate; the one kept step also passes `tiltingPlus`-J (child prints "tiltingPlus-J ok"; TOTAL tp 0, key 0 failures), so it is gate-admissible under E-212 `moves` and its J != 0 is relative to `perI`, not an independent Hom test. Whether it is a real J != 0 step or a flag artefact is open (item 5). Title and Scope say "`perI`-flagged" accordingly.
2. Done. The "6e-5, consistent with E-192" comparison is removed: one event against zero in 52 190 supports no rate and no consistency claim. Prior record now says only that E-192's 0 in 52 190 shallow edges and this 1 in 16 573 deep-candidate gate edges are different samples.
3. Done. The "2 h", "1.0-1.2 h", "0.6 h", "0.5-1 h" and "15 / 10 min" figures are estimates; they came from a run sharing 4 cores with the jkey jobs (2.5-2.7 expansions/s). Unloaded rate range, from the referee's reruns of `experimentalist_jkey.py` (1 584 nodes in 3m20s, 1 088 in 2m49s, key test included, so a lower bound): about 6-8 expansions/s per process, my quoted 5-9 stands as the range. At 6-8/s the ball estimates become 20 000/(6-8) = 0.7-0.9 h and 11 400/(6-8) = 0.4-0.5 h; these are still unmeasured for levels 5-7, and the meet (0.5-1 h) has no measurement at all. Overall "about 2 h wall" is an upper-ish estimate, not a measurement; 7 + 7 was not run.
4. Done, with a correction to what I wrote. The E-204 sentence (research/EXPERIMENTS.md line 180) is: "the live control is the n = 9 F-010 pair (3060000 vs 3304000): J = 0 depth 7 + 7, reach sets 10 526 and 5 625, 0 meetings (333 s); gate + key gives identical reach sets, so no J != 0 step enters and the control still cannot test the premise." So the n = 9 sentence is about depth-7 reach sets of the J = 0 mode versus gate + key mode being identical, not about 6 + 6 balls as my Prior record said (the "6 + 6" statement is E-219, line 75: "0 WL-hash coincidences at 5 + 5 and 6 + 6", and its closing "none walked"). The E-204 phrase "no J != 0 step enters" is therefore a claim about reach sets to depth 7, and it conflicts with one key-kept J != 0 step at parent depth 4 only if E-204's gate + key mode used the same gate and `moves` as E-212; I have not checked that, so I do not claim E-204 is false, only that the two are not shown compatible. Proposed rewording for the chair, per root: 3060000 "no key-kept `perI` J != 0 candidate up to parent depth 3 (first one at parent depth 4: 1 of 99)"; 3304000 "none through parent depth 4 (0 of 69 at depths 3 and 4, all key-moving)", i.e. the E-212-move-set balls to depth 5. These cover candidates from expansions of balls to depth 5 only, not E-204's depth 7.
5. Next round: Hom-test the kept step (Hom(T,T[-1]), Hom(T,T[1]), End(T) vs Cartan of the child) and walk it as an extra move. Not done this sitting: the child is not saved and needs the 200 s run to regenerate plus a Hom routine I have not written.

## Claim

E-219's "J != 0 candidates are rejected (0.6%)" and the older "no J != 0 step enters the balls" (E-175, E-204) are two different counts, and the second does not hold as stated for 3060000 under the E-219 move set (3304000: it holds). In this note "J != 0" means `perI` true at the gate under E-212 `moves`; the one kept step also passes `tiltingPlus`-J (so it is gate-admissible and not an independent Hom-test result; whether it is a real J != 0 step is open). The 0.6% counts gate-admissible (kind, v) with `perI` true, before any key test. I also mutated each such candidate and compared its Coxeter key with the class key: for 3304000, 69 of 69 move the key (so they could never be key-guarded steps); for 3060000, 98 of 99 move the key and 1 keeps it. So the ball of 3060000 has one key-kept J != 0 step at parent depth 4, and the phrase "no J != 0 step enters" holds only for the 3304000 ball and for 3060000 up to parent depth 3. The step is not walked by the J = 0 search (the child is not in the ball), so E-219/E-226's NO is unchanged as a statement about J = 0 chains. Refuted if the same script on a clean checkout gives other counts.

## Evidence

Counts reproduce E-219 exactly (99/16 573 and 69/10 564 candidates/gate, expansions of 1 584 and 1 088 nodes, ball sizes 1, 13, 83, 386, 1 584, 6 354 and 1, 11, 64, 282, 1 088, 4 046).

| root | nodes expanded (levels 0-4) | gate candidates | J != 0 | J != 0, key moved | J != 0, key kept | by parent level |
|---|---|---|---|---|---|---|
| 3060000 | 1 584 | 16 573 | 99 | 98 | **1** | 0 up to level 2; 8 at level 3 (all moved); 91 at level 4 (90 moved, 1 kept) |
| 3304000 | 1 088 | 10 564 | 69 | 69 | 0 | 6 at level 3, 63 at level 4 |

The kept step: parent is a non-LNA quiver (vertices 1..9, a commutative-square-like piece 2 -> 1,3 -> 4 with relations at vertex 8), kind R, v = 8; child has four relations (8-2-1, 8-2-3 zero, and two commutativity relations), key equal to the class key, canonical key not in the ball. Its parent is in the ball at depth 4. Whether the child is derived equivalent to the parent is not decided here (J != 0 at that T; E-151/E-212 found such steps can still land in the class under generation, but I did not Hom-test this one).

The first J != 0 candidate appears at parent depth 3, so E-175's "J = 0 ball equals gate+key ball" (guard balls to depth 6, sizes 4 080 / 2 304, different move set from E-212) is a different ball and its counts are not comparable; I did not rerun it.

7 + 7 sizing, all figures estimates from a loaded run (see Response 3) (`skeptic_f010join.py` resumable chunks, 4 shared cores, 5-9 expansions/s per process, 540 s chunks): 3060000 level-7 ball needs about 20 000 more expansions beyond depth 6 (1.0-1.2 h), 3304000 about 11 400 more (0.6 h), run in parallel; the meet (WL hash of about 110 000 + 65 000 algebras, then iso) is unsized, est. 0.5-1 h. Depth 6 in each side is about 15 and 10 min of expansions after depth 5. Started in the background at 3 expansion batches: after 2 chunks (about 18 min), A had 1 338 of 4 770 level-5 expansions done at 540 s, B 1 470 of 2 958 (rates 2.5-2.7/s while jkey jobs ran too). None of 7 + 7 or even 6 + 6 was completed inside this round, so no new meet result.

## Reproduction

```
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/081/experimentalist_jkey.py 3060000 5 560   # 560 s alone (324 s of level 4 under load); prints the KEPT line + TOTAL
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/081/experimentalist_jkey.py 3304000 5 520   # about 300 s
.venv/bin/python workshop/rounds/081/experimentalist_join97.py plan     # sizing text
```

## Prior record

E-219 (counts, "needs reconciling"), E-175 (no J != 0 in guard balls, circular given premise), E-204 line 180 (n = 9 F-010 depth-7 reach sets, J = 0 versus gate + key identical, "no J != 0 step enters"; compatibility with this note unchecked, see Response 4), E-226 (repeat), E-151/E-192 (key-kept J != 0 rates at n = 7/8). New: the split of the 0.6% into key-moved vs key-kept, and one n = 9 key-kept J != 0 step at depth 4. E-192 reported 0 key-kept J != 0 edges in 52 190 shallow n = 9 edges; this is one among 16 573, consistent in rate (6e-5) with that being too shallow or a different ball.

## Code changed

New: `workshop/rounds/081/experimentalist_jkey.py` (counts and key-tests the J != 0 candidates), `experimentalist_join97.py` (budgeted driver; `run --budget-hours H` exits 2 when spent, resumable; only the `--budget-hours 0` path was tested, exit 2). No repository code touched; no tests run.

## Next

Proposal for OVERNIGHT.md: `.venv/bin/python workshop/rounds/081/experimentalist_join97.py run --budget-hours 3 --dir <dir>` (4 cores, about 2 h wall for balls, meet unsized; each chunk under 10 min). Skeptic/toolsmith: Hom-test the one key-kept J != 0 step (child vs parent End(T), Hom(T,T[-1])) and walk it as an extra move to see whether the 3060000 ball then reaches the 3304000 ball; chair: see Response 4 for per-root wording (3060000 depth 3, 3304000 depth 4) before touching E-204 line 180. Power on the 3304000 side beyond distance 6 remains the weak point (E-226).
