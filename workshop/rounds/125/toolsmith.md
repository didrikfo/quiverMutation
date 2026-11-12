# A second, independent crossing exists: at n = 12, in the key group of 47 (E-314, E-317), the key-kept J != 0 step F3 at 0040033030·F3 R9 F3 changes SNF(S) tail (1,1,20) -> (1,1,12), and a 5-step guarded F/R walk joins the inequivalent LNAs 0040033030 and 0400330030; at n = 11, 4 of the 7 duality orbits of the 13-class of group 5 cross at depth 4, so E-318's "one hit in 20 roots" came from sampling one 13-class root. Every one of the 50 crossing steps found is non-tilting, changes SNF(S), has an R move on its path, and has no gate-admitted way back, so all go from the smaller SNF class to the larger

## Response to referee

1. *"Removes nothing recorded".* Done, reworded and replayed (`toolsmith_refuse19.py`, 312 s, out `toolsmith_refuse19_out.txt`). The 19 audited F-010 J != 0 steps (E-316; 20 joins, since E-233's step closes two) were replayed from the round-123 paths. `isTiltingAtVertex` is False at all 19, so the option refuses every one. Every other step of those 20 joins is J = 0, tilting and key-kept. All 20 joins start and end at the same root (3060000, 3304000, 2400230), so none of them identifies two LNAs. With the option on, the joins are lost as walks: no option-admitted F/R detour of depth <= 3 from the parent P reaches its child X (606-877 nodes per parent). E-316's "in-component" uses undirected J = 0 paths, which a directed walk need not follow. New wording below (Refusal rule, Next).
2. *F-only bound.* Done. E-317's F-only depth-5 balls from 0040033030 and from its dual (round 121, `maverick_ball47.py`) have SNF tail 20 at every node. So X12 is not within F-only distance 5 of either side, and P12 is not within F-only distance 4 (P12 -F3-> X12 is an F step). Added to claim (4).
3. *Reverse step.* Added to claim (2), citing the referee's check (not re-run by me): from none of the 50 children does a gate-admitted F or R step lead back to a node isomorphic to its parent. "Smaller -> larger" is now stated as that fact about the gate, not as a pattern in the sample.
4. *Power.* Restated. Search power comes from rediscovering E-318 (same 4 steps and paths). J-test power comes from the 9 010-step agreement. The planted edge only checks the classifier.
5. *Prior record, counts.* Done. This is the first walk-level failure of H-032 inside its named first case (the 47-group). The ball count is 72 + 13 + 4 + 3 = 92, plus the depth-2 power ball 030033030, so 93. Paired-root node counts differ (3740 / 3780, 1750 / 1782, 3194 / 3198) because `ISO=0` dedupes by the labelled `canonicalKey` (E-164). That inflates the counts and loses no coverage.
- The referee's checks 1-3 (flint SNF on brute Cartan; all 50 rows replayed; 4 of 7 orbits) are cited in Evidence. They settle the sympy-only rows.

author: toolsmith · round: 125 · kind: result + tool
thread: T10 · bears on: H-015, H-032, R-017, E-318, E-321, E-317, E-314
scope: n = 11 group 5 (72 LNAs, 13 + 59 by SNF): all 72 roots at depth 3, the 13 roots of the 13-class and 4 roots of the 59-class (030330030, 230330030, 030330230, 230330230) at depth 4. n = 12: 3 roots at depth 4, one per group of 3 of the 10 groups split by (SNF, signature) (0040033030 group of 47, 0003033030 group of 290, 0030334400 group of 264). Guarded F/R BFS with round-122 semantics (gate, illegal drop, reduce, key guard, every key-kept child expanded, root and dual). 93 balls (72 + 13 + 4 + 3 = 92, plus the depth-2 power ball 030033030); node counts are labelled (`ISO=0`, E-164), so relabelled copies count twice; 105 409 nodes, 270 396 gate-admitted steps. Not covered: 55 roots of the 59-class at depth 4; the other 7 n = 12 split groups and every other n = 12 root; groups split only by q_S at n = 12 (the sizing table skips q_S); depth >= 5. "Crossing" means INV(parent) != INV(child) for one step. That proves the step is not a derived equivalence. A crossing does not by itself mean an LNA-to-LNA join (that needs the walk on, which was checked for 2 sites).
power: search power is the rediscovery of E-318's known hit: 4 steps, 2 up to duality, the same paths. J-test power is the agreement below. The planted edges (030033030 ball -> LNA 003033030: cross = YES; a node's own J = 0 child relabelled J != 0: cross = NO) check only the classifier, one INV comparison each. The library J test `isTiltingAtVertex` agrees with `perI` on 9 010 of 9 010 gate-admitted steps of the 003033030 depth-4 F/R ball (28 J != 0, 8 982 J = 0). Both new hit sites were re-checked with round-122 code outside the library's Cartan routine (below).

## Claim

(1) **n = 12, group of 47.** Start at the LNA 0040033030 (14-class, SNF tail (1,1,20)) and walk F3 R9 F3: three J = 0, key-kept, gate-admitted steps reach P12. At P12 the right mutation F3 is gate-admitted and key-kept, but J != 0 and `tiltingPlus` / `isTiltingAtVertex` are false. The child X12 has SNF tail (1,1,12). Both are acyclic with det C = 1, S has rank 11, and d_11(S) = 20 vs 12. X12 is one J = 0 step (F7, also F9) from the LNA 0400330030, which has SNF tail (1,1,12) (33-class). So this 5-step guarded walk joins two LNAs that are not derived equivalent (E-314's split). It is the shape of E-318's walk (F3 R8 F3, F3 at n = 11) shifted by one vertex. The mirror step on the dual (F9 R3 R3, then R3) gives a child isomorphic to the dual of X12. It counts as one site.
- **Why it is independent of E-318.** It is a different n, a different key group and a different SNF change. E-317 reported this root's depth-4 ball as "all J = 0". That ball was the library's F-only walk. The crossing path uses an R move, so the two results do not conflict.

(2) **n = 11, group 5, re-counted.**
- Depth 3: none of the 72 roots has a key-kept J != 0 step in its ball.
- Depth 4, 13-class: 8 of 13 roots cross, i.e. 4 of the 7 duality orbits: {003033030, 303303000}, {203033030, 303303002}, {303305000, 500033030}, {303500030, 305003030}. That makes 48 crossing steps, 12 distinct up to isomorphism and duality, at 7 distinct parents.
- The 203033030 site is re-checked in full (brute-force Cartan, det divisors 8 -> 4). Its child is one J = 0 step (F6) from the 59-class LNA 230330030, so it is a second LNA-to-LNA join at n = 11.
- No crossing was found in the 4 59-class balls at depth 4, including 030330030 and 230330030, which sit two J = 0 steps from E-318's child X.
- **No way back (referee's check, cited).** From none of the 50 children does a gate-admitted F or R step lead to a node isomorphic to its parent. So the inverse of every crossing step is gate-refused, and a walk from the larger SNF side cannot use these edges. That is why every crossing found goes from the smaller class to the larger (13 -> 59, 14 -> 33), and why the 59-class balls show none, even though 030330030 reaches X at radius 2. It is a property of the gate on these 50 edges. It is not shown for crossing edges not yet found.
- No crossing at n = 12 from 0003033030 or 0030334400. The first of these is the n = 11 root with a 0 prepended.

(3) **Classifier.** All 50 key-kept J != 0 steps in the 93 balls are crossings. Each is classified the same way:
- tiltingPlus false, perI J != 0;
- SNF(S) changes, signature kept;
- q_S undefined, since det S = 0 in both groups, so SNF is the only separator.

So on this sample "key-kept and J != 0" and "crossing" coincide. Not one key-kept J != 0 step was harmless here. E-316's 19 n = 9 steps were all harmless, so this is not a law.

(4) **Every crossing path has an R move.** In all 50 cases, the BFS path to the step, or the step itself, contains an R move. No pure-F (library-walk) path to a crossing was seen. E-318 also found the F-only depth-6 ball of 003033030 clean. At n = 12, E-317's F-only depth-5 balls from 0040033030 and from its dual have SNF tail 20 at every node. So X12 is not within F-only distance 5 of either side, and P12 is not within F-only distance 4. This does not show that a longer F-only route to these parents cannot exist.

**Not claimed:** anything at n = 12 beyond three roots; that the 59-class side never crosses; that the library's F-only walk is safe in general; anything about H-023 (every J = 0 step here keeps INV).

**Refuted by:** a Z-congruence S_P12 ~ S_X12. That would need d_11 = 20 and 12 to be wrong in both routes. Or a wrong Cartan matrix of P12 or X12; brute force agrees with the library on both.

## Evidence

| site (up to duality) | n | group | root orbit | path, step | SNF tail parent -> child | checks |
|---|---|---|---|---|---|---|
| E-318 | 11 | 72 | 003033030 / 303303000 | F3 R8 F3, F3 and R5 | (1,2,4) -> (1,2,2) | E-318 |
| new | 11 | 72 | 203033030 / 303303002 | F3 R8 F3, F3 and R5 | same | brute Cartan, d_10 8 -> 4, joins 230330030 |
| new | 11 | 72 | 303305000 / 500033030 | 3 parents, e.g. F4 R11 R7, F11 | same | sympy SNF; referee: flint on brute Cartan |
| new | 11 | 72 | 303500030 / 305003030 | 3 parents, e.g. F1 R6 R6, F4 | same | sympy SNF; referee: flint on brute Cartan |
| **new** | **12** | **47** | **0040033030** | **F3 R9 F3, F3** | **(1,1,20) -> (1,1,12)** | **brute Cartan, d_11 20 -> 12, joins 0400330030** |

Denominator (balls, by depth):

| set | roots | depth | key-kept J != 0 steps | crossing roots |
|---|---|---|---|---|
| E-318 (round 122) | 20 (n = 9..11) | 4-5 | 6 | 1 (+ mirror) |
| n = 11 group 5, all | 72 | 3 | 0 | 0 |
| n = 11 13-class | 13 | 4 | 48 | 8 (4 of 7 orbits) |
| n = 11 59-class | 4 (+030033030 in round 122) | 4 | 0 | 0 |
| n = 12 (47, 290, 264 groups) | 3 | 4 | 2 | 1 |

**Refusal rule and what it prunes.** The rule: refuse a gate-admitted step when `mutation.isTiltingAtVertex` is false. The library change implements it as an option, default off.
- Recorded merges: F-041 balls, depth 3 per side, root and dual. At n = 8, 10 pairs, 1 600 reached keys, 0 removed, meetings 10 / 10 -> 10 / 10, cost x1.16. At n = 9, 65 pairs, 11 622 keys, 0 removed, meetings 52 -> 52 (E-321's 52 of 52), cost x1.08.
- Walks here: it prunes exactly the 50 crossing steps above, since they are all J != 0.
- The 19 audited F-010 J != 0 steps (E-316; `toolsmith_refuse19.py`): all 19 are refused, and every other step of their 20 joins is admitted. Each join starts and ends at the same root, so no LNA pair is lost. The joins themselves are not recovered within recorded depth + 2: no option-admitted detour of depth <= 3 from a parent reaches its child.
- So it removes none of the library's recorded LNA-to-LNA merges (F-041 n = 8, 9). It refuses the 19 audited J != 0 steps of the F-010 joins (E-316), whose children are in-component (undirected) and whose joins are root-to-root. It refuses every step that crossed.
- Referee's independent checks (round 125 review): flint SNF on brute-force Cartan matrices at all 50 steps (48 times (1,2,4) -> (1,2,2), 2 times (1,1,20) -> (1,1,12)), which settles the sympy-only rows; the n = 12 join rebuilt from scratch; 4 of 7 orbits recounted.

n = 12 sizing (`toolsmith_snf12_out.txt`): there are 10 key groups split by (SNF, signature), from 21 to 3 837 LNAs, 5 952 LNAs in all (the 3 837-group splits by signature only). A depth-4 ball costs 275-389 s per root at n = 12.

## Reproduction

```
ROOTS=030033030 DEPTH=2 CTRL=planted timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxclass.py 11    # 6 s, power
ISO=0 GROUP=5 SLICE=0:18 DEPTH=3 OUT=f.jsonl timeout 10m ... toolsmith_jxclass.py 11    # ~10 min per 18 roots (3 jobs parallel); slices to 72
ISO=0 ROOTS=203033030,303305000,303500030 DEPTH=4 timeout 10m ... toolsmith_jxclass.py 11   # 70-220 s per root
NOLABEL=1 ISO=0 ROOTS=0040033030 DEPTH=4 timeout 10m ... toolsmith_jxclass.py 12              # 389 s, the n = 12 hit
.venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxclass.py --tally workshop/rounds/125/toolsmith_jxclass_out.jsonl   # 2 s, tables above
R=2 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxhit12.py                        # 12 s, n = 12 checks + walk
ROOT=203033030 R=2 WALK=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxhit2.py   # ~60 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_refuse.py f041 8      # 36 s   (f041 9: 400 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_refuse.py agree 003033030 4   # 195 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_snf12.py 12           # 245 s
DMAX=3 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_refuse19.py    # 312 s, the 19 E-316 steps under the option
```

Timings were taken with 3 jobs running on 4 cores. `toolsmith_jxclass_out.jsonl` (80 KB) holds every ball head and every key-kept J != 0 step with its parent and child structure.

## Prior record

- E-318: one crossing, n = 11; reproduced here.
- E-321: no recorded n <= 10 merge has a J != 0 step; the F-041 part is reproduced via the option.
- E-317: 47-group balls all J = 0. Those were F-only walks, so there is no conflict, but its sentence "this group gives T10 no walk-level test at depth <= 5" holds only for F-only walks. Its F-only depth-5 ball from 0040033030 bounds the F-only distance to X12 (claim (4)).
- H-032: it names the 47-group as its first case. Its status cites only E-318 (n = 11) for the walk-level failure. This is the first walk-level failure of H-032 inside that first case. R-017 already holds the general statement.
- E-314: the 14 + 33 split.
- E-316: n = 9 J != 0 steps, none crossing.
- grep of `research/` for 0040033030, 0400330030, 203033030, 303305000, 303500030, `refuseNonTilting`, `isTiltingAtVertex`: nothing.

New: the n = 12 crossing and walk, the second n = 11 join, the 4-of-7 count, and the library option.

## Code changed

- `quivermutation/mutation.py`: new `isTiltingAtVertex(pathAlg, vertex)`. It is Ladkani Prop. 2.3(c), rational linear algebra on the arrow relations, about 40 lines.
- `quivermutation/search.py`: `mutationSearchDepthFirst(..., refuseNonTilting=False)` and `linesReachedFrom(..., refuseNonTilting=False)` are passed through; default off, so behaviour is unchanged. The docstring paragraph says what the option does.
- `tests/test_search_soundness.py`: 2 new tests on E-318's step. One checks that `isTiltingAtVertex` is false there while gate and key admit it. The other checks that the option removes exactly child (3) of P, and that the default equals off.
- Ran `pytest -q tests/test_search_soundness.py tests/test_mutation_procedure.py -m "not slow"`: 14 passed, 1 xfailed (pre-existing).

## Next

- Chair: F-candidate. "Key-guarded F/R walks cross derived classes at n = 11 and n = 12, in two singular-S groups. Every crossing found is J != 0 and non-tilting, and has no gate-admitted reverse step. The refusal removes none of the recorded LNA-to-LNA merges (F-041 n = 8, 9). It refuses the 19 audited F-010 J != 0 steps (E-316), whose joins are root-to-root."
- Toolsmith / overnight: the remaining 55 59-class roots at depth 4 (about 2.5 h CPU). The 47-group's other 13 14-class roots and the other n = 12 split groups, one root per class, at depth 4 (about 20 roots x 6 min). The F-041 n = 10 option check (E-321's 117 642 edges).
- Skeptic: find out whether a pure-F route reaches P12 (F-only depth 6-7 from 0040033030). The SNF re-derivation is done (review).
- Theorist: the step is "F3 after F3 R(k) F3" at the start of the line in both n. Is there a local configuration that forces J != 0 with the key kept? Why is the reverse of each crossing gate-refused?
- Chair: decide whether `refuseNonTilting` becomes the default for `linesReachedFrom` (cost 8-16 %).
