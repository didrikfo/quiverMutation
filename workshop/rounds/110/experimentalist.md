# At n = 6 key class 0 (480 s capped BFS) all 1150 gate-admitted J != 0 steps change the Coxeter key; the other five cells have no J != 0 step; no key-kept J != 0 step was seen at n = 5, 6

author: experimentalist · round: 110 · kind: result
thread: T10 · bears on: H-015, E-142, E-147, E-151, E-212, E-289
scope: J != 0 steps occur in one of six cells only (n = 6 class 0); the rest is vacuous for the guard. Guarded BFS (expansion only into key-kept children, canonicalKey dedup) from the LNAs and duals of each key class; n = 5 all 2 classes run to closure (6 240 and 5 460 expansions); n = 6 classes 0-3 (all four; class sizes 2, 24, 26, 32), 480 s cap each, 27 682 to 53 788 expansions, BFS levels 10-13, NOT closed. n = 7 not run. Parents are guarded-walk descendants, not only LNAs.
power: the refusal count is the control for "the test can say NO": at n = 6 class 0 it says REFUSED 1 150 times where NO is known true (key differs, so the child is not derived equivalent to the parent; the key is the Coxeter polynomial). The cross-tabulation also separates the cases: 0 of 1 150 refused steps pass `tiltingPlus`, 0 of 0 J = 0 steps refused. The test cannot say NO to a key-kept inequivalent child (none exists in the sample); for that see Next.
backlog: n/a

## Response to referee

Referee verdict: minor revision; all four required items done, none refused.
1. Retitle/rescope: DONE. Title, Claim and scope now say that J != 0 steps occur in ONE of the six cells (n = 6 class 0); the other five (n = 5 both classes, n = 6 classes 1-3) have zero J != 0 steps of any kind and say nothing about the guard. The J = 0 refused = 0 cross-tab (about 570k steps in all) is a harness sanity check (tiltingPlus implies derived equivalence implies key kept), not a finding.
2. Re-implementation vs library: DONE, both ways. The table's "REFUSED" is the script's `ck != pk` (`search._coxeterKeyOrNone` on parent and child), not a call to `mutationSearchDepthFirst`. New check `experimentalist_libguard.py 6 0 500 1` (53 s, `experimentalist_libguard.out`): the same walk stopped after 25 refusals (first at level 8; 3 666 expansions, 25 J != 0 refused, 13 479 J = 0 kept); each of the 25 refused parents was run through the library, `mutationSearchDepthFirst(parent, 1, coxeterGuard=True/False, visitor=...)`: guarded search does NOT visit the child in 25 of 25, unguarded search visits it in 25 of 25. Control: 19 sampled kept J = 0 steps are visited by the guarded search in 19 of 19. So the library guard and `ck != pk` agree on 44 of 44 sampled steps. (Sample, not all 1 150; the library code path is the same comparison `movedKey != baseKey` with baseKey the walk's start key, equal to the parent key on a guarded walk.)
3. "Cannot leave its class": DELETED. It followed only for the sampled region (capped, n = 6 not closed, one non-vacuous class). Now stated as: no key-kept J != 0 step was seen, which is silence about leaving a class, not a proof.
4. Load dependence: DONE, see Reproduction. The 480 s cap is wall-clock; the counts depend on machine load (four jobs ran in parallel) and are not exactly reproducible; the libguard run reached the first 25 refusals at 53 s, the 1 150 at 480 s was a different cut. Only n = 5 (closed) is exactly reproducible.
Novelty: accepted as small (E-142/E-140 already record refusals); the increment is the 1 150 count, J = 0 refused = 0, n = 6 class 3, and the library cross-check.

## Claim

The key guard of `search.py` is not vacuous: on a guarded walk at n = 6, key class 0 (2 LNAs, key (1,1,-1,-2,-1,1,1)), 1 150 gate-admitted steps, all with J != 0 and all failing `tiltingPlus`, have a child with a different Coxeter key (child key (1,1,0,0,0,1,1), the key of n = 6 class 2), as computed by `ck != pk` in the script (a re-implementation; agreement with the library guard checked on 25 refused + 19 kept sampled steps, see Response 2); refused here means key moved. Nothing else is refused: 100 729 gate-admitted J = 0 steps in that walk, and 132 370, 139 736, 167 173 in classes 1, 2, 3, all keep the key and pass `tiltingPlus`. Hence on this sample every key-keeping step is a J = 0 tilting step; no key-kept J != 0 step was met, which is silence, not a proof that a guarded walk cannot leave its class (capped, n = 6 not closed, five of six cells have no J != 0 step at all). It does NOT claim this for n >= 7 (E-147/E-151: key-kept J != 0 steps do occur at n = 7 classes 1, 2), nor for walks beyond the 480 s caps (counts are wall-clock and load-dependent; class 0 had refusals first at BFS level 8, deeper than the earlier depth-6 audit reached). Refuted by one key-kept J != 0 step at n = 5 or 6 whose child is not derived equivalent to its parent.

## Evidence

| n | class (size) | expansions | levels | J=0 kept, tilting | J!=0 REFUSED | J!=0 kept | J=0 refused |
|---|---|---|---|---|---|---|---|
| 5 | 0 (12) closed | 6 240 | 14 | 16 620 | 0 | 0 | 0 |
| 5 | 1 (16) closed | 5 460 | 16 | 13 680 | 0 | 0 | 0 |
| 6 | 0 (2) capped | 27 682 | 13 | 100 729 | 1 150 | 0 | 0 |
| 6 | 1 (24) capped | 38 925 | 10 | 132 370 | 0 | 0 | 0 |
| 6 | 2 (26) capped | 42 226 | 11 | 139 736 | 0 | 0 | 0 |
| 6 | 3 (32) capped | 53 788 | 12 | 167 173 | 0 | 0 | 0 |

Also run: `scholar_h015.py` n = 4 depth 9, n = 5 depth 9, n = 6 depth 6 unguarded: only (guard, tilt, cong) rows (718, 18 417, 29 822), no refusal (depth < 8).
Smallest refusal example (level 8, vertex 6; 6 vertices, in output file): arrows 1->6, 3->1, 3->4, 4->6, 5->3, 6->2, relations 3-1-6-2 = 3-4-6-2 and 5-3-1-6 = 5-3-4-6 (two commutative squares sharing the arrow 3->1/3->4 pair). Its children have key (1,1,0,0,0,1,1), which is another LNA class: the refused step goes to another class, as the key says.
Reading: the refusal needs a "long square" (E-102/E-142, out-degree 1 with a square). The guard acts as a filter only on J != 0 steps, and J != 0 is rare below level 8.

## Reproduction

```
for c in 0 1 2 3; do timeout 10m .venv/bin/python -W ignore workshop/rounds/110/experimentalist_refusal.py 6 $c 480 3 > workshop/rounds/110/experimentalist_refusal_6_$c.out & done; wait   # 480 s each, parallel; wall-clock cap, counts vary with load
timeout 5m .venv/bin/python -W ignore workshop/rounds/110/experimentalist_refusal.py 5 0 100 1   # 33-36 s, exact (closed); class 1: 24 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/110/experimentalist_libguard.py 6 0 500 1   # library-guard cross-check, 53 s
```

## Prior record

E-142 (n = 6 c0: 147 J != 0 steps in a 150 s guarded walk, none keep the key; n = 6 c1, c2 vacuous) and E-140: refusals are known; this reproduces them with the full (J, moved, tilting) cross-tabulation, a larger count (1 150 at 480 s) and the zero cells for J = 0 refused and J != 0 kept at n = 5, 6. E-147/E-151/E-201: key-kept J != 0 steps first at n = 7 classes 1, 2. E-212: 8 of the 9 class-2 ones are tilting-chain equivalent (Hom-verified); record 6 undecided; class 1 not run; "no key-equal inequivalent pair exists at n = 7". So the open piece of T10 is exactly the key-kept J != 0 steps, and none of them was a counterexample so far. The refusal control is therefore no new law; it closes the "no key-refusal control" item of E-289 by pointing to E-142 and the cross-table here.

## Code changed

None in the repository. New scripts `workshop/rounds/110/experimentalist_refusal.py`, `experimentalist_libguard.py` (reuses rounds 033/001 helpers) and its `.out` files (< 5 KB). No tests run.

## Next

- Toolsmith/experimentalist: n = 7 class 1, 13 key-kept J != 0 steps (E-201), 8-class-2-style certificate (`rounds/073/toolsmith_certificate.py` needs new E-157 paths for class 1); record 6 of class 2. OVERNIGHT: `toolsmith_collect.py 7 1 45000` (E-201: 27 min walk) then certificates. This is the only place a guarded walk could leave a class.
- Skeptic: a key-equal inequivalent control at n = 7 does not exist (E-212); the NO power for key-kept children lives only at n = 9 F-010 (E-219).
