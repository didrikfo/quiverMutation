# The 3304000 ball to depth 5 holds one key-kept perI J != 0 step (1 of 2 958 parents; versus 10 of 4 770 at 3060000, no rate claim); it is real and its child joins 3304000 at total 6; the 10 E-234 children fall in 8 iso classes

## Response to referee

Verdict minor revision; reproduction matched, nothing to re-run (the referee's re-runs agree; `ctl` not re-run by either of us, this sitting either).
1. S = 0 Hom values: done, stated in the Claim (Hom(T,T[-1]) = 1 at (7,1), Hom(T,T[1]) = 1 at (3,7); source lines in `toolsmith_d5hom_out.txt`). "Hom(T,T[1]) = 0" now explicitly belongs to S = -1.
2. "Not 3060000" contrast dropped (only the 3304000 ball was tested; E-219/E-226's NO is one-sided, not re-tested). "Different order" sentence replaced by the 10-events-in-8-classes statement; claim now reads "1 of 2 958 versus 10 of 4 770, no rate claim"; the percentages are removed from Claim and table.
3. E-204 sentence added: false for both roots (E-229/E-234).
Gaps 2 and 4: key equality cited as E-057/F-010; rate language now says flagged-and-kept, not all J != 0 steps.


author: toolsmith · round: 085 · kind: result
thread: T10 · bears on: E-229, E-233, E-234, E-204, E-219/E-226 (J = 0 premise)
scope: n = 9, root 3304000 (the other F-010 side; its Coxeter key equals 3060000's, asserted), E-212 `moves` ball to depth 5 (4 046 nodes), the 2 958 nodes at depth 5 scanned (depth <= 4 is E-229: 0 kept); perI-flagged gate-admissible F/R steps only, key-kept = child key equals the class key; Hom with the placements S = -1, 0; join = J = 0 chain to a depth-5 root ball with iso meet, child ball depth 3; generation by T assumed (E-168). Distinctness: quiver+relations isomorphism up to vertex relabelling, 12 children (10 of E-234, E-233's, this one). Rate comparison with 3060000 not claimed. Not covered: depth 6, multi-vertex J != 0 steps not flagged by perI, any other root.

## Claim

Scanning the 2 958 depth-5 parents of the 3304000 ball: 353 perI-flagged candidate steps (101 + 252 moved) + 1 kept; exactly ONE key-kept perI step (index 2026, path F2 R8 R9 R7 R7, step R at v = 7): 1 of 2 958 parents (1 of 354 flagged candidates), versus 10 of 4 770 (10 of 640) for 3060000; these count flagged-and-kept steps, not all J != 0 steps, and no rate difference is claimed. It is a real J != 0 step: tiltingPlus False, Hom(T,T[-1]) = 2 at S = -1 (entries (5,7), (8,7)) with Hom(T,T[1]) = 0; at S = 0 Hom(T,T[-1]) = 1 at (7,1) and Hom(T,T[1]) = 1 at (3,7) (`toolsmith_d5hom_out.txt` lines 5, 7), End(T) dims = Cartan(child)^T. Its child joins 3304000 at total 6 (child R2 || root F2 F1 F6 F5 F2), Hom-verified, meet iso. So on this side too a key-kept J != 0 step exists and does not leave the J = 0 class of its own root. The other half of E-229's "0 of 69" at depth <= 4 stands; E-204's phrase "no J != 0 step enters the balls" is false for both roots (3304000 here; 3060000 by E-234, with E-229's depth <= 4 count unchanged). Not claimed: that the rate differs between the two roots in any structural sense (1 versus 10 events; no test); anything at depth 6; the join is a consistency check only (E-228: Hom+Cartan has no NO power, the iso meet is the check). A key-kept step whose child fails to join its own root would refute "none leaves".

Distinctness of the E-234 children (answers the open item): the 10 children are NOT pairwise distinct. Isomorphism classes (up to vertex relabelling): 8, with {1915, 4722} and {1970, 2003} merged; adding E-233's child (path F1 R9 R8 R8) it joins {1915, 4722}, so the 11 steps of the 3060000 ball give 8 child classes ({1915, 4722, E-233}, {1970, 2003}, and 6 singletons). The 3304000 child is isomorphic to none of the 11 (as it must be if 3060000 and 3304000 are inequivalent; a free consistency check, since an iso would have contradicted E-057). Parents: 9 classes among the 12.

## Evidence

Ball build 116 s (sizes 1, 11, 64, 282, 1 088, 4 046; gate 10 564, J != 0 69 in the ball, matching E-229's 10 564 and 69, so the same move set). Scan in 2 shards 31 s + 36 s (0-1500: moved 101 kept 0; 1500-2958: moved 252 kept 1). Positive control `ctl` on the E-233 parent returns R 8 kept (F 6 moved). Join: root ball depth 5 111 s, child ball (543 nodes at depth 3) 10 s, meet 2 s, WL-bucket iso tests 3. Hom output: `toolsmith_d5hom_out.txt`; join: `toolsmith_d5join_out.txt`; scan: `toolsmith_d5scan_out.txt`; distinctness: `toolsmith_distinct_out.txt`.

| quantity | 3060000 (E-234) | 3304000 (this) |
|---|---|---|
| ball sizes to depth 5 | 6 354 | 4 046 |
| depth-5 parents | 4 770 | 2 958 |
| perI candidates (moved + kept) | 640 | 354 |
| key-kept | 10 (+1 at depth 4) | 1 (0 at depth <= 4) |
| per parent (descriptive only) | 10/4 770 | 1/2 958 |
| Hom(T,T[-1]) total at S = -1 | 2-4 | 2 |
| child joins own root at | 7-8 | 6 |

The 3304000 parent shares the shape of the E-233/E-234 family: perI on vertices 5, 8 (two vertices) with diamond relations (7 -> {2, 4} -> 5 and -> 8); the child quiver has 10 arrows and 4 relations.

Control for the distinctness test: it merges 1915/4722/E-233, whose child-ball sizes (1, 14, 103, 579) already agreed in round 083; it separates 9 of 12 WL hashes before iso testing; the test is the one the joins use (`skeptic_join.py` `iso`).

## Reproduction

All from the repository root; /tmp paths are regenerated by `build`:
```
mkdir -p /tmp/ts85
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/085/toolsmith_d5scan.py build 3304000 /tmp/ts85/par5.pkl    # 116 s (env DEPTH=5 default, KEYROOT=3060000 default = same key)
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/085/toolsmith_d5scan.py scan /tmp/ts85/par5.pkl 0 1500       # 31 s; then 1500 2958 (36 s)
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/085/toolsmith_d5scan.py ctl workshop/rounds/082/skeptic_kept.pkl
.venv/bin/python -W ignore workshop/rounds/085/toolsmith_d5pkl.py /tmp/ts85/par5.pkl /tmp/ts85 2026
.venv/bin/python -W ignore workshop/rounds/082/skeptic_kept.py test /tmp/ts85/k2026.pkl
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side R 3304000 5 560 /tmp/ts85/ckR        # 111 s
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side C pkl:/tmp/ts85/k2026.pkl:child 3 150 /tmp/ts85/ckC   # 10 s
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py meet /tmp/ts85/ckC /tmp/ts85/ckR
# distinctness: build 3060000 the same way into /tmp/ts85b/par5.pkl (131 s), then
.venv/bin/python -W ignore workshop/rounds/085/toolsmith_d5pkl.py /tmp/ts85b/par5.pkl /tmp/ts85b 1275 1433 1746 1818 1915 1932 1970 2003 4551 4722
.venv/bin/python -W ignore workshop/rounds/085/toolsmith_distinct.py /tmp/ts85b/k1275.pkl ... /tmp/ts85b/k4722.pkl /tmp/ts85/k2026.pkl workshop/rounds/082/skeptic_kept.pkl   # about 1 min
```
`--plan`: the build is the cost and scales as the ball (depth 3: 8 s, depth 4: 29 s, depth 5: 116 s); no run exceeded 4 min, so no sharding beyond the scan's two halves was needed. Depth 6 parents (~12 000 on this side by the growth 3.7 per level; ~19 000 for 3060000) would need a ball build of about 8 min here and 10+ min for 3060000: borderline, see Next.

## Prior record

E-229: 3304000 has 0 key-kept of 69 J != 0 candidates in balls to depth 5 (the 69 are the gate-level candidates, pre-key) and "none" at parent depth <= 4: consistent, this finds one at parent depth 5 (E-229 counted J != 0 gate candidates, mine counts perI parents at depth 5). E-234: 10 on the other root, distinctness and 3304000 listed open: both answered here. Nothing in RETRACTIONS touches these.

## Code changed

New scripts only: `toolsmith_d5scan.py` (copy of 083's with env KEYROOT and DEPTH instead of the hard-coded 3060000 base and depth 5), `toolsmith_d5pkl.py` (path to 085's scan script), `toolsmith_distinct.py`. No library change, no tests run.

## Next

Skeptic: re-run the one hit (index 2026) and check the distinctness merges by an independent isomorphism (e.g. canonical form of the relation graph). Experimentalist: the 7 + 7 NO test for the pair, now that J != 0 steps are known to be rare and J = 0-enclosed; depth-6 scan of 3304000 (about 8-9 min build, borderline, shard by checkpoint) to see whether the rate stays near 0.03%. Chair: E-234's "10 children" is "8 classes"; 3304000 gives 1 of 2 958 against 10 of 4 770; the 10 events fall in 8 classes sharing the diamond shape (R at vertex 8), so they are not independent and no rate test is possible.
