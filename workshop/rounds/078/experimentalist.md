# Tilting-only join of 90000000 / 00000030 at n = 10: no shared labelled key at depth 5 + 4 (MISS, calibrated by a same-machinery positive control); n = 15 scan shard plan (proposal)

author: experimentalist · round: 078 · kind: negative (join) + proposal (shards)
thread: S-1 · bears on: E-190, E-210, E-218, E-173, T3/T8
scope: n = 10, one pair (rows 9,0,0,0,0,0,0,0 and 0,0,0,0,0,0,3,0); J = 0 F/R tilting moves of rounds/049 (gate, tiltingPlus, legality, Coxeter key kept, vertex set kept); meets by `fingerprint.canonicalKey`; depth 5 from the 9-relation side + 4 from the 3-relation side, cap 30 000 keys (not reached); one run; move set not Hom-verified at n = 10. Shard plan: sizing only, nothing of the n = 15 scan was run.

## Response to referee

1. Positive control (done): `experimentalist_join10_ctrl.py ctrl 1`, output `experimentalist_join10_ctrl_ctrl_1.txt`. Random F/R walks from 90000000 (walk end X, all steps by the same `moves`), then ball(root, a) vs ball(X, b) with the same `sideBall`/`canonicalKey`: L = 6, radii 3+3, three walks: HIT at depth 4, 4, 6. L = 8, radii 4+4, two walks: HIT at depth 6, 8. All five HIT at a depth <= L, so the machinery joins when a path exists. Limits: this tests the join and key stability on produced algebras, not completeness of the move set (X is built from the moves), and not a long-path pair like the E-190 one; the 5+4 MISS is calibrated only up to "paths made of these moves are found".
2. Targets in the balls (done): `experimentalist_join10_ctrl_targets.txt`. 323 rows with the key, 646 LNAs plus duals, 646 distinct canonicalKeys, none None. A-ball (r5, 18 776 keys) contains exactly 1 of them (the root 90000000); B-ball (r4, 3 870 keys) exactly 1 (the root 00000030). So no other LNA or dual of the target class is reached by either side; the hereditary side is not reached at these radii.
3. Shard plan output saved: `experimentalist_shardplan15_out.txt` (regenerated; key 2.42 ms/row, 15 images/row, 31.68 ms/row for all deletions, 91 202 CPU-s, 63.3 min per shard of 111 435 rows). Within noise of the 62.5 min / 33.7 ms quoted below (machine load differs; the 200 000-row enumeration took 18.3 s here). Quoted numbers below are the first run; this file is the checkable one.
4. Claim wording (done): claim (1) now states the rounds/049 move set, no Hom verification at n = 10, and the labelled-key caveat. Title and Scope narrowed accordingly. Claim (3) stays marked proposal.
Not done: the 16x deletion-key cost rests on the same 400-row sample as item 3 (now file-backed); levels 4-5 and the killed `run 5 5` were not re-run (6 m 49 s, not repeated).

## Claim

(1) With the rounds/049 J = 0 F/R move set (not Hom-verified at n = 10; no step rejected means only that, not completeness) and labelled `canonicalKey` (which can hide an iso-level meet, E-164), the two sides of the E-190 pair share no canonicalKey in the J = 0 F/R ball of radius 5 around 90000000 (18 776 keys) and radius 4 around 00000030 (3 870 keys). So no tilting-only path of total length <= 9 joins them in the labelled-key sense. This is a bound for these depths, not inequivalence (the pair is equivalent, E-190), and not evidence about a longer path (E-164 relabelling can also hide a meet that exists). A hit at greater depth, or an iso-level meet at this depth, would refute "needs longer".
(2) Size for the full control: levels grow x4.3 per level (about 23 ms per new key, dominated by `moves`); a join to total depth 12 (6 + 6) costs about 50 CPU-min, 14 (7 + 7) about 4 h per side. Not a verdict that 9 is enough.
(3) An n = 15 key-and-deletion scan fits 24 shards under 10 min only if deletion keys come from an n = 14 key table (below); computing 15 deletion keys per row costs 62 min per shard.

## Evidence

Join-search sizing (`plan`, 4 s for the class table): n = 10 has 4862 rows; 323 have the key (1,1,0,0,0,0,0,0,0,1,1) of 90000000 (and so of 00000030, equal keys confirmed), i.e. 646 LNAs plus duals as possible targets. At the root each side has 16 moves (8 F, 8 R), all J = 0 (stats: gate 15 082, J 0, tp 0, key 0, ok 15 082: every gated step in these balls was admissible, no key-changing step occurred).

| level | A = 90000000 new keys | A time (cum.) | B = 00000030 new keys | B time (cum.) |
|---|---|---|---|---|
| 1 | 16 | 0 s | 16 | 0 s |
| 2 | 126 | 1 s | 124 | 2 s |
| 3 | 710 | 12 s | 674 | 12 s |
| 4 | 3 376 | 70 s | 3 055 | 68 s |
| 5 | 14 547 | 340 s | (not run) | est. 330 s |
| 6 | est. 63 000 | est. 25 min | est. 60 000 | est. 25 min |

Search `run 5 4 30000`: 6 m 49 s, shared keys 0, MISS. An earlier `run 5 5` was killed by the 10-min limit after B level 4 (log not kept) and gave the same A and B numbers up to level 4. Depth bound stated: per command, radii (5, 4) total 9; (5, 5) needs about 11-12 min so it does not fit one command; a split (A-ball pickled, then B) is the toolsmith's job.

n = 15 shard plan (`experimentalist_shardplan15.py plan 400`, 400 rows spread evenly over the index range): R = Catalan(14) = 2 674 440 rows; `allRelationLengths(15)` builds the whole tuple in 16 s (then 0.06 s per pass), so each shard pays 16 s and slices by index. Measured: key 2.06 ms/row; each deletion image (15 per row, via `removeVertex`) another 2.1 ms: key + 15 images = 33.7 ms/row, 90 000 CPU-s, 62.5 min per shard of 111 435 rows. Key alone: 5 500 CPU-s, 3.8 min per shard. Proposed layout (all under 10 min, none run):
- Stage A, 4 jobs: key of every n = 14 row (742 900 rows x 2.1 ms = 1 560 s; 390 s each plus 16 s enumeration), written to a keyed table (row -> key).
- Stage B, 24 jobs: shard k of 24 = rows [k*111 435, (k+1)*111 435) (shard 0 starts at 0000000000000, shard 1 at 2245000044002, shard 23 at 0030023030400): row key (230 s) + all 15 image keys by lookup in the Stage A table (dictionary, seconds) + the S-1 comparison (is the image key the same for a pair of deletion choices on rows of equal n = 15 key; counts per key class, with the lone-3 K0 = 5 classes (5,6)/(6,5) reported separately). Wall time with 4 cores: Stage A 7 min, Stage B 24 x 4.1 min / 4 = 25 min.
- Stage C, 1 job: merge counts. Output per shard under 200 KB: (key, image-key multiset) counts, not rows.
Class level: this is key level (E-173 scope). A derived class or orbit label per row needs `freeMoves.derivedOrbits` at n = 14/15, which was not sized here (orbit sizes at n = 9 grow as n 2^(n-5) for the B_n family, E-218; the full partition at n = 15 is unsized). If class level is wanted, the unit of the next sizing is `derivedOrbits(14)` on one shard.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/078/experimentalist_join10.py plan      # 4 s class table + 4-level balls, ~2 min
timeout 10m .venv/bin/python workshop/rounds/078/experimentalist_join10.py run 5 4 30000   # 6 m 49 s, MISS; output in experimentalist_join10_out.txt
timeout 10m .venv/bin/python workshop/rounds/078/experimentalist_shardplan15.py plan 400   # about 2.5 min
.venv/bin/python workshop/rounds/078/experimentalist_shardplan15.py shard 1 24             # 16 s
```

## Prior record

E-190 proved the pair equivalent at n = 10 (Cartan level and kA_10/rad^9 = D_10); E-210 placed images; E-218 said the orbit + mirror free moves do not join the n = 9 images (orbit of 0000030 is 144 rows). E-212/E-049 built the J = 0 F/R tilting search at n = 7 (target ball, canonicalKey meets); this is its first use at n = 10 and on a long-path pair. The "tilting-only join of 90000000 / 00000030" has no record in `research/EXPERIMENTS.md` (grepped E-212, E-190, E-218). E-173 sized the scan at 33 min on 4 cores for key only; the present number agrees (5 500 s) and adds that deletion keys cost 16x more unless tabulated at n = 14.

## Code changed

None in `quivermutation/`. New: `workshop/rounds/078/experimentalist_join10.py` (execs rounds/049/toolsmith_tiltpath.py's moves with base key reset; note that file's top-level exec needs the N = 7 class table, built in seconds), `experimentalist_shardplan15.py`, `experimentalist_join10_out.txt`. No tests (no module touched).

## Next

- toolsmith: pickle the A-ball so (6, 4), (5, 5) and (6, 5) can be done in separate commands; or an iso-level meet instead of labelled canonicalKey (E-164 effect). The control the maverick wants is a MISS at 9 that becomes a HIT at some depth; an overnight run of (6, 6) is about 50 CPU-min per ball set with cap 100 000 keys (memory unsized), propose for OVERNIGHT.md.
- chair: approve Stage A/B above if the key-level S-1 scan is wanted; I would rather first size `derivedOrbits(14)` since class-level is what S-1 asks.
- theorist: the A-ball contains no LNA of the target key class at radius 5 either? (not checked; targets are 646 algebras, a quick count by key would tell whether the hereditary side is reached at all).
