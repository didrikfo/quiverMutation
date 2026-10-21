# Experimentalist notebook (rewritten round 078)

## What I now believe
- Row convention: an n-vertex LNA is an n - 2 entry row; `nk.LinearNakayamaAlgebra(n, list(row))` builds it (`fromRelationString` wants the table's 'a;b|c' form, not digits). Check digit count against n.
- 075: A_n = rad^(n-1) vs B_n (3-relation on last four vertices), n = 6..12: same key, A_n has no move, B_n orbit n*2^(n-5), never joined by orbit + mirror except n = 5. Free-stretch rule (D >= 3) mostly vacuous; E-190 pair is outside it.
- 078: tilting-only J = 0 F/R join of 90000000 / 00000030 at n = 10: MISS at radius 5 + 4 (18 776 + 3 870 keys, 6m49s). Ball growth x4.3 per level, 23 ms/key; 16 moves at each root; no key-changing or J != 0 step ever occurred. 6 + 6 ~ 50 CPU-min. A miss is a bound (labelled canonicalKey, E-164 effect).
- 078: n = 15 scan: key 2.06 ms/row (3.8 min per 24th), each deletion key 2.1 ms (62 min per shard if 15 taken per row). Layout: n = 14 key table in 4 jobs, then 24 shards with lookups. `allRelationLengths(15)` tuple takes 16 s, then free. Class level (orbit labels) unsized.
- Coxeter key == class at n = 8, 9; no fixed i or class+i+coverage rule predicts the image (E-194, E-200). Exact classifyLength(10): seed 61 % in 311 s, main pass unsized.

## What I tried
- 078: `rounds/078/experimentalist_join10.py` (plan / run DA DB CAP; execs rounds/049 tiltpath moves), `experimentalist_shardplan15.py` (plan / shard k K). Gotcha: the exec'd tiltpath source overwrites a variable named ARGV; and a 10-min `run 5 5` is killed (needs ~11 min).
- 075: `experimentalist_join9.py`, `_delpair.py`. 073: `_pair.py`, `_plan10.py`. 070: `_free10*.py`. 067/065: deletion tables.
- Lessons: grep EXPERIMENTS.md first; mirror needs len(row)+2; forward orbit absent != not equivalent; piping through tail hides progress, tee to a file.

## What I would do next
1. Check whether the radius-5 A-ball reaches any LNA of the key class (646 targets): is the hereditary side even approached? Cheap.
2. Pickle the balls and do (6,4) / (5,5) in separate commands; overnight (6,6) with iso meets.
3. Size `derivedOrbits(14)` on one shard for class-level S-1; then run Stage A/B if approved.
4. S-1 item 3 (long core equivalent to quipu LNA) and item 2 (H-020 head/tail at n = 12 -> 11).
5. Exact classifyLength(10) overnight, then redo free10 on exact labels.
