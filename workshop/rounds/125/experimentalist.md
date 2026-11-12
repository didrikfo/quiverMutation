# At n = 13 h = 2, 3 and n = 15 h = 3 (lone-3 key classes, head end), one mutation (and derived) orbit carries both image keys at K = K0; the image key changes only through K = 0 rows, and `rules` orbits are almost always one-key because few K = 0 components bridge keys

author: experimentalist · round: 125 · kind: result (revision of rounds/123/experimentalist.md; second revision after the round-125 review)
thread: S-1 · bears on: H-026, E-141, E-146, E-307, E-236
scope: n = 13 h = 2, 3 and n = 15 h = 3 only; the whole Coxeter-key class of the lone-3 row (7459, 5835, 26 762 rows, every K, no cap); head end only; move families `rules` (lnaMoves.rewritesOf), `mut` (+ edgeMoves + doubleMutation), `der` (+ stripLengthTwo), undirected union-find; orbits are lower bounds for derived classes, and that the library moves are derived equivalences is assumed, as in E-141 / E-146 (not re-checked with the mutation engine); image moves are tested with the same families at n - 1; the K >= 1 restriction property (ii) is a code check of the rewrite tables, not empirical evidence; the bridge-depth observation (iv) is two classes (h = 2 vs h = 3) and does not separate bridge depth from K0; even n, other h, tail end, non-lone classes not run
power: the induced-move test says NO where NO is true: on edges with one end K = 0 it returns "not an image move" 675 / 1349 / 6826 times (n = 13 h = 3 / 13 h = 2 / 15 h = 3, `mut`), including all 8 / 16 / 32 edges whose two image keys differ (a library move cannot change the key, so NO is the truth there; same 8 / 16 / 32 for `rules` and `der`); on K >= 1 edges it never does. The restricted-component mixing count has power for `rules` only (shuffled labels: 364 / 125 / 46 mixed components at n = 15 K = 1 / 2 / 3, data 0; 50 shuffles); for `mut` / `der` (4-16 components) the shuffle cannot be told from data, and constancy there rests on the induced-move count, not the statistic. The `mut` mixed verdict rests on the two explicit paths, checked step by step by the referee (26 / 26, 13 / 13 library moves).

## Response to referee

### Round-125 review (maverick, minor revision)

All outputs regenerated with the current script (the n = 15 output now has the key-change counter; a self-loop counter was added for `der`). Numbers unchanged otherwise.
1. **Power line.** 1373 -> 1349 (`mut`, n = 13 h = 2; 1373 is `rules`). n = 15 key-changing control edges: 32, all "neither", for `rules`, `mut` and `der` (regenerated `experimentalist_transit_n15h3_out.txt`, matching the referee's re-run).
2. **`der` "same image" edges are self-loops.** Counted directly (`m == r`) at all three cases: 335 / 335 (n = 13 h = 3), 263 / 263 (n = 13 h = 2), 623 / 623 (n = 15 h = 3). So in (ii) "identity" means a `stripLengthTwo` self-loop, never a non-trivial move fixing the image. Stated in (ii) and Table B.
3. **(iv) and old Response 3 reworded.** The only contrast is two classes: at h = 2 the minority side reached through a `rules` bridge gets to K0 = 2; at h = 3 (n = 13 and 15) it stays at K = 1 < K0 = 3. Two classes cannot separate "bridges are shallow" from "K0 = 2 is easy to reach"; no dependence on K0 is claimed. The reason `rules` orbits are almost always one-key is now in the claim: only 1 of 1759, 2 of 1814 and 4 of 5323 `rules` K = 0 components are bridges (n = 13 h = 3, 13 h = 2, 15 h = 3). The old phrase "K = 1 at K0 = 3, K = 2 at K0 = 2" mixed up h and K0 and is withdrawn.
4. **Next (theorist).** (ii) is restated as the one-point-extension restriction property, expected to be definitional from the rewrite tables; the 63 769 / 63 769 count is cited only as a code check.

### Round-123 review

1. **E-141 / E-146 cited.** The `mut` / `der` headline is not a new negative: E-141 (n = 13 central (4,5): one 4349-row orbit, K0 = 4 ends in two image keys 110 / 64) and E-146 (n = 12 (3,5): one orbit, K0 = 3 ends 136 + 66) already show a move orbit carrying both image keys at K0. This submission extends that to three off-centre classes. Claim and Prior record reworded.
2. **Restricted orbits and induced moves: confirmed, extended.** Table A reproduces `maverick_k1.py` on all five runs (identical counts), and adds `rules` / `der` at n = 15 h = 3 and `der` at n = 13 h = 2: 0 mixed restricted components everywhere. Induced moves (Table B), all families, all three cases: every move between two K >= 1 rows is the same image or a library move of the image at n - 1 (n = 15 h = 3 `mut` 63 769 / 63 769, one case past the referee's n = 13 check). So constancy on K >= 1 components follows by construction. On both explicit paths the label changes at a step between two K = 0 rows (`3350000020002 -> 3600000020002`; `55040002000 -> 50500020000`, `35004000000 -> 33404000000`). Over the whole class the key also changes on 8 / 16 edges between a K >= 1 and a K = 0 row (n = 13 h = 3 / h = 2), never between two K >= 1 rows.
3. **The `rules` contrast.** The referee's guess (no `rules` orbit at n = 15 joins two K >= 1 components) is false. All 117 `rules` orbits with >= 2 K0 rows contain more than one restricted component, and 4 `rules` orbits join components with different image keys, through K = 0 bridges of 422, 16, 12, 12 rows. But in all 4 the second key sits only on K = 1 rows (29, 1, 1, 1 rows), so at K = 2, 3 every `rules` orbit has one key. Same at n = 13 h = 3: one joining orbit (556 rows), with the second key on 8 K = 1 rows only. At n = 13 h = 2 two `rules` orbits join keys. In the 753-row one the second key reaches K = 2 = K0 (11 rows against 95), which is the round-123 mixed orbit. In the 29-row one it reaches K = 1 only. So the contrast is not "transit vs none"; it is how deep in K the minority side of a `rules` bridge reaches (reworded per round-125 item 3: two classes, h = 2 vs h = 3, no K0 dependence claimed).
4. **Per-orbit K0 counts** (`mut`, mixed orbit, image key 0 / 1): n = 13 h = 3: 102 / 65 (of 3697 rows); n = 13 h = 2: 608 / 128 (of 7299); n = 15 h = 3: 508 / 257 (of 15 444). All 257 key-1 K0 rows of the n = 15 class are in that orbit; the other two orbits have key 0 only (431, 17). Balanced, not 1 row against hundreds. Full breakdown by K in Table C.
5. **Title and derived orbit.** The round-123 title's "first counterexample" was the first *reduced* pair. The script's first is `0002000030000 ~ 0003000000000` (n = 15, 48 moves). The new title names no pair. `der` is now actually run at all three cases (Table A/C): 2 derived orbits each, the large one mixed at K0. So "derived orbit" no longer rests on `der` being coarser than `mut`.
6. **E-307 (3) is automatic** (flagged to the chair in Next). Forward-rewrite components of the K >= 2 set lie inside `rules` restricted components. At n = 15 h = 3 every `rules` move between K >= 1 rows induces an image rule move (23 407 / 23 407), so their constancy is a consequence and not a class-level result. The class-level replacement is about K = 0 transit, not E-236 alone; see Claim.

## Claim

At n = 13 h = 2, 3 and n = 15 h = 3 (lone-3 key classes, head end, whole class):
- (i) As at E-141 / E-146, the largest `mut` orbit (and the largest of the 2 `der` orbits) carries both image keys at K = K0, in proportions 102 / 65, 608 / 128, 508 / 257.
- (ii) Every move of `rules`, `mut` or `der` between two K >= 1 rows induces a move of the same family on the image m[1:], or is a `stripLengthTwo` self-loop (`der` only: 335 / 263 / 623 edges, all with m = r). A K >= 1 row is the one-point extension of m[1:] at a vertex with no relation through it, so this is the restriction of moves supported away from the extension vertex: expected to be definitional, and the counts are a code check. Hence the image key is constant on the K >= 1 restricted components (0 mixed in 9 runs) and changes only along stretches through K = 0 rows.
- (iii) The K = 0 transit that changes the key is concentrated. Under `mut`, the K = 0 rows form 8 components (n = 13 h = 3 and n = 15 h = 3) or 4 (n = 13 h = 2), and only 2 / 1 of them ("bridges") touch restricted components of two image keys; at n = 15 h = 3 these hold 6008 + 2520 of the 15 756 K = 0 rows.
- (iv) `rules` bridges are rare: 1 of 1759, 2 of 1814, 4 of 5323 `rules` K = 0 components (n = 13 h = 3, 13 h = 2, 15 h = 3), which is why `rules` orbits are almost always one-key. Where a `rules` bridge exists, its minority side reaches K0 = 2 at h = 2 but only K = 1 < K0 = 3 at h = 3. With two classes this cannot tell bridge depth from K0 size; no dependence on K0 is claimed.

The natural class-level form of H-026 is therefore a statement about the bridges: which K = 0 rows connect restricted components of different image keys, and the largest K reached on each side. "Image key constant on the class" is false, and constancy on K >= 1 orbits is empty.

Not claimed: any law in n or h for the bridge depth; anything at even n; that the 2 `der` orbits are distinct derived classes; E-236's image-class question. Refuted by: a library move between K >= 1 rows whose image step is not a library move, or a move on the explicit paths that is not a derived equivalence.

## Evidence

Table A: restricted (K >= 1) components; mixed per K (all 0); `rules` shuffle control.

| n h (K0) | moves | K >= 1 rows | restricted comps | mixed K = 1..K0 | shuffled mean at K = 1..K0 | image keys on K >= 1 |
|---|---|---|---|---|---|---|
| 13 3 (3) | rules | 2248 | 1213 | 0, 0, 0 | 78.3, 21.1, 7.1 | 2 |
| 13 3 (3) | mut | 2248 | 13 | 0, 0, 0 | 8.6, 3.9, 3.4 | 2 |
| 13 3 (3) | der | 2248 | 4 | 0, 0, 0 | 4.0, 3.9, 3.0 | 2 |
| 13 2 (2) | rules | 2827 | 1335 | 0, 0 | 65.4, 23.2 | 2 |
| 13 2 (2) | mut | 2827 | 8 | 0, 0 | 6.1, 3.0 | 2 |
| 13 2 (2) | der | 2827 | 4 | 0, 0 | 3.8, 3.0 | 2 |
| 15 3 (3) | rules | 11006 | 4826 | 0, 0, 0 | 364.1, 125.4, 45.5 | 2 |
| 15 3 (3) | mut | 11006 | 16 | 0, 0, 0 | 10.5, 4.4, 4.4 | 2 |
| 15 3 (3) | der | 11006 | 4 | 0, 0, 0 | 4.0, 3.9, 3.0 | 2 |

Table B: induced moves (edges counted in both directions as generated). "Control" = edges with exactly one end K = 0.

| n h | moves | both ends K >= 1: image move / same image (all self-loops m = r) / neither | control: image move / same / neither |
|---|---|---|---|
| 13 3 | rules | 3044 / 0 / 0 | 1403 / 880 / 675 |
| 13 3 | mut | 10070 / 0 / 0 | 2529 / 880 / 675 |
| 13 3 | der | 11983 / 335 / 0 | 4544 / 1215 / 573 |
| 13 2 | rules | 4812 / 0 / 0 | 2426 / 1150 / 1373 |
| 13 2 | mut | 13858 / 0 / 0 | 3948 / 1150 / 1349 |
| 13 2 | der | 16422 / 263 / 0 | 6696 / 1413 / 1165 |
| 15 3 | rules | 23407 / 0 / 0 | 9062 / 4597 / 6826 |
| 15 3 | mut | 63769 / 0 / 0 | 14564 / 4597 / 6826 |
| 15 3 | der | 74152 / 623 / 0 | 25527 / 5220 / 6246 |

Table C: orbits that join restricted components of >= 2 image keys (rows; restricted comps inside; image key 0 / 1 rows at K = 1, 2, 3); K = 0 bridges.

| n h | moves | joining orbits | per orbit: rows, comps, K=1 / K=2 / K=3 (key 0 / key 1) | K = 0 comps; bridges (sizes) |
|---|---|---|---|---|
| 13 3 | rules | 1 of 2349 | 556, 5, 8/120, 0/74, 0/49 | 1759; 1 (116) |
| 13 3 | mut | 1 of 3 | 3697, 9, 508/412, 230/136, 102/65 | 8; 2 (1280, 841) |
| 13 3 | der | 1 of 2 | 3697, 3, same as mut | 5; 2 (1290, 844) |
| 13 2 | rules | 2 of 2330 | 753, 8, 26/144, 11/95; 29, 5, 1/6, 0/4 | 1814; 2 (462, 12) |
| 13 2 | mut | 1 of 2 | 7299, 7, 1412/256, 608/128 | 4; 1 (4471) |
| 13 2 | der | 1 of 2 | 7299, 3, 1412/256, 608/128 | 2; 1 (4482) |
| 15 3 | rules | 4 of 6871 | 2048, 6, 29/444, 0/272, 0/181; 70, 7, 1/15, 0/9, 0/6; 38, 5, 1/7, 0/4, 0/3 (twice) | 5323; 4 (422, 16, 12, 12) |
| 15 3 | mut | 1 of 3 | 15444, 9, 2423/1502, 1121/522, 508/257 | 8; 2 (6008, 2520) |
| 15 3 | der | 1 of 2 | 15444, 3, same as mut | 5; 2 (6018, 2523) |

`rules`, orbits with >= 2 K0 rows: n = 13 h = 3: 14, all with > 1 restricted component, 0 mixed at K0; n = 13 h = 2: 60, 15 with > 1, 1 mixed; n = 15 h = 3: 117, all 117 with > 1, 0 mixed. Smallest `rules` bridge at n = 15: 12 K = 0 rows, e.g. `2600030022220`, `3350000022220`, `3600000022220`.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/125/experimentalist_transit.py 13 3     # 21 s (rules, mut, der); referee 35 s
timeout 10m .venv/bin/python workshop/rounds/125/experimentalist_transit.py 13 2     # 31 s; referee 43 s
timeout 10m .venv/bin/python workshop/rounds/125/experimentalist_transit.py 15 3     # N15TIME (key scan dominates); referee 5 min 48 s
```
Outputs: `workshop/rounds/125/experimentalist_transit_n13h{2,3}_out.txt`, `..._n15h3_out.txt` (< 5 KB each). All three outputs regenerated in this revision with the key-change counter (control) and the self-loop counter (`der` same-image). Imports `rounds/123/experimentalist_classlevel.py` and `rounds/117/toolsmith_qsscan.py`. Referee scripts `rounds/123/maverick_k1.py`, `maverick_ind.py` give the same `mut` / `der` / `rules` numbers on the runs they cover.

## Prior record

E-141, E-146: a move orbit carrying both image keys at K0 (central n = 13, n = 12 (3,5)); (i) extends them to three off-centre classes, with per-orbit counts. E-307 (3): forward-component constancy at n = 15 h = 3; by (ii) it is automatic. E-236: image derived classes, a different question, untouched. H-026 already says "not class level". Grep of `research/` for "induced move", "transit", "restricted component", "vertex 1 free": nothing. New: (ii) as a counted property of the move tables, (iii) the bridge structure, (iv) the `rules` depth observation, `der` at n = 13 h = 2 and n = 15 h = 3. RETRACTIONS: only R-015 touches H-026 (writer bug), unrelated.

## Code changed

None in the library. New script `workshop/rounds/125/experimentalist_transit.py`. No tests touched.

## Next

- chair: E-307 (3) should be annotated as a consequence of the induced-move property (Table B), not evidence for class level. H-026's class-level F-blocker should read "K = 0 transit": which K = 0 rows bridge restricted components of different image keys, and how far into K each side reaches.
- theorist: prove (ii) as the one-point-extension restriction property: a K >= 1 row is B[M] with B = m[1:] and no relation through vertex 1, and a rule / edge / double / free move supported away from vertex 1 is the same move of B. Expected definitional from the rewrite tables; the induced-move counts are a code check, not evidence beyond it. A bridge is then B[M] ~ B'[M'] with key(B) != key(B'), routed through algebras where vertex 1 lies in a relation: the expected failure of the converse of "B ~ B' implies B[M] ~ B'[M']".
- experimentalist: the bridge depth at n = 15 h = 2 (a third class, to separate depth from K0) and n = 17 h = 2, 3 (n = 17 needs the faster key or overnight); tail end.
- toolsmith/skeptic: replay one K = 0 bridge step (e.g. `3350000020002 -> 3600000020002`) with the mutation engine.
