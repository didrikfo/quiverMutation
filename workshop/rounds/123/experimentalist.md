# At n = 13 (h = 2, 3) and n = 15 (h = 3) the H-026 image key is NOT constant on a mutation orbit, a derived orbit or the q_S class at K = K0: one orbit carries both image keys (first counterexample 0003000000000 ~ 0003335000000 at n = 15); it is constant only on orbits of the rewrite table alone, and even that fails at n = 13 h = 2

author: experimentalist · round: 123 · kind: negative
thread: S-1 · bears on: H-026, E-307, E-236
scope: n = 13 h = 2, 3 and n = 15 h = 3 only; the whole Coxeter-key class of the lone-3 row (every row of length n with that key, every K: 7459, 5835, 26 762 rows; no row cap); head end only; orbits = undirected union-find closure under the library's move tables (`rules` = lnaMoves.rewritesOf; `mut` = rules + edgeMoves + doubleMutation; `der` = mut + free move), a lower bound for the derived class; q_S = (signature, |det S|, q_S histogram) from rounds/117/toolsmith_qsscan.inv; moves along the counterexample paths not re-checked with the mutation engine; even n, other h, tail end, non-lone classes not run
power: the mixing count does say "mixed" where mixing is true -- shuffled image labels on the same orbits give 35-66 mixed orbits at n = 15 h = 3 rules K = 3 (data 0) and 2-3 at n = 13/15 mut (data 1), 200 shuffles each; it says "constant" where constancy is true -- a planted orbit function as label gives 0 mixed orbits in all 8 runs, and every K > K0 (single image key) gives 0. The counterexample rows' image keys were re-computed with `piecewiseHereditary.removeVertex` + `coxeterTables.lnaCoxeterKey` (agree).

## Claim

The question was: is the image key (Coxeter key of the vertex-1 deletion) of H-026 a function of the class of the source row, not just of its forward-rewrite component (E-307)? At the three cases run, **no, at K = K0**: (i) the key class is one q_S class (q_S does not split it: n = 13 h = 3, h = 2, n = 15 h = 3 each give a single (sig, |det S|, hist)), so at q_S level the class-level question is the key-level split itself; (ii) the key class falls into 3 / 2 / 3 mutation orbits (`mut`; 2 derived orbits at n = 13 h = 3 with `der`), and the largest orbit carries **both** image keys at every K <= K0, also when restricted to reduced rows (no two-arrow relation, so K is not moved by a free move). First counterexample at the E-307 case n = 15 h = 3, K = 3, reduced: the lone row `0003000000000` (image key 1) and `0003335000000` (image key 0), joined by a 26-move path of rule and double-mutation steps (no edge move). (iii) Under the rewrite table alone (`rules`, symmetric, whole class) the image key IS constant at K = 2, 3 for n = 13 h = 3 and n = 15 h = 3 (0 mixed of 601 / 1225 orbits at n = 15, shuffle mean 51 / 159), so E-307's sub-orbit result extends to symmetric rules orbits on the whole class there -- but at n = 13 h = 2 even `rules` orbits mix at K = K0 = 2 (`00300000000` img 1 ~ `00350000000` img 0, 13 rule moves). Consequence: the K-threshold of H-026 is a statement about keys, not a class invariant; the "class level" F-blocker cannot be met in the form "image key constant on the class". Not claimed: anything at even n, at other h, about whether the 3 mutation orbits are one derived class (q_S cannot separate them), or that the image *derived classes* are distinct (E-236's question). Refuted by: a move on the paths below that is not a derived equivalence (then the orbits are wrong).

## Evidence

Per run: key-class rows; orbits; at each K <= K0, orbits containing K-rows / orbits with >= 2 image keys; reduced rows at K0.

| n h (K0) | moves | class rows | orbits | K=K0-1: hit/mixed | K=K0: hit/mixed | K0, reduced only: rows, mixed | K > K0 mixed | q_S classes | time |
|---|---|---|---|---|---|---|---|---|---|
| 13 3 (3) | rules | 5835 | 2349 | 301 / 0 | 129 / 0 | 50, 0 | 0 | 1 | 80 s |
| 13 3 (3) | mut | 5835 | 3 (3697, 2115, 23) | 2 / 1 | 3 / 1 | 50, 1 | 0 | 1 | 11 s |
| 13 3 (3) | der | 5835 | 2 (3697, 2138) | 2 / 1 | 2 / 1 | 50, 1 | 0 | 1 | 11 s |
| 13 2 (2) | rules | 7459 | 2330 | 745 / 2 | 356 / 1 | 64, 1 | 0 | 1 | 78 s |
| 13 2 (2) | mut | 7459 | 2 (7299, 160) | 2 / 1 | 2 / 1 | 64, 1 | 0 | 1 | 19 s |
| 15 3 (3) | rules | 26762 | 6871 | 1225 / 0 | 601 / 0 | 86, 0 | 0 | 1 | 9 min 17 s |
| 15 3 (3) | mut | 26762 | 3 (15444, 11154, 164) | 3 / 1 | 3 / 1 | 86, 1 | 0 | 1 | 2 min 23 s |

Image keys by K on the whole n = 15 h = 3 class (label: rows): K = 1 {0: 4907, 1: 1502}, K = 2 {0: 2183, 1: 522}, K = 3 {0: 956, 1: 257}, K = 4..8 {0} only -- the key-level K0 = 3 of E-304 reproduces on the full class (E-307 used K >= 2 only). Whole orbits (any K) always mix, since K = 0 rows carry up to 7 image keys.

Counterexamples (first found by the script: smallest mixed orbit, lexicographically first rows):
- n = 15 h = 3, K = 3, mut: `0002000030000` (i0) ~ `0003000000000` (i1), 48 moves; reduced pair `0003000000000` (i1) ~ `0003335000000` (i0), 26 moves; the label changes at the rule step `3350000020002 -> 3600000020002` (a K = 0 row), the rest of the K <= 1 stretch is 15 double mutations at a repeated vertex.
- n = 13 h = 3, K = 3: `00020030000` ~ `00030000000` (mut 35 moves, der 22); reduced `00030000000` ~ `00033350000`.
- n = 13 h = 2, K = 2: `00200003000` ~ `00300000000` (mut 21); rules only: `00300000000` ~ `00350000000` (13 rule moves, via K = 0 rows).
Image keys (n - 1 = 12): `0030000000` -> (1,1,0,-1,-2,-2,-2,-2,-2,-1,0,1,1); `0033350000`, `0020030000`, `0000030000` -> (1,1,0,-1,-2,-3,-3,-3,-2,-1,0,1,1) (library recomputation).

Controls (200 shuffles of the image labels over the same K-rows and orbits; planted label = orbit index mod 2):

| run | K | real mixed | shuffled mean (min-max) | planted |
|---|---|---|---|---|
| 15 3 rules | 3 / 2 | 0 / 0 | 51.3 (35-66) / 159.3 (136-179) | 0 / 0 |
| 15 3 mut | 3 / 2 | 1 / 1 | 3.0 (2-3) / 2.6 (2-3) | 0 / 0 |
| 13 3 rules | 3 / 2 | 0 / 0 | 8.3 (5-12) / 21.8 (13-29) | 0 / 0 |
| 13 2 rules | 2 / 1 | 1 / 2 | 23.8 (14-33) / 70.3 (55-86) | 0 / 0 |

With 2-3 mut orbits the shuffle can hardly be told from the data; the mut verdict rests on the explicit paths, not on the statistic. q_S was checked on 3 extra members of each orbit: 0 disagreements (a YES-side consistency check: moves preserve q_S).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/123/experimentalist_classlevel.py plan 13 3 15 3      # 15 s: 208 012 / 2 674 440 rows, class scan 7 s / 165 s on 4 cores
timeout 10m .venv/bin/python workshop/rounds/123/experimentalist_classlevel.py run 13 3 mut        # 11 s (also: rules 80 s, der 11 s; run 13 2 mut 19 s, rules 78 s)
timeout 10m .venv/bin/python workshop/rounds/123/experimentalist_classlevel.py run 15 3 mut        # 2 min 23 s
timeout 10m .venv/bin/python workshop/rounds/123/experimentalist_classlevel.py run 15 3 rules      # 9 min 17 s (the 400 shuffles dominate; close to the limit)
```
Outputs: `workshop/rounds/123/experimentalist_classlevel_n13h{2,3}_{mut,rules,der}_out.txt`, `..._n15h3_{mut,rules}_out.txt` (each < 6 KB). Uses 4 cores for the key scan.

## Prior record

E-307 (3): image key constant on forward-rewrite components of the K >= 2 set at n = 15 h = 3, "not class level"; this extends it to symmetric `rules` orbits on the whole class (same case) and shows it fails one move family wider and at n = 13 h = 2. E-236 (round 083) asked a different class-level question (how many `derivedOrbits(n-1)` orbits the *images* fall into) and left it undecided; nothing in `research/` tests the source-orbit dependence. H-026 already says "not class level". Grep of RETRACTIONS for image key / H-026: only R-015 (writer bug), unrelated. Not a rediscovery as far as grep shows.

## Code changed

None in the library. New script `workshop/rounds/123/experimentalist_classlevel.py` (imports `workshop/rounds/117/toolsmith_qsscan.inv`). No tests touched.

## Next

- chair/theorist: reword H-026's open item: "class level" in the sense "image key constant on the class" is refuted at n = 13, 15; a class-level form must be different, e.g. "the image *set* of keys at K = K0 is the same for every orbit" or a statement about which derived class the image lands in (E-236's question).
- theorist: why do rules-only orbits keep the image key at n = 15 h = 3 but not n = 13 h = 2? The n = 13 h = 2 rule path goes through K = 0 rows `55040002000 -> 50500020000`; is h = 2 = K0 the boundary where a head rule can reach the deleted vertex?
- toolsmith/skeptic: re-verify the 26-move n = 15 path (and the 13-move n = 13 h = 2 rule path) with the mutation engine, step by step, since the counterexample rests on the move tables.
- experimentalist: same test at n = 17 h = 5 (class scan needs a faster key or overnight), at even n non-central classes, and with the K counted on the reduced row.
