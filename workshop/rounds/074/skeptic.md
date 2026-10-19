# n = 7, key classes 1 and 2, parent depth <= 8: the 70 key-kept failing edges have 60 distinct canonicalKeys but only 12 distinct children up to relabelling, and 4 of the 12 contain no E-151 window child (also none is the opposite of any of the 12)

author: skeptic · round: 074 · kind: result
thread: T10 · bears on: H-015 (key guard), E-204, E-157, E-163, E-212
scope: n = 7, key classes 1 and 2 only; parent depth <= 8 (complete, no shard stopped by its budget); children compared as presentations (quiver with relations), with and without taking opposite algebras, not as derived-equivalence classes; "joined" = joined under the J = 0 premise (E-157/E-160/E-163); the E-151 window is reconstructed by position; class 1 was NOT re-run through `toolsmith_collect.py` (exceeds 10 min), so the class-1 window membership rests on positional reconstruction plus agreement of counts (16, 10 keys) with E-157.

## Response to referee

1. Op/dual (required 1): done, `skeptic_ctrl.py` (output `skeptic_ctrl_out.txt`). For each of the 4 new classes I took the opposite presentation (arrows reversed, paths reversed, same coefficients) and tested it against the representatives of all 12 classes with the structure-graph isomorphism (sign-blind, then sign-aware where blind matched): 0 matches, so none of the 4 is the op of a window class (nor of another new class, nor self-dual). Sanity: op(op(X)) = X for all 12 reps. As a side fact no window rep is the op of another window rep either. The 4 / 14 counts stand; the claim now says so.
2. Iso-test controls (required 2): done, same script. Positive: 22 children (up to 3 per window class, 8 classes) each given 3 random relabellings (vertex permutation, permutation of parallel arrow keys, random +-1 sign gauge on arrows, shuffled relation and term order): 66 of 66 recognised in their own class by the blind test and by the sign-aware test; and 0 matches of any relabelled class rep against the other 11 classes. Negative: the 13 pairs of classes with equal (arrow, relation) counts are all separated by exact isomorphism (0 isomorphic, 0 even with equal WL hash). A weaker control I tried, flipping one term sign of a multi-term relation, separated 0 of 12 reps under the sign-aware test: that is not evidence either way (arrow rescaling can absorb such a flip), so I do not use it. Limit: the control exercises the test on presentations as stored, not the BFS pipeline, whose merging of 31 copies is the other evidence.
3. Title: now carries "n = 7, key classes 1 and 2, parent depth <= 8"; and the wording says "contain no E-151 window child (as isomorphism classes of presentations)".
4. Class 1 through `toolsmith_collect.py`: not done this sitting (exceeds 10 min, killed twice). Stated in the Scope line: class 1 window membership rests on positional reconstruction.
5. Gap 4 (loose wording): "over-counts by a factor of 5" is the average 60 / 12, not per class; class sizes (16, 10, 10, 6, 6, 6, 4, 4, 2, 2, 2, 2 edges) are as in the table. Gap 3 and the BFS rebuild: unchanged, the reviewer trusted the pickles.

## Claim

E-204's 70 key-kept `tiltingPlus`-false edges (34 class 1, 36 class 2; rebuilt here, same 34 and 36 as E-204) give 70 children with 60 distinct `canonicalKey`s (24 in class 1, 36 in class 2, none `None`) but only **12 distinct children when vertex and arrow labels are forgotten** (7 in class 1, 5 in class 2, none shared between classes). Every one of the 12 classes holds several labelled keys, so `canonicalKey` over-counts by a factor of 5; the labelled count (60) is not a count of different algebras. The 25 children of E-151/E-154 (the ones E-157, E-160, E-163 joined; 19 labelled keys, 10 + 9) fall into 8 of the 12 classes. **4 classes (14 of the 70 edges: 4 + 2 + 6 + 2, all with parent depth 8) contain no E-151 child**, and none of the 4 is the opposite algebra of any of the 12 classes (so none is joined by duality either), so no recorded join covers them. Under the labelled key it reads 41 of 60 keys outside the E-151 keys, but 27 of those (14 belong to the new classes, 14 distinct keys) are relabelled copies of joined classes; the iso count is the honest one. Not claimed: that the 4 classes are, or are not, derived equivalent to their parents or to an LNA (no join was searched); "isomorphic" below is isomorphism of presentations, which makes the joined classes' joins transfer by relabelling but proves nothing about the 4. Refuted if a structure-preserving isomorphism exists between a new class and a window child that my test missed (my test would only err toward splitting if it had a bug, see Evidence).

## Evidence

1. Rebuild: `skeptic_failkids.py` repeats E-204's walk (gate, `tiltingPlus`, illegal-relation drop, reduce, key guard, `canonicalKey` dedupe; logic copied from `rounds/070/toolsmith_cyccount.py`) to parent depth 8, stores every failing edge with parent and child. Classes 1, 2: 16 404 + 3 x 4 370 and 22 195 + 3 x 5 411 expansions, matching E-204's totals (34 + 36 failing edges); 34 + 36 = 70 failing edges, distinct as (parent key, v), no parent key `None`, no child key `None`. The depth-8 frontier is checkpointed once per class and four shards expand its quarters (each under 10 min).
2. Iso test (`skeptic_failkids_an.py`): a graph with vertex nodes, one node per arrow (parallel arrows kept distinct), relation nodes, term nodes labelled by |coefficient| and length, position-labelled edges to the arrows; networkx isomorphism after a Weisfeiler-Lehman bucket. This forgets vertex labels and arrow keys. Signs are not in the graph, so I then checked them: within each class every member has a graph isomorphism whose relabelled relations have the same sign-gauge-canonical form (`fingerprint._relationsUnder` with an explicit arrow renaming) as the class's first member: 0 failures, 0 undecided (cap 3000 isomorphisms per pair). So the 12 classes are classes of presentations including relation signs up to arrow rescaling. Parallel-arrow children: 12 of 70 (arrow counts 7 in 58, 8 in 12).
3. Window: the E-151/E-154 set is the first 20 000 expansions of the same BFS. A prefix parent (depth < 8) is always inside (at most 12 033 / 16 784 expansions); a depth-8 parent is inside iff prefix + its frontier position < 20 000. This gives **16 class-1 and 9 class-2 edges, exactly E-157's 16 and 9, each distinct (parent key, v)**, 10 + 9 = 19 distinct child keys as E-157 prints ("10 + 9 distinct canonical keys"). Class 2 cross-check: `rounds/049/toolsmith_collect.py 7 2 20000` rerun, 9 failing records, its 9 child keys equal my window keys. Class 1 through `toolsmith_collect.py` did not finish in 10 min (killed twice, once alone), so class 1 rests on the positional reconstruction plus the count and key agreement above.
4. Overlap table (iso classes; size = failing edges):

| class (cls, edges) | in the E-151 window | child arrows / relations (rep.) |
|---|---|---|
| 8 classes | yes (25 window edges + 31 more edges that are relabelled copies, 56 in all) | various |
| cls 1, 4 edges | NO | 7 / 5 |
| cls 2, 2 edges | NO | 7 / 4 |
| cls 2, 6 edges | NO | 7 / 3 |
| cls 2, 2 edges | NO | 7 / 5 |

   Class sizes: 16, 10, 10, 6, 6, 6, 4, 4, 2, 2, 2, 2. Parents are not compared (the question was children); the 4 new classes have parents only at depth 8.
5. **Hand-checked child** (class 2, new class of 6 edges, v = 6, parent depth 8): parent arrows 1->7, 3->1, 3->4, 4->5, 5->6, 6->2, 7->6 with the single relation (3,4,5,6,2) - (3,1,7,6,2) (signs as stored). Mutating at 6 (one outgoing arrow 6->2): the arrows 5->6, 7->6 become composites 5->2, 7->2, the arrow 6->2 reverses to 2->6, and the relation becomes (3,4,5,2) - (3,1,7,2) = 0 plus the two new zero relations (5,2,6) = (7,2,6) = 0. This is the printed child (7 arrows, 3 relations), 25 -> 19 in dimension. I checked it by my own path count (`skeptic_failkids_hand.py`, using only `ap.idealBasis` through `rounds/050/skeptic_tilt.py`): parent Cartan row 3 = [1 1 1 1 1 2 1], child [1 1 1 1 1 0 1]; both Coxeter polynomials (computed here from -C^{-T}C) are (x-1)^2 (x+1) (x^2+x+1)^2, so the key is kept. Why `tiltingPlus` is false, by hand: with T_6 = (P_6 -> P_2) and the other summands P_j, the only candidate for Hom(T, T[-1]) is a degree-0 map P_3 -> P_6, i.e. a combination of the two paths 3 ~> 6 ((3,4,5,6) and (3,1,7,6)) which the arrow 6->2 sends to zero. The two paths times 6->2 are exactly the two terms of the relation, so the combination with the relation's signs is killed: one nonzero chain map, not null-homotopic (T_6[-1] has nothing in degree -1). The program agrees: Hom(T,T[-1]) is 1-dimensional, only at (i, j) = (3, 6), Hom(T,T[+1]) = 0, and the m = 0 table equals the child's Cartan matrix (End(T) has the child's dimension data). So this step is a key-kept J != 0 step with a 1-dimensional obstruction, as in E-174/E-181; whether this child is derived equivalent to its parent stays open (the T-based argument cannot decide it).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/074/skeptic_failkids.py 7 C 9 590 0/4 DIR/cC_s0.pkl   # C = 1, 2; writes DIR/cC_s0.pkl.frontier (about 320 s class 1, 390 s class 2)
timeout 10m .venv/bin/python -W ignore workshop/rounds/074/skeptic_failkids.py 7 C 9 590 S/4 DIR/cC_sS.pkl   # S = 1, 2, 3 after the frontier exists (100 to 200 s each; run 6 together on 4 cores)
timeout 10m .venv/bin/python -W ignore workshop/rounds/049/toolsmith_collect.py 7 2 20000 DIR/e151_c2.pkl 100   # class 2 only; class 1 exceeds 10 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/074/skeptic_failkids_an.py DIR   # about 2 min; output skeptic_failkids_an_out.txt
timeout 10m .venv/bin/python -W ignore workshop/rounds/074/skeptic_ctrl.py DIR   # op comparison and iso controls, under a minute; needs an_summary.pkl; output skeptic_ctrl_out.txt
.venv/bin/python -W ignore workshop/rounds/074/skeptic_failkids_hand.py DIR    # seconds; output skeptic_failkids_hand_out.txt
```
DIR = /tmp/claude-0/sk74 here (pickles not committed; the checkpoint pickles are tens of MB).

## Prior record

E-204 counted 70 edges and states they are not distinct children; E-192 counted 25 at a cap; E-157/E-160/E-163 joined the 25 (under J = 0), E-212 class 2 8 of 9 with isomorphism at the meets; STATE asks for "distinct canonicalKey count of the 70 failing children and the J = 0 join of those not in E-157". Not in `research/`: the count 60 / 12, the 8 + 4 split, the positional reproduction of E-157's 16 / 9 inside the full ball (confirms E-201: 16 = parent depth 7 + 8 splits 4 + 12 as in E-151). Consistent with E-201's remark that counts depend on the cap: the full ball gives 70 edges where the 20 000-expansion window gives 25, and the 45 extra edges add only 4 new presentations (14 edges; the other 31 are copies of window classes).

## Code changed

None in `quivermutation/`. New scripts in `workshop/rounds/074/`: `skeptic_failkids.py`, `skeptic_failkids_an.py`, `skeptic_failkids_hand.py`, `skeptic_ctrl.py` (+ the `_out.txt` files). No tests touched.

## Next

- toolsmith: the J = 0 join (F/R chains with the E-212 Hom test and isomorphism at the meets) of the 4 new classes and of their parents (depth 8); class 1 first (rep.: child arrows 2->4, 3->2, 4->5, 4->6, 1->7, 5->1, 6->1, v = 7 of the parent printed in `skeptic_failkids_an_out.txt`). A miss is a bound, not an exit.
- The 4 classes are not duals (op) of any of the 12 classes (round 074 response, `skeptic_ctrl.py`).
- The same count at n = 8 is unavailable (E-192: 0 key-kept J != 0 edges there).
