# A loop of one Oppermann mutation of a graded quiver at i exists iff the quiver has one at i; for parents with relation arrows degree >= 1 loops appear at other vertices after one step, and the 25 failing n = 7 parents (acyclic with relations, checked) have no loop at i of any degree

author: theorist · round: 090 (revision of 089) · kind: proof (one-step lemma) + negative (the 089 headline was too strong)
thread: T10/T5 · bears on: E-189, E-235, E-240, E-248, H-015
scope: quiver level only (arrows, endpoints, degrees; no differential). Lemma: all graded parents. Machine: 2 574 random ungraded parents (089), plus 2 000 random acyclic parents with degree-1 (and some degree-2) relation arrows n = 3..6, 9 020 one-step mutations; iterated depth <= 3 on 100 such parents (5 707 nodes). The 25 failing parents of class 1 (16) and class 2 (9) at n = 7, rebuilt by collection (E-151 walk, 20 000 expansions). v1 of arXiv:1504.02617 only. Not covered: whether the iterated loops survive in a minimal model; the cancellation proof with such a loop at i.

## Response to referee

1. Retitle, rewrite (2)-(4) to the narrowed scope. Done. "Ungraded implies no loop of degree >= 1" is withdrawn. It holds for Q with all arrows in degree 0 and no relation arrows (a hereditary algebra, or a quiver used as is), not for a parent given as kQ0/(R); Construction `const.usual_algebra` (l. 291) and the A_3/rad^2 example (l. 296-307) add an arrow omega of degree 1 for each minimal relation. Corrected consequences below.
2. Counterexample and relation-arrow variant of Parts A/B. Done. Referee's example reproduced: `0->1, 1->2, 0->2, omega: 0->2 (deg 1)`, mutate at 0 gives a degree-1 loop at 2, type (c). Without the parallel arrow c (A_3/rad^2 itself) no mutation at any vertex gives a loop. Parts A/B rerun with relation arrows: below.
3. The 25 parents are acyclic with relations. Done, with a named command; below. A relation arrow omega runs from the first to the last vertex of the relation, so it is a loop at i only for a relation that is an oriented cycle.
4. Title parenthetical and E-216. Done. The garbled clause is removed. E-216's `CyclicDescendantsError` refuses on the Gabriel quiver of H_0 (a cyclic parent is rejected), not on the dg quiver; Thm 1.1's hypothesis "no degree-0 loop at i" is on the dg quiver and is a different condition. 089 conflated them; they coincide in effect only for acyclic parents.
5. Chair's sub-bullet: may be closed for loops at i, for the 25 parents only, on the checks below. It is not closed for other vertices or for iterated mutation.

## Claim

Let Q be a graded quiver as in Oppermann Sect. 2, i a vertex with no degree-0 loop, A the degree-0 arrows out of i, M the mutated quiver (Steps 1-4, `mutation.tex` l. 320-375). A loop of M at v is exactly one of: (a) a loop of Q at v; (b) alpha.phi, phi: v -> i, alpha: i -> v in A, degree |phi| - [v=i]; (c) phi.alpha^-1, phi: i -> v not in A (a degree-0 arrow out of i is in A, so phi has degree >= 1 or is a loop at i), alpha in A, degree |phi|; (d) alpha.phi.beta^-1, phi a loop at i of degree >= 1, degree |phi|. Consequences, with scope:
(1) For every Q: M has a loop at i iff Q has one at i, of the same degree.
(2) If every arrow of Q has degree 0, loops of M are (a), (b) only, all of degree 0 for (b): no loop of degree >= 1 after one step. This covers hereditary parents; it does NOT cover a parent given by relations, where omega arrows have degree 1.
(3) If Q is acyclic as a quiver with arrows of any degree (every arrow, omega included, runs between distinct vertices), (a), (b), (d) are empty and the loops of M are exactly the (c) loops: pairs (phi out of i of degree >= 1, alpha in A) with the same target. Such a loop sits at a vertex v != i, and has degree |phi| >= 1. So an acyclic parent with relation arrows can have M with loops of degree >= 1 after one step, never at i.
(4) The 25 failing parents (n = 7, classes 1, 2) are acyclic with relations (command below). Hence by (1) none has a loop at i of any degree, E-235's hypothesis holds for them, and the failure of the 25 children is not a cancellation failure.
What this does not claim: that the (c) loops at other vertices are harmless for later steps; the repo passes to H_0 each step so (see below) they do not feed the next mutation, but this is the repo procedure, not a theorem here. Refutation: a loop of M at i for a Q with none at i, or a failing parent with a cyclic relation.

## Evidence

Proof of the classification: as in 089 (arrows of M are Q_1 \ A, alpha^*, and the composites; a loop needs equal endpoints; alpha is never a loop). Degrees from Steps 1-4 (referee checked l. 326-358). For (3): acyclic Q has no 2-cycle (kills (b)) and no loop at i (kills (d)); (c) needs phi, alpha parallel i -> v. Weakest step: "acyclic" for Q with relation arrows of higher degree means every arrow of any degree joins distinct vertices; I use that Ext^k(S_s, S_t) != 0 needs a path s -> t in the Gabriel quiver of an acyclic algebra, so Ext-arrows of any degree are non-loops. This is the standard reading, not checked in the source.

Part A with relation arrows (`theorist_relarrows.py rel`): 2 000 random acyclic parents (degree-0 arrows incl. parallel; degree-1 omega s->t where a degree-0 path of length >= 2 joins s, t; some degree-2 arrows), all vertices i: 9 020 steps; loops at i: 0; steps with a loop anywhere: 681 (1 330 loops; degree 1: 1 018, degree 2: 312); loop set equals the predicted (c) set in all 9 020 steps (0 mismatches).
Part B (dg mutation iterated, no H_0, no consecutive repeat), 100 acyclic parents with relation arrows:

| depth <= | nodes | nodes with a loop | by depth | blocked by degree-0 loop at i | legal steps at a vertex with a loop of degree >= 1 |
|---|---|---|---|---|---|
| 2 | 1 706 | 196 | 1: 25, 2: 171 | 0 | 27 |
| 3 | 5 707 | 1 048 | 1: 25, 2: 171, 3: 852 | 27 | 187 |

So with relation arrows the "positive-degree loop at i" case of E-235 already occurs at depth 2 of the dg mutation, one step earlier than for hereditary starts (089: depth 3). It is reached only if one mutates dg quivers without passing to H_0.

The 25 parents: `.venv/bin/python -W ignore workshop/rounds/090/theorist_relarrows.py parents /tmp/tsm/c1.pkl 1` (and `c2.pkl 2`). Pickles built by the E-151 walk collector (`workshop/rounds/050/toolsmith_collect.py 7 <cls> 20000 /tmp/tsm/c<cls>.pkl 100`; 442 s and 335 s, 16 and 9 failing steps, total 25). Result: class 1, 16 parents: 0 with an oriented cycle or loop in the Gabriel quiver, 57 relations, 0 with equal endpoints; class 2, 9 parents: 0, 21 relations, 0 equal endpoints. All 25 parents have relations (none hereditary), so the omega arrows are present in the dg model and all join distinct vertices.

## Reproduction

```
.venv/bin/python workshop/rounds/090/theorist_relarrows.py ex                          # counterexample, <1 s
.venv/bin/python workshop/rounds/090/theorist_relarrows.py rel 2000 3                  # Parts A/B, <1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/tsm/c1.pkl 100   # 442 s; cls 2: 335 s (mkdir -p /tmp/tsm)
.venv/bin/python -W ignore workshop/rounds/090/theorist_relarrows.py parents /tmp/tsm/c1.pkl 1                   # seconds; cls 2 likewise
```
Outputs: `workshop/rounds/090/theorist_relarrows_out.txt`, `theorist_relarrows_parents_out.txt`. The script imports `workshop/rounds/089/theorist_loopreach.py`.

## Prior record

E-248 records the 089 submission and the referee's counterexample. E-198 / E-168: the printed-path parents are acyclic (Gabriel quiver; E-168 via `theorist_loops.py`, loops and cycles on edges). New here: the check on the 25 failing parents themselves, relation endpoints included, and the relation-arrow Parts A/B. E-189 (hypothesis "no degree-0 loop"), E-235/E-240 (cancellation, no loop at i), E-216. No RETRACTIONS entry found.

## Code changed

None in the repository. New: `workshop/rounds/090/theorist_relarrows.py`.

## Next

- Skeptic: confirm that Ext-arrows of any degree of an acyclic algebra are non-loops, and that omega endpoints are first and last vertex of the relation.
- Theorist: the (c) loops at v != i after one step: does a minimal model of kM have them essential, and does H_0 of kM still equal End(T) (E-199)? Not recommended before the NO-capable test (T10).
