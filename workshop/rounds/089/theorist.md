# Loops of degree >= 1 at i are unreachable in the repo's single-step setting: one Oppermann mutation of an ungraded parent yields a quiver with no loop of degree >= 1 at all, so E-235's open case never occurs for the 25 failing children (they are reachable only by iterating the dg mutation)

author: theorist · round: 089 · kind: proof (one-step lemma) + negative (open case is vacuous for the repo)
thread: T10/T5 · bears on: E-189, E-235, E-240, E-157, H-015
scope: quiver level only (arrows, sources, targets, degrees; the differential is not tracked). Lemma: all parents. Machine: 2 574 random one-step parents n = 2..5 (cycles, parallel arrows, degrees 0..2), and iterated mutation from random acyclic ungraded parents n = 2..5, depth <= 3 (21 915 nodes) and <= 4 (400 parents, 16 698 nodes). v1 of arXiv:1504.02617 only. Not covered: whether the iterated-case loops survive in a minimal model; the cancellation proof when such a loop is present.

## Claim

Let Q be a graded quiver as in Oppermann Sect. 2 and i a vertex without degree-0 loop (Thm 1.1's only hypothesis, E-189), A the degree-0 arrows out of i, M the mutated quiver (Steps 1-4, `mutation.tex` l. 320-375). A loop of M at a vertex v is exactly one of:
(a) a loop of Q at v (degree |phi| - [v=i] + [v=i] = |phi|, unchanged);
(b) alpha.phi with phi: v -> i, alpha: i -> v in A (a 2-cycle of Q through i and v), degree |phi| - [v=i] (Step 2; for v != i this is degree |phi| and = 0 when phi has degree 0);
(c) phi.alpha^{-1} with phi: i -> v in Q_1 \ A and alpha: i -> v in A (parallel arrows i -> v, phi of degree >= 1 because it is not in A), degree |phi|;
(d) alpha.phi.beta^{-1} (Step 4) with phi a loop of Q at i (degree >= 1), alpha, beta: i -> v, degree |phi|.
Consequences. (1) M has a loop at i iff Q has one at i, of the same degree: (b), (c), (d) with v = i would need alpha: i -> i, excluded since A has no loops. (2) If Q is concentrated in degree 0 (every ordinary algebra, acyclic or not), (c) and (d) are empty and (b) gives degree-0 loops only: M has no loop of degree >= 1 anywhere. (3) If Q is acyclic and ungraded, M has no loop at all. (4) The repo takes H_0(kM) = End(T) as the next parent (quiver of an ordinary algebra, all arrows degree 0), and refuses to mutate at a vertex with a loop (`CyclicDescendantsError`, E-216). So the parent of every repo step has all arrows in degree 0 and the premise of Thm 1.1 and of E-235 is "no loop at i", which `noCycles` gives. The case "loop of degree >= 1 at i" cannot occur for any of the 25 failing n = 7 children or their 25 parents (all acyclic, E-168 / E-198). This answers the open item behind the J = 0 premise: Thm 1.1 needs only "no degree-0 loop" and the repo meets the stronger condition; cancellation (E-235, E-240) applies, and the failure of the 25 children is not a cancellation failure. It remains H_1 != 0, i.e. J != 0, Hom(T,T[-1]) != 0 (E-189, E-161), unchanged.

What this does not claim: that positive-degree loops at i are harmless (E-235's control shows d_1 stops being monomial there, 6 968 of 10 206 loop letters), or that they cannot occur if one iterates the dg mutation without passing to H_0 (below).

## Evidence

Proof of the classification: M's arrows are the Q arrows not in A (re-graded), alpha^* (source target(alpha) != i, target i), and the composite arrows of Steps 2-4. A loop needs equal source and target. Q-arrows: loops of Q (a). alpha^*: source j != i since alpha is not a loop, so never a loop. alpha.phi: source = source(phi), target = target(alpha); equal iff phi: v -> i, alpha: i -> v, case (b) (v = i would make alpha a loop). phi.alpha^{-1}: source = target(alpha), target = target(phi); phi leaves i, so loop iff alpha, phi parallel i -> v (or phi a loop at i and target(alpha) = i, impossible); phi not in A and not a degree-0 non-loop, hence degree >= 1 or phi is a loop at i, covered by (d). Step 4 likewise. Degrees from Step 1-4 formulas. Weakest step: using the printed "Q_1 \ A" to mean arrows of positive degree; a degree-0 loop at i is excluded by hypothesis, a degree-0 non-loop out of i is in A by definition (l. 318).

Machine (`theorist_loopreach.py`, 3 s): Part A, 2 574 steps from random parents with cycles, parallel arrows and degrees 0..2: loops of M by cause (a) 1 434, (b) 1 223, (c) 146, (d) 138; 0 loops outside the predicted list (prediction compared per loop). Ungraded case: only (b), degree 0.

Part B, iterated dg mutation (quiver level, any vertex at every step, no repetition filter), from random ACYCLIC ungraded parents:

| depth | nodes with a loop | loop degrees seen |
|---|---|---|
| 0, 1, 2 | 0 | none |
| 3 | 676 of 21 915 (depth <= 3 run) | all degree 1 |
| 4 | 1 404 of 16 698 (depth <= 4 run) | degree 0 and 1 |

In no node was the next mutation blocked by a degree-0 loop at i in the depth <= 3 run, so every depth-3 node with a loop has legal next steps at the looped vertex: the "degree >= 1 loop at i" case of E-235 does occur, but only after three dg mutations. Smallest example (exhaustive search, n <= 4, arrows <= 4, no consecutive repeat): Q = 1 <- 0 -> 2 (A_3, ungraded), mutate at 1, 0, 2, 0; the last step creates a degree-1 loop at vertex 1 of type (c) (the composite phi.alpha^{-1}, phi of degree 1 coming from Step 1 at an earlier step). Whether this loop is essential in a minimal model I do not know (differential not tracked).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/089/theorist_loopreach.py 3000 3     # 3 s
timeout 10m .venv/bin/python -c "import sys; sys.path.insert(0,'workshop/rounds/089'); import theorist_loopreach as t; print(t.minimal_example(4, True))"
```

## Prior record

E-189 (Thm 1.1 hypothesis: no degree-0 loop), E-235 / E-240 (cancellation proved for no loop at i of any degree; positive-degree loops listed as not covered), STATE T5 and Agenda item 1 ("loops of degree >= 1 at i"), E-216 (CyclicDescendantsError), E-168 / E-198 (printed-path parents acyclic). Not in `research/` before: the classification (a)-(d) of loops of M, and the reachability statement. grep of RETRACTIONS found nothing on this. The classification is elementary and a reader of Steps 1-4 may consider it known; I found no record of it.

## Code changed

None in the repository. New file `workshop/rounds/089/theorist_loopreach.py` (stand-alone; no tests needed).

## Next

- Scholar/skeptic: check the classification against the LaTeX Steps 1-4 independently (especially the degree of phi.alpha^{-1} and whether Q_1 \ A means positive degree), and whether H_0 of a minimal model of the iterated dg mutation can contain the depth-3 loop.
- Theorist: only if someone mutates dg quivers without passing to H_0 (not the repo's procedure) does the d_1-non-monomial case matter; then the A_3 example 1 <- 0 -> 2 with sequence 1, 0, 2, 0 needs a differential to be tracked. Not recommended before the NO-capable test (T10).
- Chair: item 1 sub-bullet "loops of degree >= 1 at i" can be closed as not reachable by the repo's walks.
