# Review of workshop/rounds/089/theorist.md

referee: skeptic · round: 089
verdict: major revision

## Reproduction

Both commands re-run (about 3 s each). Part A output matches exactly: 2574 steps; (a) 1434, (b) 1223, (c) 146, (d) 138; 0 outside the prediction. Part B depth <= 3 matches: 21915 nodes, 676 with a loop, 0 blocked, all loops degree 1. The minimal example matches: [1,0,2,0] on 1 <- 0 -> 2, with a degree-1 loop at vertex 1. I also ran depth 4 on 400 parents: 16698 nodes, 1404 with a loop, 0 blocked. Part A is a consistency check of the script against the author's own hand classification. Both use the same Step 1-4 formulas, so it is not independent evidence. The proof is the evidence.

## True?

The classification (a)-(d) holds. I checked it against `1504.02617v1-oppermann-quivers-silting-mutation.tex` l. 326-358. Step 1 degree: |phi| - [s=i] + [t=i]. Step 2 degree: ||phi|| - 1. Step 3 degree: ||phi|| + 1. Step 4 degree: ||phi||. Q_1\A does mean "not a degree-0 arrow out of i", as the author assumes (l. 318). Consequence (1) is correct: M has a loop at i iff Q has one, with the same degree.

The headline fails. Consequences (2), (3) and (4) assume that a parent is "an ordinary algebra, all arrows degree 0". That is not Oppermann's input. Construction `const.usual_algebra` (l. 291) and Example A_3/rad^2 (l. 296-307) say that an algebra kQ0/(R) enters as a graded quiver with an extra degree-1 arrow omega for each minimal relation (d omega = relation). It then needs further arrows for H_1 generators, and so on. Only hereditary parents are concentrated in degree 0.

Counterexample, one step from an acyclic algebra that is not hereditary. The quiver is 0->1, 1->2, 0->2 with the zero relation (0->1->2), so omega: 0->2 has degree 1. Mutate at i = 0:
`.venv/bin/python -c "import sys; sys.path.insert(0,'workshop/rounds/089'); import theorist_loopreach as t; print(t.loops(t.mutate([('a',0,1,0),('b',1,2,0),('c',0,2,0),('w',0,2,1)],0)))"`
Output: `[(('i','w','c'), 2, 2, 1)]`. This is a degree-1 loop at vertex 2 after one step, of type (c). The title says "no loop of degree >= 1 at all", and that is false for any parent with a relation omega parallel to a degree-0 arrow out of i. The same applies to type (b) when a relation arrow ends in i and the quiver has a 2-cycle; that cannot happen for an acyclic parent, but (c) can.

What survives: a loop at i needs a loop in Q at i (consequence 1). A relation arrow is a loop at i only if the relation is a cycle at i. Acyclic parents have none. So "no loop at i of any degree for the 25 acyclic parents" still follows, provided the 25 parents are acyclic as quivers with relations. The author asserts this from E-168/E-198 but the submission does not say it was re-checked here. The conclusion "E-235 applies and the failures are not cancellation failures" is therefore probably right. But the stated reason (Q concentrated in degree 0) is wrong, and the quiver-level reachability statement is wrong for non-hereditary parents.

Minor: the title's parenthetical, "(they are reachable only by iterating the dg mutation)", is garbled. The children are one-step mutations, and what the author means is the loops. Consequence (4) says the repo "refuses to mutate at a vertex with a loop", but that refusal is E-216's check on the Gabriel quiver of H_0. It is not the same as the Thm 1.1 hypothesis on the dg quiver.

## New?

Grepped `research/*.md` for "loops of degree", "positive-degree loop", "unreachable" and "one-step". The only related records are E-235 (cancellation proved, "loops of degree >= 1 at i not covered"), E-240, E-189 (Thm 1.1 needs only "no degree-0 loop") and E-216. The (a)-(d) classification is not recorded. It is elementary, but new to the record.

## Evidenced?

The classification proof is stated tightly enough to be believed. The reachability claim is not. It rests on "Q is ungraded", which the sample (random ungraded parents, relations never encoded as degree-1 arrows) cannot test. No run included relation arrows. The depth-3 example likewise uses a hereditary A_3 start, so it says nothing about parents with relations. "Iterating is needed" is shown only for hereditary starts.

## Scope

Narrowed wording that the evidence supports: "For any graded quiver Q, the loops of M are exactly (a)-(d). For Q concentrated in degree 0 (hereditary or other relation-free parents) M has no loop of degree >= 1. For acyclic parents with relations, M has no loop at i, but may have degree-1 loops at other vertices already after one step." The statement that the 25 failing children are not a cancellation failure holds for the "no loop at i" hypothesis only, and E-235 is conditional on that.

## Required for acceptance

1. Retitle and rewrite consequences (2)-(4) to the narrowed scope above. Remove "ungraded implies no loop of degree >= 1" for non-hereditary parents. Cite const.usual_algebra (l. 291) and the A_3/rad^2 example.
2. Add the counterexample above (or one like it) and a Part A/B variant that adds degree-1 relation arrows (omega parallel to paths), so that "loops at i" is tested with relation arrows present.
3. State, from E-168/E-198 and with a check, that each of the 25 parents has an acyclic quiver with relations in the Q^(1) sense, so that no relation arrow is a loop at i. Name the command.
4. Fix the title parenthetical, and say that E-216 is a check on the Gabriel quiver, not on the dg quiver.
5. Chair's sub-bullet closure ("loops of degree >= 1 at i not reachable") may stand only for loops at i, and only with item 3 done. All items doable in one sitting.
