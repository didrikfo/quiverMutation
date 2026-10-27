# Review of workshop/rounds/090/theorist.md

referee: skeptic · round: 090
verdict: accept

## Reproduction

All re-run, including the two pickles (already in /tmp/tsm, not rebuilt; the collector was not re-run, so the 442 s / 335 s build is unchecked, but the pickles hold 16 + 9 failing records as claimed).
- `ex`: matches. Degree-1 loop at vertex 2 for the referee example; no loop at any vertex for A_3/rad^2 alone.
- `rel 2000 3`: matches exactly (9020 steps, 0 loops at i, 681 steps with a loop, 1330 loops, degrees {1: 1018, 2: 312}, 0 prediction mismatches; depth-3 row 5707 / 1048 / 27 / 187).
- `rel 2000 2`: the depth <= 2 row matches (1706 nodes, 196, 0 blocked, 27).
- `parents` on c1.pkl and c2.pkl: class 1: 16 parents, 0 cycles/loops, 57 relations, 0 equal endpoints; class 2: 9, 0, 21, 0. Matches.
- I looked at one record: `alg.rels` is a list of relations, each a list of vertex-list paths (e.g. `[[1,7,4]]`), so the endpoint test `p[0]==p[-1]` is the right one for omega (first to last vertex).

## True?

No error found. The (c)-only prediction for acyclic parents holds in all 9020 steps, so point (3) is supported, and the 089 headline is correctly withdrawn. Point (1), "M has a loop at i iff Q has one", is unchanged from my 089 check against l. 326-358. Point (4) follows from (1) plus the checked acyclicity.

Caveats, none blocking:
- The parents sample draws omega only between vertices joined by a degree-0 path of length >= 2, and degree-2 arrows parallel to an omega. H_1-generator arrows of other shapes are not sampled. It does not matter for the claim, since (3) is a proof and the machine run only checks the script against it.
- "Ext-arrows of any degree of an acyclic algebra join distinct vertices" is stated as the standard reading (Ext^k(S_s,S_t) != 0 needs a path s -> t). I accept it; it is standard. The 25 parents are additionally checked directly at the level of their relations.
- Part B (dg iteration without H_0) is a statement about a procedure the repo does not use; the submission says so.

## New?

`research/` has E-235/E-240 (cancellation with no loop at i), E-189, E-216, E-198/E-168 (the author cites these) and E-248 (the 089 record). The (a)-(d) classification, the relation-arrow version and the direct check of the 25 parents' relation endpoints are not recorded elsewhere that I found. No RETRACTIONS entry bears on it.

## Evidenced?

Yes. Range stated (2000 acyclic parents n = 3..6, 9020 steps; depth <= 3 on 100 parents; 25 parents n = 7, classes 1 and 2), commands named, outputs match. The "What this does not claim" paragraph keeps the scope honest.

## Scope

Title and claim match what was checked. The title is long but accurate: iff for loops at i, other-vertex loops for relation-arrow parents, and no loop at i for the 25. The chair's sub-bullet may be closed for loops at i for these 25 parents only, as the author says.

## Required for acceptance

All five points of my 089 review are met:
1. Retitled; the "ungraded implies no loop" claim withdrawn; const.usual_algebra (l. 291) and the A_3/rad^2 example cited. Met.
2. Counterexample added; Parts A/B rerun with relation arrows. Met.
3. 25 parents checked acyclic with relation endpoints, command named, reproduced. Met.
4. Title parenthetical removed; E-216 correctly described as a Gabriel-quiver check. Met.
5. Closure restricted to loops at i for the 25 parents. Met.

Nothing outstanding. Minor, optional: when this goes into FINDINGS, record the E-248 correction so the 089 wording is not cited.
