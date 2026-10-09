# The Hom(T,T[m]) test accepts all 33 edges of the three E-158 paths: under the J = 0 premise each edge's T is self-orthogonal with Cartan(End T) equal to the child's Cartan matrix (not End T ~ child); all 25 key-keeping `tiltingPlus`-failing n = 7 children have such a path, 23 by their own path and 2 (c1 children 13 and 15) only by key equality

author: skeptic · round: 053 · kind: result
thread: T10 · bears on: H-015, E-155, E-158, E-159, E-161
scope: n = 7, key classes 1 and 2; Cartan-matrix level only (not End(T) ~ child), under the premise J = 0; the three paths printed in E-158 (c2 child 6, c1 child 5, c1 child 12; 12 + 12 + 9 = 33 edges); E-149 walk rebuilt at 20 000 expansions (c1: 16 failing children, c2: 9); Cartan-level End(T), generation assumed. Not covered: see "What the test does not check".

## Response to referee

1. Done. Title, Claim and Scope now say "Cartan(End T) equals the child's Cartan matrix" and "under the J = 0 premise"; no derived-equivalence claim. The abstract wording "independent tilting test" is read as "same verdict as the gate, independently coded".
2. Done. For a two-term complex T, Hom(T,T[m]) can be nonzero only for m in {-1,0,1}, so testing m = +-1 is enough; the "not argued" remark is replaced by this (it was an overstatement of the gap). Generation is still not tested.
3. Done. Per-edge output committed: `skeptic_replay3_c2_6.txt`, `skeptic_replay3_c1_5.txt`, `skeptic_replay3_c1_12.txt` (about 1 KB each; each ends `RESULT edges N all-accept True meet True`), plus `skeptic_key13_out.txt` and `skeptic_key15_out.txt` (new script `skeptic_key15.py`).
4. Done. Two children are covered by key equality only: c1 child 13 (= child 5, key f7abe9) and c1 child 15 (= child 14, key b32eca; keys equal and non-None, re-run now, parent depth 8 and v = 7 for both). So the 25 split as 23 with their own edge-tested path (E-155, E-158, E-161 paths) and 2 by key equality. The title was wrong to name only child 13; fixed.
5. Next round: End(T) ~ child with relations. Not done.

## Claim

Replaying each printed E-158 path edge by edge (F v: step at v; R v: step at v of the dual, carried back) with the E-159 test (`skeptic_tilt.py`, written without `tiltingPlus`, `perI` or the gate), every one of the 33 edges has Hom_K(T,T[-1]) = 0 and Hom_K(T,T[+1]) = 0 for T = T_v(a), and the dimension matrix of Hom(T_i,T_j) equals the Cartan matrix of the stated child (transposed, labelled). Both ends of each path have equal, non-None `canonicalKey`. Together with E-159 (E-155 paths) and E-161 (c1 14/15), every one of the 25 E-149 failing children (23 by own path, children 13 and 15 by key equality only) is now joined to an LNA/dual of its class by a path whose edges were each passed by an independent tilting test. Refuted by any edge of these paths with a nonzero Hom(T,T[+-1]), which did not occur.

It does not claim the 25 children are derived equivalent to the LNA: see below.

## Evidence

| child | key | LNA/dual | edges (child + LNA side) | tilt | cartan | ends meet |
|---|---|---|---|---|---|---|
| c2 child 6 | 0d3498 | #0 | 7 + 5 = 12 | 12/12 | 12/12 | yes |
| c1 child 5 | f7abe9 | #9 | 7 + 5 = 12 | 12/12 | 12/12 | yes |
| c1 child 12 | dc45e9 | #6 | 4 + 5 = 9 | 9/9 | 9/9 | yes |

(Each printed line also shows the gate context: J = False and tiltingPlus = True at every edge, as the paths were generated; and every Hom(T,T[-1]), Hom(T,T[1]) dimension is 0.) Moves replayed, in order, child side / LNA side: c2 6: `F1 F7 R3 R5 R5 F4 R6` / `F1 F1 F2 R6 R5`; c1 5: `F1 F3 F5 F4 R7 F3 R1` / `F7 R2 R1 F2 F7`; c1 12: `R7 R4 R7 R1` / `F1 F5 F4 R7 F3`. Per-edge output is committed as `skeptic_replay3_<class>_<child>.txt` in this folder. c1 child 13: `canonicalKey` equals that of child 5, non-None, same labelled quiver, parent depth 8 and v = 7 in both (`skeptic_key13.py`); so child 13 (and, likewise, child 15 with child 14: `skeptic_key15.py`) is covered only by being the same algebra as child 5 up to the key, not by its own replay. Power of the test was shown in E-159 (rejects all 25 failing J != 0 steps, 200 J = 0 controls accepted); I did not re-run those here, only reused the same functions. All 33 edges having J = False means these edges are also ones the gate-side test already accepted; in that sense the test agreeing here is the expected outcome (E-159 saw gate = tiltingPlus = mine on 1081 step tests), and what is added is code independence, not information.

### What the test does not check (exactly)

1. Generation. T is checked to be self-orthogonal (Hom(T,T[m]) = 0, m = -1, 1; m beyond +-1 is not needed: T is a two-term complex, so Hom(T,T[m]) = 0 automatically for |m| >= 2) but never to generate K^b(proj A). Okuyama-Rickard is assumed for "replace one summand by its left approximation".
2. End(T) ~ child as algebras. Only the dimension matrix Hom(T_i,T_j) is compared with the Cartan matrix of the repo's child. The child's relations (e.g. which parallel-arrow relations hold) are not compared; two algebras with equal Cartan matrices can differ. This is the main gap: the claim "edge x -> y is a derived equivalence" would need End(T) ~ y.
3. The input algebra's relations are trusted (`procedure.relationsFrom`, repo `reducePathAlgebra`, `dualPathAlgebra`); mutation output is the source of both T and the child.
4. R edges use op-duality: the test is run on dual(x) and the child is dualized back; the duality D: mod A -> mod A^op is assumed to carry derived equivalences, which is standard but not tested.
5. Path closure is by `canonicalKey` equality (isomorphism of quivers with relations up to relabelling, as the repo defines it); the key's cap (None above 720 relabelings) was not hit at the two ends or at any intermediate node of these paths (all printed node tags are non-None), but I did not check that the key distinguishes non-isomorphic algebras at these nodes.
6. Whether a vertex v whose step stays inside the same key class is a mutation at a vertex without loops is the gate's job; I trust that `sorted(c.vertices()) == V` (vertex-preserving) is the only condition I check.
7. Child 13 as above. The parent -> child edge of these three children (parent depth 8) was not replayed here; E-159 reports its 25 parent paths separately and I did not check which parents they cover.
8. Nothing here says the premise "J = 0 and tiltingPlus implies tilting" holds in general; it holds on 33 + 324 + 13 tested edges.

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 1 20000 /tmp/sk53/c1.pkl 100   # 500 s, checkpoints; rerun resumes (13 s); c2: "7 2" 381 s
S=workshop/rounds/053/skeptic_replay3.py
timeout 10m .venv/bin/python -u $S /tmp/sk53/c2.pkl 2 6 0 'F1 F7 R3 R5 R5 F4 R6' 'F1 F1 F2 R6 R5'      # ~1 min
timeout 10m .venv/bin/python -u $S /tmp/sk53/c1.pkl 1 5 9 'F1 F3 F5 F4 R7 F3 R1' 'F7 R2 R1 F2 F7'      # ~1 min
timeout 10m .venv/bin/python -u $S /tmp/sk53/c1.pkl 1 12 6 'R7 R4 R7 R1' 'F1 F5 F4 R7 F3'              # ~1 min
.venv/bin/python workshop/rounds/053/skeptic_key13.py /tmp/sk53/c1.pkl
```
Each replay ends with `RESULT edges N all-accept True meet True`.

## Prior record

E-158 printed these paths and says they were not replayed by the skeptic's test; E-159 and E-161 did the E-155 and c1 14/15 paths. This closes that gap; no new mathematics, and no entry in RETRACTIONS.md concerns it. It discharges the "Hom replay of the 3 E-158 paths" item of STATE.

## Code changed

None in the library. New: `workshop/rounds/053/skeptic_replay3.py` (replay_13 generalised: class, indices and move lists as arguments; also prints Hom(T,T[+1]) and node key tags), `skeptic_key13.py`. No tests apply.

## Next

- skeptic (next): quiver-level End(T) on the 13 E-161 and 33 E-158 edges (the real gap, item 2 above), e.g. compute End(T) as quiver with relations from the chain-map basis and compare with the child's, at least on the Cartan-equal-but-relations-differ cases.
- toolsmith: the three parent -> child edges (parent depth 8) of these children under the same test, if the parents are in the pickle.
- theorist: which hypothesis of AI 2.31/2.32 a J = 0 step needs (generation), as in the agenda.
