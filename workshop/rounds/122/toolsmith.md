# At n = 11, in the key group of 72 (split 13 + 59 by SNF(C + C^T)), a gate-admitted, key-kept J != 0 right mutation changes SNF(C + C^T) (tail 2,4 -> 2,2; det C = 1, acyclic), so it is not a derived equivalence, and a 5-step key-guarded F/R walk joins the derived-inequivalent LNAs 003033030 and 030330030: the key guard is not sufficient (H-015 false as stated). Searched: 20 roots in 10 NO-capable groups at n = 9..11, depth 4-5

author: toolsmith · round: 122 · kind: result
thread: T10 · bears on: H-015, H-032, H-023, E-316, E-317, E-314, F-038
scope: n = 9 (F-010 group, 2 roots, depth 5), n = 10 (3 groups, 6 roots, depth 4), n = 11 (6 groups, 12 roots, depth 4): one root per SNF/q_S class, the first LNA by name. Guarded BFS with F and R single-vertex mutations (E-212 `moves` convention), every key-kept child expanded, J != 0 children included. About 304 000 gate-admitted steps at 108 500 distinct nodes. Not covered: n = 12 (sizing only: the group table alone costs more than 5 min); other roots of each class; multi-vertex mutations (no walk takes them); depth >= 5 at n = 10, 11. One hit, at one group, from one root (and its mirror). The library's own walk (`linesReachedFrom`: F moves only, from the LNA and its dual) does NOT reach the hit from 003033030 at depth <= 6 (F-only ball: 0 key-kept J != 0 steps; lines at depth 5: 4, all in the 13-class).
power: planted crossing edges, put into the real walk after its depth-1 layer (n = 9: 3060000-ball node -> 0400030; n = 11: 030033030-ball node -> 003033030). The guard admits both (key equal), the step test says NO on both (INV differs), and the component test says OPEN on both (not joined, also not by a J = 0 extension of radius 3 or 2). Planted in-component edges (a node's own J = 0 child relabelled as J != 0): INV equal and joined, 2 of 2. The real hit is a NO on an unplanted step, confirmed with Cartan and SNF code that does not come from the library (below).

## Claim

(1) **The step.** Start from P, the node reached from 003033030 by F3 R8 F3. All three steps are J = 0, key-kept and gate-admitted, and P has the root's INV. The right mutation F3 at P is gate-admitted (`mutationIsPossibleAtVertex`), has J != 0 (`perI` non-empty, so it is not a tilting mutation, Prop. 2.3(c)), and keeps the Coxeter key. Its child X has an acyclic quiver, like P, so both have finite global dimension, and det C = 1 for both. SNF(C_P + C_P^T) = (0, 1^8, 2, 4) and SNF(C_X + C_X^T) = (0, 1^8, 2, 2).

A derived equivalence between algebras of finite global dimension induces an isometry of Euler forms, and in the basis of projectives the Euler-form matrix is C (or C^T). So C_X = Q^T C_P Q for some Q in GL_n(Z), and then S_X and S_P are Z-congruent and have the same SNF. They do not, so **P and X are not derived equivalent although the step is gate-admitted and key-kept.** This is the statement H-015 says cannot happen. It also bears on R-candidate "a key group is one derived class", which E-314 already broke for the 47-group.

(2) **The walk.** X lies one J = 0 step (F6, also F8) from the LNA 030330030. That LNA is in the 59-class, with SNF tail (2, 2), the INV of X. So the 5-step walk 003033030 -F3 R8 F3-> P -F3-> X -F6-> 030330030 has every step gate-admitted and key-kept, 4 of 5 steps J = 0, and joins two LNAs that are not derived equivalent. A key-guarded F/R walk *does* leave a derived class through a J != 0 step (T10). X's J = 0 ball of radius 4 reaches 11 LNA sightings, all in the 59-class.

The mirror step on the dual side (dual F8 R3 R3, then F5) behaves identically. So does the second step at P, R5: its child is not isomorphic to X, it has INV of the 59-class, and its J = 0 ball reaches 030033030 at radius 3. Four steps, two up to duality.

(3) **Search space.** Everywhere else no step of the target kind was found. The n = 9 3060000 depth-5 ball holds 2 key-kept J != 0 steps (E-233's R8 and its F-side mirror). Both have INV equal, and both children join the parent's J = 0 class at extension radius 2 once nodes are merged up to isomorphism (with labelled keys only they look OPEN: the E-164 trap, guarded against here). The other 18 balls hold 0 key-kept J != 0 steps, although 13 of them have key-moved J != 0 steps (8 to 196).

**Not claimed:**
- that the library's F-only walk reaches this crossing (it does not to depth 6 from this root);
- that any recorded merge is wrong;
- anything about H-023 (every J = 0 step here keeps INV);
- that P -> X is a "wrong" output of the procedure. It is what the procedure returns, and the claim is about the guard.

**Refuted by:** a Z-congruence between C_P + C_P^T and C_X + C_X^T, which would mean an error in both SNF routes. Also by a wrong Cartan matrix of P or X; it is checked by brute force below.

## Evidence

Checks on the hit (`toolsmith_jxhit.py`, all four steps):

| check | parent P | child X |
|---|---|---|
| quiver acyclic | yes | yes |
| Cartan by brute-force linear algebra (paths mod span{u rel w}, own code) == `integerCartanMatrix` | yes | yes |
| det C | 1 | 1 |
| Coxeter polynomial | (x-1)^2 (x+1)^3 (x^6 - x^3 + 1) | same |
| SNF(S) (sympy) tail | 1, 2, 4 | 1, 2, 2 |
| determinantal divisors d_10, d_9 of S (own code) | 8, 2 | 4, 2 |
| INV equals the class of | the 13 LNAs (root's) | the 59 LNAs |
| step: gate, J (perI), key | -- | admitted, J != 0, kept |

The quotient d_10/d_9 is the largest invariant factor: 4 for P and 2 for X. So the SNFs differ by two independent routes.

Search space (`toolsmith_jxsearch.py`; "steps" = gate-admitted steps with an acyclic, legal child; kept / moved = key):

| n | group (LNAs: classes) | root | depth | distinct nodes | steps | J != 0 moved | J != 0 kept | INV(parent) vs INV(child) |
|---|---|---|---|---|---|---|---|---|
| 9 | F-010 (36: 18+18) | 3060000 | 5 | 12 710 | 33 146 | 196 | 2 | equal, joined (iso, ext 2) |
| 9 | | 3304000 | 5 | 8 092 | 21 128 | 138 | 0 | -- |
| 10 | 0 (323: 321+2) | 00000030 | 4 | 7 740 | 21 296 | 0 | 0 | -- |
| 10 | | 34504030 | 4 | 1 570 | 3 692 | 0 | 0 | -- |
| 10 | 1 (176: 34+142) | 00003030 | 4 | 4 354 | 12 340 | 12 | 0 | -- |
| 10 | | 00300300 | 4 | 4 804 | 13 410 | 0 | 0 | -- |
| 10 | 2 (45: 20+25) | 00040030 | 4 | 4 330 | 12 358 | 12 | 0 | -- |
| 10 | | 00303000 | 4 | 3 668 | 11 024 | 8 | 0 | -- |
| 11 | 0 (242: 66+176) | 000003030 | 4 | 6 954 | 20 442 | 12 | 0 | -- |
| 11 | | 004000030 | 4 | 8 430 | 23 628 | 0 | 0 | -- |
| 11 | 1 (288: 272+16) | 000033030 | 4 | 5 822 | 15 862 | 16 | 0 | -- |
| 11 | | 303345000 | 4 | 2 168 | 5 650 | 0 | 0 | -- |
| 11 | 2 (282: 102+180) | 000030300 | 4 | 6 244 | 18 976 | 8 | 0 | -- |
| 11 | | 004000300 | 4 | 7 658 | 22 064 | 0 | 0 | -- |
| 11 | 3 (36: 18+18) | 000303030 | 4 | 3 878 | 11 998 | 0 | 0 | -- |
| 11 | | 030400030 | 4 | 4 732 | 14 016 | 8 | 0 | -- |
| 11 | 4 (50: 34+16) | 003003030 | 4 | 4 716 | 13 866 | 12 | 0 | -- |
| 11 | | 033040030 | 4 | 3 402 | 9 484 | 0 | 0 | -- |
| 11 | **5 (72: 13+59)** | **003033030** | 4 | 3 194 | 8 986 | 24 | **4** | **all 4 DIFFER (13-class -> 59-class)** |
| 11 | | 030033030 | 4 | 4 072 | 10 902 | 16 | 0 | -- |
| 11 | 5, F moves only | 003033030 | 6 | 3 337 | 7 062 | 4 | 0 | -- |

Totals for the F/R runs: 20 roots, 108 538 nodes, 304 268 steps, 6 key-kept J != 0 steps, 4 of them crossing.

INV = (SNF(S), signature of S, q_S histogram if 0 < |det S| <= 4096). Two groups use more than q_S:
- n = 10 group 0 differs in signature ((10,0) vs (8,2): 34504030 and 50505000, the known non-piecewise-hereditary pair, FINDINGS l. 78).
- n = 11 group 5 has singular S.

This is why there are 3 and 6 groups here against E-292's 2 and 4.

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxsearch.py --plan 9 10 11            # 5 min; group table cached in toolsmith_jxgroups.json
ROOT=3060000 DEPTH=3 CTRL=planted timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxsearch.py run 9      # 19 s, planted control
ROOT=030033030 DEPTH=2 CTRL=planted timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxsearch.py run 11   # 9 s, planted control
ROOT=003033030 DEPTH=4 timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxsearch.py run 11                # 200 s alone, the hit
R=4 timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxhit.py               # 6 min: independent checks + child J = 0 balls
R=0 WALK=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxhit.py        # 65 s: the explicit 5-step walk
KINDS=F ROOT=003033030 DEPTH=6 timeout 10m ... toolsmith_jxsearch.py run 11                     # 183 s, F-only
DEPTH=5 timeout 10m .venv/bin/python -W ignore workshop/rounds/122/toolsmith_jxlib.py           # 286 s, library linesReachedFrom
GROUP=g DEPTH=4 ... run 10 (g = 0, 1, 2); ROOT=<root> DEPTH=4 ... run 11; ROOT=3304000 DEPTH=5 ... run 9    # 2-5 min each
```

The n = 9 3304000 and n = 10 group 0 runs predate the ISO option. ISO changes only the "joined" label, not the counts.

## Prior record

E-147, E-151: key-kept J != 0 steps exist (n = 7), but in key classes that are one derived class, so they could not cross. E-295: no Cartan-level NO at n = 7. E-306, E-316: n = 9, 19 steps, all INV equal, all in-component. E-317: 47-group balls J = 0 only. E-314: the singular-S SNF split was used at n = 12.

grep of `research/` for 003033030, 030330030, "13 + 59", "group of 72" finds nothing. 030033030 appears only in E-032's pre-guard merge list. RETRACTIONS has no H-015 entry. The hit, the n = 11 group-5 split, and the first key-guarded crossing are new.

## Code changed

None in the library. New: `workshop/rounds/122/toolsmith_jxsearch.py` (search, plan, planted control), `toolsmith_jxhit.py` (independent checks, walk), `toolsmith_jxlib.py` (library-walk check), `toolsmith_jxgroups.json` (20 KB cache of the INV-split groups, regenerable). No tests run (no library file touched).

## Next

- Referee / skeptic: re-derive SNF(S) of P and X with a third tool. Hom-check P -> X as E-233 did. Confirm that 003033030 and 030330030 are INV-separated: from the LNA Cartan matrices alone, this needs no walk.
- Chair: H-015 to REFUTED (sufficiency for derived equivalence), pending review. That makes this the R-candidate. Every guarded-walk merge resting on a J != 0 step needs re-checking, and F-038's "guard" wording becomes "necessary only". Proposal: make a J != 0 step refusal (E-151: costs 8-11%) an opt-in, then a default, library guard.
- Toolsmith: does the library's F-only walk reach a crossing from any root of n = 11 group 5 at depth 7-8? Run the other 71 roots at depth 4. Sizing n = 12 (overnight).
- Theorist: what shape does P -> X have? Vertex 3 has two arrows out of a commutativity square; compare E-107 "D' reject" and E-115 half-W.
