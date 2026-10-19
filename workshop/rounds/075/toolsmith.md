# At n = 7 (key classes 1 and 2), the printed representatives of the 4 unjoined presentation classes of E-213 and their parents each reach an LNA by J = 0 F/R paths of at most 3 + 4 edges when met up to isomorphism (tiltingPlus-level, not Hom-verified; parent and child may reach different seeds); a labelled-`canonicalKey` search (single run, round-049 move set) missed the class-1 child at 6 + 5

author: toolsmith · round: 075 · kind: result
thread: T10 · bears on: H-015, E-157, E-163, E-164, E-212, E-213
scope: n = 7; one representative per each of the 4 classes (class 1 once, class 2 three times), child and parent of it: 8 algebras, not the 14 edges (the other members of each class are relabelled copies, E-213); J = 0 premise (edges are gate + `tiltingPlus` + key-kept F/R moves, NOT re-verified by the Hom test); meets are sign-aware isomorphisms with a target ball to depth 4 (complete, uncapped); path edge lists were not recorded (depths only); generation of K^b(proj) assumed (E-168).

## Response to referee

1. Search-1 stdout committed: re-ran `toolsmith_join.py R 5 6 400000 DIR` for R = 0..3 (four in parallel, under 10 min, uncapped) into `toolsmith_join_r0_out.txt` .. `r3_out.txt`. Output reproduces the table below exactly (rep 0 child MISS at 6 + 5, levels 1, 6, 29, 109, 405, 1586, 6335; all other entries HIT with the same totals); only timings differ. The 6 + 5 sentence stays in the title, now marked as a single run with the round-049 move set.
2. Title and claim now say "tiltingPlus-level, not Hom-verified" and that parent and child may be joined to different seeds (so nothing about parent ~ child).
3. Reproduction block now says all scripts must be run from the repo root.
4. [next round] F/R move lists and the 073 verifier: not cheap (needs both scripts changed to store paths and a verifier harness); not done, stays under Next.

## Claim

Under the J = 0 premise, each of the 8 algebras (child and parent of the 4 E-213 representatives, names rep 0..3 = class 1 v 7, class 2 v 3, class 2 v 6, class 2 v 7) is joined to a seed LNA (or its dual) by a mixed F/R chain of `tiltingPlus`-true, key-kept edges: the child of rep 0, 1, 2 is isomorphic to a node at distance 2 from a seed (total 2), the child of rep 3 at distance 3, the parents at total 6, 6, 5, 6. So the 4 classes are not exceptions to "every failing child is joined to an LNA", under the same J = 0 premise as E-157/E-163, and inherit its limitation: the premise (Hom(T,T[-1]) = 0 for the steps on the path) is checked by `tiltingPlus`, not by an independent Hom test. It does NOT claim: that parent and child are derived equivalent (that needs the two seeds linked and the J != 0 step is the one that is not a tilting step; E-212 does that for 8 of 9 window steps, not here), nor anything about the other members of the classes beyond that they are relabelled copies (E-213), nor that a labelled search is wrong, only that it can miss. Refuted by: a Hom(T,T[-1]) != 0 edge on one of these chains, or an isomorphism accepted by the sign-aware test that is not one.

## Evidence

Step 0, the inputs: the 4 parent/child pairs are parsed from the printed reps in `rounds/074/skeptic_failkids_an_out.txt` (`NEW class rep` blocks) and rebuilt with `procedure.toPathAlgebra`. For each: parent and child have the Coxeter key of class 1 / 2; the printed step is gate-admitted, `perI` nonzero (J != 0) and `tiltingPlus` False (printed by `toolsmith_join.py`).

Search 1, labelled `canonicalKey` meets (`toolsmith_join.py`, tiltpath moves F/R, target ball = 12 (class 1) / 14 (class 2) seeds incl. duals, depth 5 = 15 862 / 25 452 keys, uncapped; algebra ball depth 6, cap 400 000 never reached):

| rep | class | child | parent |
|---|---|---|---|
| 0 | 1 | MISS (levels 1, 6, 29, 109, 405, 1586, 6335; 159 s) | HIT 3 + 5 = 8 |
| 1 | 2 | HIT 4 + 5 = 9 | HIT 3 + 5 = 8 |
| 2 | 2 | HIT 5 + 5 = 10 | HIT 2 + 5 = 7 |
| 3 | 2 | HIT 4 + 5 = 9 | HIT 3 + 5 = 8 |

Search 2, same moves, meets up to isomorphism (`toolsmith_isojoin.py`: WL hash of the structure graph with arrows as nodes, then explicit networkx isomorphism whose relabelled relations have the same gauge-canonical form, the rounds/074 `signedIso`; target ball depth 4, complete, 970 / 682 hash buckets; algebra ball depth 5, cap 400 000, nothing stopped):

| rep | child (algebra side + target side) | parent |
|---|---|---|
| 0 | 0 + 2 = 2 | 2 + 4 = 6 |
| 1 | 0 + 2 = 2 | 2 + 4 = 6 |
| 2 | 0 + 2 = 2 | 1 + 4 = 5 |
| 3 | 0 + 3 = 3 | 2 + 4 = 6 |

The class-1 child of rep 0 is the point: the labelled search found nothing in 6 + 5 (and took 159 s) but the child is a vertex-relabelled copy of a node two F/R steps from a seed. Totals in search 2 are the first level at which a meet occurs, so they are upper bounds on the shortest path; "alg side 0" means the child itself is isomorphic to a ball node. Nothing in either search is a no-hit for the Hom reading: no Hom test was run.

What the search can and cannot rule out: there was no miss in search 2, so there is no "none found" to bound. In search 1 the single miss (rep 0 child, labelled, 6 + 5, uncapped by CAP) is explained by labelled keys (the E-164 point, `canonicalKey` is not relabelling invariant, E-213 finds 5 labelled keys per class). Search 2 cannot rule out an error in the premise: it says nothing about whether the J != 0 step itself preserves the derived class (parent ~ child), and the 8 joins are to LNAs, possibly different ones.

## Reproduction

Run every command from the repo root (the scripts open `workshop/rounds/075/toolsmith_join.py` by relative path).

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/075/toolsmith_join.py R 5 6 400000 DIR     # R = 0..3, ~5-10 min each when four run together (target ball 230-330 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/075/toolsmith_isojoin.py R 4 5 400000 500   # R = 0..3, 60-90 s each (four together); prints to stdout
```
Outputs kept: `toolsmith_isojoin_r0_out.txt` .. `r3_out.txt`. Search-1 stdout is kept in `toolsmith_join_r0_out.txt` .. `r3_out.txt` (R = 0..3; add `> file`). DIR holds the target-ball pickle (default /tmp/claude-0/ts75; create it first).

## Prior record

E-157/E-160/E-163 joined the 25 window children (labelled keys); E-212 joined parent and child of 8 of 9 class-2 window steps with Hom-verified chains and iso meets; E-213 left the 4 classes as the open set; E-164 records that `canonicalKey` is labelled. New: joins for the 4 classes, and that iso meets reach them at total <= 6 where a labelled search needed more (the class-1 child never at 6 + 5).

## Code changed

New: `rounds/075/toolsmith_join.py`, `toolsmith_isojoin.py` (they exec the round-049 `toolsmith_tiltpath.py` preamble and copy the round-074 `struct` / `signedIso`). No library file touched; no tests run (none apply).

## Next

- toolsmith: record the F/R move lists of these 8 chains and push them through the 073 Hom verifier (`toolsmith_certificate.py`); then the seed-to-seed link between each parent's and child's LNA (needed for parent ~ child, the real question).
- toolsmith/skeptic: rerun E-157's record-6 class-2 MISS with iso meets (it may be a relabelling miss); make iso-meet the default in `meetingPoints` (E-164 backlog).
- Only reps of the 4 classes were joined; since the other members are relabelled copies (E-213), nothing else is needed for the J = 0 join of the 70.
