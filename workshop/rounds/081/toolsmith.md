# At n = 7, class 2, seeds 1 and 8 of the class are joined by a 4-edge Hom-verified J = 0 chain, and the E-222 rep 2 chain child -> seed 1 -> seed 8 -> parent (11 edges) is replayed as one Hom-checked run under the J = 0 premise; record 6 of E-212 gets 8 edges, pieces verified separately (meets are the necessary E-212 isomorphism test)

author: toolsmith · round: 081 · kind: result
thread: T10 · bears on: E-222, E-216, E-212, E-157, E-213
scope: n = 7, key class 2 only; meets by the necessary-only iso test; Hom+Cartan demonstrated no NO power on 26 key-kept corruptions (the meet iso rejected all 26) (rep 2 of E-222, record 6 of E-212); balls of depth 3 from each of seed 1 and seed 8 of `classes[base]` (not the whole class), record 6 searched to target depth 4 / algebra depth 6 with iso meets; J = 0 premise as for the F/R moves (gate, `tiltingPlus`, key kept); generation of K^b(proj) assumed (E-168); class 1 and the other 5 seed pairs not run.

## Response to referee

1. Title and Claim narrowed (done). The meet test is the sign-aware structure-graph isomorphism plus the E-212 `iso`; both collapse the arrow identity inside relation terms and count parallel arrows only, so they are necessary conditions, not a proof that M and M' are the same algebra. The Claim now says so. "Closes" is kept only as "under the J = 0 premise plus E-168, with meets by the necessary test".
2. Rep 2 chain replayed as one run (done): `toolsmith_chain11.py`, one process, same seed objects S[1], S[8], 11 edges (2 + 4 + 5), 3 meets, output `toolsmith_chain11_out.txt` (a few minutes with the preamble). All 11 edges pass, 3 meets pass. Limit of the replay: each edge is checked in the direction seed -> end (the E-212 `horn` is directional), so the "child -> seed 1" and "parent" legs are the seed legs read backwards, glued by the meet isomorphisms; I do not compose the isomorphisms into one explicit vertex map. "Single run" means one process and one set of objects, not one directed walk. Record 6 (12 edges) is not part of this run; its two pieces stay separately verified (referee re-ran them).
3. Per-edge Hom lines (done): see the Evidence section; each edge has sum Hom(T,T[-1]), sum Hom(T,T[+1]), sum dim Hom(T,T) against the Cartan sum. All the Hom(T,T[+-1]) sums are 0 and every dim-Hom sum equals the Cartan sum (19, 20, 21, 19, 23, 24, 22, 25, 25 along the edges).
4. Corruption control (done, cheap): replacing one move of the 4 seed-link edges by each other (F/R, v = 1..7) gives 52 corruptions: 26 NOSTEP (key-moving), 0 Hom/Cartan rejections, 26 with all edges passing and then 26 of 26 meets fail the signed iso. So the Hom+Cartan test alone showed NO power here (0 of 26 key-kept swaps rejected; consistent with E-212's 1 in 195), and the discriminating step is the meet isomorphism, which is the necessary-only test. This weakens "Hom-verified" as a statement of evidence: the evidence that these edges are tilting steps is the premise J = 0 plus E-168, with Hom as a consistency check.
5. Point 4 of the review (rec 6 and rep 2 share the seed pair): accepted, the title no longer says it closes two results.

## Claim

Seeds 1 and 8 of `classes[base]` (class 2, n = 7) are not isomorphic and are joined by the chain seed1 -F6-> M <-R5 R5 R6- seed8 (1 + 3 = 4 edges), every edge passing the E-212 Hom+Cartan test (Hom(T,T[+-1]) = 0, dim Hom = Cartan of the start) and the meet M being isomorphic by the sign-aware structure-graph test (and the E-212 `iso`); both are NECESSARY conditions (arrow identity inside relation terms collapsed, parallel arrows only counted), not a proof of algebra isomorphism. Hence the E-222 rep 2 child (ends at seed 1) and parent (ends at seed 8) are linked: child -> seed 1 -> seed 8 -> parent is a chain of 2 + 4 + 5 = 11 edges replayed in one process (`toolsmith_chain11.py`, each edge Hom-checked seed -> end, three meets by the necessary iso), so under the J = 0 premise plus E-168 all four reps of E-222 now have parent ~ child by explicit chains. Record 6 of E-212 (undecided there: "no child path") is also joined: child -> seed 1 by F7 F6 (2 edges, iso meet), parent -> seed 8 by F1 F2 | F1 F1 R6 R7 (6 edges, iso meet), both Hom-verified; with the 4-edge seed link this is 12 edges. The old miss was the labelled `canonicalKey` meeting (E-164 effect), not an inequivalence. Not claimed: that the premise holds (these are all J = 0 steps, so the chain cannot test it), that depth is minimal (first hit by total), or that the pieces compose with the isomorphism relabellings checked as one run (each piece is verified separately; meets are isomorphisms, so the composite is consistent and the replay is one process but not one directed walk with a composed vertex map). The claim is refuted by a `!` (nonzero Hom) or NOSTEP in the output, or a failing meet iso; neither occurs.

## Evidence

| item | side | moves (to meet) | edges | Hom+Cartan | meet |
|---|---|---|---|---|---|
| seed link | seed 1 | F6 | 1 | pass | signed iso, E-212 iso |
| | seed 8 | R5 R5 R6 | 3 | pass | |
| record 6 | child -> seed 1 | (none) \| seed side F7 F6 | 2 | pass | signed iso |
| | parent -> seed 8 | algebra side F1 F2 \| seed side F1 F1 R6 R7 | 6 | pass | signed iso |

 Rep 2 chains are those of E-222 (child: seed side F6 F7, 2 edges; parent: algebra side F5, seed side F1 F1 F2 R7, 5 edges), now replayed in `toolsmith_chain11_out.txt` with per-edge lines, e.g. `seed8 edge 2 R5: sum Hom(T,T[-1]) = 0, sum Hom(T,T[+1]) = 0, sum dim Hom(T,T) = 21, sum Cartan(child) = 21 -> pass` (seed link: F6 19/19; R5 20/20, R5 21/21, R6 19/19; all in that file). Seed ends of record 6 (child 1, parent 8 of 14) are the same pair as rep 2, so one seed link closes both. Ball sizes: seed 1 ball depth 3 has 129 keys, seed 8 ball 192; target ball (all 14 seeds, depth 4) 7208 keys. Because the seed-to-seed search is a separate ball per seed (E-222's pooled ball assigns each key to one seed), this is new information: pooling hid the link. Seeds 1 and 8 have different labelled keys and the minimal total is 4 (a total 0 would have been an isomorphism); other seed pairs not tested. Positive control for the verifier is E-212's (power thin: 1 Hom rejection in 195 corruptions); negative control now added (response 4): Hom+Cartan rejected 0 of 26 key-kept single-move swaps of the seed link; the meet iso rejected all 26.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/081/toolsmith_seedlink.py 2 3 3 > workshop/rounds/081/toolsmith_seedlink_out.txt   # 3 s of search, ~60 s with preamble
timeout 10m .venv/bin/python -W ignore workshop/rounds/049/toolsmith_collect.py 7 2 20000 /tmp/ts73/c2.pkl 100    # 405 s here (311 s in E-212)
timeout 10m .venv/bin/python -W ignore workshop/rounds/081/toolsmith_rec6.py /tmp/ts73/c2.pkl 6 4 6 > workshop/rounds/081/toolsmith_rec6_out.txt   # ~100 s
```
```
timeout 10m .venv/bin/python -W ignore workshop/rounds/081/toolsmith_chain11.py > workshop/rounds/081/toolsmith_chain11_out.txt   # item 2-4: 11-edge replay, Hom lines, corruption control
```
Outputs are small and committed.

## Prior record

E-222 left rep 2's seed link open (1 vs 8); E-212 left record 6 undecided; E-157 missed it at 6 + 5 labelled. grep of EXPERIMENTS/RETRACTIONS for a seed-1/seed-8 link finds none. New: both closed; the E-212 "seed-to-seed link" remark (parent and child can reach different seeds) is now answered by an explicit link for the one case that occurs.

## Code changed

None in `quivermutation/`. New scripts `toolsmith_seedlink.py`, `toolsmith_rec6.py` (the latter is `078/toolsmith_chains.py` with the record's pair; a global named `P` is the prime in the exec'd preamble, so the pair is called PAR, CHI). No tests run (no library code touched).

## Next

Chair: record this as an E-222 amendment (rep 2 chain, 11 edges, premise + E-168, necessary iso at meets) and E-212 record 6 pieces verified (class 2 9 of 9 only under the premise and the necessary-only iso). Experimentalist/skeptic: all 14 seeds of class 2 pairwise (ball of depth 3 each) to see whether the seed graph is connected by J = 0 chains; a key-equal inequivalent pair with a power-bearing NO test is still the open item (F-010 5 + 5 overnight). Toolsmith: class 1 record equivalent (E-212 class 1 not run).
