# At n = 8, 9 the only key-equal derived-inequivalent LNA pair is the F-010 pair at n = 9 (q_S separates it); q_S is constant on all steps of three balls, as guaranteed for J = 0 (derived-equivalence) steps, so those steps are a code check and only the 2 J != 0 steps could test it

author: toolsmith · round: 117 · kind: backlog
thread: T10 · bears on: H-015, R-candidate "a key-guarded walk stays in one class", H-024, E-292, E-300, E-238, F-010
scope: item 1 (re-verified by referee at n = 8, 9) and walk (depth 3 re-run by referee; depths 4, 5 not): every LNA at n = 8 (429) and n = 9 (1430), all 27 key groups with >= 2 LNAs, invariant (signature, |det S|, q_S histogram of S = C + C^T). Walk test: n = 9, balls of 3060000 (depth 4 + J != 0 scan of the depth-4 parents; 1584 nodes), 3304000 (depth 5 + J != 0 scan of depth-5 parents; 4046 nodes) and 2400230 (depth 3; 130 nodes) under the repo's key-kept J = 0 moves and key-kept perI J != 0 steps. Item 2: all 44 remaining single-edge A-side corruptions of the step-8 chain. NOT covered: n >= 10, walk depth beyond 4/5 (3060000 depth 5 and the 3304000 depth-6 parents not run), steps that leave the ball, HH, Z[Phi]-module class.
power: q_S says NO as an algebra invariant where NO is true: 3060000 vs 3304000/2400230 (certified inequivalent, F-010, quipu theorem E-238) differ on all cross pairs of the balls (386 x 282 nodes at depth 3: invariant sets disjoint); A-side step-8 corruptions: 44 of 44 rejected (see item 2). Walk side: the test could not fail for the J = 0 steps (derived equivalences keep q_S), so it has no NO power there; only the 2 J != 0 steps could have fired, and they were already joined. The cross control shows q_S separates the halves, not that any step could cross.
backlog: item 1 settled (pair exists at n = 9, none at n = 8; instrument exists, T10 premise still untested because no crossing step was found) · item 2 settled [next round: independent re-verification] (78 of 78 A-side corruptions rejected; E-300's "50 left" was 44)

## Response to referee

Verdict minor revision; the three required items:
1. Done. Claim (c) and the title now say that q_S invariance on the J = 0 F/R steps follows from their being tilting mutations (derived equivalences preserve the Cartan form up to Z-congruence, hence signature, |det S|, q_S). The 14 800 (14 815 J = 0) steps check the code (Cartan and q_S computation), not T10. The NO-capable effective sample is the 2 J != 0 steps (E-233, E-239), both already joined by chains, outcome not informative. Accepted in full.
2. Done. "never-triggered" is replaced by "cannot trigger on J = 0 steps"; the Power line now says the walk test could not fail for J = 0 and that the cross control shows only that q_S separates the two halves as algebras, not that any step could cross.
3. Item 2 marked `[next round]` for independent re-verification: it needs `/tmp/skeptic_c1_13_ck1.pkl` (about 18 min of rebuild in two calls) and four 6-8 min slices on 4 cores, about 40 min in total, which I did not run this sitting; the pickle is gone again. The 78 of 78 count stays as author-reported, unreproduced by the referee.
Also noted: item 1(a) is conditional on F-032's placement and F-010 (no cospectral quipu pair at n <= 8), neither re-checked by the referee. The "walk side has no NO case" statement stands, but now with the reason: J = 0 steps cannot say NO.

## Claim

Item 1. (a) n = 8: all 10 key groups with >= 2 LNAs have a single (signature, det, q_S) class, and no inequivalent key-equal LNA pair exists at all: all 429 LNAs are placed in quipu classes (`lnaStatus`: 429 QUIPU) and there is no cospectral quipu pair at order <= 8 (F-010), so equal key means equal quipu class (conditional on F-032's placement). (b) n = 9: exactly one group of 36 LNAs (the F-010 group) splits, 18 + 18 by q_S (same signature (8,1), same |det S| = 12, q_S histograms differ); 3060000 is in one half, 3304000 and 2400230 in the other (matching E-238's naming); the other 16 multi-LNA groups at n = 9 are single-class. So the pair asked for exists at n = 9 (and is the only one), and an instrument that says NO on it exists: q_S of the Cartan matrix. (c) Walk test (a code check, not a T10 test): along 14 817 key-kept steps of the three balls (J = 0: 4017, 10 495 and 303 in the three balls; J != 0: the one step of E-233 in 3060000, the one of E-239 in 3304000) parent and child have equal (signature, det, q_S): 0 NO. This is guaranteed, not found: J = 0 F/R steps are derived equivalences (tilting), which preserve the Cartan form up to Z-congruence, so q_S, signature and |det S| cannot change; these 14 815 steps test only the Cartan/q_S code. The only steps that could trigger are the 2 J != 0 steps, whose outcome is not informative since both were already joined by chains. Every node of each ball has the root's invariant (1 class per ball), as expected. Not claimed: the J = 0 premise (T5) or "a guarded walk stays in a class" (R-candidate) is supported beyond what E-233/E-239/E-244 already joined; both J != 0 steps were already joined by chains (E-233, E-239).

Item 2 [next round: independent re-verification, not re-run by referee]. All 44 untested A-side single-edge corruptions (78 candidates in all, not 84: 6 edges x (2 x 7 labels) minus the 6 true moves; the 34 done in round 114 plus these 44) are rejected, by the meet isomorphism every time (21 with an additional "nostep" failure, 23 meet only; Hom alone never decides). With round 114: A-side 78 of 78 rejected, B-side 77 of 78 (the 1 accepted is a valid alternative edge), end swaps 12 of 12. Same limit as before: the control tests the meet iso, not the Hom conditions.

## Evidence

| quantity | n = 8 | n = 9 |
|---|---|---|
| LNAs / key groups with >= 2 LNAs | 429 / 10 | 1430 / 17 |
| groups split by (sig, det, q_S) | 0 | 1 (36 LNAs: 18 + 18) |
| cospectral quipu groups (F-010) | 0 | 1 |
| LNAs placed in a quipu class | 429 | 1421 (9 not-quipu, own keys) |

Walk balls (n = 9, key = F-010 key): 3060000 depth 4: 1584 nodes, 4017 J = 0 steps + 1 J != 0 step, all invariant-equal. 3304000 depth 5: 4046 nodes, 10 495 J = 0 + 1 J != 0, all equal. 2400230 depth 3 (130 nodes, 303 steps) equal and has the 3304000 invariant. Cross control: the 3060000 ball vs the 3304000 ball have disjoint invariant sets (depth 3: 386 vs 282 nodes). Raw: `toolsmith_qsscan.txt`, `toolsmith_qsball.txt`, `toolsmith_s8ctl_Arest.txt`.

Verdict on "is there a NO-capable instrument for T10 at n <= 9": YES at the Cartan level for the one n = 9 key group (q_S, 2 classes); NO instrument at n = 8 (nothing to separate: all key-equal LNAs are equivalent) and at n = 7 (E-298, E-305). But the instrument only decides class membership of an algebra; a T10 NO needs a key-kept step from one half to the other, and none of 14 817 explored steps is one. Whether the walk can reach the other half at all is the question the instrument is waiting for; it needs deeper balls (depth 6 for 3060000 is ~19 000 parents, E-250 sizing), or a search for the other-half LNAs inside the walk of one half (a "hit" test), not a better invariant.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/117/toolsmith_qsscan.py 8 9      # about 1 min
DEPTH=4 timeout 10m .venv/bin/python workshop/rounds/117/toolsmith_qsball.py 3060000   # about 2 min
DEPTH=5 timeout 10m .venv/bin/python workshop/rounds/117/toolsmith_qsball.py 3304000   # about 7 min
DEPTH=3 timeout 10m .venv/bin/python workshop/rounds/117/toolsmith_qsball.py 3060000 3304000 2400230   # about 1 min, includes the cross control
# item 2: needs /tmp/skeptic_c1_13_ck1.pkl (two calls, 9 min each: workshop/rounds/069/skeptic_c1_13.py 1 540 23000), then
SLICE=34:47 timeout 10m .venv/bin/python workshop/rounds/117/toolsmith_s8ctl_Arest.py   # also 47:60 60:73 73:84; 4 cores, 6-8 min each
```

## Prior record

F-010, E-238 (naming of 3060000 / 3304000 / 2400230 by the quipu theorem; n = 8: 0 conflicts), E-292 and E-265 (q_S separation of E-259 pairs, n = 9 group found there with 1 pair: this entry adds the full n = 8, 9 scan and the three-ball step test), E-233/E-239/E-244 (the J != 0 steps, joined by chains), E-300 and the round-114 controls (this finishes them). q_S on the F-010 pair is therefore not new in itself; new: n = 8 has no pair, the split is exactly 18 + 18 and matches E-238's naming, q_S is invariant along all 14.8 k explored steps (expected for J = 0, not evidence).

## Code changed

None in `quivermutation/`. New scripts only, listed above; no tests needed.

## Next

H-entry proposed (converted): "from a 3060000-half algebra, does any key-kept (J = 0 or J != 0) step ever produce a 3304000-half Cartan invariant?" Instrument: `toolsmith_qsball.py` at depth 6 for 3060000 (~19 000 parents; split via the `moves` scan across 4 cores, needs a checkpointed build as 087), pre-registered refutation = one step with different q_S. Experimentalist: run it overnight; theorist: it is only informative for T10 if such a step exists, so also search for it directly by Cartan-class hits (q_S of the other half) at the end of long random walks.
