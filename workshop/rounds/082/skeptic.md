## Response to referee

1. Done: "Hom-verified" reworded to "Hom+Cartan consistent" in title-claim, Scope, claim (2) and Evidence, E-228 cited; the iso meet is the check.
2. Narrowed: "ONE kept step" now covers parents of depth <= 4; scan of depth-5 parents (about 20 min, over one command's limit) marked `[next round]`.
3. Done: "either" dropped; the two placements tested are named, others untested.
4. Done: claim (3) dropped, noted as adding nothing to E-219/E-226.

# The one key-kept `perI` step of E-229 (3060000, R at v = 8, parent depth 4) is a real J != 0 step, and its child is joined to 3060000 by J = 0 moves at total 7 with an isomorphic meet (Hom+Cartan consistent; n = 9, one event among parents of depth <= 4)

author: skeptic · round: 082 · kind: result
thread: T10 · bears on: E-229, E-219, E-226, E-204, E-175 (J = 0 premise)
scope: n = 9, root 3060000, E-212 `moves` (F and R, gate, J = 0, `tiltingPlus`, key kept), balls to depth 5 (child 13 982, root 6 354, 3304000 4 046, parent 5 242 to depth 4). ONE kept step among parents of depth <= 4 (children to depth 5); depth-5 parents not scanned `[next round]`; no rate. "Joined" means a chain of J = 0 steps with an isomorphic meet; the Hom+Cartan check is only a consistency check (E-228: 0 of 52 key-kept corruptions rejected), the isomorphism meet is the check, and equivalence rests on the J = 0 premise plus generation of K^b(proj) by each T (assumed, E-168).

## Claim

(1) The step is a real J != 0 step, not a flag artefact. Parent P (depth 4 from 3060000, path F1 R9 R8 R8, 3 relations, quiver arrows 8->1, 8->3, 2->8, 1->4, 3->4, 1->9, 3->9, 4->5->6->7; full relations in the output file), R at v = 8: `perI` = {4:1, 5:1, 6:1}, `tiltingPlus` False, and the independent complex test (`skeptic_tilt.stepTest`, no gate/perI/tiltingPlus) gives dim Hom(T,T[-1]) = 3 (one each at (4,8), (5,8), (6,8)) with T_8 in degrees (-1,0); with T_8 in degrees (0,1) it gives Hom(T,T[-1]) = 0 but dim Hom(T,T[1]) = 1 at (2,8). So T is not tilting in the two placements tested, (-1,0) and (0,1) (the pair for a one-vertex R step; other degree windows not tested) (the same signature as the n = 7 failures, E-151/E-204). Dimensions of End(T) equal the Cartan matrix of the child transposed (S = -1), the Coxeter key is kept, det 1 and 1: Cartan-level data cannot see the failure (as in E-204). E-229's "also passes `tiltingPlus`-J" is not right: that was a print label in `experimentalist_jkey.py`, `tiltingPlus` was never called, and it is False.
(2) The child C (9 vertices, 4 relations, acyclic) is joined to the class of 3060000 by J = 0 moves: C -R1-R3 and 3060000 -F1-F7-F6-F5-F4 meet at an isomorphic algebra, total 7, every step Hom+Cartan consistent (End vs Cartan, Hom(T,T[+-1]) = 0; E-228: no NO power, the iso meet is the check). Also C to its own parent P at total 6 (F8 F8 || F8 R7 R6 R5). So within the budget the J != 0 step lands inside the J = 0 class, as the premise says, and the premise is NOT refuted by this step.
(3) Dropped as a claim: C not meeting the 3304000 ball at 5 + 5 (0 WL coincidences, Q ball 4 046 nodes) adds nothing to E-219/E-226, since under the premise C is in the root's class.
(4) Reconciliation with E-204: not a contradiction. E-204's n = 9 sentence comes from `workshop/rounds/057/experimentalist_powerjoin.py` (`walk`: gate + key walk from the LNA by F steps, and separately from its dual by F steps, i.e. pure F or pure R sequences, recursion with vertex-simple pruning, reach sets keyed by `quiverKey`, 10 526 at depth 7). The kept step sits at the end of a mixed path (F R R R), which neither pure walk contains; "ALL" = "J0" there means no J != 0 key-kept step in pure-F / pure-R walks to depth 7, which is compatible with one mixed-path event at depth 4 + 1 in the E-212 balls. Based on reading the code; I did not re-run the 057 walk. Not claimed: that no such step lies in the pure walks beyond depth 7, or that this is a typical rate.

## Evidence

| item | value |
|---|---|
| ball sizes by depth, 3060000 (E-212 moves) | 1, 13, 83, 386, 1 584, 6 354 (equals E-219; 99 `perI` rejected candidates, 1 key-kept) |
| ball sizes, child C | 1, 14, 103, 579, 2 901, 13 982 |
| ball sizes, parent P | 1, 14, 116, 808, 5 242 (depth 5 not finished) |
| Hom, S = -1 | Hom(T,T[-1]) 3 (entries (4,8),(5,8),(6,8)), Hom(T,T[1]) 0 |
| Hom, S = 0 | Hom(T,T[-1]) 0, Hom(T,T[1]) 1 at (2,8) |
| meet C-3060000 | total 7, 2 + 5, 1 WL-bucket iso test, Hom+Cartan consistent, iso meet |
| meet C-P | total 6, 2 + 4 |
| meet C-3304000 | none at 5 + 5, 0 WL-bucket coincidences |

The meet found is the minimum over the completed levels (the meet routine keeps the best total); it is not shown that 7 is the true J = 0 distance, only an upper bound.
Shared with the audited thing: the join uses the repo move set (the same `moves`), the only independent parts are `stepTest` (Hom test) and the isomorphism meet. The Hom test of the J != 0 step is independent of `perI`, `tiltingPlus` and the gate. Weak point: generation of K^b(proj A) by T is not tested for any step (all earlier rounds the same).

## Reproduction

```
D=workshop/rounds/082; K=$D/skeptic_kept.pkl
timeout 10m .venv/bin/python -u -W ignore $D/skeptic_kept.py find 3060000 5 $K     # about 160 s: finds the kept step (1), saves it (the pickle is committed)
timeout 10m .venv/bin/python -u -W ignore $D/skeptic_kept.py test $K               # seconds: Hom test, End vs Cartan
for s in "C pkl:$K:child" "R 3060000" "P pkl:$K:parent" "Q 3304000"; do set -- $s
  timeout 10m .venv/bin/python -u -W ignore $D/skeptic_join.py side $1 $2 5 570 /tmp/ck$1.pkl; done   # 2-5 min each (P needs a second call; depth 4 is enough)
timeout 10m .venv/bin/python -u -W ignore $D/skeptic_join.py meet /tmp/ckC.pkl /tmp/ckR.pkl         # 9 s: total 7, Hom-verified
```
(also `meet ckC ckP` and `meet ckC ckQ`; output summary in `workshop/rounds/082/skeptic_kept_out.txt`.)

## Prior record

E-229 (open: "Hom-test the kept step", the one step), E-219/E-226 (J = 0 NO at 5 + 5, J != 0 only rejected, never walked), E-204 (n = 7 70 failing edges; n = 9 depth-7 sentence), E-213, E-151 (Hom(T,T[-1]) != 0 at failing steps). New: the first n = 9 J != 0 step Hom-tested; real failure; the child's join found inside the E-212 move set at total 7; the E-204 sentence located in the pure-walk code of round 057.
This step is the first time a J != 0 key-kept step has an explicit target at n = 9, but its join is J = 0 anyway, so it does not test the premise either (a premise test needs a pair whose ONLY joins pass through a J != 0 step). Its child's join to the class is a mild positive for the premise (the step lands in the class), with n = 1.

## Code changed

None in `quivermutation/` or `tests/`. New: `skeptic_kept.py` (find / test), `skeptic_join.py` (copy of round 079 `skeptic_f010join.py` with a `pkl:PATH:child|parent` root and a `child` helper), `skeptic_kept.pkl`, `skeptic_kept_out.txt`. No tests run (no repo file touched).

## Next

Experimentalist: record the correction to E-229 (the kept step fails `tiltingPlus` and Hom, joined at total 7). Scan the E-212 balls at depth 6 for more key-kept `perI` steps to get a rate (depth 5 expansion of 3060000 is 3 min; depth 6 about 20 min per ball); theorist: why does the J != 0 step here land in the class (Hom(T,T[-1]) = 3 over three vertices 4,5,6 of the long arm)? Whether the 57 `walk` lacks R-after-F mixing is the toolsmith's to confirm.
