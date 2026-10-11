# The n = 8 class 2 loose end is a bug in step 7 of the mutation rewrite (`reduceAgainstPivots` is not a normal form), not a gap in `tiltingPlus`; the rejections of E-086 are confirmed by the Cartan congruence

author: theorist · round: 015 · kind: result
thread: T5 · bears on: H-015, E-086, E-080, E-059, F-038, R-005

## Claim

(1) For the replayed non-tilting parents of E-086 (n = 7 class 0, n = 8 class 2, n = 9 class 0; 5 + 1 + 5 recorded
rejections; the n = 8 and n = 9 files list the others the same way) the Cartan congruence, computed without
`tiltingPlus`, fails in every one: the predicted matrix R C R^T has an entry -1 (impossible for a Cartan matrix),
and the child's key differs. So `tiltingPlus` False and "Cartan congruence fails" agree on all 11 rejecting parents
replayed; nothing new about E-086's main claim.
(2) The n = 8 class 2 loose end (10 steps: gate True, `tiltingPlus` True, key moves) is explained. In all 10 the
Cartan congruence fails too (the child has dimension 2 where R C R^T predicts 1 at one or two entries), but `tiltingPlus`
is right: it is vacuous or true, and R C R^T is the Cartan matrix of End(T). The wrong thing is the **child**: step 7
of `procedure.mutateAtVertex` misses a relation out of the mutated vertex because `arrowPaths.reduceAgainstPivots`
reduces only while the leading term is a pivot, so its "residue" is not a normal form (two congruent elements get
different residues), and `_kernelOverIdeal` solves its linear system on those residues. With `reduceAgainstPivots`
replaced by a full reduction (monkeypatch, `theorist_fix.py`), all 10 steps give a child with congruent Cartan matrix
and the parent's key (parents identical to the recorded ones); the 3 rejection sets stay rejected.
Parallel arrows are incidental: 6 of the 10 lines have one pair, 4 have none (the 10 lines are 7 distinct parents).
It does **not** claim: that this is the only defect of the rewrite, how often the bug fires (only the 7 distinct
parents of this one file were examined; the other 13 classes of E-086 recorded no such step), or that E-086's counts
(rejections, 1.29e6 guard-admitted steps) are unchanged: the BFS ran with the buggy rewrite, so child algebras
that depended on a missing relation were possibly wrong and are not re-run. No library change was made.

## Evidence

Parent 1 of the 10 (n = 8, class 2, start index 6, path 1,1,2,2,4,5,6; mutate at 3): arrows 1>2, 1>5, 2>7, 3>1, 4>5,
5>6, 5>7, 6>7, 7>8; relations 3>1>2 = 0, 3>1>5>6 = 0, 1>5>7 = 0, 1>2>7 + 1>5>6>7 = 0 (mod 1>5>7), 4>5>6>7 = 0 mod
4>5>7. Vertex 3 has no arrow in, so `tiltingPlus` is vacuous and the Okuyama-Rickard complex is tilting (not
independently verified beyond Ladkani 2.3(c)); R C R^T is then End(T)'s Cartan matrix: (7,3) = 1. The rewrite gives
2: it has arrows 3*>2, 3*>6 (from the relations 3>1>2, 3>1>5>6) and no relation `3*>2>7 + 3*>6>7`. Step 7's
kernel at target 7 is spanned by that relation: the shadows are 1>2>7 and -1>5>6>7, and 1>2>7 = -1>5>6>7 in A.
`reduceAgainstPivots` gives 1>2>7 the residue `-1567 + 157'` (157' = 1>5>7, itself a pivot) and -1567 the residue
`-1567`, so the system has no solution and the kernel is empty (`theorist_step7.py` prints the trace). The echelon
basis from `idealBasis` is not reduced (the 1>2>7 row contains the pivot 157'), and the loop stops at the first
non-pivot head. After the fix the extra relation appears and everything matches.

| set | lines | gate | tiltingPlus | cong (reduced and unreduced child) | key equal | after full-reduction patch |
|---|---|---|---|---|---|---|
| rejections n = 7 c0, 8 c2, 9 c0 | 5 + 1 + 5 | True | False | False, predicted entry -1 | no | unchanged (False; n = 8 c0 also run: 4 lines) |
| n = 8 c2 'M' lines | 10 | True | True | False, child 2 vs predicted 1 | no | cong True, key equal, 10 of 10 |

Not a retraction of E-059 / E-086's `tiltingPlus` result: this is the first case where `tiltingPlus` and the congruence
disagree, and the congruence is the one that is right about End(T) while the rewritten child is the faulty object.
The same non-canonical residue is used by `_isForcedByNearer` (decided by `reduceAgainstPivots`; zero test, so safe)
and by `tiltingPlus`'s `rank` of residues of p·beta (not safe when a path p lies in the ideal; not tested here).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/015/theorist_cartan.py 8 2   # 30 s; also "7 0", "9 0" (outputs saved)
timeout 10m .venv/bin/python workshop/rounds/015/theorist_detail.py       # 30 s, where the congruence fails
timeout 10m .venv/bin/python workshop/rounds/015/theorist_step7.py        # 5 s, trace of step 7 for parent 1
timeout 10m .venv/bin/python workshop/rounds/015/theorist_fix.py 8 2      # 30 s, with full reduction; also "7 0", "8 0"
.venv/bin/python workshop/rounds/015/theorist_child.py 1                  # 5 s, parent, raw and reduced child
```
Outputs: `theorist_cartan_n{7c0,8c2,9c0}.txt`, `theorist_fix_n8c2.txt`, `theorist_detail_out.txt`, `theorist_child_1.txt`.

## Prior record

E-080 (the congruence as a check for the length-4 square), E-086 (the loose end, "not pursued"), E-059 ("0 disagreements"
between `tiltingPlus` and the congruence, depth-limited). `grep` of `research/` for `reduceAgainstPivots`, residue
and normal form finds nothing about this defect. So the defect is new, as far as the record goes.

## Code changed

None in the library. New scripts `theorist_{cartan,detail,step7,fix,child}.py` in this folder. `theorist_fix.py` monkeypatches only.
Suggested library change (for the toolsmith, not made): in `arrowPaths.reduceAgainstPivots`, keep reducing past a
non-pivot leading term (my `fullReduce`), or make `idealBasis` fully reduced. Ran
`pytest -q tests -k "arrowPaths or procedure" -m "not slow"` on the unchanged tree (24 passed, 1 xfailed).

## Next

- Toolsmith: apply the fix, add a unit test on this parent (child Cartan = R C R^T; 6 vertices are enough: the 3>1 / 1>2 / 1>5>6 pattern), run the fast mutation tests, then ask whether anything recorded depended on the old behaviour (E-086 n = 8 class 2 BFS size 54,326 will move).
- Experimentalist: re-run the n = 8 class 2 walk with the fix; look for any step with congruence False and `tiltingPlus` True (after the fix there should be none; if some remain, `tiltingPlus` really is incomplete).
- Chair: R-005 gains a possible second source of "rewrite is not a derived equivalence": a missing relation, not only an inadmissible vertex. Worth a record line once the fix is checked.
