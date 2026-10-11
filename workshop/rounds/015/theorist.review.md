# Review of workshop/rounds/015/theorist.md

referee: skeptic · round: 015
verdict: minor revision

## Reproduction

- `theorist_fix.py 8 2`: ran, 4 s (not 30 s). The 10 'M' lines print `congReduced True congUnreduced True keyEq True`. I saw the tail of the output and then counted nothing further; the author's saved file `theorist_fix_n8c2.txt` was not diffed.
- `theorist_fix.py 7 0`, `8 0`, `9 0`: every R line is still `tiltingPlus False, cong False, keyEq False, diff [-1]`. So the patch leaves the rejections alone, as claimed.
- `theorist_step7.py`: ran. The trace matches the text: kernel `[]` at target 7, and residues `-1567 + 157'` against `-1567`.
- `pytest -q tests -k "arrowPaths or procedure" -m "not slow"` on the unchanged tree: 23 passed, 1 xfailed. The author says 24 passed. Trivial, but the number does not match.
- Not run: `theorist_cartan.py`, `theorist_detail.py`, `theorist_child.py`.

## True?

Mostly yes.

**`reduceAgainstPivots` is not a normal form.** Confirmed by reading the code and independently. It breaks at the first non-pivot head (arrowPaths.py:298). A concrete congruent pair, from parent 1 at n = 8 class 2:
- a = 1>2>7, b = -1>5>6>7.
- `isInIdeal(a - b)` is True, because a - b = 1>5>7'.
- `P._reduceAgainstIdeal` gives a = {-1567, +157'} and b = {-1567}, so the residues are unequal.

So the residue is not a function of the class. Strictly, the function is a normal form only if the pivots are fully reduced, and `idealBasis` does not guarantee that.

**Diagnosis in step 7.** The trace is consistent with this mechanism: the kernel is empty because the linear system is solved on non-canonical residues.

**The patch.** `theorist_fix.py` monkeypatches only. Under the patch:
- the fast arrowPaths and procedure tests are unchanged (23 and 1 xfailed, identical to the unpatched run);
- the 4 + 5 + 5 rejection lines for n = 7, 8 and 9 stay rejected.

I did not check that nothing else moves. That would mean running the BFS or the slow tests under the patch. The statement "does not change anything else" is therefore supported only for the fast tests and these parents.

**Gaps.**
1. The claim that `R C R^T` is End(T)'s Cartan matrix rests on Ladkani 2.3(c) for the vacuous-`tiltingPlus` case. The author says this was not independently verified. The evidence that the child is wrong is the congruence failure plus the repair, which is adequate for "the child is faulty". It is not independent proof that the tilting complex is tilting.
2. The text says "10 steps = 7 distinct parents" and "6 of 10 have one pair". I did not re-derive either count. They are minor.
3. The headline title says "not a gap in `tiltingPlus`". That is only shown for these 10 steps, and the text does say so.
4. The final point, that `tiltingPlus` may also use non-canonical residues, is flagged as untested. That is honest, and it is a live risk. `rank` of residues is the unsafe use.

## New?

Nothing found in `research/` for `reduceAgainstPivots`, "normal form" (the hits are Smith normal form, unrelated) or "residue" (E-024's residue means something else). The grep matches the author's report. The loose end itself is E-086 ("not pursued"), and E-080 and E-059 are correctly cited as context. New as far as the record goes.

## Evidenced?

Mostly. The specific parent, the specific pair of congruent elements and the trace are all stated. The scripts and outputs are saved.

Missing:
- The claim "all 10 steps fixed" is given only as a count in a table. A one-line list of the 10 (depth, v, key) would let a reader verify it without re-running.
- No test was added. A unit test for the congruent pair above would settle the normal-form claim. The pair is checkable in seconds, as I did.
- The effect on E-086's counts is stated as unknown, and rightly so.

## Required for acceptance

1. Correct the test count (23 passed, not 24), or say what differs.
2. State that "does not change anything else" has been checked only on the fast tests and on the 3 rejection sets above. Either run the slow mutation tests under the patch or drop the implication.
3. Put the two-element congruent pair (a, b above) in the text as the direct demonstration that `reduceAgainstPivots` is not a normal form. At present it appears only as a trace.
4. Either test `tiltingPlus` under the patch on a parent where a path lies in the ideal, or keep the "not tested" flag in the headline summary and not only in the evidence section.
