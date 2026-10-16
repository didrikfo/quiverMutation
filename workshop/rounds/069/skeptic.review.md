# Review of workshop/rounds/069/skeptic.md

referee: theorist · round: 069
verdict: minor revision

## Reproduction

Not re-run: the walk (3 segments, ~27 min total; over the 10 min limit). Smaller check: ran the author's `skeptic_c1_13_prefix.py` (~30 s) and `skeptic_c1_13_rebuild.py /tmp/skeptic_c1_13_cls1.pkl` (~20 s) on the author's checkpoint in /tmp. Output matched the table exactly: first-seen expansions 7369 ... 22450 and parallel-parent flags for all 13; prefix counts 64/10/16 at 20000, 67/13/19 at 22450 and 24032, 68/14/20 at 24033; hand rebuild: 8 non-parallel parents give gate True, tiltingPlus False, key same, child == hand child, J as listed; 5 parallel ones skipped. So the post-processing is verified, but the walk itself is trusted from the checkpoint (I did not check the BFS loop independently). Also read `rounds/043/skeptic_c2.py` (the E-147 script): level-by-level BFS, canonicalKey dedup, key guard, time budget; same walk shape as claimed.

## True?

No error found in what was checked. Weak points:
1. The identification "E-147's 13 = the first 13 D = 0 steps" rests on count coincidence (67 and 13 at once, in the window 22450-24032 expansions). E-147 recorded no expansion count and no step list, so step-by-step identity cannot be shown; only "the counts agree at one point of the same walk". 67 alone is a weak coincidence; 67 together with 13 and a window of 1583 expansions is decent but is one point. Class 2 (E-147: 64 and 9) was not matched the same way, though the title/claim does not say it was.
2. "12 strict H1/H2" in E-147 versus the 13 here: not compared.
3. "Round 43 ran at 500 s" being 9-10 min here is inferred from this machine's speed; not measurable on E-147's.

## New?

Partly. Containment is E-179's open item ("class-1 containment open", E-147 13 vs E-151 16); this closes it conditionally. The cap-dependence of counts is already flagged in E-151's scope ("Step counts, not distinct steps (E-147: 13 / 9 distinct)"), but not the 10 + 3 split or the prefix identification. The parallel-parent fact for class 1 is in E-147 ("parallel arrows") and E-154 (8 of 16 parents parallel; here 5 of 13 Cartan-distinct steps, consistent, different counting unit). Hom(T,T[-1]) != 0 for all key-keepers is E-163/E-128, as the author says.

## Evidenced?

Mostly. The table is specific (expansion, depth, v, J, parent). Missing: the exact E-147 convention for "distinct" is reconciled only by assertion ((C_A, C_B, v) vs (canonicalKey, v)); the 16-vs-13 mapping (13 Cartan steps = 19 pairs, 10 steps = 16 pairs) is given and checks out in my run. Checkpoint and `_cls1.pkl` are in /tmp, uncommitted, so the 20 s/30 s scripts are not re-runnable by a third party without 27 min; the `.out` files are kept, which helps.

## Scope

Title says "E-147's 13 ... are the first 13 ... D = 0 steps of the E-151 BFS" as fact. Narrowed wording: "the first 13 D = 0 steps of one rebuilt walk, which reaches E-147's counts (67 J != 0, 13 D = 0) at 22450-24032 expansions; E-147's own step list is not preserved so identity is by counts". "Cap artefacts" for 13/16/9 is supported for class 1 only; class 2 not rerun (stated).

## Required for acceptance

1. Reword title/claim as above (identity by count agreement at one window, not step-for-step).
2. State that the walk was not independently re-run by the referee and that the checkpoint is not committed; add a one-line sha or size/mtime so the .out files can be tied to it.
3. State the E-147 "12 strict H1/H2" figure and whether the 13 here contain 12 such steps (if cheap, in the existing scripts) or mark [next round].
4. Label the "Theorist" pattern remark (|supp J| = 1 in 12 of 13, Hom(T,T[-1]) dim 1 or 2) as an observation at one sample, not a regularity: the 111-pair tally (94/14/3) is the only support.
