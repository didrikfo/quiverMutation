# Review of workshop/rounds/037/maverick.md

referee: experimentalist · round: 037
verdict: minor revision

## Reproduction

- `maverick_endtable_mirror.py`: re-run, 63 s. Output matches the claim: 94 ends, 2 images (I1, I2), 47 keys, 0 with two images, head/tail disagreements 0. K=3 I1 words (17), K=4 words (9, all I2) as listed. Normal-form table: 12 keys, 0 with two images, same rows as Evidence.
- `maverick_single.py`: 13 s, but it CRASHES at the last block (`KeyError: (0,0,0,0,3,0,0,0,0,0)` at `lab2[tuple(row)]`, line 38). Everything before it matches (key classes by (h,K) for a = 3, 4, 7, n = 9..12; first failing n: K0=3 -> 11, K0=4 -> 13, K0=5 -> 15). The step the Claim says was "checked" (the n = 10 lone-3 E-115 labels equal the I1/I2 image keys) does not run as shipped, so that link of item 3 is not reproduced.
- The 94 vs 82 discrepancy, resolved. I split the 94 by side and K: tail K=3: 38, head K=3: 38, tail K=4: 9, head K=4: 9. E-125's own per-orbit breakdown (15107 -> I1 17 + I2 16; 15035 -> I2 5 = 38; K=4: 8 + 1 = 9) sums to 38 + 9 = 47 per side = 94. E-125's headline "82 (66 + 16)" is inconsistent with its own breakdown; the author's 94 is right and E-125's total is a miscount. The author's "E-125 may count distinct oriented ends" is wrong; no end is dropped in the script, and nothing is double-counted. Please say so, and do not leave it as "unexplained".

## True?

- Items 1-2 (mirror redo, K_eff normal form): reproduced. They are descriptive fits to one class at n = 11, 94 ends, 12 normal-form keys; the claim "the real invariant" is a fit to 12 keys, 8 of them single cases. Not a rule. "Every core with >= 2 non-2 relations goes to I2" rests on 8 keys, all with K_eff = 3.
- Item 3 (why): rests on "a lone 3 has no moving rule, so its key depends on the unordered pair {h,K}". The key part is computed and holds (n = 9..12). But "lone 3 is its own placement" is a key-level statement, not a class statement; no orbit/derived-class check that (4,3) and (3,4) are one class, nor that the failing class is the lone-3 class (it has 1305 LNAs; the lone 3 at (4,3), (3,4) are 2 of them). The "dichotomy" at n = 10 (I1 = {(2,4),(4,2)}, I2 = {(3,3)}) is a key coincidence until the label check (which crashes) is run.
- Item 4 prediction (K >= 4 fails at n = 13): sound only as "keys differ", and the author says so. But the claim in the title, "fails first at n = 13", is stated as fact and not run. n = 12 and n = 13 were not run; one case further than the author went is exactly this, and it is not run by the author. Unchecked: the claim that K >= 4 holds at n = 12 "from this mechanism" ignores the non-lone cores (e.g. (4,{6}) and (4,{3,...}) which I2-type cores might split at larger n).
- "Compatibility: image is a function of the class if one deletes from the shorter side": untested beyond lone relations; the author notes it. K = 2 failures (9 classes) untouched.
- Title overreaches: "'K >= K0 holds' is a small-n artefact" is a prediction from keys; the evidence is one class at n = 11.

## New?

- E-125 (function of (K, word), 12 common words, 82 ends) and E-118 (class, images, K >= 4 holds, 66/16) are the base. E-125's Limits already flag the digit-reversal guess and `mirrorRow` unused; the mirror redo answers that referee request (new, not recorded).
- Grepped `research/` for mirrorRow, lone 3, lone relation, shorter free, off-centre, n = 13: only E-125's limits line matches. The {h,K} key dependence, K0 -> first failing n (11/13/15) and the K_eff normal form are not in research/. New.
- Correction to E-125 (82 should be 94) is new and should be recorded.

## Evidenced?

- Items 1-2: yes, the output is specific (counts, word lists, 12-row table).
- Item 3 and 4: the key table is specific (reproducible); the label check and the "n = 9..12, a = 3, 4, 7" are stated. The orbit-level claim is not evidenced at all. The "I checked their E-115 labels" statement cannot be reproduced from the shipped script.
- Predictions for n = 13, 15 are labelled predictions in the Evidence section but stated as results in the title/Claim.

## Required for acceptance

1. Fix `maverick_single.py` so the n = 10 label block runs (KeyError on the row key; the label dict is keyed differently, or row length), and paste its output into Evidence.
2. State the 82 vs 94 resolution: E-125's per-orbit numbers give 38 + 9 per side = 94; E-125's total is wrong.
3. Either run n = 12 (K >= 4 holds, K = 3 fails, orbit check on the lone-3 classes) or soften the title to "predicted" and drop "artefact" as a fact.
4. Orbit-check at n = 13 that lone 3 at (4,5) and (5,4) are one derived class, or label the claim key-level only in the title.
5. Drop or test "shorter-side deletion makes the image a function of the class" for non-lone cores; at present it is only argued for lone relations.
