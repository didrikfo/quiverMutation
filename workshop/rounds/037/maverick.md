# The I1/I2 split of E-125 is a single relation sitting off-centre: deleting from its shorter free side and from its longer side land in different classes, so "K >= K0 holds" is a small-n artefact and fails first at n = 13 (K >= 4), 15 (K >= 5)

author: maverick · round: 037 · kind: result
thread: S-1 · bears on: E-118, E-125, H-020 (not derived from the rule table; see Claim)

## Claim

Speculation level: tested on small cases, plus a computed (key) lower bound.

1. **Mirror redo.** Reading every K >= 3 end of the failing n = 11 class (key 1,1,0,-1,-2,-3,-3,-2,-1,0,1,1, 1305 LNAs) as a *tail* word, with head ends passed through `freeMoves.mirrorRow` (not digit reversal): the image is still a function of (K, word), 47 keys, 0 with two images, head and tail never disagree. E-125's "image is a function" survives. Its 12 common words shrink to 8 (3, 23, 203, 223, 2003, 2023, 2203, 2223); the other four of E-125 were artefacts of digit reversal. The "last letter 3" split is a fit, not a rule: at K = 3, 16 of 21 words ending in 3 go to I1, 5 go to I2 (20333, 22333, 2333, 333, 2334...), and 1 word not ending in 3 (`7`) goes to I1.
2. **The real invariant.** Length-2 relations are free (derived-class level), so drop them and measure the free run K_eff beyond the last non-2 relation. Over all 94 ends there are only 12 keys (K_eff, stripped relation list), 0 with two images:
   - I1 exactly when the stripped core is ONE relation, `3` or `7`, with K_eff = 3.
   - Every other case (the same single `3` or `6` at K_eff = 4; every core with two or more non-2 relations) goes to I2.
   Words like `302`, `2302` are the single 3 with K_eff = 4, which is why they sit with I2 though they end in 2. "Last letter 3" is a proxy for "stripped core = lone 3 with K_eff = 3".
3. **Why (derivation for the single relation, not from the rule table).** The VERIFIED_MOVES table has no rule that moves a lone 3 (only 2-relation, 4+4-pairs, and 3-with-neighbour windows), so a lone 3 at (h free at head, K free at tail) is its own placement, and the key depends on the unordered pair {h, K} (computed n = 9..12, a = 3, 4, 7). At n = 11 the lone-3 class is {(4,3), (3,4)} (a mirror pair). Deleting a free vertex from the tail of (4,3) gives (4,2); from the head gives (3,3). At n = 10 these are different classes: I1 = {(2,4),(4,2)} and I2 = {(3,3)}; I checked that their E-115 labels are exactly the two image keys of E-118/E-125. K = 3 is the shorter free side (delete -> off-centre, I1); K = 4 is the longer side (delete -> central, I2). Not "room to move": the same dichotomy arises for any h != K, and for h = K both deletions agree.
4. **Consequences (key-level, so sound for "differs").** A class containing a lone 3 at (h,K), h != K, both >= K0, has two ends with K >= K0 whose images {h,K-1} and {h-1,K} have different keys. First such n: K0 = 3 at n = 11 (this is E-118's failure), K0 = 4 at n = 13 (h,K) = (4,5), K0 = 5 at n = 15. So "K >= 4 holds" (10/10 at n = 11) is expected to fail at n = 13; the threshold is n-dependent, not a property of the core. Compatibility (S-1 Q1): the image is a function of the class if one deletes from the *shorter* free side at both ends of the class (mirror-aware choice), at least for lone relations.
Not claimed: that other cores obey this; that n = 12 was run (it was not); that K = 2 failures (9 classes at n = 11) are the same mechanism (untested); that the 94 end count equals E-125's 82 (see Prior record).

## Evidence

- Mirror redo: 94 end deletions, 47 (K, mirrored tail word) keys, 0 conflicts. K = 3 I1 words (17): 3, 23, 203, 223, 2003, 2023, 2203, 2223, 20003, 20023, 20203, 20223, 22003, 22023, 22203, 22223, 7 (the `7` is the lone 7). K = 4: 9 words, all I2.
- Normal form table (K_eff, stripped relations as (offset, length)): (3,{3}) I1; (3,{7}) I1; (4,{3}) I2; (4,{6}) I2; (3, {3,6}), (3,{6,6}), (3,{6,5 at 2}), (3,{6, 4 at 3}), (3,{6, 3 at 4}), (3,{3,3,3}), (3,{3,3,4}), (3,{3,3,5}) all I2.
- Lone 3 keys at n = 10 by (h,K): (0,6) (1,5) (2,4)=(4,2) (3,3) all distinct; I1 key = (2,4), I2 key = (3,3). Lone 7 at n = 11 (0,3): deleting the tail gives (0,2), key equals I1 (key coincidence, not class identity).
- Predicted first failing n by threshold (computed from keys): K0 = 3: 11; K0 = 4: 13; K0 = 5: 15.
- n = 12 sizing: not done. By the same keys the n = 12 lone-3 classes are {(3,5),(5,3)} (fails at K >= 3) and (4,4) (compatible); so I predict K >= 4 still holds at n = 12 from this mechanism and K = 3 fails. That is a prediction, not a run; an orbit check would still be needed to say the source class is one derived class (key only).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/037/maverick_endtable_mirror.py   # about 90 s
.venv/bin/python workshop/rounds/037/maverick_single.py                        # a few s
```

## Prior record

E-125 gives the (K, word) function and the 12 same-core words (head by digit reversal, flagged by its referee); E-118 the class, the images, K >= 4 holds. Not in `research/`: the mirrorRow redo (E-125's referee asked for it), the free-move normal form, the lone-relation reading, the predicted failures at n = 13, 15. Unexplained: E-125 says 82 ends (66 + 16), my loop over the same class counts 94 ends; I did not find why (E-125 may count distinct oriented ends). The 47 vs 77 key count differs for the same reason, so do not compare the tables line by line.

## Code changed

None in the library. New scripts `workshop/rounds/037/maverick_endtable_mirror.py`, `maverick_single.py`; they import r027/r030 helpers. No tests touched.

## Next

- experimentalist/toolsmith: test the n = 13 prediction directly: the lone 3 at (4,5) and (5,4), K >= 4 ends (orbit check that both are one derived class; their images at n = 12 have different keys by `maverick_single.py`).
- theorist: derive the {h,K} dependence of a lone relation's key from H-020 (head/tail for a core with no moves); the rule-table derivation proper is not done.
- Reframing for S-1 Q1: the class of a lone relation is the mirror-orbit {h,K}; compatible deletion = shorter side.
