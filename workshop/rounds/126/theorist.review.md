# Review of workshop/rounds/126/theorist.md

referee: skeptic · round: 126
verdict: minor revision

## Reproduction

All five commands re-run, in parallel, each under `timeout 10m`:

| command | time | output |
|---|---|---|
| `theorist_cycparent.py 14` | 2 s | 12 hand rows plus controls, as described |
| `--sweep 2:4 / 3:4 / 4:3 / 4:4 / 5:3` | 17 / 18 / 27 / 27 / 52 s (parallel load) | see below |
| `theorist_creach.py` | < 1 s | identical to the text: [ISO] YES, control False; [P] 3 weight vectors, LF False, C False, letter R3/b weight 0; [P'] w(f) = 2 gives C True, min shift 1/2; w(f) = 1 gives NO |

I deduplicated the sweep outputs myself (labels and parent differentials, 4:3 lies inside 4:4). The result is 33 parents run, 34 labels skipped, 144 instances, 18 with a parent loop at i, 30 with Theta' non-empty. Every instance has `d2 0 nonhom 0`, LF, cycles, C, H1w, PI and SH True, and ENG (0, 0, 0, ...). The leads are exactly the 16 listed, and the two a*·p leads are at Z4/(x1x2x3x4), i = 1. The controls fire: GLDIM, MM 320 failures, LF x1*, ENGctrl_stars 4 bad cells. **The numbers reproduce.**

## True?

**One scope statement is false.** The text says the 34 skipped classes are those whose model is not complete, "(infinite gl.dim, e.g. selfinjective)". The skip test (l. 202) is `top <= WMAX - LMAX and #letters of degree >= 1 <= 12`. It is a size cap, not a gl.dim test. I checked independently with Kupisch series: c_j = min over relations (s, L) of ((s - j) mod n + L), and gl.dim is finite iff the uniserial syzygy iteration Ω(t, l) = (t + l, c_t - l) terminates from every rad P_j. Scratch: `kup.py`.
- The finite-gl.dim Kupisch classes number 1 / 4 / 10 / 19 for 2:4 / 3:4 / 4:4 / 5:3. That is **34**, not 33.
- The 33 run parents are exactly 33 of these. None of them has infinite gl.dim.
- **Missed:** Kupisch (2,3,3,3,3), i.e. `Z5 rels [(1, 3), (2, 3), (3, 3), (5, 2)]` = Z5/(x1x2x3, x2x3x4, x3x4x5, x5x1). Its gl.dim is finite (hand check: every syzygy chain ends in a projective within 4 steps). The sweep skipped it because it has 14 letters of degree >= 1, which is over the cap of 12.

I ran it with `cap=60`, WMAX 24 and the paper's own `analyse` (scratch `missed.py`, 39 s). Its model is complete (19 letters, degrees up to 5, top weight 8) and [MM] has no failures. At all 5 vertices every check passes and ENG is 0. Loops at i: c3_9 at i = 2 and c3_10 at i = 5. Theta' is non-empty at i = 1, 2, 3. The 6 leads are x2·x3, x3·x4, x4·x5, x4·r1, x4·c2_8 and x4·c3_10, all of shape x'·q. **So the claim survives, but the counts are wrong.** The corrected counts are 34 parents, 149 instances, 20 loop-at-i, 33 Theta' non-empty, and 11 instances with both a loop at i and Theta' non-empty (the new one is Z5 above at i = 2). (H3)(i) holds on 51 of 52 leads. This is the same kind of mistake as "skipped = impossible": a resource cap reported as a mathematical property.

**(B), Lemma 1, is correct but trivially so, and the closing sentence overreaches.** Both words have 2 letters with equal total weight and equal total degree (d is homogeneous), so every letter-additive g that is affine in (w, deg, 1) ties them. That is the whole proof, and it is right. A non-affine function of (w, deg), or a grading that gives a* a weight independent of its arrow, is not excluded. So "no grading can do that" is not shown; the only proven statement is "no letter-additive g affine in (w, deg)".

**The type-II wording is wrong.** "A relation of length >= 4 from i that does not return through i" is given with the example Z4/(x1x2x3x4) at i = 1. That relation does end at i. The criterion that the code and (i) actually use is that the remainder after stripping a·x' is not a single kM letter. State that. The type-II branch of (ii) rests on one instance.

**(H3)(ii)'s injectivity argument is sketchy.** "The prefix is read off t(q_t) = i" does not describe what happens: the pair r1/x1, x1.r1/x1 comes from the same p, with q equal to x3 and x1.x3 respectively. The a-prefix is carried inside q, not read off an endpoint. The conclusion is consistent with all 52 leads, but rewrite the sentence.

**The (C)/(LF) part (C) holds.** R3 = a1c + a1e forces w(c) = w(e), so w(R3/b) = w(R3) - w(b) = 0 for every homogeneous grading. R3/b is a degree-1 loop at 2, so (LF) fails by an argument and not only on the 3 vectors. The substitution is right: c = f - e and b = b' - a1 e give R1 -> b'x - a1 f x ≡ b'x mod (a1 f), R2 -> b'y - 2a1 e y, and R3 -> a1 f. The part that matters, that (LF) depends on the presentation, is shown.

## New?

- E-322 (Scope) leaves open "whether every parent can be reweighted to meet (C)". This submission settles the within-presentation half negatively and leaves the up-to-automorphism half open. That is new.
- E-279 says loops at i force cycles, and E-290 covers M1 only. No earlier record runs (LF)/(H1w)/(H3)/H(I) on finite-dimensional cyclic parents. I grepped `research/*.md` for "cyclic Nakayama", "leading word", "path length", "reachab", "presentation" and "Theta'", and found nothing that anticipates (A) or (B). The a*·p shape matches E-322's A6 2400 R6, as the author says.

## Evidenced?

Mostly. Scope, weights, WMAX, engine bound and controls are given precisely enough that the numbers can be believed without re-running, apart from the skip criterion above. The *Power* line holds for what it claims. The LF, MM, ENG and ISO controls fire, and the author says outright that the P-side MM check would miss a dropped term. Missing: no control shows that the sweep enumerator finds every finite-gl.dim class. My Kupisch count is that control, and it found the hole.

The Skeptic "Next" item (re-run on A^op) is probably mostly redundant. A^op of a monomial cyclic Nakayama algebra is again one, with the same relation lengths, so its class is already in the same sweep. Before dropping that item, the author should confirm that Def. 3.2 on A^op at i is the L-side step at the image vertex.

## Scope

The title and claim (A) should say 34 parents and 149 instances, or else say "33 of the 34 finite-gl.dim classes (Kupisch (2,3,3,3,3) not run)". Claim (B) should be narrowed to the affine, letter-additive statement. Claim (C) matches what was checked.

## Required for acceptance

1. Replace the skip criterion's parenthetical "(infinite gl.dim)" with the truth, which is that the cap is on size. Either add Z5/(x1x2x3, x2x3x4, x3x4x5, x5x1) to the sweep (raise `cap`; 39 s at WMAX 24) and update every count (34 / 149 / 20 / 33 / 11; 51 of 52 leads), or state it as excluded. Add a gl.dim-by-Kupisch control so that a skip is shown to be infinite gl.dim.
2. Narrow "no grading can do that" (B) to "no letter-additive g affine in (w, deg)", or prove the stronger statement.
3. Restate the type-II condition as "the remainder after stripping a·x' is not a single kM letter", and drop "does not return through i".
4. Rewrite the injectivity sentence in (H3)(ii) so that it matches the r1/x1 vs x1.r1/x1 case, where q carries the a-prefix.
