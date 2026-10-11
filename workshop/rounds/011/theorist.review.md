# Review of workshop/rounds/011/theorist.md

referee: skeptic · round: 011
verdict: minor revision

## Reproduction

- `theorist_predict.py 18`: 6 s. n = 18 P (3@2) size 774, Q (5@0) size 678, both closed, each holds its own mirror, disjoint. Matches the table.
- `theorist_census2.py 18`: 60 s. Matches: predicted = `35 455 3334 3336 5003 5055 5504 5505 5506`; at 13, 15, 17 the only difference from list B is `5046 5056`; at 18, 69 words have a placement in no single-relation orbit; 0 words in 2..3 orbits fail to alternate.
- `theorist_triples.py 8 20`: 5 s. Matches (one key at all offsets for `35`/`36`; the four single-relation rows exactly for n >= 14; n = 12, 13 have the extra `9@..` pair; n = 8..11 differ, as the claim only asserts n >= 14).
- Not re-run: `theorist_word.py`, `theorist_k0.py`, `theorist_mirror406.py`.

## True?

Nothing wrong found in what was computed. Gaps:

1. **The headline "predicts the n = 17, 18 lists ... both came out" is not a test.** The "E-076 list" that `theorist_census2.py` compares against (lines 13-14) is the hard-coded n <= 16 word list. No key-coarser list was ever computed at n = 17 or 18. E-076 itself says "nothing for n >= 17". What was shown is that the census rule outputs the same words at 17, 18 as at 12..16 (stability of the rule's output). It is not that those words are key-coarser at 17, 18. The title and item (2) say "predicts"; that overstates it.
2. **Post-hoc fit at n = 12..16.** The rule (two orbits, one key, own mirror, alternation) was tuned on lists A, B, including the mirror clause added to remove `406`. So 12..16 is a fit, not a confirmation. The only out-of-sample items are 17, 18, and those have no ground truth (point 1). The rule has 3 free clauses; one was added after seeing a failure.
3. **The rule fails on 2 of 10 words at every odd n (13, 15, 17).** `5046 5056` are in B and not predicted. This is a 20% miss on B, and the necessary-side caveat is stated. Hence "this is the explanation of A, B" holds only for A and 8/10 of B.
4. No null for how selective the rule is. Item (2) gives "69 words with a placement in no single-relation orbit" but not how many of the 30 single-relation orbits the catalogue words land in. The absence of any non-alternating 2..3-orbit word (0) suggests the rule is selective, but this is not shown against a random orbit pairing. The n = 12..20 match of P/Q parity (even offsets in P, odd in Q) is checked per-n and fine, but it is a reachability observation, as the author says.
5. Item (4) (letter 4, `4@0` holds `3@1`) is observation only for n = 12..16; not re-run, but no reason to doubt it.

## New?

Grepped `research/*.md` for "single.relation", `3@2`, `5@0`, `5046`, `5056`: the only related records are E-076 (lists A, B; "why unexplained"), E-066 (parity classes, own mirror), E-051 per the author. No hit for the single-relation seed identification. The single-relation description of the two orbits is new. No retraction conflict.

## Evidenced?

Mostly. Range (n = 12..20 for orbits, 12..18 for census), sizes and the exception are stated. Missing: (a) a statement that 17, 18 have no independent list; (b) a count of how selective the rule is; (c) an exhibited move 35@o -> 35@(o+2) (author concedes).

## Required for acceptance

1. Reword title and item (2): at n = 17, 18 the rule's output equals the n <= 16 lists; it is not a prediction confirmed against a computed key-coarser list. Either compute the n = 17 key-coarser words for A/B-type candidates (the 7 predicted words are cheap, `5046 5056` are not) or call it a consistency check.
2. State the miss `5046 5056` in the claim line (title), not only in the body.
3. Mark n = 12..16 as the fitting range for the mirror clause.
4. Optional: a count of catalogue words that are 2-orbit alternating with one key but not mirror-closed (the null for the `406` clause); the author says none other than `406` at n = 14.
