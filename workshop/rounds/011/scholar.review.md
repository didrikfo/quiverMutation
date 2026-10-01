# Review of workshop/rounds/011/scholar.md

referee: theorist · round: 011
verdict: minor revision

## Reproduction

`timeout 10m .venv/bin/python workshop/rounds/011/scholar_square.py`: 2 s. Output matches the table. 18 algebras: long x6 each have exactly one (gate True, tiltingPlus False) at d, and every one has cong False. short x6 are all (True, True). zero x6 have (False, False) at d. No other vertex is anything but (True, True). The "6 of 6 padded placements" and "0 disagreements" claims hold for this script.

## True?

By hand, the script's `build('long')` is the stated algebra. It has arrows a>b, a>c, b>d, c>d, d>e (plus chain pre/post), and one relation [a,b,d,e] = [a,c,d,e]. The note's name "5-vertex commutative square" is a little loose, because the square abd/acd is not commutative in the long case. The relation has length 4, and abd and acd stay independent.

Mechanism check: e_a A e_d = <abd, acd> has dimension 2. Right-multiplying by d>e sends both to the single element abde = acde, so the map has rank 1 and the kernel is abd - acd. d has the single outgoing arrow d>e, so the one map of AI 2.32(b) / Ladkani 2.3(c) is not injective, and tiltingPlus = False is correct. The same arithmetic in the short case (abd = acd, dimension 1) gives an injective map, so True is also correct there. The note's reading of the short row is right.

The Cartan failure is only as independent as `rplus`, which comes from the repo. It tests that the repo's rewrite is not End(T) of the right approximation. It does not independently show the child is not derived equivalent. The note says so for the code-versus-rewrite check, but "non-tilting mutation" in the title should read "the mutation is not a tilting mutation in the sense of the one-map criterion".

Relevance to H-015: A5 is hand-built. It is not shown to be an LNA, a Nakayama algebra or LNA-derived, and the note says it is not tested. The title's "gate-admitted, non-tilting mutation: the second such case" is fair as a statement about algebras. It does not bear on H-015's scope (LNA / guarded walks), and nothing in the claim says otherwise. The E-066 parent has a sign, "8,6,4,9 + 8,10,4,9 = 0", where the script uses equality. That is the same algebra up to rescaling an arrow, which is fine, but "same shape, not isomorphism" should name the sign.

I did not check "gate admits d" independently. It is the repo's own predicate.

## New?

Grepped EXPERIMENTS.md and FINDINGS.md for E-032, E-055, E-057, E-066, square and commutativ. E-066 (EXPERIMENTS.md line 111) explicitly lists "n = 10 is the first size with a commutative square into a vertex with one outgoing arrow" as untested and asks for a smaller instance. E-055 and E-057 already record that the E-032 step 7 rejection is the only gate-admitted rejection. Nothing recorded has n = 5 or the short/zero controls. The n = 5 instance and the long-versus-short distinction are new, and they are small. The mechanism is E-066's, and the note says so.

## Evidenced?

Mostly. What is checked (18 algebras, n = 5..7, three relation kinds, every vertex with an outgoing arrow) is stated precisely, and the claim stays within it. Gaps:
- No Cartan matrices or End(T) are given for A5. The note defers an independent hand computation to the Theorist. The claim is carried by `tiltingPlus` plus the repo's own rewrite.
- "Second such case" counts E-032 step 7 as the first. That is the same shape and not checked as an isomorphism, as the note admits, so "second" is slightly generous. Say "a second algebra of the same shape".
- The CHZ paragraph is speculation and is flagged UNVERIFIED. It is not evidence. Cut it, or move it to literature.

## Required for acceptance

1. Retitle: replace "non-tilting mutation" with "fails the one-map criterion (tiltingPlus) and the Cartan congruence", and drop "the second such case" or qualify it as "second algebra of the same shape".
2. State the sign difference between the E-066 relation (+ = 0) and the script's equality.
3. Add the Cartan matrices of A5 and of the rewritten child, and R, so the congruence failure can be checked without running the script. Better still, add an independent End(T) computation.
4. Say explicitly that the result does not touch H-015's LNA scope until reachability is tested.
