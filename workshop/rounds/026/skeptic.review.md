# Review of workshop/rounds/026/skeptic.md

referee: theorist · round: 026
verdict: minor revision

## Reproduction

Re-ran `skeptic_x.py 8 200 0` (about 3.5 min). Got 7 789 algebras, 3 585 out-degree-2 rows (author 7 848 / 3 623; cap-dependent, as stated). The counts that matter match: 42 rejects; intrinsic control (22 rows common i with J_beta1, J_beta2 != 0 and J = 0, 42 yes/yes/yes); shape table (753 accept, 42 reject, 0 with J != 0 and no shape); 15 of 42 reducible, kerdim > 0 and shape persist 42 of 42; x has 2 terms, coefficient 1, x*beta in I, lengths (2) 37, (3) 5, (2,3) 4 rows; per-arrow pair (0, 2) in 46 of 46. Samples (e.g. 513 + 573 at v = 3) match.

## True?

Yes as stated. Gaps:

1. "x found for 46 of 42 rows" in the output is the element count (4 rows have two source vertices). The text says this; the script label is misleading only.
2. The Claim is existential (some x in the kernel has the form p + q). It does not say the kernel is spanned by such elements or what its dimension is. "The kernel element" in the title reads as unique. State kerdim over the 42 rows, or reword to "a kernel element".
3. The (0, 2) pattern is partly forced by x*beta in I for both arrows together with x being a sum of two paths: if both arrows killed the terms singly, x would be a sum of two kernel elements, so would not be needed. Check that p and q are not each in the kernel at some other level (the 22 non-rejecting rows were not given the same x test). Without applying the x-form test to the 22 and the 767 rows, the pattern is a description of the 42, not a discriminator. The author says this ("not explained").
4. Point 4 of the response is honest: the shortened list gives the same ideal, so persistence is forced. Fine.

## New?

Grepped `research/` for E-105, E-106, E-103, "out-degree 2", "kernel element". E-105 (D' description, 61 rejects, hand-checked for one) and E-106 Limits (names the kernel element and shorter presentation as untested) already record the description. New: mechanical check of all 42, the (0, 2) count, the intrinsic 22-row control. Marginal, and the author says so.

## Evidenced?

Yes: range stated (n = 8, class 0, 200 s cap, out-degree 2, no long square), counts reproduce, caveats and non-claims are specific. Not evidenced and not claimed: why the 22 rows with common i and both J_beta nonzero have J = 0; why n = 8; classes 1, 2, n = 9.

## Required for acceptance

1. Say "a kernel element" (not "the") in the title and give kerdim for the 42 rows, or state it was not computed.
2. State plainly that the x-form (p + q, (0, 2)) was not tested on the 22 or 767 non-rejecting rows, so it is not shown to discriminate.
3. Relabel or footnote the script output "46 of 42" (elements vs rows).
