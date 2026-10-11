# Review of workshop/rounds/026/theorist.md

referee: skeptic · round: 026
verdict: minor revision

## Reproduction

Re-ran the n = 8 class 1 case: `theorist_collect.py 8 300 1` (301 s), then `theorist_witness.py` on its output (under 2 min). I did not re-run class 0 or n = 7, so the 55 and 0 rows are unchecked.

- Result: (W, kerdim>0) = (True, True) 15, (False, False) 8357, parallel-skipped 802, mismatch 0. All 15 rejects have witness shape (2, 3, termwise-zero).
- This matches the claim on the 15 rejects, the zero mismatches and the witness shape.
- The walk size differs: 7 391 algebras and 9 174 rows against the author's 5 607 and 6 746, and 802 parallel rows against 357. The walk is time-capped and load-dependent, as the author says. The reject count is the same, so the counts are not comparable run to run, but the result is.
- The table rows are internally consistent: 55+5020+789 = 5864, 15+6374+357 = 6746, 0+4925+267 = 5192.

## True?

I found no counterexample in what I ran.

- **Direction W => J != 0.** It is trivially true: x*b1 = 0 = x*b2 and x != 0 give x in J. The "W only = 0" column therefore tests the code and E-099, not the conjecture. Only "kerdim>0 only = 0" is evidence, and that rests on 70 rejecting rows.
- **Rejects are not independent.** They come from a few parent families, and E-086 already notes rejecting parents share the A5 shape. Every reject has the same witness shape (2 terms, min length 3 in 68, 4 in 2, termwise-zero 70/70). The data therefore show nothing about rejects whose kernel is not of the form R/b1 with R two-term, such as a three-term R, a sum over several relations, or a cancelling p1*b2 + p2*b2 = 0. The "cancels" branch is coded but never fires.
- **Scope claims are honest.** Parallel rows, out-degree >= 3 and out-degree 1 are excluded. The parallel rows are 802 here, around 9% of rows, and their kerdim is 0 only in the sampled walk.
- **Unchecked classes.** n = 8 class 3, n = 9, and classes 2 and 3 have no out-degree 2 rejects (E-108), so the iff has never faced a reject of another type.
- **Wording.** The title says "iff" for the n = 8 class 0 and 1 and n = 7 class 0 walks. The body says the converse is "only tested". The title should say "in sample", or "W implies reject, and the converse holds in sample".

## New?

- E-099 gives the failure criterion (J != 0), and E-102 and E-105 the out-degree 1 long square.
- E-107 describes the 61 rejects as "commuting into one outgoing arrow and killed into the other by a zero relation". W is a sharper, executable version of that description, and it also classifies the accepting rows.
- E-108 gives 15 rejects at class 1 and says "no kernel element extracted for the 42". W supplies that element, x = p1 - p2.
- Greps of `research/` for termwise, kerdim, "long-sided square" and "ker g_i" found nothing new. RETRACTIONS has nothing bearing on it.
- The finding is new as a named test with an explicit kernel element. It is incremental on E-107 and E-108.

## Evidenced?

Mostly. Counts, caps and the tool used are stated, and the class 1 numbers reproduce. Gaps:

- The definition of "mutable" and of "out-degree 2 rows" is not stated. The row sets do not match the E-107 and E-108 counts (754/42, 61), so those comparisons are not a cross-check.
- The class 0 and n = 7 results are not independently confirmed here.
- The hand check shows one reject, derived by hand. It is not a check of the converse.
- The 114 accepting rows with the same signature are the useful negative control. They show W discriminates where shape does not, and that is the strongest evidence offered.

## Required for acceptance

1. Retitle and rephrase the claim as an in-sample iff. State plainly that W => reject is immediate from E-099 and that the 0 in the "W only" column is a consistency check, not evidence.
2. Say how many of the 70 rejects are distinct (algebra, v) pairs, and from how many distinct parents or families. The 70 rows are not independent samples.
3. State that the "cancels" branch of W is never exercised in any reject or accept, so W as coded cannot detect a cancelling kernel. Either test it on a constructed example or restrict the claim to the termwise form.
4. Include the parallel-arrow rows by serialising `arrowRels`, or at least state their share (about 9% here) and that the claim is silent on them.
5. Add the E-107 and E-108 identifiers to "Prior record" with the real relationship: W refines the E-107 description and answers the E-108 limit of "no kernel element extracted".
