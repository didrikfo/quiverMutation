# Review of workshop/rounds/063/scholar.md

referee: skeptic · round: 063
verdict: minor revision

## Reproduction

Fetched arXiv:1504.02617 e-print into an empty scratchpad dir: gzip of one file, 1173 lines, 78 817 bytes, sha256 prefix 7cdf284fbc3dcd97 -- matches. Read Thm 1.1, Setup, Sect. 3 (Steps 1-4, Def. of d~), Sect. 7 (895-1043). No computation to re-run.

## True?

Substance holds. Hypothesis is exactly "no loops of degree 0 at i" (lines 264, 315, 321; setup adds only: non-negative, free on degree-wise finite generators over k^n, k a field). Non-negativity of M re-derived by me from the four arrow families: matches, so End(T) = H_0(kM) and J = H_1 are fair readings. Q-arrows of degree >= 3 give M-arrows of degree >= 2 in all four families: confirmed.

Errors in the write-up:
- Line numbers are off by one or more: theorem hypothesis is line 264 (not 263), "none of these arrows is a loop" is 321 (not 322), the Cone is 900-901 (not 905-908). Sect. 7 "Theorem" is 991, Lemma 963, Proposition 1027: these match.
- "cycles through i are handled by Step 4 (phi : i -> i)". Step 4 (line 354) handles only loops phi : i -> i of degree >= 1 (the degree-0 ones are excluded). Longer cycles through i are handled by Steps 2/3 (compositions into/out of i). Wording wrong; conclusion (no cycle exclusion) is right.
- "Proof contains no step that breaks for parallel arrows or cycles" is a statement about a sketched proof (the author concedes item 5 is one paragraph, "similar" for the second family). It is "I saw no hypothesis used", not "no step breaks". Narrow it.
- Item 1 of the proof reading (left approximation, loop used) is the author's own inference, flagged as such. Fine.

## New?

Mostly not new. `research/literature/1504.02617-quivers-for-silting-mutation.md` lines 36-40 already says "One hypothesis: no degree-0 loop at i. That is all." EXPERIMENTS.md E-181, E-183 and the paragraph near line 140 already have: conditional on Thm 1.1, J lives in H_1, G2 is the live item. Genuinely new: (i) the proof was read and non-negativity of M checked; (ii) Q^(<=2) suffices for H_0. Also E-168 holds the loop check; I did not recheck its 370/25 counts.

## Evidenced?

Adequate for a reading: line ranges named, derivations short. Missing: Sect. 8 (well-definedness, d~^2 = 0) unread, so "the theorem as proved" rests on the paper. The Q^(<=2) claim is derived by hand with no check against E-183's generator count.

## Scope

Title says "gives End(T) = H_0(kM) with no hypothesis beyond no degree-0 loop". Hypothesis accurate. "does not give J = 0" is correct. Narrow "proof contains no step that breaks" to "the hypotheses invoked in Sect. 7 are only the loop condition; step 5 sketched".

## Required for acceptance

1. Correct line numbers (264, 321, 900-901) and the Step 4 description (loops of degree >= 1, not cycles).
2. Say that the hypothesis list duplicates literature/1504.02617-quivers-for-silting-mutation.md:36-40 and E-181; state the new content as (i) and (ii) only.
3. Soften "no step breaks" per Scope.
