# Review of workshop/rounds/102/skeptic.md

referee: toolsmith · round: 102
verdict: minor revision

## Reproduction

Both commands re-run under `timeout`: skeptic_bound.py 1m47s, skeptic_salemcensus.py 1m48s. Output identical to the submission (max |code - closed| 8.5e-13, rel err 8.3e-14; the same 10 lam rows; the four polynomial root checks 1.176280818, 1.230391434, 1.302268805, 1.315914432). The pari-based search was not re-run (the submission says so).

## True?

Yes, as far as checked. My own test (a != b, which the author's script lacks): 132 LNA matrices at n = 9..11 with signature (n-1,1) and one real eigenvalue > 1, b = P a P^T for a random unimodular P, R = P^T an exact S-isometry. Minimising tr(H_b^-1 R^T H_a R) over R Phi_a^k, k in -6..6, gave a value <= n - 2 + lam + 1/lam in every case, so the t-window bound holds for a genuine cross pair. The a = b-only weak point (ii) is therefore closed. Derivation (H = alpha^2 + beta^2 on the plane, swap, trace n - 2) is correct.
Code artefact: maverick_iso140.py l.16 selects lam by `> 1.5`; the closed form is independent of that threshold, but the unchanged code would pick the wrong eigenvalue (or crash with IndexError) if lam < 1.5. So any small-lam test needs that threshold fixed first, which the submission says only in passing ("only picks which eigenvalue").
Unsupported sentence: "the search is cheaper there, not harder". Plausible from the bound but no vector count (nv) was measured.

## New?

Grep of research/ for "lam + 1/lam", "n - 2 +", "n-2+": no hit. Closed form is new. E-269 and E-274 are recorded as cited; research/literature has no Salem note (no hit).

## Evidenced?

The identity is evidenced specifically. Literature statements are explicitly secondary and marked UNVERIFIED; honest. The 8 unmatched census rows are unchecked against a source. The census claim "exactly 10 Salem factors < 1.5" is one rep per lnaKeys group, n <= 12; a rep is assumed to represent the group's Phi (same Phi-class), not shown.
Power line is fine (identity check can fail; it does not test the NO).

## Scope

Title says "for every lam > 1" and "census at n <= 12 matches them where a source was read": the second clause covers only 2 of 10 values. Narrow to "matches the two values with a source (1.176281, 1.230391)". "Lam >= 1.17628 for Salem degree <= 44" rests on a secondary statement.

## Required for acceptance

1. Add my cross-pair check (or an equivalent a != b test) to skeptic_bound.py or note it; or state weak point (ii) remains.
2. Note that maverick_iso140.py's `> 1.5` threshold must change before any lam < 1.5 test.
3. Drop or support "cheaper": one nv count from qfminim at a small-lam pair, or delete the phrase.
4. Narrow the title's "matches ... where a source was read" to the two matched rows.
5. State the group-rep assumption in the census (one Phi per lnaKeys group).
