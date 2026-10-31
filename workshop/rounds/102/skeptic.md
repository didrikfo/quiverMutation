# E-269's search bound is the closed form n - 2 + lam + 1/lam (+ slack) for every lam > 1, so no published Salem bound is needed for it; the published lists (secondary, UNVERIFIED) only say which lam can occur (lam >= 1.17628 for Salem degree <= 44), and the repo's census at n <= 12 matches the two values with a source read (1.176281, 1.230391)

author: skeptic · round: 102 · kind: literature
thread: T3/T8 · bears on: H-024, E-265, E-269, E-274
scope: (a) Re-derivation by hand + numerical check of the bound on 450 LNA (a = b) and 160 cross pairs a != b (n = 9..12, b = P a P^T) Cartan matrices (n = 10, 11, 12, S of signature (n-1,1), det S != 0, exactly one real eigenvalue > 1 of Phi; <= 150 per n, random sample of lnaKeys strings). (b) Literature: two arXiv e-prints read in LaTeX (2409.11159 v3, 2008.10965 v2). Mossinghoff's web lists and the primary papers (Lehmer 1933, Boyd 1977/78, Flammang-Grandcolas-Rhin 1999, Mossinghoff-Rhin-Wu 2008) were NOT read: the wayback/cecm hosts do not resolve here. The bound's search itself (qfminim on the commutant lattice) was not rerun (cypari absent). Not a proof of E-269's NO.
power: n/a (literature and an algebraic identity). The identity check can fail (relative error 8e-14 observed, i.e. it held); it does not test E-269's NO.

## Response to referee

Referee (toolsmith, round 102, minor revision). Point by point:

1. **Cross pair (a != b): done.** New `skeptic_cross.py` (160 pairs, n = 9..12, b = P a P^T with random unimodular P, R = P^T exact S-isometry with R Phi_b = Phi_a R, minimise over R Phi_b^k). Result: 0 violations of n - 2 + lam + 1/lam, smallest lam 1.1763. Note: with k limited to -6..6 as in the referee's range, 9 of 160 failed, because t then lies outside the lam^(+-k) reach for lam near 1.2; k = -40..40 gives 0 violations. So the window argument needs k large enough to bring t into [lam^-1/2, lam^1/2] (always possible; k is unbounded in the theory). Weak point (ii) is closed. Reproduce: `timeout 10m .venv/bin/python -W ignore workshop/rounds/102/skeptic_cross.py`.
2. **1.5 threshold: noted.** `maverick_iso140.py` l.16 (`w[k].real > 1.5`) must be lowered (e.g. > 1.0001 and exactly one such eigenvalue) before any lam < 1.5 test, else IndexError or the wrong eigenvalue. My scripts use > 1.0001 and require uniqueness. The closed form itself does not depend on it.
3. **"Cheaper": dropped.** No nv count was made (cypari absent, cannot run qfminim here). The bound being smaller is shown; whether the search is cheaper is not. The sentence is removed from Claim 2.
4. **Title narrowed** to the two matched values; "lam >= 1.17628 for Salem degree <= 44" is now marked secondary/UNVERIFIED in the title.
5. **Group-rep assumption stated** in Claim 4 below: the census takes one representative per `lnaKeys` group and assumes it represents the group's Phi (same Phi-class); this was not verified.

## Claim

1. **E-269's bound is exactly n - 2 + lam + 1/lam, plus slack.** In E-269's notation the code's BOUND (`maverick_iso140.py` l. 28-30) equals tr(H_b^-1 S_b) + max_t [t^2 + t^-2] + slack; the first term is n - 2 and the max is lam + 1/lam, attained at the ends of the t-window. No lam-specific table enters. Derivation (checked by me, by hand and numerically): with p, q the Phi-eigenvectors for lam, 1/lam, S(p,q) = 1, S(p,p) = S(q,q) = 0, d = q - p, H = S + (Sd)(Sd)^T: on the plane span(p,q) in coordinates (alpha, beta), S(x,x) = 2 alpha beta and S(d,x) = alpha - beta, so H = alpha^2 + beta^2 (positive definite) and H = S on the S-orthogonal complement. Hence H^-1 S is the swap on the plane (trace 0) and the identity on the complement (trace n - 2). For an isometry R, R^T H_a R = S_b + (S_b z)(S_b z)^T with z = R^-1 d_a = t q_b - p_b/t, and z^T S_b H_b^-1 S_b z = t^2 + t^-2 (S swaps coordinates, H^-1 = I on the plane). Multiplying R by Phi_a^k rescales t by lam^k, so t can be taken in [lam^-1/2, lam^1/2]; t^2 + t^-2 is convex in log t, maximal at the ends: lam + 1/lam. Nothing used lam > 1.5: the code's 1.5 threshold only picks which real eigenvalue is lam (agrees with the 101 referee's check). Sign of t is immaterial.
2. **Consequence for lam in (1.1, 1.5):** lam + 1/lam is in (2.0, 2.17), the bound without slack is in (n + 0.009, n + 0.17) for n = 12, just above the norm n of the identity (which sits at t = 1: bound minus norm = lam + 1/lam - 2 >= 0). The bound is smaller there (search cost not measured: no qfminim nv count, cypari absent). The "NO for lam < 1.5 untestable at n <= 12" of STATE/E-274 is an input gap (no certified-inequivalent pair with a small-lam Salem factor exists in the census), not a bound gap. The bound is not an integer: n - 2 + lam + 1/lam is irrational; lam + 1/lam is an algebraic integer of degree <= n/2, so an "exact integer bound" does not exist in this form.
3. **Literature, which lam can occur (all secondary for the primary lists):**
   - Lehmer's number 1.176280818 (poly x^10+x^9-x^7-x^6-x^5-x^4-x^3+x+1) is the smallest KNOWN Salem number: supported by 2409.11159, Example in Sec. 1 ("Salem number tau_0 ... smallest known Salem number to date") and 2008.10965, Sec. 5 text before Prop. 1 (alpha_L = 1.176281). Whether anything smaller exists is Lehmer's problem, open (not addressed by any source I accepted).
   - Complete lists: 2409.11159 states Mossinghoff's list of 47 Salem numbers < 1.3 was "certified complete up to degree 40 in [Flammang-Grandcolas-Rhin 1999] and later up to degree 44 in [Mossinghoff-Rhin-Wu 2008]" (subsection "Small Salem numbers"; no theorem number, text only; the primary papers UNVERIFIED). 2008.10965 says (L) is "the complete list of all Salem numbers of degree <= 44, due to Boyd, Mossinghoff and others". So for Salem factors of degree <= 12 <= 44: lam >= 1.17628, and the interval (1.1, 1.17628) is empty at n <= 12. This is a table lookup of a published (secondary) statement; it is not a proof I checked.
   - 2008.10965 Prop. 1 (first numbered proposition by count, label `238`): tau_2 = 1.230391, min poly t^10 - t^7 - t^5 - t^3 + 1, "seventh smallest Salem number in the list (L)". Supported; matches the repo value for n = 11.
   - 2409.11159, longtable (unnumbered), rows 2d = 12: 1.302268805 (coeffs 1 -1 0 0 0 -1 1) and 1.315914432 (1 0 0 0 -1 -1 -1); I reconstructed both polynomials and recomputed the roots: they match to 9 digits. Its table starts at 1.3, so it says nothing about degrees <= 12 below 1.3.
   - A published statement turning the gap into a lookup "which lam occur at n <= 12 below 1.5": NOT FOUND in anything read; Mossinghoff's fixed-degree lists would be the source (cited as [MossinghoffList2] in 2409.11159, not retrievable).
4. **Census vs publications (computation, one rep per `lnaKeys` group, n = 4..12, all signatures; ASSUMPTION, not shown: a rep's Phi represents its whole group):** Salem factors with lam < 1.5 are exactly 10 (table below). 1.176281 (deg 10, n = 10) and 1.230391 (deg 10, n = 11) match the polynomials in item 3. The other 8 are not checked against a source (Mossinghoff's list unreachable). The two degree-12 values in (1.3, 1.3243) of 2409.11159 (1.3023, 1.3159) do not occur as Phi-factors of LNAs up to n = 12: not a contradiction (not claimed to be realisable).
5. **Not addressed / not found:** a published exact bound for Salem numbers in (1.1, 1.5) as a function of degree beyond the lists above; Lehmer-type lower bounds such as Dobrowolski's bound (from memory, UNVERIFIED, not searched in a source) are far too weak at degree <= 12 to matter (they would give lam > 1 + c (log log d / log d)^3), and I do not use them. Two arXiv papers claim a proof of Lehmer's conjecture (1709.03771, 1911.10590): found by search, NOT read, not relied on.

## Evidence

Queries (export.arxiv.org API): `abs:Salem AND abs:Lehmer AND abs:smallest`; `ti:Salem AND ti:numbers AND ti:small`; `abs:Salem AND abs:degree AND abs:list AND abs:Mossinghoff`; `ti:Lehmer AND ti:conjecture`; `ti:Salem AND ti:Coxeter`; `ti:"Salem numbers" AND abs:smallest AND abs:degree`; two others returned nothing. Identifiers seen: 0905.2396, 2404.03288, 2409.11159, 1709.03771, 1911.10590, 1401.3128, 2008.10965, 1306.3443, 1905.04023. E-prints fetched (LaTeX, scratchpad `lit/`): 2409.11159 (v3 per API: 2025-01-31 source), 2008.10965 (v2), 0905.2396, 1401.3128, 1306.3443, 1905.04023, 2404.03288 (the last five grepped only for the values 1.2612/1.2303: no hits except 2008.10965). Wayback/cecm hosts: DNS failure.

Identity check (skeptic_bound.py): 450 matrices, max |code bound - (n - 2 + lam + 1/lam)| = 8.5e-13; identity tr(H^-1 R^T H R) = n - 2 + lam^{2k} + lam^{-2k} for R = Phi^k, k = -1, 0, 1, 2 (exact isometries, asserted), max relative error 8.3e-14. Smallest lam in the sample: 1.1763 (n=10), 1.2304 (n=11), 1.2612 (n=12), 1.3600, 1.3836, 1.4013.

| lam | deg | n (group reps) | source match |
|---|---|---|---|
| 1.176281 | 10 | 10 | Lehmer poly, 2409.11159 / 2008.10965 |
| 1.230391 | 10 | 11 | 2008.10965 Prop. 1 |
| 1.261231 | 10 | 12 | none read |
| 1.280638 | 8 | 9 | none read |
| 1.360000 | 8 | 10 | none read |
| 1.383637 | 10 | 12 | none read |
| 1.401268 | 6 | 8, 10, 11, 12 | none read |
| 1.425005 | 8 | 10, 12 | none read |
| 1.457987 | 8 | 9, 11 | none read |
| 1.488565 | 12 | 12 | none read |

Weak points: (i) my derivation shares the E-269 set-up (R Phi_b = Phi_a R, no swap of lam and 1/lam: forced, since R maps the lam-eigenline to the lam-eigenline); (ii) CLOSED by skeptic_cross.py (160 cross pairs a != b, 0 violations); (iii) float, not interval arithmetic: unchanged from E-269, but with the closed form the bound can be evaluated to any precision, so E-269's float majorant could be replaced by exact lam + 1/lam from the Salem polynomial.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/102/skeptic_bound.py        # ~1 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/102/skeptic_salemcensus.py  # ~3 min
```

## Prior record

Grep of `research/` : E-269 (the search; "float bound" in the Scope), E-265, E-274 (Scope, Result), STATE T3/T8 and the skeptic request. E-269's Reduction/Evidence derives the same inequality but gives the number 13.11 = 10 + 2.611 + 0.5 without the closed form; the closed form n - 2 + lam + 1/lam is new here (the 101 maverick text also says only "depends on lam through the t-window"). Small-lam NO control is in E-274 (census, no pair below 1.661). Grep of research/literature/ for Lehmer/Salem: no hits, so no note on small Salem numbers exists; the two e-prints above are new to the record.

## Code changed

None in the repo. New: `workshop/rounds/102/skeptic_bound.py`, `skeptic_salemcensus.py`, `skeptic_cross.py` (numpy/sympy only; no cypari).

## Next

- Maverick/chair: E-269 can quote "bound = n - 2 + lam + 1/lam + slack" (exact, from the Salem polynomial) and drop the float majorant; the remaining float is the qfminim Gram rounding (1e-6).
- Scholar: fetch Mossinghoff-Rhin-Wu 2008 / Flammang-Grandcolas-Rhin 1999 (not on arXiv, unreachable here) to turn the secondary "complete to degree 44" into a primary quotation; and Mossinghoff's fixed-degree lists to check the 8 unmatched rows above.
- Skeptic/maverick: a small-lam NO control needs pairs at n >= 13 or non-quipu certificates; the bound is no obstacle (it shrinks to n + 0.1 as lam -> 1.2).
