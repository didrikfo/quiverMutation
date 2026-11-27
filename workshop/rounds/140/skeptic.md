# Skeptic, round 140 (closing conference, cycle 2)

**Most promising question (next cycle).** Does a *proved* FCY control exist that is neither Dynkin nor hypersurface, so F-058's "two shifts => not FCY" is tested where the YES side is not just N_n(2)/N_n(3)? Candidate: N_15(7) (nu^7 ~ [8] on objects, r138 side find, still unproved). Reason: every F-058 control so far is a hypersurface or A_n type; the certificate's power on the genuinely new groups (N_14(4), N_14(6), N_15(3)) rests on the NO side being elementary, which is fine, but a false NO would only show on a non-trivial YES. Pre-register c = m b / a at m = 7k, both codes, then ask scholar for the source.

**Outside the main thread.** T10 in-class loss at depth >= 5 on E-316's n = 9 group (E-350 B). It is the only measured mode where the option costs nodes (93-180 per depth-4 walk), (A) at n = 11 cannot see it, and it is cheap (depth 5 sharded) compared with depth 8. Also dormant: H-029 widths 8-9 (retire unless real cases exist).

**Weakest claim.** "The J = 0 premise / child = End(T)" (T5, H-023) behind every "25 of 25 joined" and every Cartan-level NO chain. E-345 is n <= 6, acyclic, monomial; E-357 re-proves J = 0 steps only at n = 11 and by three computations of one quantity (so it tests code, not the premise). Second: E-269's float bound with q_S-only NO side (no q_S-independent NO control).

**Closest to a finding.** E-353/E-355 (F-058 size control + second prime): shifts predicted as numbers before the run, both codes agree on 4176/3888 pairs, non-integer predictions give 0 returns. One check still needed: a full return at m = 132 and a non-uniserial orbit check (nu2 compares dimension vectors in one degree; matters for N_14(4), N_14(6), whose certificate is uniserial-only), or the proved non-hypersurface control above.

**Needs.** Scholar: a source (or refutation) for N_15(7) being FCY, and a q_S-independent NO pair for F-055. Toolsmith: nu2 comparing isomorphism classes, depth-5 E-316 run. Theorist: Anick (H3)(ii) general case, or a statement that H-023 for the library does not need it.

## Closing ledger (my threads)
- T5 planted controls (E-348): SETTLED for monomial same-Cartan wrong children (341/341 NO, 978/978 ISO). STILL OPEN: modulus plant through the mutated vertex; non-monomial wrong ideals with an arrow parallel to a path; the 5 undecided steps; cyclic children, n >= 7.
- H-023 / H-031 premise: STILL OPEN (scoped at n <= 6); do not retire H-031 on the strength of E-279 alone.
- T10 option default: SETTLED as "off" for this cycle; depth 8 and E-316 depth 5 STILL OPEN.
- F-058 independent rebuild: second prime done, own-code controls done; SETTLED at its scope; m = 132 control and char 0 (argued, not run) still open.
- F-057 rebuild: SETTLED (E-354 sameIdeal + E-357 Hom, 0 disagreements); note E-357 is code-testing, not new evidence for J = 0 premise.
- F-055 rebuild (E-356): YES side SETTLED; NO side RETIRED as independent (rests on q_S); lam < 1.5 open.
- F-056 label: SETTLED to "Scherotzke 5.7 + Fosse 4.6, three primaries unread".
- H-029: RETIRE unless width-8 cases are exhibited. H-017b: min defect >= 1 census n = 10 still the only falsifier, OPEN, low priority.

## T10 recommendation: refuseNonTilting default
OFF. The option's only measured cost is in-class (93-180 nodes per depth-4 walk at n = 9, 19 of 19 steps lose their own child); its measured benefit (refusing E-318/E-323 crossings, seeds 514/3277) is already delivered by the key guard plus cross-class counts, and the one instance where it matters (J != 0 key-kept step) is a documented limitation, not a default. Threshold for flipping to ON, pre-registered now: lost in-class nodes <= 1% of the in-class ball AND >= 1 key-kept cross-class refusal at depth <= 8 on some n >= 11 root that the key guard otherwise admits. Neither holds (0 cross-class admits, ~10-20% loss). Keep as opt-in with regression = F-041 counts unchanged AND E-318/E-323 refused AND seeds 514/3277 refused.
