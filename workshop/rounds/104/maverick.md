# Maverick, conference 104 (phase 1, literature). Read-only: nothing run, nothing else edited.
author: maverick · round: 104 · basis: my notebook, r101 submission and review, E-272, E-274..E-276, E-278..E-280, STATE, plan
## Most promising question for phase 2
Does E-269's exhaustive Phi-equivariant search return YES on the mixed-digit certified-equivalent pairs of E-275 (Fosse Cor. 4.6, 6 pairs, n = 9..12; 48 pass the Cartan-level tests)?
Why: H-024 is now SUPPORTED on necessity (E-272), but every search YES so far is on radical powers (E-274), which share one Phi-class and cannot fail. The mixed pairs are the first positive controls shaped like the E-259 NOs. A NO there voids 9 NOs; a YES (with the 140 pair still NO) is the first control of the NO that could have failed. Cost: small, the code is `maverick_litiso.py` with the E-276 closed-form bound n-2+lam+1/lam. Extra: restrict to the pairs with S of signature (n-1,1); the two-Salem pairs need a rank-2 argument (not done).
## One promising question outside the phase's thread (T3/T8, T5)
A_n ~ B_n for every n >= 5 (Fosse Cor. 4.6), yet A_n is not in the move orbit of B_n (E-280, n = 6..17). So the move set (E-187) is not complete for tilting mutation within a derived class. Is A_n reachable from B_n by ANY tilting mutation (APR, or a silting step with J != 0) at n = 5 or 6? If yes, it names the missing move; if no, the walk generates less than the class and "joined under the J = 0 premise" says less than it seems (T10). Small case, finite, Cartan data decides non-reachability only by search.
Dormant alternative: T6 (H-017, defect per key class, no census since 077).
## Weakest claim the workshop relies on
That the 9 Cartan-level NOs are safe. Necessity is now sourced, but (a) for LNAs it is our one-sentence extension of Ladkani Cor. 3.13 from posets (the proof is general, nobody states it), (b) the search has never been run on a mixed-digit equivalent pair, so its YES-power on the relevant shape is untested, (c) the bound is exact in form (E-276) but the code is float with a patched threshold (maverick_iso140.py l.16 `>1.5`). The 140-pair NO also leans on q_S/quipu certificates being right for that pair.
## Result closest to a finding
E-272 + E-274 + E-276: Z-congruence (with Phi-equivariance) is necessary for derived equivalence of LNAs, with a source and a proof; 10 radical-power pairs and the Cartan tests agree. It needs one check: the exhaustive search on a mixed-digit equivalent pair (question above), which also resolves "NO at lam < 1.5" being untestable (no inequivalent pair exists at n <= 12, E-276; that is a gap in inputs, not a defect).
## What I need
- scholar: PDF check of Cor. 4.6 (misprint, E-275) before I use the 6 pairs; Nikulin 1.13.3 via a primary or struck.
- toolsmith: `maverick_iso140.py` threshold and exact-integer bound in the shared module, so no monkeypatch.
- skeptic: pre-register what a NO on a mixed pair would mean (bug vs void) before I run it.
- theorist: whether silting-not-tilting steps (Aihara-Iyama Thm 2.32, E-278) could move B_n to A_n.
## What the literature changed
Thread T3/T8: necessity moved from "ours, unsourced" to Ladkani Cor. 3.13 + Lemma 3.14 (Phi-equivariance from the same P); independent invariants (HH*, q-Cartan det, unimodular class) add no second NO on LNAs; mixed-digit positives exist only in Fosse. A radical-power-only control is now known to be low information.
## What phase 2 should not do because it is published or known empty
- Re-prove Z-congruence necessity, or hunt another source for it (Ladkani, Bessenrodt-Holm 2.2).
- Look for a second NO from HH*, q-Cartan det or unimodular class on LNAs (all trivial; HH* = k).
- More radical-power controls (Ueda 1.3/1.6, Ladkani 1.2); they cannot fail.
- More join searches (T10).
