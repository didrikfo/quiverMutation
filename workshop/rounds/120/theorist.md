# Theorist, round 120 (cycle 1 closing conference)

Basis: my notebook, rounds/118 (M4/M5 write-up, E-309), STATE, plan, backlog. I wrote nothing in 117 or 119 (no theorist files there).

**Most promising question for the next phase.** Is H-031 (general pivot-reducible loop cancellation) true, with Hom/End(T) closed? Every proved model (M1, M2, M4, M5) has the same shape: kM = B * T(E), E contractible, B of Lemma-A form, and a g-filtration reduces B to Lemma A. That makes a proof of the general statement look achievable, and it is the only route to deciding H-023 (is every gate-admitted J = 0 step a derived equivalence?), which underlies the "25 of 25 joined" claim. Falsifier: a loop model with H(I) != 0 or no reachable pivot shape. Next concrete step: M3/M3b (two loops, junction) by hand, then a finite-dimensional variant with an added relation so Hom/End(T) is computable.

**Promising question outside that thread.** T7/H-010: does the end-only argument already exclude gap-2 lowering? E-313 is a null with no power, so a short proof (or a counterexample shape) is the only way to retire it honestly. Second choice: the self-mirror exception in H-026 (E-307), where a derivation would turn a fit into a law.

**Weakest claim the workshop relies on.** The J = 0 premise (H-023), together with the signs and printed differentials my proofs use: `last = -1` is taken from E-202, and the proofs rest on the engine's printed d, not on Oppermann's journal version. Hom/End(T) for M4/M5 (infinite-dimensional) is not done, so "cancellation" is shown at the dg level only. Runner-up: Nikulin Cor. 1.13.3 still UNVERIFIED, behind every Cartan-level NO (though E-308 now gives Ladkani + Rickard + Aihara-Iyama as an independent source for necessity).

**Result closest to a finding.** E-309: M4 and M5 cancel for all weights (pivot change of generators + g-filtration to Lemma A), with the structure checks in theorist_m45_struct.py (generator bijection, word counts per cell to length 5, B and E separately d-closed, explicit inverse substitution). One check it still needs: an independent second implementation of the d identities, plus Lemma A with base k^3 (M5) and finiteness of the cells for the g-filtration (skeptic). Until then it is a checked proof, not an F-entry.

**What I need.** Skeptic: the independent d-identity check and the Lemma A / k^3 / cell-finiteness check above. Toolsmith: nothing new; the walk checkpoint (about 11 min rebuild) is lost with the session. Scholar: the primary for Anick Thm 3.1/3.2 and Oppermann's journal version, to confirm the sign conventions.

**Backlog items.**
- Settled (my part): 9 (M4/M5 proved, converted to H-031 for cycle 2); 13 (E-305 points: the 2 NO children are not trees, HH^2(ab = 0 square) = 1 by Bardzell, corruption counts, the 78 vs 84 bookkeeping), as far as I own them.
- Converted to cycle 2: 9 remainder (M3/M3b, Hom/End, finite-dimensional step), 19 (H-017 to H-017b), 14 if still "narrowed only" (convert, do not carry).
- Retired: 18 stays UNVERIFIED by declaration (human PDFs unavailable); 17 not approved (already decided); 15 non-minimal part (already R-016).
- Still open: 2 (author-only), 8 (H-025 search), 10 residual (my retirement argument for gap-2), 1 caveat (no key-kept step crossing the F-010 halves found, H-015 untested), H-029 widths 8-9, E-080 shape of the 8 anomalies.
