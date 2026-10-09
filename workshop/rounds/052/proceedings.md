# Round 052 -- conference proceedings

All six personas wrote position statements (`rounds/052/<id>.md`); nothing was run. Five of six name the same weakest claim: that a J = 0 `tiltingPlus` step is a derived equivalence (the Hom test E-161/E-163 is Cartan-level, generation assumed; AI 2.31/2.32 and CHZ 3.6 cited from memory, UNVERIFIED), and that the key-guard law is class 0 only. The theorist's weakest claim is the same plus the caveat that the two should not be stated together.

## Ledger
E-entries since the last conference (E-150 .. E-165) touched these:
- **H-015** (guard sufficient): stays OPEN. E-150/E-151/E-157/E-160/E-163 say the 25 key-keeping `tiltingPlus` failures look like class members *under the J = 0 premise*; the premise is Cartan-level only. Open point named by the scholar: vertex-level versus whole-T tilting condition.
- **H-010** (overlap reducible only at an end): stays SUPPORTED. E-152/E-153 refine F-022 (run-of-three, L1 locality n <= 10); no proof; k-step window untested.
- **H-017** (quipu has more relations than cords): stays OPEN. E-162 shows the Euler signature is a Cartan function, so cannot test it; no relation-sensitive search exists.
- **H-020** (distance-to-ends): stays SUPPORTED. E-159, E-165 extend the floating-rule check (widths <= 8, one length each for 6..8). The theorist notes round 051's ablation shows the table changes no verdict in the reduced walk: the law may be about the move set, not the table; to be restated before any F-entry.
- **H-021**: stays OPEN, status line already reflects round 048 (only `s = n - k(c)` pairing survives); T1 dormant.
- No new finding (F-) and no retraction (R-): every result since 048 is conditional on the J = 0 premise or on a capped sample. No thread closed.

## Proposed agenda (ranked)
1. **T10 premise, non-Cartan test.** Skeptic: Hom replay of the 3 E-160 paths and quiver-level End(T) on the 13 E-163 edges (or why impossible). Theorist: which hypothesis of AI 2.31/2.32 is needed for a J = 0 step and whether generation is automatic for LNAs; scholar: record the statements as UNVERIFIED. (Skeptic, theorist, toolsmith.)
2. **Group-A witness path** 05040330 -> 33460000 by a tilting-only, relabelling-aware search, then `tiltingPlus` replay. (Experimentalist, toolsmith.) No overnight.
3. **Outside T10: T3/T8, P vs Q** (named by experimentalist, skeptic, scholar, maverick): is there a finer, non-Cartan derived invariant separating the two mirror-closed key-sharing orbits at n = 10; first a power control (equal-Cartan, different-class pairs at n <= 9, or a statement that none exist). (Maverick, theorist, toolsmith.)
4. **T4 / H-020**: cores with >= 3 relations at n = 13, 14 where rules = [] changes a verdict; else restate H-020 for the move set. Plus H6 in a rules-only walk. (Theorist.)
5. **Breadth: T7 / H-010**: run-of-three bystander at k = 4 two arrows away, one cell with a positive control; and a kernel/Cartan reason for exactly 2 LNA-to-LNA mutations (E-153). (Toolsmith/theorist.)
Parked: S-1 n = 15 K0 = 5 sizing (maverick; only after item 1 settles), literature (PDFs).

Disagreement: none on priorities; personas differ only on who runs item 1 (skeptic vs toolsmith vs maverick's invariant route). The scholar prefers a whole-T tilting test; the maverick a non-Cartan invariant with a power control.

## Questions for the steering committee
1. **Agenda**: approve the list above or change it in `STEERING.md`. Recommend: approve.
2. **PDFs of arXiv:1009.3370 and 2509.12983** (AI 2.31/2.32, CHZ 3.6): would settle the weakest claim directly. Recommend: yes, if you can supply them.
3. **`canonicalKey` DEFAULT_CAP 5040** (toolsmith): raise it together with the docstring rewords? Recommend: not yet; do it in a toolsmith round alongside the docstring fixes with the E-160 wrapper test as a regression.

## Decisions taken for the steering committee
- Round 051 q1 (overnight `merges.py 10 --depths 5 6 7 --witness`): no; group-A witness by a relabelling-aware search first.
- Round 051 q2 (PDFs): wanted if the human supplies them; literature stays parked (arxiv.org 403).
