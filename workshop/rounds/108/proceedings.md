# Round 108 -- proceedings (conference, cycle 1, phase 2: research)

Step 0.5: no open questions from round 107 (it asked none); nothing to decide. Merge of origin/main: already up to date. All six personas wrote positions (nothing run).

## Positions (one line each: most promising question)
- experimentalist: does A_n stay reachable from B_n by tilting exchanges at n >= 7, first hit at depth n - 2?
- theorist: can Prop. 7.6's loop case be proved by a two-stage elimination (filter by Theta letters, then unitriangular)?
- skeptic: along the E-287 paths, is every intermediate End(T) an LNA, and where does the gate or the key refuse or change?
- scholar: is the E-287 path a genuine tilting/silting mutation sequence, with LNA intermediates (a move-set-independent certificate for T10/H-027)?
- toolsmith: same as the skeptic, plus the depth law at n = 7..9.
- maverick: can a planted-error control make E-269's NO side work at lam < 1.5?

Agreement: four of six (experimentalist, skeptic, scholar, toolsmith; theorist as its outside question) put the E-287 path at the centre. Disagreement: the skeptic calls "H-027 SUPPORTED" the weakest claim (A_n ~ B_n is Fosse Cor. 4.6, so hits were guaranteed; the hits say nothing about gate or guard) while the experimentalist and toolsmith call E-287 closest to a finding. Weakest claim named: J = 0 premise (experimentalist, toolsmith, theorist); Cartan-level necessity step / E-269 NO power (scholar, maverick); H-027 (skeptic). Chair's reading: both are right -- E-287 is a sound existence result for the path, and its worth is only in the intermediate-End question.

## Ledger
- **H-024 (SUPPORTED):** no status change. **Strongest F-entry candidate.** Evidence: necessity published (E-272), E-281 YES on 6 of 6 mixed pairs and NO on 4 of 4 cross pairs. **Blocking check:** a planted-error NO on one certified-equivalent mixed pair (perturb one Cartan entry, expect NO), pre-registered by the skeptic, together with the maverick's revision of round 107 (E-276 restated, group count 15 vs 12 reconciled). Assigned to round 109.
- **Strongest R-entry candidate:** the claim in the `mutationSearchDepthFirst` docstring (l. 322) and the guard's reputation that a key-guarded walk stays in one derived class. Unsupported since round 045; E-282 shows 7 of 8 key-kept J != 0 edges are silting, not tilting, but none is shown to leave the class. **Blocking check:** for the E-287 paths, End(T) of every intermediate step with its Gabriel quiver, and the gate/key verdict per step (theorist + toolsmith). No R-entry yet; the docstring fix is a backlog item.
- **H-027 stays SUPPORTED** with the skeptic's caveat recorded here: existence of the path was guaranteed by Fosse Cor. 4.6; what is new is the explicit short tilting path. Generation (Bruning) and silting-order check still open.
- **H-023 (T5):** stays OPEN. E-284 is a computation in ten models, not a proof. H-010, H-017, H-020, H-025, H-026: no change.
- Threads: none closed. T1, T2 already closed; T9 parked.
- No F- or R-entry written this conference; both candidates have a named blocking check.

## Proposed agenda for phase 3 (research, rounds 109-112; conference 112)
1. **H-024 power (maverick, skeptic).** Round 109: the maverick's revision (due) with a planted-error control; skeptic pre-registers the perturbations. Then a certified-inequivalent pair at lam < 1.5 (n = 13..15, toolsmith enumerates).
2. **T10 / H-027 path audit (theorist, toolsmith; skeptic referee).** Intermediate End(T) of the E-287 paths, LNA or not, gate verdict and key per step; then n = 7 depth to first hit with a per-depth log; generation chain for the n = 5 hits. Scholar: Bruning and Aihara-Iyama theorem numbers from the primary.
3. **T5 / H-023 loop case (theorist, toolsmith, skeptic).** Two-stage elimination proof attempt; Hom-computed End(T) for M1 against dim H(kM); a degree-0 loop control that can say NO (planted failure, pre-registered by the skeptic).
4. **Outside the main thread: S-1 item 3 / H-026 (experimentalist, skeptic).** Pre-registered prediction K0 = 6 at n = 17; does the law survive with "Cartan Z-congruent" replacing "same key". Alternate slot: H-025 converse, literature only (scholar), or the singular-S groups of E-269 (maverick).
5. Backlog seeds stay for phase 5. T7 gets no slot: E-288 has no control and the overnight run is not approved. T6 waits for a statistic with power.
Disagreement to note: the skeptic would not spend more rounds on H-027 beyond item 2; experimentalist wants n >= 7 first. Item 2 sequences them: audit the n = 5, 6 intermediates before n = 7.

## Questions for the steering committee
None.

## Decisions taken for the steering committee
Step 0.5: nothing to decide. No overnight run approved.
