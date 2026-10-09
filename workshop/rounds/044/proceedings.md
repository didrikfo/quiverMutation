# Round 044 -- proceedings (conference)

All six personas wrote position statements (`rounds/044/<id>.md`); nothing was run, nothing promoted.

## Positions
| persona | most promising question | weakest claim |
|---|---|---|
| experimentalist | do key-preserving J != 0 steps exist at n = 8, 9 c1, c2 and does the x^2 law hold on them | E-145 class-1 count rests on one script |
| theorist | why e_i = F^s e_w with c_2 = 0 on every gate-admitted J != 0 LNA-walk step (module derivation) | "J != 0 leaves the key" (class 0 only); Q_2 = 1 proved for s = 1 only |
| skeptic | is a key-preserving D = 0 child at n = 7 c1/c2 reached by tiltingPlus + Cartan-congruent steps, and derived-equivalent by an independent test | H-015 as the derived-equivalence test; key guard is also the BFS filter |
| scholar | why such steps exist; does it falsify H-015 on walks or show the key is no class invariant | H-015 cited from memory (AI 2.31/2.32, CHZ 3.6 unread) |
| toolsmith | does a key-preserving walk step ever fail tiltingPlus (is the BFS filter a larger set than tilting) | "key preserved implies tilting"; reverse control shallow, 10.3% edges lost |
| maverick | does the S-1 K-threshold law hold at n = 15 (K0 = 5) | key-guard law; S-1 law rests on two lengths |

Agreement: five of six name the E-145 key-preserving J != 0 steps (or their consequence for H-015) as the live question, and four ask the skeptic for an independent hand rebuild of the class-1 steps. Disagreement: maverick alone prefers S-1 n = 15; theorist wants a module-level derivation where toolsmith/skeptic want empirical soundness checks first. Scholar asks the human for PDFs.

## Proposed agenda (proposed, round 044), replacing the round-040 agenda
1. **T5 key guard, independent check.** Skeptic rebuilds the 13 class-1 E-145 steps at n = 7 by hand (or a parallel-arrow-safe loader); experimentalist then runs a guard-off n = 7 census (classes 0-2, `--plan` first) with per-step `tiltingPlus` and Cartan congruence, saving the n = 8 J != 0 steps. First question: is any key-preserving step tilting-compatible and Cartan-congruent?
2. **T5 theory.** Theorist: orbit data giving D = 0 (B = 0 = 1 - c_s - c_{-s}, s = 10 in c2); module-level reason for c_2 = 0; the s = -2 family. Says which property the meet needs (tilting, silting, or key).
3. **T5 soundness of the tools.** Toolsmith: does a key-preserving walk step ever fail `tiltingPlus` (per-step tally); loss-by-depth tally for the reverse search and a deep control. No overnight.
4. **Literature (scholar).** Retry the arXiv fetch of 1009.3370 / 2509.12983 under the standing permission; if still blocked, state what H-015 needs from them and park.
5. **S-1 n = 15 K0 = 5** (maverick, toolsmith): size by orbit enumeration with `--plan`; which of min(h, K) or the core key governs first failure.

## Decisions taken for the steering committee
- round 043, question 1 (agenda): keep the round-040 agenda; item 1 = why key-preserving J != 0 steps exist at n = 7 c1, c2 -- decided by the chair of round 044; no answer from the human. Superseded by the proposed agenda above.
- round 043, question 2 (overnight): none; the 3 h reverse job waits for a deep control -- decided by the chair of round 044; no answer from the human.
- round 043, question 3 (literature): no PDFs supplied; item stays parked, scholar retries the arXiv fetch under the standing permission -- decided by the chair of round 044; no answer from the human.

## Questions for the steering committee
1. **Agenda:** approve the proposed agenda above. Recommend yes.
2. **Overnight:** none proposed. Recommend none.
3. **Literature:** PDFs of arXiv:1009.3370 and 2509.12983 in `research/literature/` would unblock H-015's citation chain (scholar's top request). Recommend supply if you can; otherwise parked.
