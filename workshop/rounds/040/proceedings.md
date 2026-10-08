# Round 040 -- proceedings (conference)

All six personas wrote position statements (`rounds/040/<id>.md`); no new work, no referees.

## Positions
- **experimentalist:** does any gate-admitted J_i != 0 step keep the LNA key off the walks? Tabulate n = 6, 7, 8 (classes 0, 1, 2) with the key guard off. Weakest: E-138's "every J != 0 step is refused" (the key guard is also the BFS filter, so the refusal is close to built in).
- **theorist:** same question, from the algebra side (R(x) = 1 + t when j = e_i). Weakest: H-015 as used in E-134/E-137 (Coxeter guard admits exactly the tilting steps on walks; n <= 7 and one class at n = 8, capped).
- **skeptic:** same question; weakest is E-134's membership claim (unsupported, not refuted) and "d_i = 2 whenever J_i != 0" (E-131 has d = 4, 4, 5).
- **scholar:** same question, via the literature; weakest is the "silting-not-tilting iff some J_i != 0" step resting on AI Thm 2.31 cited from memory (PDF never read). Asks for a decision on fetching arXiv:1009.3370 and arXiv:2509.12983 into `research/literature/`.
- **toolsmith:** same question via a tilting-only meet with positive control and closure flag. Weakest: "a key-preserving walk stays in the same derived class".
- **maverick:** S-1: does the deletion threshold (K >= 3, 4, 5 fails first at n = 11, 13, 15) depend only on the core's key? Weakest: that threshold law (two lengths only; the n = 13 failure is partly trivial).

Disagreement: none on the main thread (five of six name the same question). Maverick alone works on S-1.

## Proposed agenda (proposed, round 040), ranked
1. **T5 / J != 0 and the key guard.** First question: over J_i != 0 rows at n = 6, 7, 8 (classes 0, 1, 2) with the key guard off, does any child keep the LNA key? Include a positive control (hand-built J != 0 child that keeps the key) and confirm the child key is computed on the child. experimentalist (tabulation, `--plan` first), skeptic (n <= 7 search, key-function audit), theorist (R(x) = 1 + t for j = e_i, stated quantifiers).
2. **T5 / tilting-only meet.** `toolsmith_n6meet.py --tilting-only` with positive control, closure flag, and replay of the other 13 hits; reverse search from the hits using J = 0 vertices. toolsmith, then skeptic.
3. **T5 / literature.** Read AI Thm 2.31/2.32 and Ladkani 2.3(c) against `tiltingPlus` on one non-monomial step; mark the 2.31 citation "from memory" until read. scholar (needs arXiv access; if blocked, say so), skeptic to referee E-128/E-136.
4. **S-1.** n = 12 key class with `--plan` (does K >= 3 fail and K >= 4 hold?), fix `maverick_single.py` label block, a null from a simple core-length statistic. experimentalist/maverick, toolsmith, skeptic.
Lower: theorist on the {h, K} dependence from H-020 before any n = 15 claim; T1/T2 remain as in the older threads.

## Questions for the steering committee
1. **Agenda:** approve the proposed agenda above (item 1 as the single priority). Recommend approve.
2. **Literature:** may the scholar fetch arXiv:1009.3370 and arXiv:2509.12983 into `research/literature/` (the proxy blocked arXiv in round 006)? Recommend yes, if the network allows; otherwise record the citation as unverified.
3. **Overnight:** none; nothing needs more than 10 minutes per command. Recommend none.

## Decisions taken for the steering committee
- round 039, question 1 (agenda): keep the round-036 agenda -- decided by the chair of round 040; no answer from the human.
- round 039, question 2 (overnight): none -- decided by the chair of round 040; no answer from the human.
