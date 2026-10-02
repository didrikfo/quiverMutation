# Round 025 -- proceedings

Worked: scholar (revision), skeptic, experimentalist. Referees: experimentalist (scholar), theorist (skeptic), skeptic (experimentalist).

## scholar (revision of 023) -- accept
Claim: at n = 8 class 0 the walk-level "reject iff long square" fails; 61 distinct out-degree 2 rejects, each commuting into one arrow and killed into the other by a zero relation, all failing the Cartan test through the real rewrite. Referee reproduced the 61 and checked all of them by an independent per-arrow kernel computation. Decision: accept; wording fixes (all-61 check, two separate runs, n = 9 coverage) are small and carried into E-105. Promoted: **E-105**, glossary "D' reject", limit note added to E-100.

## skeptic -- revise
Claim: the 42 reproduce, are walk-reached and are "the E-103 kind, minimal presentations". Referee (theorist): numbers reproduce, but "key equals base" is vacuous (BFS), "same mechanism" is a shape match (control: 754 same-shape steps accept vs 42 reject), and the title's "minimal" contradicts the author's own caveat. Decision: revise (due round 026: extract x for the 42, test the shorter presentation). The runs and the control are recorded in **E-106**.

## experimentalist -- note (runs recorded; claim rejected as stated)
Claim: out-degree 2 rejects recur at n = 8 class 1, not at class 3 or n = 9 class 0. Referee (skeptic): counts are load-dependent (15 vs 38 at class 1), "a second way to break the iff" is wrong (the two out-degree 1 rejects are parallel-arrow long squares missed by `longSquare`), "out-degree 2" means >= 2, the zeros support nothing, E-103 already records most of it. I agree and do not ask for a revision: the only new content, the class 1 recurrence, is in **E-106**.

## Questions for the steering committee
1. **Agenda:** keep the round-024 agenda; item 1 now turns on why D' appears at n = 8 and a test of walk reachability beyond the BFS prefix (recommend yes).
2. **Overnight:** a checkpointed n = 9 class 0 walk with the reject classifier would exceed 10 min; not proposed until the toolsmith adds a checkpoint (recommend none yet).

## Decisions taken for the steering committee
- Round 024, question 1 (agenda): approved unchanged -- decided by the chair of round 025; no answer from the human.
- Round 024, question 2 (overnight): none -- decided by the chair of round 025; no answer from the human.
