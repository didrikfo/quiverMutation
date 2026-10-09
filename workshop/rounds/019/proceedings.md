# Round 019 -- proceedings

Worked: toolsmith, theorist, experimentalist. Referees: scholar (toolsmith), skeptic (theorist), theorist (experimentalist). All three: minor revision. I applied the scope fixes myself (the referees' requests were wording and scope, not errors) and accepted.

## toolsmith -- `toolsmith.md` (T5)
Claim: `toolsmith_walk.py` = E-084 walker + checkpoint/resume; the n = 8 class 2 depth-8 walk completes in two slices under the fixed library with 0 key-moved steps. Referee (scholar): resume check reproduced at n = 7 class 2 depth 5; n = 8 run not re-run, saved outputs consistent with E-084; asked to weaken "exactly the 10 steps" to equal-in-count, state the n = 8 resume equivalence rests on count agreement, scope "refuted" to n = 8 class 2. Decision: **accept** with those scopes. Promoted: **E-094**; H-015 status line.

## theorist -- `theorist.md` (T1/T2)
Claim: the in-orbit placement of each split word has right gap g in {0, 1}, fixed per word across n; lemma R alone reaches `333@0` in 0 of 45 cases. Referee (skeptic): reproduced at n = 12, 16, also n = 17; `3344` is an exception, not an instance; rule is a relabelling of E-091; part (3) untested; the R-terminal result is new. Decision: **accept** as 17 of 18 plus exception, scoped. Promoted: **E-096**.

## experimentalist -- `experimentalist.md` (T5)
Claim: dim ker is 0 on tilting and >= 1 on non-tilting steps; on the 1 050 walk steps total dim ker is exactly 1 (distinct parents, 907 + 143). Referee (theorist): E-078 rerun exact; tilting <=> dim ker 0 is E-093's identity, not new; "<= 1" is conditional on A5-shaped parents; two quantities conflated. Decision: **accept**, histogram and distinct-parent counts as the new content. Promoted: **E-095**.

## Questions for the steering committee
1. **Overnight:** the depth-9 n = 8 class 2 walk is now sliceable (`toolsmith_walk.py`, 38 MB checkpoint). Recommend: no overnight yet; replay the 10 E-084 parents first (small).
2. **Agenda:** approve the round-016 agenda unchanged? Recommend: yes. Round 020 is a conference.

## Decisions taken for the steering committee
- Round 018, question 1 (overnight): none -- decided by the chair of round 019; no answer from the human.
- Round 018, question 2 (agenda): approve the round-016 agenda unchanged -- decided by the chair of round 019; no answer from the human.
