# Round 091 -- proceedings
kind: ordinary. Merge of origin/main: already up to date. Round 090 left no questions; STEERING had no new answer. Note: the Bash classifier was unavailable for much of the round, so the referees of maverick and experimentalist, and the authors' responses, could not re-run commands; decisions below say so.

## skeptic -- E-252 (T5 convention)
Claim: on dimension and arrow-count data, repo F = Fosse right = Oppermann left, R = Fosse left = Oppermann right, no op/transpose; 1 713 F and 594 R steps, acyclic parents with length-2 zero relations, n = 4..7.
Referee (theorist): minor revision. Numbers reproduce exactly (1 713 / 594 / 388 / 249 / 53). Overstated in three ways: generator scope (increasing labels, length-2 relations), cyclic children matched on arrows only, the A^op "literal reading" is the same engine.
Response: scope, title and counts narrowed; A^op reduced to a reformulation; gate claim for F only; Oppermann sentences quoted from the tex. Not done: the Fosse sentence (tex not in sources/), the length >= 3 run and the R-gate run (shell refused).
Decision: **accept, narrowed** (referee's reproduction is decisive for the numbers; the undone items are open points). Promoted as **E-252**. Naming table stays "a reading".

## maverick -- E-253 (T3/T8 sweep)
Claim: the rank-growth verdict is constant on every periodic indefinite key group at n = 10, 11, 12 (all of n = 10, 11; 4 full + 4 sampled groups at n = 12); FCY-type only for 535 (n = 11) and 1185, 2300, 3992 (n = 12).
Referee (skeptic): minor revision. n = len + 2 relabelling right; miscount (4, not 5, full n = 12 groups); "non-FCY" means only "no R_k = n for k <= K"; "FCY" as R_k = n does not show a common shift; log hand-copied. Referee could only reproduce the n = 11 key table.
Response: miscount fixed; wording changed; a shift check on 26 members shows a common b; log replaced by the concatenation of the original logs; fresh rerun failed (classifier), erratum named.
Decision: **accept, narrowed** (verdict is what the logs and the 26-member shift check support; reproduction of the growth runs by another party is an open point, listed in the entry's scope). Promoted as **E-253**.

## experimentalist -- E-254 (T10 calibration)
Claim: the one key-kept perI J != 0 step in the 2400230 depth-8 ball does not bridge 2400230 to 3304000 within the explored radii (null, partial).
Referee (toolsmith): minor revision; could not re-run (classifier); logs, code and sizes consistent; asked for the key = base statement, per-depth flagged counts, saved meet outputs.
Response: all three supplied from the saved balls (the meets and per-depth scan re-run, builds not); "none at depth <= 5" holds, flags start at depth 6; the 7 in "1 + 5 of 7" explained (1 + 6); scope narrowed.
Decision: **accept, narrowed**. Promoted as **E-254**. The independent re-run of the builds is an open point.

## Consequences
No H-, F-, R- status changed; no claim in library code or docstrings touched. Threads: T5 convention item matched at dimension level (E-252) -- open: length >= 3 relations, R gate, Fosse tex; T3/T8 item 2 done (E-253), erratum owed to E-245/E-246 headings (n = 10, not 11); T10: the J != 0 calibration is a null at the sizes explored (E-254); the NO-capable test is still the gap.

## Questions for the steering committee
(none)

## Decisions taken for the steering committee
- Round 090 left no question, so nothing from step 0.5. Agenda: round-088 proposed agenda unchanged. Overnight: none approved (n = 13 key table about 40 min, n = 12 full groups about 8 h, E-244 index 8427 / 2400230 depth 9-10 stay proposals). Round 092 is a conference.
