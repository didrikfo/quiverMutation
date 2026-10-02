# Round 026 -- proceedings

Ordinary. Called: skeptic (revision), theorist, toolsmith. Referees: theorist (skeptic), skeptic (theorist), experimentalist (toolsmith). Questions from round 025 were unanswered; the chair decided them (below).

## skeptic (revision of 025) -- accept, with the referee's corrections
Claim: in all 42 n = 8 class-0 rejects x is a two-term sum killed termwise by one out-arrow and only as a sum by the other; the reject survives the shorter presentation (15 of 42 reducible). Referee (theorist): minor revision, reproduces; marginally new (E-105/E-106 give the cases); fix "a kernel element", say the x-form was not tested on accepting rows, relabel "46 of 42". Decision: accept; the fixes are wording and I applied them in the entry. Promoted as E-108.

## theorist -- accept as a conjecture, not a result
Claim: on 17 802 out-degree 2 rows (n = 8 c0, c1, n = 7 c0) reject iff witness W (two-term relation through v, x = p1 - p2 != 0, killed by the other arrow), 0 mismatches. Referee (skeptic): minor revision; class 1 reproduces (15, 0 mismatches), class 0 and n = 7 not re-run; W => reject is trivial; the converse rests on 70 non-independent rejects of one shape and the "cancels" branch is never exercised. Decision: accept as E-107 headed "Conjecture W", limits as the referee states. Not promoted to HYPOTHESES (empirical only).

## toolsmith -- accept
Claim: `longSquare` fixed for parallel arrows (4 of 16 out-degree 1 J != 0 steps at n = 8 c1 were missed), resumable reject walk with `--ckpt`/`--budget-hours`, n = 9 class 0 sized. Referee (experimentalist): minor revision; tests pass, resume check identical, `--budget-sec` exits 2; n = 6, 7 sweeps vacuous; 4 positives unverified for irredundancy; penultimate-arrow rule is a definition choice. Decision: accept; caveats recorded. Promoted as E-109. No library change (the test was only in workshop scripts).

## Questions for the steering committee
1. **Overnight, n = 9 class 0 reject walk** (`toolsmith_rejwalk.py 9 --class 0 --budget-hours 7 --ckpt ...`, multi-night, may never close: new algebras grow by about 2.5 per level). Recommend no yet; first test W on a short n = 9 prefix (theorist/experimentalist).
2. **Agenda:** keep the round-024 agenda; item 1 now has a candidate rule (W) to test on classes 2-3, n = 9, parallel arrows, and by an independent proof of the converse. Recommend yes.

## Decisions taken for the steering committee
- Round 025, question 1 (agenda): kept unchanged -- decided by the chair of round 026; no answer from the human.
- Round 025, question 2 (overnight): none; toolsmith built the checkpoint first -- decided by the chair of round 026; no answer from the human.
