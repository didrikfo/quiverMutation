# Round 010 -- call
kind: ordinary
## Assignments
- toolsmith: Add a candidate-index argument and `--budget-hours` (exit 2 when spent) to `workshop/rounds/009/maverick_verify.py` (or its committed successor), keeping its behaviour otherwise; run it for one n = 9 candidate at depth 5 to confirm it reproduces E-074, and size depth 6 for one shard. (thread T6; answers STEERING round 009 q1)
- skeptic: Design and run a neighbour-aware null for E-073's "the seed `aaa` collapses to `34` only for a = 4": ask whether a = 4 is special among a = 2..9 compared with what other seeds/neighbour pairs do, so that the claim is either supported against a null or withdrawn. (thread T2; answers STEERING round 009 q2)
- experimentalist: Print the key-coarser core lists (cores where orbit+mirror is strictly finer than the key) at n = 14, 15, 16 for `--max-word 4`; does the even-n list repeat across 14 and 16, the odd across 13, 15, and is size 20300 of `348/349` at 16 the `4056` orbit? (thread T1/T3)
## Revisions due
(none)
## Referees
- toolsmith's submission: theorist
- skeptic's submission: experimentalist
- experimentalist's submission: skeptic
