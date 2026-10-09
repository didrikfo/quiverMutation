# Round 045 -- call
kind: ordinary
Step 0.5: the round-044 questions were unanswered; decided by this chair (round-044 agenda approved with the guard audit ahead of it; no overnight; literature parked). Special request of STEERING (retrospective 1) followed: STATE.md rewritten from scratch, ledger pass for H-015 only, thread T10 "guard audit" opened as agenda item 1. No suggested question taken (S-1 had no slot; the previous call took none either).
## Assignments
- toolsmith: T10 guard audit (a). What does an added J = 0 (or `tiltingPlus`) check cost per step, measured on one n = 7 and one n = 8 guarded walk (same walk, guard as is vs. with the added check; wall time per step and total, number of steps the check refuses)? Also: does any step of those walks that keeps the key fail `tiltingPlus` (tally)? Propose where in `mutationSearchDepthFirst` the check would go; do not change the library default.
- experimentalist: T10 guard audit (b). Does any class merge recorded in `research/FINDINGS.md` at n <= 8 depend on a J != 0 step? Re-walk a recorded merge with the tilting-only guard (J = 0 / `tiltingPlus` only, key guard off) and compare the classes found with the recorded ones. Pick the merges whose size allows `--plan` within `max_command_minutes`; say which merges were and were not covered.
- theorist: T7 (breadth slot; dormant since before round 006, H-010 never assigned in a call). Either give the proof of H-010 (overlap reducible only at an end) from the procedure, using step 7 of arXiv:2112.08129 as the non-monomial test case (do not run `probe.py --steps 7`), or write a short note that proposes closing T7 with the reason (what the proof needs that is missing).
## Revisions due
(none)
## Referees
- toolsmith: skeptic
- experimentalist: skeptic
- theorist: maverick
