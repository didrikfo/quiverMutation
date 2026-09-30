# Steering

Owned by the human. The chair reads this first every round, and it outranks
everything else, including `config.yaml`. Edit freely; keep it short.

status: active
<!-- set to `paused` and the round stops after reading this file -->

## Direction

Work on the open hypotheses in `research/HYPOTHESES.md`, preferring the ones
the most recent findings bear on (H-020, H-021 at the time of writing). Aim for
results that would survive into `research/FINDINGS.md`: specific, reproducible
by a command, and checked against what is already recorded.

## Priorities

1. Settle or sharpen an open hypothesis rather than start a new line.
2. Before claiming anything, check it is not already in `research/` (the
   project's main cost is rediscovery).
3. Proofs or mechanisms for things that are already well supported by data.

## Avoid

- Runs longer than `max_command_minutes`. Propose them for `OVERNIGHT.md`
  instead.
- Refactoring code that no experiment needs.

## Special requests for the next round

<!-- e.g. "conference", "only skeptic and theorist", "everyone on H-021" -->
(none)

## Answers to the chair

<!-- answer questions from the latest proceedings here; the chair reads them.
     Unanswered questions are decided by the next chair and recorded here. -->
- round 001, question 1 (H-021) (applied, round 002): settle the **strict** reading of "mirror"
  first (the mirror of a placement at a different offset from those the orbit
  holds). Restate H-021 only if its failure persists under that reading.
  The rebuilt census script, `workshop/rounds/001/experimentalist_census.py`,
  records both readings.
- round 001, question 2 (overnight): yes. Both are in `OVERNIGHT.md`, Menu 4:
  the `n = 14` census (4 shards) and the Ladkani audit at `n = 9, 10`. The
  human runs them; results come back as E-entries.
- round 001, question 3 (`tiltingPlus`) (applied, round 002): wait for a second non-monomial
  negative example before promoting it to the library.
- round 002, question 1 (`isTilting`) (applied, round 003: not promoted): hold off on promoting it to the library until the chair judges the evidence solid enough to decide (e.g. a gate-admitted rejection from the overnight audit). Answered by the human.
- round 002, question 2 (H-021) (applied, round 003): restate H-021 without the mirror clause; theorist restates it in round 003, experimentalist runs the 7 survivors at n = 14. Answered by the human.
- round 003, question 1 (censuses n = 12, 14): yes, added to `OVERNIGHT.md` Menu 4 -- decided by the chair of round 004; no answer from the human (applied, round 004)
- round 004, question 1 (H-017 overnight): approve the n = 9 depth-7 run only, not n = 10 until a positive control exists -- decided by the chair of round 005; no answer from the human (applied, round 005)
- round 004, question 2: keep round 005 a conference -- decided by the chair of round 005; no answer from the human (applied, round 005)
- round 006, question 1 (agenda): approve the round-005 agenda unchanged -- decided by the chair of round 007; no answer from the human (applied, round 007)
- round 006, question 2 (overnight): no overnight run proposed yet; toolsmith sizes `--max-word 5` at n = 14 first -- decided by the chair of round 007; no answer from the human (applied, round 007)
