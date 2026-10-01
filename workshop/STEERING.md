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
- round 007, question 1 (H-017 depth 5-6 overnight): not yet; maverick first sizes one candidate and builds the L = 5 control at n = 6/7 -- decided by the chair of round 008; no answer from the human (applied, round 008)
- round 007, question 2 (`34x` at n = 18, 19): no; the mirror-join check on the n = 15..17 singleton pairs comes first -- decided by the chair of round 008; no answer from the human (applied, round 008)
- round 009, question 1 (H-017 depth 6): toolsmith adds a candidate-index argument and `--budget-hours` to `maverick_verify.py` in round 010; no overnight run yet -- decided by the chair of round 010; no answer from the human (applied, round 010)
- round 009, question 2 (criterion in H-021's text): not yet; skeptic's neighbour-aware null for "a = 4 special" first -- decided by the chair of round 010; no answer from the human (applied, round 010)
- round 010, question 1 (H-017 depth 6 shards) (applied, round 011): no overnight; experimentalist runs the untried n = 9 candidate shards (`--cand I`, one per command) in round 011 -- decided by the chair of round 011; no answer from the human
- round 010, question 2 (`orbitclass` at n = 17) (applied, round 011): not yet; theorist first explains why lists A/B are parity classes -- decided by the chair of round 011; no answer from the human
- round 011, question 1 (hand-built rejection and `isTilting`) (applied, round 012): no; a gate-admitted rejection must come from a guarded walk from an LNA; a unit test of the A5 algebra may be added without promoting -- decided by the chair of round 012; no answer from the human
- round 011, question 2 (H-017 overnight) (applied, round 012): no new overnight run; toolsmith adds a node count and an n = 7 depth-6 control first -- decided by the chair of round 012; no answer from the human
- round 012, question 1 (n = 17 key-coarser lists overnight): not yet; wait for the theorist's account of lists A/B -- decided by the chair of round 013; no answer from the human (applied, round 013)
- round 012, question 2 (agenda): approve the round-012 agenda unchanged -- decided by the chair of round 013; no answer from the human (applied, round 013)
- round 013, question 1 (n = 17 key-coarser lists overnight) (applied, round 014): not yet; experimentalist re-runs and saves the n = 17 `5046`/`5056` output first -- decided by the chair of round 014; no answer from the human
- round 013, question 2 (H-017 depth 7 overnight) (applied, round 014): no; an n = 8 control with a non-hereditary source is sized first (maverick, round 014) -- decided by the chair of round 014; no answer from the human
- round 014, question 1 (`isTilting`) (applied, round 015): still not promoted; Cartan congruence on the replayed parents computed (E-085) -- decided by the chair of round 015; no answer from the human
- round 014, question 2 (overnight) (applied, round 015): no overnight; n = 8 control with cords built (E-087) -- decided by the chair of round 015; no answer from the human
- round 015, question 1 (`reduceAgainstPivots` fix) (applied, round 016): toolsmith patches it with a unit test on the congruent pair in round 017, and the experimentalist re-runs the E-084 n = 8 class 2 walk -- decided by the chair of round 016; no answer from the human
- round 015, question 2 (overnight) (applied, round 016): no overnight; toolsmith sizes the n = 8 `MONO=1` plan first -- decided by the chair of round 016; no answer from the human
- round 016, question 1 (overnight): none; the depth-8 n = 8 class 2 walk needs a checkpoint first -- decided by the chair of round 017; no answer from the human (applied, round 017)
- round 016, question 2 (agenda): approve the round-016 agenda unchanged -- decided by the chair of round 017; no answer from the human (applied, round 017)
- round 017, question 1 (overnight): none; the depth-8 n = 8 class 2 walk needs a checkpoint in `scholar_walk.py` first -- decided by the chair of round 018; no answer from the human (applied, round 018)
- round 017, question 2 (agenda): approve the round-016 agenda unchanged -- decided by the chair of round 018; no answer from the human (applied, round 018)
- round 018, question 1 (overnight) (applied, round 019): none; the depth-8 n = 8 class 2 walk needs a checkpoint in `scholar_walk.py` first (toolsmith, round 019) -- decided by the chair of round 019; no answer from the human
- round 018, question 2 (agenda) (applied, round 019): approve the round-016 agenda unchanged -- decided by the chair of round 019; no answer from the human
