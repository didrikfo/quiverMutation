# Retrospective 4: cycle 2 (rounds 121-140)

*2026-11-28 · written after the round-140 closing conference by the session
that runs the workshop's schedule, as `CYCLE.md` asks. Cycle 2 was the first
with no predefined phase focus, researchers and referees on opus, and one
round per day. A copy is in `workshop/retrospectives/004-cycle-002.md`.*

## The numbers

| | cycle 1, rounds 101-120 | cycle 2, rounds 121-140 |
|---|---|---|
| rounds | 20 (15 ordinary, 5 conferences) + planning, 8 days | 20 (15 ordinary, 5 conferences) + planning, 2026-11-07 to 11-27 (21 days) |
| models | sonnet throughout | researchers and referees opus; chair and conferences sonnet |
| reviews | 45: 42 minor, 3 major, 0 accept | 45: 38 minor, **4 accept**, 3 major, 0 reject |
| E-entries | 43 (E-272..E-314) | 44 (E-315..E-358) |
| H-entries written | 8 | 1 (H-032, at planning) |
| status changes | 3 (all to SUPPORTED) | **5**: H-015 -> REFUTED (122, R-017; F-057); H-027 -> CONFIRMED (128, F-056); H-029 -> SUPPORTED (130; widths 8-9 retired at 140); H-025 -> CONFIRMED (136, F-058); H-030 -> SUPPORTED and closed (137) |
| F-entries | 0 | **4**: F-055 (124), F-056 (128), F-057 (132), F-058 (136), each scoped in its title |
| R-entries | 3 (R-014..R-016) | 2: R-017 (a key group is not one derived class), R-018 (the gate does not decide tilting on cyclic monomial rows) |
| library changes | docstrings, 5 regression tests | `mutationSearchDepthFirst` docstring; opt-in `refuseNonTilting` (default off, decided at 140); a warning on cycle-upstream gate calls |
| questions to the human | 1 | 0 (PDFs still wanted, not needed) |

**Cost (session usage, real figures).**

| | per session | sessions | total |
|---|---|---|---|
| ordinary rounds, phases 1-4 | $13.03-$20.72, mean $17.25 | 12 | $207.03 |
| ordinary rounds, phase 5 (replication) | $15.44, $12.99, $10.69 | 3 | $39.12 |
| conferences (sonnet) | $2.18-$2.81, mean $2.43 | 5 | $12.14 |
| planning step | $4.90 | 1 | $4.90 |
| **cycle** | | 21 | **$263.19** |

All 15 ordinary rounds: mean $16.41 (range $10.69-$20.72), against about $5 on
sonnet in cycle 1. The cycle cost about three times cycle 1 ($86). The
ordinary rounds are nearly all of it: conferences and planning are 6%.

## The plan's progress outcomes

| outcome (`cycles/002/plan.md`) | verdict |
|---|---|
| T10: a recorded count deciding whether n = 9 can say NO; the 47-group answered | **met, and more.** E-316: n = 9 cannot test T10; E-317: the 47-group is 9 census orbits, none mixed. Then a cheaper instrument than the depth-6 ball found the crossing at n = 11 (E-318, E-323, E-332, E-333): H-015, open since round 045, is refuted, and F-057 says so at its scope. The depth-6 ball the cycle-1 retrospective asked the human to approve turned out unnecessary. |
| H-023: an F-entry scoped to library steps, or the blocking lemma stated and tested | **not met.** The child = End(T) test holds at 587 of 587 J = 0 steps at n <= 6 (acyclic monomial; E-345), with power shown on planted wrong children (E-348); the (LF)/(H1w)/(H3) lemmas are proved on their models (E-319, E-322). But no F-entry, and the plan's failure line, "conditional on the omega-reduction still standing at 140", is true. H-023 has now been the weakest claim for two cycles. |
| H-026: an F-entry at class level, or a counterexample; H-032 written with a verdict | **not met.** E-325 and E-349 describe how keys change (only through K = 0 bridges), but there is no class-level F-entry and no counterexample. H-032 was written, but it is still OPEN. |
| H-030 closed; every F-entry rebuilt by a non-author | **met.** E-351 shows the injectivity test equals Hom, and H-030 is closed. Phase 5 rebuilt all four F-entries independently (E-353..E-357), with no retraction. |

That makes two of four met, and the cycle's biggest result (F-057) was the
one the plan expected least of.

## Did the format do its job?

**Opus for researchers and referees: yes, at about three times the cost.**
- **Output.** The cycle wrote four F-entries; the first 120 rounds wrote
  none. It also refuted the workshop's oldest open hypothesis.
- **Review.** Review accepted outright for the first time (4 of 45). It
  still returned 3 majors, so it discriminates rather than rubber-stamps.
- **Computation.** Submissions carried their own Hom-level checks and
  planted controls without being asked: E-333, E-334, E-348.
- **Caveat on attribution.** This cannot be separated cleanly from the
  scoped-findings rule (below), which came in at the same time. But the
  rule only says *when* to write an F-entry; the results under F-055 to
  F-058 were new work.

**Scoped findings, with a blocker named once (cycle-1 change 4): works.**
- Every conference wrote an F-entry whose blocking check had been done
  (F-055 to F-058), with the next blocker as its open point.
- At 140, the candidate "option is harmless in-class" was retired rather
  than renamed, as the rule asks.

**No predefined focus, with a partial literature focus in phase 1: fine.**
- **The literature part (121-123).** It was small and cheap. It settled a
  real question: no published splitting lemma covers H-031 (E-315), and the
  ENS cancellation needs (T)+(P) (E-320). It did not crowd out research:
  phase 1 also found the first crossing (E-318).
- **Phases 2-4.** These were plain research and moved T10, H-025 and H-027
  to findings.

**Backlog and replication in phase 5 (137-139): yes. This is the cycle's
format lesson.**
- **What it did.** Each F-entry was rebuilt by someone who had not written
  it, with that person's own code: F-055 by the theorist, F-057 by the
  experimentalist and the toolsmith, F-058 by the skeptic and the maverick.
  F-056's literature route was re-sourced.
- **What it found.**
  - No finding fell.
  - Every rebuild sharpened a finding's stated scope:
    - F-055's NO side holds only through q_S (E-356);
    - nu2's uniserial returns compare dimension vectors, not modules
      (E-355);
    - F-057's NO still rests on E-308/E-272 (E-354).
- **Cost.** These were the cheapest ordinary rounds of the cycle ($10.69 to
  $15.44). Rebuilding is cheaper than discovering.
- **The backlog half.** It settled `refuseNonTilting` (default off, with a
  pre-registered flip rule), retired H-029's widths and four non-load-bearing
  citations, and wrote the owed 1310.1907 file.

**The closing conference (cycle-1 change 5): works.** `cycles/002/backlog.md`
marks every open item settled, retired, converted or carried, with reasons.

## What did not work

- **Code consequences still carry over.**
  - Two small library changes have been owed since round 127 and are
    "carried" again at 140:
    - the `isTiltingAtVertex` illegal-relation refusal (E-351);
    - the `isMutable` refusal on an upstream cycle (R-018).
  - The closing backlog makes this visible, but nothing forces it. The
    rule for the cycle-1 backlog phase ("done in this phase or declined
    with a reason") was not applied, because phase 5 was a combined
    backlog-and-replication phase and the replication won the slots.
- **H-023 has no proof push that lands.**
  - The plan gave the theorist a proof push in phases 2-3. It produced
    lemmas and tests but not the End(T) step.
  - A thread that is everyone's weakest claim for 40 rounds may need a
    phase of its own (a *problem session*: every slot on one question,
    with the conference judging whether a proof exists), not one persona's
    third priority.
- **Long runs are never decided.**
  - The cycle declined the depth-8 T10 run and the E-306 78-of-78 rebuild
    (~40 min), and S-1 at n = 17 needs an int64-key overnight run. The
    reasons were sound: a null could not change the default.
  - But no chair can approve a run longer than 10 minutes, so these items
    can only be carried or retired. This is a standing limit, not a cycle-2
    problem. The human approved such a run once (the 7 + 7 join, round
    090).
- **Few new hypotheses.** One H-entry, against eight in cycle 1. Partly this
  is healthy: the cycle closed what it opened. But the cycle-3 agenda
  (H-023 at algebra level, F-055 without q_S, F-058's two tau-orbits, T10 at
  n = 12) is entirely follow-up.

## Changes for cycle 3 (made now)

1. **Replication is a named focus.**
   - `CYCLE.md` describes it: each assignment rebuilds one result by a
     non-author with independent code and states what scope survived.
   - The planning step is told that a replication part in the last phase
     worked in cycle 2, and should consider it whenever the cycle wrote
     F-entries.
   - `README.md` lists the focus.
2. **Owed code goes first.**
   - A library change that a closing backlog carries goes to the first
     toolsmith slot of the next cycle's first round. The planning step
     must schedule it in `plan.md`.
   - A change carried twice is either made by the conference chair of the
     next phase or declined in `research/` with a reason. It is never
     carried a third time.
   - This covers the two changes owed since round 127.
3. **The planning step may give a phase a problem-session focus.**
   - `CYCLE.md` describes it: every slot works on one question named by
     its H-entry. Each round's call says which step of the argument or
     which experiment each persona takes. The phase's conference decides
     proved, refuted, or the exact missing step.
   - The plan decides whether to use it. H-023 is the obvious candidate,
     but this retrospective does not impose it.

Not changed:
- models;
- one round per day;
- five phases of four rounds;
- three researchers and one referee;
- the 10-minute command limit.

The cost per round is the human's decision (below).

## For the human

- **Cost.**
  - The figures:
    - Opus personas cost about $16.40 per ordinary round, against about
      $5 on sonnet.
    - A conference costs about $2.40.
    - A cycle costs about $263 over three weeks.
  - What it bought: four findings and a refuted hypothesis, against none
    in the 120 sonnet rounds before.
  - If the rate is too high, there are two levers:
    - fewer ordinary rounds per phase (3 + conference becomes 2 + conference);
    - referees back on sonnet. Referees are a minority of each round's
      cost; this lever is untested.
  - The format can stay as it is either way.
- **Long runs.** Three items wait on runs longer than a round allows:
  - the E-306 rebuild (~40 min on 4 cores);
  - S-1 at n = 17 (overnight);
  - T10 at n = 12 depth 8 (not yet sized).

  None is urgent. Approving one in `STEERING.md` is the only way they move.
- **PDFs.** Nikulin 1979 (Cor. 1.13.3) is the one that matters: it
  underwrites every Cartan-level NO (F-055, F-057). Anick 1982 is the
  second.
- **Checkpoint to `main`.** It is still outstanding. Since checkpoint 1
  the record has added about 200 E-entries, 4 F-entries and 5 R-entries.
  Library changes to review are the `CyclicDescendantsError` refusal, the
  docstrings, the regression tests and the opt-in `refuseNonTilting`.
