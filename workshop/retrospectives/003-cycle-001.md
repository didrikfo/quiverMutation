# Retrospective 3: cycle 1 (rounds 101-120)

*2026-11-06 · written after the round-120 conference by the session that
runs the workshop's schedule, as `CYCLE.md` asks. It looks at cycle 1 as a
research process; the format changes at the end carry out the human's plan
of 2026-11-03 (`STEERING.md`, "Plan after round 120") and are already in
`CYCLE.md`, `ROUND.md`, `config.yaml` and `README.md`. A copy is in
`workshop/retrospectives/003-cycle-001.md`.*

## The numbers

| | rounds 045-100 (retro 2) | cycle 1, rounds 101-120 |
|---|---|---|
| rounds | 56 (42 ordinary, 14 conferences) | 20 (15 ordinary, 5 conferences) + planning step, 2026-10-30 to 11-06 |
| reviews | 126: 113 minor, 10 accept, 3 major, 0 reject | 45: 42 minor, 3 major, 0 accept, 0 reject |
| E-entries | 122 | 43 (E-272..E-314) |
| H-entries written by the workshop | 0 | **8** (H-023..H-027, H-029..H-031; H-028 retired as a citation, not a hypothesis) |
| status changes | 1 | **3**: H-024 OPEN -> SUPPORTED (101), H-027 -> SUPPORTED (107), H-026 -> SUPPORTED (113, narrowed to odd n in 117); H-025 written SUPPORTED |
| R-entries | 0 | **3**: R-014 (E-218's "only at n = 9, 10"), R-015 (a writer bug in the round-097 scan), R-016 (a harness bug in round 118) |
| F-entries | 0 | **0** |
| library changes | 3 | docstrings of `search.py`, `fingerprint.py`, `procedure.py` (E-311, round 120), 5 regression tests (E-312) |
| questions to the human | ~25 | 1 (PDFs, not blocking), asked three times |
| cost (session usage) | ~$4.3 / ordinary, ~$0.7 / conference | ~$5.0 / ordinary (3.3-7.6), ~$1.9 / conference (sonnet instead of haiku), $1.7 planning; ~$86 for the cycle |

## What the cycle did, by line

- **Cartan-level separation (H-024).** The cycle's clearest progress.
  - The necessity step was sourced in round 101 (Ladkani, E-272) and
    completed with Rickard and Aihara-Iyama in round 117 (E-308).
  - Controls that could fail were run: YES on Fosse's mixed-digit pairs
    (E-281) and NO on 19 of 19 q_S-separated pairs (E-292). A planted test
    gave NO on 21 of 21 in one n = 11 group (E-314).
  - The open gap is the absence of any planted NO at λ < 1.5, and power
    shown at only one group.
  - A new R-candidate, "a key group is one derived class", came out of
    E-314: the n = 12 group of 47 splits 14 + 33.
- **A_n ~ B_n (H-027, R-014).** This came from literature work: Fosse
  Cor. 4.6 covers every n ≥ 5. A_n was then reached by tilting outside the
  workshop's move orbit (E-287, E-289), and its generation was sourced
  (E-308).
- **The J = 0 premise (H-023, H-031).**
  - Literature reframed the premise: Oppermann needs only "no loop of
    degree 0 at i" (E-273). A loop of degree ≥ 1 at i forces a cycle
    (E-279).
  - The remaining gap is the workshop's own Lemma A for loop letters. Four
    formal models are proved (E-290, E-294, E-309), but the library's
    Hom/End(T) is not.
- **T10, does a guarded walk stay in one class? (H-015).** Still stuck,
  but now for a stated reason.
  - Every test at n ≤ 7 says YES by construction: the Hom chains (E-297,
    E-300), the Cartan test (E-295) and the HH test (E-305).
  - The first instrument that can say NO exists only at n = 9 (q_S splits
    F-010 18 + 18, E-306). Using it needs a depth-6 ball, which is longer
    than a round allows.
- **K-threshold law (H-026).**
  - Pre-registered predictions held at n = 21 (E-299) and at three more
    cells (E-304).
  - One failed at n = 14 (E-307): the law holds at key level for odd n.
  - The rule was corrected off-centre (E-301), and a writer bug was
    retracted (R-015).
- **Library.**
  - The `search.py` docstring overclaim, open since round 045, is fixed
    (E-311).
  - An apparent library bug was a harness artefact (R-016). The real
    residue is 8 over-admissions in 936 minimal presentations (H-030).

## Did each phase's focus do its job?

**Phase 1, literature (101-104): yes, clearly, and more than any other phase.**
- **What it produced.**
  - H-024's status change.
  - The Oppermann reframing of H-023.
  - Fosse Cor. 4.6, which led to H-027 and R-014.
  - Closed-form E-269 bounds.
  - The A_n/B_n correction that phase 2 built on.
- **What it changed downstream.** Phases 2-4 ran no join searches, which
  every plan rejected after phase 1. Phase 2's first control (E-281) used
  the pairs phase 1 found.
- **Caveats.**
  - Part of the value came from a backlog of literature neglected for
    100 rounds; a second literature phase would find less.
  - The scholar went on doing literature work in research phases (E-285,
    E-308, the round-110 recheck). That pattern is cheaper than a whole
    phase.
- **Unreachable.** Four citations stayed UNVERIFIED for want of PDFs
  (Nikulin 1979, Happel–Seidel, BGS, Anick 1982).

**Phases 2-4, research (105-116): mixed.**
- **Under the new rules.**
  - Each conference named what blocked a finding, and the next round
    worked on it.
  - Pre-registered predictions and planted controls appeared, and twice
    they failed honestly: the n = 14 prediction (E-307) and the planted
    control with no power (E-291).
  - The *Power* line produced "(no power)" titles instead of silent
    overclaims: E-288, E-291, E-313.
- **Moving targets.**
  - The blocker rule turned into moving goalposts. Each blocking check,
    once done, revealed the next one, so no F-entry was written:
    - H-024 was blocked in turn on mixed pairs, then a planted NO, then
      a second group.
    - H-026 was blocked on a re-run, then class level, then a rival
      formula.
  - Results are scoped carefully, but none is ever declared done at its
    scope.
- **Threads.** T10 took a slot in most rounds and stayed YES-side for
  12 rounds before the round-116 conference said so plainly.

**Phase 5, backlog (117-120): yes.**
- **The list.** Of 19 items, 13 were settled, 2 retired (with reasons) and
  2 converted into hypotheses; 5 carried over by name.
- **Long-pending items.** It closed the consequence items that had been
  pending since rounds 045-050: the guard and `canonicalKey` docstrings,
  the `allowCyclic` KeyError, and the reruns after the
  `CyclicDescendantsError` change.
- **Errors caught.** It found and retracted the workshop's own errors
  (R-016) and retired four overnight proposals that had sat undecided for
  up to 50 rounds.
- **Why a phase.** None of this happened in 50 earlier rounds that had the
  same rule ("consequences for the record") but no phase that had to act
  on it.

**Conferences.**
- On sonnet, positions still converge: four or five of six name the same
  question.
- They also argue now:
  - The skeptic called H-027's support guaranteed by Fosse (round 108).
  - Personas named the result closest to a finding each time.
- The ledger made 3 status changes and 3 retractions in 5 conferences,
  against 1 in 14 before.

## What worked, what did not

| | verdict |
|---|---|
| live questions as H-entries (retro 2, change 1) | **works**: 8 written, 3 moved; the ledger has something to move |
| ledger names candidates and blocking checks (change 2) | **works for direction, not for closure**: no F-entry, because the check after the blocker is always another blocker |
| *Power* line (change 3) | **works**: "(no power)" is now said in titles; planted controls became routine |
| `STATE.md` byte limit (change 4) | **works**: under 12 KB every round, no rediscovery of board items this cycle |
| sonnet conferences (change 5) | **works modestly**: disagreement and "closest to a finding" lines, at ~$1.2 more per conference |
| no agenda question (change 6) | **works**: one question to the human in 20 rounds, and it was a real one |
| backlog phase for pending consequences (change 7) | **works**: see phase 5 |
| literature phase | **worked this time**, largely as a catch-up |

Not working:
- **No F-entry in 120 rounds.** The workshop's results are by now
  well-scoped and controlled, for example:
  - H-024 at λ ≥ 1.661 on the groups tested;
  - the H-026 law at key level for odd n;
  - the E-187/E-205 proofs of retrospective 2.
  The rule treats a finding as something with no open check, and every
  result here has one.
- **T10 is out of reach with 10-minute commands.** Its only instrument
  needs a depth-6 ball at n = 9, and the dispatcher ran such a job once
  (7 + 7, 70 min). The chair can only propose such runs, and none was
  approved this cycle.
- **Review never accepts outright.** None of 45 reviews was an accept;
  every submission needed a same-round response. That is affordable, but
  it means review rarely distinguishes a solid result from a weak one.

## Changes for cycle 2 (made now)

From the human's plan (`STEERING.md`, 2026-11-03):

1. **No predefined phase focus.**
   - `config.yaml` now lists five `research` phases.
   - The planning step may give a phase, part of one, or some personas a
     *literature* or *backlog* focus; it writes the choice into `plan.md`
     under *Phase focus*, and the round chair follows it.
   - The descriptions of both focuses stay in `CYCLE.md` for when the plan
     uses them.
   - The evaluation above suggests the planning step should consider one
     backlog phase (or a backlog half of phase 5) per cycle; it is told
     so, but not bound.
2. **Opus for researchers and referees.** Conferences and the chair stay on
   sonnet. This is now set in `config.yaml`.
3. **One round per day**, at 00:43 UTC; the scheduled task is updated.
   Cycle 2 (rounds 121-140) takes about three weeks.

From this retrospective:

4. **Scoped findings.**
   - At a conference, a candidate whose blocking check was done is
     written as an F-entry *at the scope it reached*, with its title
     carrying the scope. The next blocker becomes the finding's open
     point, not a reason to wait.
   - A blocker may be named only once per candidate: the second time, the
     chair writes the F-entry at its current scope or retires the
     candidate.
5. **The cycle's last conference always marks the backlog.** It covers
   every open item in `STATE.md`: settled, retired, converted, or carried
   over with a reason. This holds even when no phase had a backlog focus,
   so pending consequences cannot sit silently again.

Not changed:
- personas, `researchers_per_round` (3), one referee;
- the same-round response, decide-and-record;
- the 10-minute limit;
- phase length (4) and the number of phases (5).

## For the human

- **Cost to watch.** Opus for researchers and referees will raise the cost
  per ordinary round from about $5. Cycle 2's first phase summary
  (`cycles/002/phase-1.md`) should give the real figure. Check it after a
  few rounds, then adjust the daily rate.
- **T10 needs one longer run.** The round-120 chair asked for the depth-6
  ball of 3060000 to be sized with `--plan` first. If the estimate is a
  few hours, approving it, as you did for the 7 + 7 join, is the only way
  this thread moves.
- **Four PDFs.** Nikulin 1979, Happel–Seidel, BGS and Anick 1982 would
  settle the four citations still marked UNVERIFIED.
- **Checkpoint to `main`.** The record has grown by about 150 E-entries,
  8 H-entries and 3 R-entries since checkpoint 1. Library changes to
  review include the `CyclicDescendantsError` refusal, the docstrings
  (E-311) and the regression tests (E-312).
