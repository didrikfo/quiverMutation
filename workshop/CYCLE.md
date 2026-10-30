# Cycles and phases

The workshop's rounds are grouped into **cycles**. A cycle is a planning step,
five **phases** of four rounds each (three ordinary rounds and a conference),
and a cycle retrospective. Each phase has a **focus** that changes what the
chair assigns and what the conference is about. This file is the procedure
for the parts that are not a plain round; `ROUND.md` says where a round uses
it. Introduced by retrospective 2 (`retrospectives/002-rounds-045-100.md`);
every part of it is meant to be changed by a cycle retrospective.

```
planning step ── phase 1 ── phase 2 ── phase 3 ── phase 4 ── phase 5 ── cycle retrospective
 (not a round)   literature  research   research   research   backlog    (not a round)
                 3 + conf    3 + conf   3 + conf   3 + conf   3 + conf
```

## Where a round is in its cycle

`config.yaml` has a `cycle:` block: `first_round` of cycle 1, `rounds_per_phase`
and the list of `phases` (their focus, in order). For round `r`:

* cycle `c = (r - first_round) // (rounds_per_phase * len(phases)) + 1`;
* position `k = (r - first_round) % (rounds_per_phase * len(phases))`;
* phase `p = k // rounds_per_phase + 1`, focus `phases[p - 1]`;
* the round is the phase's conference when `(k + 1) % rounds_per_phase == 0`
  (with the default numbers this is the same as "a multiple of
  `conference_every`").

`STATE.md` says the same in its header (`cycle`, `phase`, `focus`), so a persona
need not compute it. The cycle's files live in `workshop/cycles/CCC/`
(three digits): `plan.md`, one `plan/<persona>.md` per persona,
`phase-P.md` (the chair's short summary written at each phase's conference),
and `retrospective.md`.

## Planning step (before the first round of a cycle)

A separate session, started by the dispatcher when `STATE.md` says
`next_round_kind: planning`. It does no research and runs nothing longer
than a grep.

1. Set up as in `ROUND.md` step 0 (branch, merge `main`, `STEERING.md`, pause
   switch). Make `workshop/cycles/CCC/plan/`.
2. Brief every persona in the roster at once, in parallel,
   `model: <models.researcher>`:

   > You are **<Name>**; your persona is in `workshop/personas/<id>.md`. The
   > workshop is planning cycle CCC: rounds AAA-BBB, five phases of four rounds
   > (literature, research, research, research, backlog; see
   > `workshop/CYCLE.md`). Read your notebook, `workshop/STATE.md`,
   > `workshop/STEERING.md`, the last cycle's `retrospective.md` (or
   > `workshop/retrospectives/` for the latest), and grep `research/` for what
   > you need. Write `workshop/cycles/CCC/plan/<id>.md` in the shape of
   > `workshop/templates/plan.md`, at most 40 lines. It is a loose plan, not a
   > commitment: what you would like to find out this cycle and why; for
   > **phase 1**, the papers or kinds of paper you would look for and which
   > open claim each would check; for **phases 2-4**, the two or three
   > questions you would work on and the hypothesis each bears on; for
   > **phase 5**, what you expect to be left open and want closed; and what
   > result by the end of the cycle you would call progress. Do not run
   > anything and edit no other file. Reply with your main question in one
   > line.

3. The chair writes `workshop/cycles/CCC/plan.md` (at most 80 lines):
   * **Questions of the cycle.** Every question the plans propose to work on
     is tied to an H-entry. If one has none, the chair writes it now in
     `research/HYPOTHESES.md` (next free number, falsifiable statement, status
     `OPEN` or `SUPPORTED` with pointers, `· *workshop cycle CCC planning*` in the
     date line). A thread in `STATE.md` must name its H-entry.
   * **Phase 1 reading list**: papers and searches, each with the claim or
     inconclusive result it should check, and who reads it.
   * **Loose agenda for phases 2-4**: ranked threads, owners, first questions.
     Not binding: phase conferences change it.
   * **Expected backlog**: what the plans expect to need closing in phase 5.
   * **What would count as progress** by the end of the cycle (two to four
     concrete outcomes, e.g. "F-entry for X written, or the blocking check
     named").
   * Where the personas disagree.
4. Rewrite `STATE.md` for the cycle (header `cycle`, `phase: 1`,
   `focus: literature`, `next_round_kind: ordinary`), with the threads of the
   plan. Add a DIGEST entry `## Cycle CCC planning -- <date>` (at most 15
   lines), ending with one bold line inviting the human to edit `plan.md` or
   `STEERING.md` if they want the cycle to go differently. Commit and push as
   in `ROUND.md` step 5.

## Phase focus: what changes in a round

### Literature (phase 1)

* **Call.** Every assignment is literature work on an open thread: find
  papers (arXiv search and reading, `research/literature/`), and check
  against them an open hypothesis, an inconclusive or "no power" experiment,
  a claim marked UNVERIFIED, or a term the project coined. The scholar works
  every round of the phase; the other two slots go to personas whose threads
  have the most unchecked claims. No breadth slot and no suggested-question
  slot: the phase already ranges over all threads.
* **Work.** Computation only to check a published statement against the
  repo (reproduce a paper's example, test a stated theorem on small cases).
  A submission of kind `literature` lists what was searched (queries,
  identifiers), what was read and in which version, and for each claim
  checked: *supported by*, *contradicted by*, *not addressed*, with the exact
  theorem number matched in the source. A new paper summary goes in the
  submission for the chair to promote to `research/literature/`.
* **Review.** The referee checks the cited statements against the source
  (the LaTeX, if fetched) and its hypotheses against the case here.
* **Conference.** Positions answer: what did the literature change for my
  thread; what is now known to be known; what should phase 2 not do
  because it is published. The ledger also covers literature results.

### Research (phases 2-4)

Rounds run as `ROUND.md` describes, with the breadth slot and the
suggested-question slot. The cycle plan's loose agenda is the starting
agenda of phase 2; each conference re-proposes it.

### Backlog (phase 5)

The research backlog, not the code backlog.

* **First round of the phase.** Before the call, the chair writes the
  backlog list in `workshop/cycles/CCC/backlog.md`: every open item from the
  `STATE.md` threads, *Awaiting revision*, *Requests between personas*,
  "open:" points of the cycle's E-entries, claims marked UNVERIFIED,
  overnight proposals never decided, and pending consequences for library
  code or docstrings. One line each: item, source identifier, what would
  close it. Rank by what it would settle.
* **Call.** Each assignment takes one or more backlog items and must end
  with each item **settled** (answered, with evidence), **retired** (not
  worth doing, with the reason), or **converted** (it is really a new
  question: written as a proposed H-entry for the next cycle). A backlog
  assignment that only narrows an item is a failure to report as such.
* **Consequences.** A pending consequence for library code or a docstring
  is done in this phase or explicitly declined with a reason; it does not
  carry over silently.
* **Conference (the cycle's last round).** The ledger covers the whole
  cycle, not just the phase: every hypothesis the cycle worked on gets a
  decision, and the chair marks each backlog item settled, retired,
  converted, or carried over (with the reason).

## Phase summary (at each conference)

The conference chair writes `workshop/cycles/CCC/phase-P.md` (at most 40
lines): what the phase did, by thread; what entered `research/` (with
identifiers); what the ledger decided. This is where history goes; it is
removed from `STATE.md` at the same time.

## Cycle retrospective (after the last round of a cycle)

Done by the dispatcher session, not a round session, when `STATE.md` says
`next_round_kind: retrospective`. It writes
`workshop/cycles/CCC/retrospective.md` and a copy in
`workshop/retrospectives/` (numbered after the last one), looking at the
cycle as a research process:

* the cycle plan's progress outcomes: met, not met, why;
* each phase: did its focus do what it was for (did phase 1 find things that
  changed phase 2; did phase 5 close what it listed);
* the numbers (rounds, reviews, entries, statuses changed, cost);
* **the format itself**: number of phases, their order and focus, whether to
  add another focus (examples: a *consolidation* phase that only writes
  F-entries; a *proof* phase; a *replication* phase where results are
  rebuilt independently; a *problem-session* phase on one question), the
  phase length, the conference cadence.

It then makes the changes it decides on (this file, `ROUND.md`, `config.yaml`,
templates, `README.md`), says what it changed in the retrospective, sets
`next_round_kind: planning` for the next cycle, and pushes.
