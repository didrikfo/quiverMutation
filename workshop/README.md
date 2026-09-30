# The workshop

A small multi-agent research group that works on this repo unattended, in Claude
Code on the web, and reports back. Several researcher personas work
independently, then read and referee each other's work through fixed channels,
and a chair writes up what happened. You come back, read one file, and steer.

Nothing here is code. The harness is a procedure ([`ROUND.md`](ROUND.md)) that a
cloud session follows, plus the files it reads and writes. Every part of it can
be changed by editing Markdown.

---

## How a round works

One **round** is one cloud session. It is started by a Routine on a schedule, or
by hand. The session plays the **chair** and runs four phases:

```
  STEERING.md (you)        STATE.md (chair's board)
          \                     /
           v                   v
   1. CALL        chair picks 2-4 personas and gives each a question
           |
   2. WORK        each persona is a separate subagent, in parallel; it sees
           |      the call, the board and its own notebook, never the
           |      others' work this round; it writes one submission
           v
   3. REVIEW      each submission gets a referee from a different persona,
           |      who re-runs the reproduction command and gives a verdict
           v
   4. PROCEEDINGS chair decides, updates the board, promotes accepted results
                  into research/, writes DIGEST.md, commits and pushes
```

A submission sent back for revision comes up again at the next round: its
author answers the review (a **rebuttal**), and the referee looks again. Every
fourth round is a **conference** instead: no new work, each persona gives a
short position ("what I now think is most promising, and why"), and the chair
proposes an agenda for you to approve in `STEERING.md`.

The personas are in [`personas/`](personas/). Each keeps a short
**notebook** in [`notebooks/`](notebooks/) that it rewrites every round it works;
that is its memory between rounds, and nothing else carries over.

## Where things are

| file | owner | holds |
|---|---|---|
| [`STEERING.md`](STEERING.md) | **you** | direction, priorities, things to avoid, the pause switch |
| [`DIGEST.md`](DIGEST.md) | chair | one short entry per round, newest first -- **read this when you come back** |
| [`STATE.md`](STATE.md) | chair | the board: round counter, open threads, pending revisions, the rota |
| [`config.yaml`](config.yaml) | you | roster, how many work per round, models, time limits |
| [`ROUND.md`](ROUND.md) | you | the procedure a round follows |
| `personas/*.md` | you | who the researchers are |
| `notebooks/*.md` | each persona | its own working memory |
| `rounds/NNN/` | the round | `call.md`, one `<persona>.md` per submission, `<persona>.review.md` per review, `proceedings.md` |
| [`templates/`](templates/) | you | the shape of a submission and a review |

Results that pass review are written into [`../research/`](../research/) by the
chair, in that directory's conventions, marked as coming from the workshop.
Everything happens on the **`workshop` branch**. Nothing reaches `main` until
you open a pull request from `workshop` and merge it, so that PR is where you
check what the group claims before it becomes the record.

## Running it

**By hand, once.** Start a Claude Code session on this repo and send the prompt
in [Round prompt](#round-prompt) below. Watch the first round or two before
putting it on a schedule.

**On a schedule.** A Routine that starts a fresh session on each firing does
*not* get the repository attached, and could not push (found on the first
attempt, 2026-09-29). So the schedule is two-stage: a Routine fires a short
turn into a **dispatcher** session every 2 hours, and the dispatcher starts
the round with `create_session`, giving it the repository as its source and
`workshop` as its outcome branch, which is what lets it push:

```
create_session(
  source_url      = "https://github.com/didrikfo/quiverMutation",
  source_revision = "workshop",
  outcome_branch  = "workshop",
  model           = "claude-sonnet-5-5",
  title           = "Workshop round NNN",
  prompt          = <the round prompt below>)
```

The dispatcher also checks that the previous round pushed, and does not start
a round while one is still running. A round is one session, so the schedule is
the token budget. Pause it by disabling the Routine, or more cheaply by
setting `status: paused` in `STEERING.md` (the round then reads one file and
stops).

### Round prompt

```
Run one round of the research workshop in this repository. Make sure you are
on the workshop branch (git checkout -B workshop origin/workshop). Then read
workshop/ROUND.md and follow it exactly, as the chair. At the end push to the
workshop branch (git push -u origin HEAD:workshop). If a push is refused, end
your turn with the exact error. Do not open a pull request. End with the
DIGEST.md entry you wrote.
```

## Steering

Edit `STEERING.md` on the `workshop` branch -- in the GitHub web editor, or by
asking any Claude session to do it. The chair reads it first, every round,
and it outranks everything else it reads. Useful things to put there:

* a question you want worked on, or one to drop;
* "persona X sits out" or "only the skeptic and the experimentalist this week";
* "the next round is a conference";
* `status: paused`.

You can also answer the chair directly: its proceedings end with **questions
for the steering committee**, and whatever you write under them in
`STEERING.md` is read as your answer. The chair says which option it recommends; if you
have not answered by the next round, that round's chair decides, records the
decision under *Answers to the chair* marked as its own, and carries on. You
can overturn it there whenever you like.

## Keeping the cost down

The knobs, most effective first:

1. **How often rounds run** (the Routine's schedule).
2. **`researchers_per_round`** in `config.yaml`. Two is a real round; four is
   the most that fits comfortably.
3. **Models.** Researchers and referees default to a cheaper model than the
   chair; see `config.yaml`.
4. **Reading discipline.** `research/*.md` is ~9000 lines. Personas are told to
   grep it for identifiers, never to read it whole, and to start from
   `STATE.md` and their notebook. If a round's cost jumps, this is the first
   thing to check in its transcript.
5. **Compute limits.** One command may run at most `max_command_minutes`.
   Long censuses belong in `OVERNIGHT.md` on your own machine; a persona can
   *propose* one there (as a submission of kind `proposal`) rather than run it.

## Changing the shape

This is a first version and meant to be changed. Things that are one edit
away:

* **Add or change a persona**: a file in `personas/` and a line in
  `config.yaml`.
* **Two referees instead of one**: `reviewers_per_submission` in `config.yaml`.
* **A different cadence of conferences**: `conference_every`.
* **A different venue altogether** (a reading group on one paper, a
  "problem session" where everyone attacks the same question, a debate between
  two personas): add it under *Special rounds* in `ROUND.md` and ask for it in
  `STEERING.md`.

Things that would need more than an edit, for a later version: personas in
separate long-lived cloud sessions rather than subagents of one session
(more independence and longer work, at the price of git coordination and
several times the tokens), and rounds that pick up long computations started
by an earlier round.
