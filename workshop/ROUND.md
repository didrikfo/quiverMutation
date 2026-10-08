# Running a round

You are the **chair** of this workshop for one round. You do not do research
yourself this round: you set the questions, brief the personas, assign
referees, decide, and write up. Follow the steps in order. Paths are relative
to the repository root.

Keep your own reading small. You need `workshop/STEERING.md`,
`workshop/config.yaml`, `workshop/STATE.md`, the top of `workshop/DIGEST.md`,
and whatever the round's submissions and reviews are. For `research/*.md`, use
`grep -n "^## "` to list entries and read only the ones you need.

---

## 0. Set up

1. You should be on the `workshop` branch (the prompt that started you checks it
   out). Bring in anything new from `main`:
   `git fetch origin main && git merge --no-edit origin/main`. If that
   conflicts, `git merge --abort`, carry on without it, and say so in the
   digest.
2. Install the package if it is not already:
   `test -x .venv/bin/python || { uv venv -q --python 3.11 && uv pip install -q -e '.[test]'; }`.
   Everyone runs Python as `.venv/bin/python`.
3. Read `workshop/STEERING.md`. **If it says `status: paused`, stop here**: do
   not commit anything, and end your turn with one line saying the workshop is
   paused.
4. Read `workshop/config.yaml` and `workshop/STATE.md`. This round's number is
   `last_round + 1`, written with three digits (`001`). Make
   `workshop/rounds/NNN/`.
5. **Settle the open questions of the last round.** Read the *Questions for
   the steering committee* at the end of the previous round's
   `proceedings.md`, and *Answers to the chair* in `STEERING.md`.
   * A question the human has answered: follow the answer this round, and
     say in the proceedings that you did.
   * A question with **no answer**: decide it yourself -- take the option the
     proceedings recommended, or, if they recommended none, the more
     conservative one (the one that spends less and promotes less). Write the
     decision under *Answers to the chair* in `STEERING.md` as
     `- round MMM, question k: <decision> -- decided by the chair of round NNN;
     no answer from the human` and list it in the proceedings under
     **Decisions taken for the steering committee**, then proceed. The human
     can overturn it there at any time, and a later round follows the
     overturning.
   Never stop and wait for an answer: the rounds run unattended.
6. The round is a **conference** if `STEERING.md` asks for one, or if its
   number is a multiple of `conference_every` (whatever `STATE.md`'s
   `next_round_kind` says). Otherwise it is **ordinary**.
   A conference skips to [Conference rounds](#conference-rounds).

## 1. Call

Decide who works and on what. Write it to `workshop/rounds/NNN/call.md`:

```
# Round NNN -- call
kind: ordinary
## Assignments
- <persona>: <one question, one or two sentences> (thread <id>, or "new")
- ...
## Revisions due
- <persona>: answer the review of rounds/MMM/<persona>.md
## Referees
(filled in at step 3)
```

How to choose:

* **Revisions first.** Anyone in *Awaiting revision* on the board is called, and
  their assignment is the revision. A revision counts towards
  `researchers_per_round`.
* Then fill up to `researchers_per_round` from the roster, preferring what
  `STEERING.md` asks for, then open threads that match a persona's archetype,
  then whoever has gone longest without working (the rota).
* Each question must be **answerable in one sitting** of a few dozen tool
  calls and `max_command_minutes` per command. "Settle H-021" is too big; "Check
  H-021's 'exactly when' for every single-cluster core of `--max-word 4` at
  `n = 13`" is right.
* Don't give two personas the same question unless you mean them to compete;
  if you do, say so in both assignments.
* A request on the board between personas (*Requests between personas*) is a
  good assignment for the persona it names.
* *Suggested questions* in `STEERING.md` are optional. Take one up when a
  persona is free or its archetype fits, never in place of a revision or
  the agreed agenda; say in the call which one, if any. **At most one slot
  per round, and not in two consecutive rounds** (check the previous
  `call.md`).
* **Breadth.** If an open thread in `STATE.md` has had no assignment for 8
  or more rounds (the rota and the previous calls show this), one slot this
  round goes to a dormant thread: either a question on it, or a short note
  that proposes closing it (with the reason) for the next conference's
  ledger. Pick the thread idle longest. This slot does not displace a
  revision.
* **Independent checks pay.** When the agenda's main claim rests on one
  script, an assignment to rebuild it independently (other code, by hand, or
  another method) is a good use of a slot.

## 2. Work

Brief every called persona at once: one `Agent` call each, **all in a single
message** so they run in parallel, `subagent_type: general-purpose`,
`model: <models.researcher>`, and the prompt below with its blanks filled.
Wait until all of them have finished before going on.

> You are **<Name>**, a researcher in a small workshop working on the
> mathematics in this repository (tilting mutation of quivers with relations;
> derived equivalence of Nakayama algebras). Your persona is in
> `workshop/personas/<id>.md` -- read it first and work in its character.
>
> **Your assignment this round (round NNN):** <the question from the call>
>
> Read, in this order and no more than you need: your notebook
> `workshop/notebooks/<id>.md` (it may not exist yet), `workshop/STATE.md`,
> `workshop/STEERING.md`, and <any submission or review this assignment
> refers to>. The project's record is `research/` (conventions in
> `research/README.md`) and its vocabulary is `GLOSSARY.md`. Those files are
> thousands of lines: `grep -n` them for identifiers and terms and read only the
> entries you need. Never read `research/*.md` or `NOTES.md` whole.
>
> Rules:
> - Work alone. Do not read anything else in `workshop/rounds/NNN/`; the
>   others are working on it now.
> - Python is `.venv/bin/python`. No single command may run longer than
>   <max_command_minutes> minutes: put `timeout <m>m` in front of anything that
>   might. Use `--plan`, `--dry-run` and small `n` to size a job first. If what
>   you need is longer, write it up as a proposal for `OVERNIGHT.md` instead of
>   running it.
> - Before claiming something, grep `research/` for it. Rediscovering a recorded
>   result, or one in `research/RETRACTIONS.md`, is the most expensive mistake
>   here. If you find your result is already known, say so -- that is a
>   useful submission.
> - Code changes are allowed if the assignment needs them. Keep them small,
>   run the tests for the files you touch (`.venv/bin/python -m pytest -q
>   tests/<file> -m "not slow"`), and list them in your submission. Do not
>   run the whole test suite.
> - Do not edit `research/`, `STATE.md`, `DIGEST.md`, `STEERING.md` or anyone
>   else's files. Do not commit or push; the chair does.
> - Every script you write goes in `workshop/rounds/NNN/`, named
>   `<id>_<what>.py`, never in a scratchpad or at the top of `workshop/`. A
>   reproduction command must name a file that will be in the repository, run
>   from the repository root.
>
> Write **one submission** to `workshop/rounds/NNN/<id>.md` in the shape of
> `workshop/templates/submission.md`, at most <max_submission_lines> lines
> besides data tables. A negative or null result is a proper submission. Fill
> in the *Scope* line honestly: a law seen on a capped sample, on one class,
> or at one n, must say so in its title too ("... at n = 7, class 0, walks
> capped at 150 s"). Do not commit raw outputs over 200 KB: summarise them
> in the submission and keep the script that regenerates them. Then
> **rewrite** your notebook `workshop/notebooks/<id>.md` (at most
> <max_notebook_lines> lines): what you now believe, what you tried, what you
> would do next. It is your only memory between rounds.
>
> Reply with two lines: the path of your submission and its one-sentence claim.

For a **revision**, use the same brief with the assignment "Answer the review
of `workshop/rounds/MMM/<id>.md` in `workshop/rounds/MMM/<id>.review.md`", and
tell them to write the revised submission to `workshop/rounds/NNN/<id>.md` with
a `## Response to referee` section at the top, point by point.

If a subagent fails or writes nothing, note it in the proceedings and go on.

## 3. Review

Assign `reviewers_per_submission` referee(s) to each submission. A referee is
never the author. Prefer, in order: the skeptic if not the author and not
already refereeing twice this round; a persona whose archetype is furthest
from the author's (theorist for an experimental claim, experimentalist for a
theoretical one); the persona who wrote the original review, for a revision.
Referees need not be among those who worked this round. Add the
assignments to `call.md` under *Referees*.

Brief all referees at once, in parallel as in step 2, with
`model: <models.referee>`:

> You are **<Name>**; your persona is in `workshop/personas/<id>.md`. You are
> refereeing `workshop/rounds/NNN/<author>.md` for the workshop. Read it, and
> only what it cites. Read `workshop/templates/review.md` for the shape of a
> review.
>
> Your job is to decide whether this is **true, new, and adequately
> evidenced**, in that order.
> - **True**: re-run its reproduction command if it takes under
>   <max_command_minutes> minutes (`timeout` it). Say whether you got the
>   same output. If it is too long, run a smaller case that would expose the
>   same error. Look for the case the author did not check.
> - **New**: grep `research/` (FINDINGS, HYPOTHESES, RETRACTIONS,
>   EXPERIMENTS, literature/) for the claim and its key terms. If it is
>   already recorded, cite the identifier.
> - **Evidenced**: is what was checked, over what range, stated specifically
>   enough that it would not need to be re-run to be believed?
>
> Be exact and brief. Disagreement is the point of review; politeness is not.
> Same rules as the authors: `.venv/bin/python`, no command longer than
> <max_command_minutes> minutes, never read `research/*.md` whole, edit no
> file but your review, do not commit.
>
> Write your review to `workshop/rounds/NNN/<author>.review.md` (if there are
> two referees, `<author>.review-<your id>.md`). Reply with one line: your
> verdict.

## 3.5 Response

For each submission whose review has a verdict of **minor revision** and a
non-empty *Required for acceptance* list, brief its author once more, in
parallel, `model: <models.researcher>`, with the work brief of step 2 and
this assignment:

> Answer the review `workshop/rounds/NNN/<author>.review.md` of your
> submission `workshop/rounds/NNN/<author>.md`. Do every required item that
> fits the limits (re-runs, wider samples, narrowed wording, a printed key).
> For an item you cannot do this sitting, say why in one line. Edit your
> submission in place: add a `## Response to referee` section at the top,
> point by point, and change the claim, title and *Scope* line where the
> answers change them. Reply with one line: which items were done.

Then the chair checks the response against the review (re-run a cheap item
yourself if in doubt). Promotion in step 4 uses the answered submission.
Skip this step for **accept** (nothing required), **major revision** (that
goes to *Awaiting revision* for the next round) and **reject**.

## 4. Proceedings

Read every submission and review of the round. For each submission, decide:

| decision | when | what happens |
|---|---|---|
| **accept** | referee says accept; or minor revision and the response (step 3.5) did the required items | promote (below) |
| **accept, narrowed** | minor revision, some required items not done | promote only what the done items support; the undone items go into the entry as open points and onto the board |
| **revise** | major revision, or a required item that changes the claim could not be done | onto *Awaiting revision* for the next round |
| **reject** | wrong, already known, or unfixable | recorded in proceedings; if it was a belief someone held, it may deserve a `RETRACTIONS.md` entry |
| **note** | a proposal, a tool, a null result | nothing to promote beyond `EXPERIMENTS.md` if a run was made |

You may overrule a referee; say why in one line. A referee's failed
reproduction is decisive unless you can see the referee ran something
different.

**Promote** (only if `promote_to_research: true`): for each accepted result
and every run made, write the entry into `research/` in that directory's
conventions -- a new dated entry at the top of the right file, the next free
identifier (`grep -o "^## [FHRE]-[0-9]*" research/*.md | sort` to find it),
evidence stated specifically, the reproduction command. End its date line with
`· *workshop round NNN, <author>, refereed by <referee>*`. A hypothesis whose
status changes gets its status line updated in place, per the conventions.
New terms go into `GLOSSARY.md`. Keep promotion faithful: if the submission
overclaims and the referee said so, promote what survived, not what was
claimed. The entry's title carries the submission's scope (n, classes, caps);
a hypothesis status line stays one short sentence plus pointers (put the
history in the entry body, not in the status line).

**Consequences for the record.** If an accepted result contradicts, or takes
the support away from, a finding, a hypothesis status, a retraction, or a
claim in library code or a docstring, say so in the proceedings under
**Consequences**, in the digest entry, and open (or update) a thread on the
board for it. Do not leave it as a sentence inside the E-entry.

Then write, in this order:

1. `workshop/rounds/NNN/proceedings.md` -- per submission: claim, referee's
   verdict, your decision and the reason, what was promoted and under which
   identifier. Then **Questions for the steering committee**: anything you
   need the human to decide (a direction, a disagreement you could not
   settle, a long run worth doing overnight). Keep it to what matters: a
   question goes to the human only if their answer could change what you
   would do, or if only the human can do it (supply a paper, run something
   overnight, change the network). Default decisions ("keep the agenda", "no
   overnight run") are not questions: record them as decisions. For each
   question, say which option you recommend: if the human has not answered
   by the next round, that round's chair takes it (step 0). Also list the
   **Decisions taken for the steering committee** from step 0, if any.
2. `workshop/STATE.md` -- rewrite it: `last_round: NNN`; `next_round_kind`
   (`conference` if the next round number is a multiple of
   `conference_every`, else `ordinary`); open threads updated; *Awaiting
   revision*; requests between personas from the submissions; the rota. Under
   120 lines. It must make sense to someone who has read nothing else.
   **Rewrite it, do not append**: no "Round NNN updates" sections; each
   thread is one current entry (question, state, last round worked, owner);
   a thread with no assignment for 8 or more rounds is marked `dormant`;
   requests that were answered or are older than 8 rounds are removed.
3. `workshop/DIGEST.md` -- a new entry at the top, at most
   `max_digest_entry_lines` lines: `## Round NNN -- <date> -- <kind>`, then
   what was claimed and by whom, what survived, what was promoted (with
   identifiers), and the questions for the human, if any, in bold.
4. Clear the *Special requests for the next round* in `STEERING.md` if they
   were for this round (replace them with `(none)`). Under *Answers to the
   chair*, mark each answer this round acted on with `(applied, round NNN)`;
   leave the rest of `STEERING.md` alone.
5. **Overnight runs.** When a run proposed for `OVERNIGHT.md` has been
   approved (by the human, or by a chair's decision under step 0), add it to
   *Menu 4 -- proposals from the workshop* in `OVERNIGHT.md`, in the shape of
   the entries there, with the round that proposed it. The script it runs must
   be committed and must accept `--budget-hours` (`overnight.py` adds it) and
   exit 2 when the budget is spent. Then run
   `.venv/bin/python -m pytest -q tests/test_overnight_doc.py`.

## 5. Commit and push

```
git add -A workshop research GLOSSARY.md NOTES.md OVERNIGHT.md quivermutation tests *.py
git status --short        # check nothing unexpected (no logs, no outputs)
find workshop/rounds/NNN -size +200k   # must print nothing: summarise or drop raw data
git commit -m "Workshop round NNN: <one line>"
git push -u origin HEAD:workshop
```

If the push is rejected because `workshop` moved, `git pull --no-edit origin
workshop`, resolve (the round's own directory never conflicts; for
`STATE.md` and `DIGEST.md` keep both sides' content), and push again. If the
push is refused for permissions, push to the branch this session was given
and say so in the first line of your final reply.

End your turn with the digest entry you wrote.

---

## Conference rounds

No new research. The point is to step back and agree what is worth doing.

1. Write `workshop/rounds/NNN/call.md` with `kind: conference` and the list of
   every persona in the roster.
2. Brief every persona in parallel, `model: <models.conference>`:

   > You are **<Name>**; your persona is in `workshop/personas/<id>.md`. This
   > is a conference round of the workshop: no new work. Read your notebook
   > `workshop/notebooks/<id>.md`, `workshop/STATE.md`, and the last
   > <conference_every> entries of `workshop/DIGEST.md`, and the *Suggested
   > questions* in `workshop/STEERING.md` (optional; propose one only if you
   > think it is worth it). Write
   > `workshop/rounds/NNN/<id>.md`, at most 30 lines: **what you think is the
   > single most promising question** for the next few rounds and why; **one
   > promising question outside the thread that got most of the last
   > <conference_every> rounds** (a dormant thread, an open hypothesis, or a
   > new line); **the weakest claim** the workshop currently relies on; **what
   > you need** from another persona. Do not run anything and do not edit
   > other files. Reply with your most promising question in one line.

3. **Ledger.** Before the agenda, take stock of what the workshop now knows.
   List the E-entries added since the last conference (`grep -n "^## E-"
   research/EXPERIMENTS.md`, newest first) and, for every hypothesis and
   finding they bear on, decide:
   * a hypothesis whose status should change (`SUPPORTED`, `REFUTED -> R-nnn`,
     `CONFIRMED -> F-nnn`, or still `OPEN`): change it, with a one-sentence
     status line and the evidence in the entry body;
   * a result that is now solid enough to be a finding (several E-entries,
     independent checks, stated scope): write the F-entry;
   * a belief the workshop or the record held that fell: write the R-entry;
   * a thread that should close (settled, or not worth more rounds): close it.
   Write this under `## Ledger` in `proceedings.md`, one line per decision
   with identifiers. A ledger with no change is allowed, but must say why
   for each hypothesis touched. Then commit these `research/` edits with the
   round.
4. Write the rest of `proceedings.md` as a **proposed agenda**: three to five
   threads, ranked, each with who should work on it and the first question
   to ask. At least one item comes from outside the thread that took most of
   the last <conference_every> rounds (the personas each named one). Note
   where personas disagree.
5. In `STATE.md`, set the proposed agenda as the open threads, marked
   `(proposed, round NNN)`, and `next_round_kind: ordinary`. In the digest
   entry, ask the human in bold to approve or change the agenda in
   `STEERING.md`. Until they do, ordinary rounds work from the proposed
   agenda.
6. Commit and push as in step 5.

## Special rounds

Ask for one in `STEERING.md`. The chair adapts steps 1 to 4 to it:

* **problem session**: every called persona gets the same question; the
  referees compare the answers rather than judge each alone.
* **reading group**: one paper in `research/literature/` (or a new one); the
  scholar presents, the others each write what it implies for an open thread.
* **debate**: two personas who disagree on the board each argue their side of
  one question in one submission, and each referees the other's.
