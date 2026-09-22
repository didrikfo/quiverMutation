# Research harness

A written record of what we have tried, what we believe, what we got wrong, and
why — so that none of it has to be rediscovered.

This exists because rediscovery is the main cost in this project. Several
results found in the course of the work were things already known and forgotten,
and more than one line of work was pursued on an assumption that had already been
shown false. The point of these files is to make that cheap to check.

Terms used across these files are defined in [`../GLOSSARY.md`](../GLOSSARY.md).
When an entry coins a word, add it there.

## The files

| file | holds |
|---|---|
| [`FINDINGS.md`](FINDINGS.md) | things established, with the evidence for them |
| [`HYPOTHESES.md`](HYPOTHESES.md) | things suspected, open or resolved, with what would settle them |
| [`RETRACTIONS.md`](RETRACTIONS.md) | things believed and then found false, kept with what corrected them |
| [`EXPERIMENTS.md`](EXPERIMENTS.md) | runs made, with parameters and outcome, so they are not repeated |
| [`literature/`](literature/) | one summary per paper: results, lemmas, caveats, with the arXiv reference |

## Conventions

**Every entry is dated and has an identifier.** `F-003`, `H-007`, `R-002`,
`E-011`. Identifiers are permanent; cross-reference by them.

**Nothing is deleted.** A finding that turns out to be wrong is *not* removed or
edited into correctness. It gets a `**RETRACTED yyyy-mm-dd**` line at the top
saying what superseded it, stays where it is, and the correction goes in
`RETRACTIONS.md` with its own identifier. A wrong belief that has been recorded
is worth more than one that has been tidied away, because it stops the same
reasoning being repeated.

**Status is explicit.** Every hypothesis carries one of

- `OPEN` — nobody has tested it
- `SUPPORTED` — evidence for it, not settled
- `CONFIRMED → F-nnn` — established; the finding is where the evidence lives
- `REFUTED → R-nnn` — shown false
- `PARKED` — deliberately not being pursued, with the reason

**Evidence is specific.** "Verified" on its own is not a record. Say what was
checked, over what range, and how many cases: *"checked at every window position
of every LNA of lengths 7 to 10, 22 confirmations, no failures"*. A claim whose
evidence is not written down cannot be trusted later, and will be re-run.

**Reproduction is a command.** Where a result came from code, give the command
or the function that produces it, so it can be re-run rather than re-derived.

## Adding an entry

Copy the shape of the entries already there. New entries go at the **top** of
their file, so the most recent work is first. Keep each to what a reader needs:
the claim, the evidence, the date, and anything that would mislead someone who
only read the claim.

## When to write

- Before a long run: write the hypothesis it is meant to test.
- After a run, whatever the outcome: record it in `EXPERIMENTS.md`. A run that
  found nothing is worth recording precisely so it is not repeated.
- The moment something believed turns out to be false: `RETRACTIONS.md`, before
  fixing the code, while the reasoning is still fresh.
- When reading a paper: a summary in `literature/`, not the paper itself.
