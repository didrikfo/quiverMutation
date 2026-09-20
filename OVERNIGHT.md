# Overnight runs

How to leave this machine working for a night, and what to ask it. Every command
below is complete: copy one, run it, go to bed.

Nothing here needs the code changed. Pick a length, pick a width, run it. The
ledgers make every run resumable and every rerun cheap, so a night that is cut
off is a night's data and not a wasted night.

---

## The command

One shape, always:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'JOB' --run 'JOB'"
```

`overnight.py` keeps the Windows host awake, restarts a job that dies, gives
every job the same `--budget-hours`, and writes each one's output to
`logs/<name>-<stamp>.log`. Each `--run` takes one job as a quoted string; there
can be as many as there are cores to spare.

**The venv is a WSL venv.** Windows Python has no `networkx` and will not run
this. `wsl -e bash -lc "..."` is not optional.

**Leave the laptop lid open.** The keep-awake stops the machine idling to sleep;
it has no say over the lid.

Check before committing a night to it — `--dry-run` prints the commands and
starts nothing:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --dry-run --hours 9 --run 'batch.py cores 16 --max-word 4 --jobs 7'"
```

**Sixteen cores.** Keep the `--jobs` across all the night's jobs at **14 or
fewer**, or they fight each other and the machine gets slower, not faster.

**Size a job before you pick it.** `--plan` prints the ledger it would write,
how many units it is, and how many of them are already done. It runs nothing:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python batch.py cores 17 --max-word 4 --plan"
```

---

## Tonight, if you have no particular question

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 14 --max-word 4 --jobs 2' --run 'batch.py cores 16 --max-word 4 --jobs 3' --run 'batch.py cores 17 --max-word 4 --core-limit 250 --jobs 3' --run 'batch.py cores 15 --max-word 2 --pair-word 2 --gaps 1,2,3,4 --jobs 3' --run 'batch.py sample 13 --count 50000 --jobs 2' --run 'batch.py sample 17 --count 20000 --jobs 1'"
```

Three new lengths of the core census — `n = 17` deliberately cut to its first
250 cores rather than cut by the clock — the two-cluster shapes that have never
been run at all, and two sampling lengths.

**Give every job more work than the night can finish.** A job that runs out of
units exits and leaves its cores idle until morning; a job that runs out of
budget exits 2 with everything it did on disk and continues next time. So set
`--count` high, take the longer length, and let the budget do the cutting —
except for a census, where the cut should be `--core-limit` and not the clock.

---

## Menu 1 — `cores`: slide a configuration along a line

Put one overlapping configuration in an otherwise empty line, at every offset,
and ask at each whether the moves carry it to an almost separate LNA. Three
verdicts and never two: `inside` (a certificate), `outside` (the forward orbit
**closed** without one), `undecided` (the walk hit its cap — the budget was
measured, not the moves).

The summary prints each core's slide as `i`/`o`/`?` left to right, with `head`
and `tail` — how many offsets at the source end and at the sink end are inside.
**Those two columns are the thing to compare between lengths.**

### What to vary

| flag | what it changes | values worth running |
|---|---|---|
| `length` | the line | `11 12 13 14 15 16 17 18` |
| `--max-word` | how many vertices one core spans | `2` `3` `4` `5` |
| `--max-arrows` | the longest relation in a core | `6` (default), `8`, `10` |
| `--pair-word` | each half of a two-cluster core | `2` `3` |
| `--gaps` | free arrows between the two halves | `1,2,3` (default), `1,2,3,4,5`, empty for none |
| `--cores` | run only these words | `45,504` — `245,2045,2245,2555,2556,3344` |
| `--core-limit` | only the first N words of the catalogue | `40` `60` `120` `400` |
| `--orbit-limit` | rows before a placement is undecided | `20000` (default), `60000`, `200000` |
| `--join-limit` | rows per side for the two-ended join | `6000` (default), `40000` |

`--cores` and `--core-limit` are **filters**: they narrow what a run does
without changing what an answer means, so they share the ledger with the full
census and two nights can split one census between them. Everything else is in
the ledger's name, because it changes what a verdict says.

**The same `--core-limit` asks about the same cores at every length** — the
catalogue is built the same way regardless of `n`. That is what makes two
part-finished censuses comparable, and it is what the last run got wrong.

### What a census costs, on one core

| command | placements | one core | 6 workers |
|---|---|---|---|
| `cores 13 --max-word 4` | 1705 | 1.1 h | ~11 min |
| `cores 15 --max-word 4` | 2591 | 3.8 h | ~38 min |
| `cores 16 --max-word 4` | 3034 | 8.7 h | ~1.5 h |
| `cores 17 --max-word 4` | 3477 | 37.5 h | ~6 h |

The catalogue grows slowly and the cost per placement does not: 13 to 16 is
about a doubling a length, and 16 to 17 is more than four times. `n = 18` at
this width has not been timed — expect a long night at best, so start it with a
`--core-limit` and let it resume. `--max-word 3` is about a third of
`--max-word 4`; `--max-word 5` is several times more.

**Always `--plan` first at a new length.** The number that matters is
`left`, and the cost per unit is what the table above is for.

### Nights worth running

**A run of lengths, for H-020.** The head and the tail should not move with `n`.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 11 --max-word 4 --jobs 2' --run 'batch.py cores 12 --max-word 4 --jobs 2' --run 'batch.py cores 14 --max-word 4 --jobs 3' --run 'batch.py cores 16 --max-word 4 --jobs 5'"
```

**The two-cluster shapes, for F-040 and H-018.** These have never been run: the
catalogue puts them after every single core, and no run has ever reached them.
A length of 15 or more is needed for a barricade to fit at all.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 15 --max-word 2 --pair-word 2 --gaps 1,2,3,4 --jobs 5' --run 'batch.py cores 17 --max-word 2 --pair-word 2 --gaps 1,2,3,4 --jobs 5' --run 'batch.py cores 16 --max-word 2 --pair-word 3 --gaps 1,2,3 --core-limit 900 --jobs 4'"
```

`--pair-word 3` without the `--core-limit` is 21072 placements, which is several
nights; 900 makes it one, and the rest resume into the same ledger later.

**The six undecideds, sharpened.** They sit exactly where the verdict changes,
which is where a wrong head or tail would come from. This writes its own ledger,
because a different limit is a different question.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 6 --run 'batch.py cores 15 --max-word 4 --orbit-limit 60000 --join-limit 40000 --cores 245,2045,2245,2555,2556,3344 --jobs 6'"
```

**Wider cores, one length.** Everything so far is words of at most four
vertices. Five is the first width nothing has looked at.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 14 --max-word 5 --gaps , --jobs 7' --run 'batch.py cores 16 --max-word 5 --gaps , --core-limit 400 --jobs 7'"
```

(`--gaps ,` leaves the pairs out, so the night is single cores only.)

**Longer relations.** `--max-arrows` has never been above 6, and a relation of
eight arrows needs a long line before it means anything.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 17 --max-word 3 --max-arrows 8 --gaps , --jobs 7' --run 'batch.py cores 18 --max-word 3 --max-arrows 8 --gaps , --jobs 7'"
```

**A handful of cores, every length.** Cheap enough to run in the foreground, and
the fastest way to see whether a head or a tail moves. One length per `--run`,
so it is the same command shape as everything else:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 4 --run 'batch.py cores 11 --max-word 4 --cores 45,504,455,605,3344 --jobs 2' --run 'batch.py cores 13 --max-word 4 --cores 45,504,455,605,3344 --jobs 2' --run 'batch.py cores 15 --max-word 4 --cores 45,504,455,605,3344 --jobs 3' --run 'batch.py cores 17 --max-word 4 --cores 45,504,455,605,3344 --jobs 3' --run 'batch.py cores 18 --max-word 4 --cores 45,504,455,605,3344 --jobs 4'"
```

Because `--cores` is a filter, those rows land in the full census's ledger and
a later full run of the length does not redo them.

---

## Menu 2 — `sample`: draw LNAs at a length too big to enumerate

Uniform over all Catalan(n-1) LNAs of the length. Each draw is named by the
theorem, carried by the moves, or is a **leftover**. The summary splits the
leftovers into the ones whose orbit **closed** and the ones whose walk hit the
cap, because those are two different facts and only the first is about the
moves.

### What to vary

| flag | what it changes | values worth running |
|---|---|---|
| `length` | the line | `12 13 14 16 17 18 20 22` |
| `--count` | how many draws | `4000` `20000` `50000` |
| `--seed` | an independent replicate | `0` `1` `2` |
| `--depth` | mutation search out of every leftover | `0` (cheap), `4`, `5` |
| `--orbit-limit` | rows before a draw is a leftover | `20000` (default), `100000` |

`--seed`, `--depth` and `--orbit-limit` are all in the ledger's name. A second
seed is a genuinely independent sample and the honest way to get an error bar.

**A sample cut off by the budget is still a sample.** Each draw is generated
from its own index, and the workers take the indices in order, so stopping
early leaves a prefix — which is as uniform as the whole would have been. Set
`--count` high and let the budget decide. A **census** cut off by the budget is
*not* a sample, because the catalogue is ordered: that is what `--core-limit`
is for.

**`--depth 0` and `--depth 4` are two different jobs, not one.** At `n = 15` a
depth-4 draw cost 257 s of one core on the night of E-045 — the search runs on
every leftover, and at that length most draws are leftovers. Depth 0 drops the
search and leaves the orbit walk, which is itself most of the remaining cost at
a long length and very little of it below about `n = 13`. Run depth 0 at a
large count for the *rate*, and depth 4 at a small count for the *structure*.

Neither is worth guessing at: start a new length with a short foreground run and
read the per-unit rate off the progress lines.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python batch.py sample 17 --count 40 --jobs 4"
```

### Nights worth running

**The rate, at every length, cheaply.** This is the series 0.7, 5.4, 16, 28, …
measured by one instrument instead of three.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py sample 12 --count 50000 --jobs 2' --run 'batch.py sample 13 --count 50000 --jobs 2' --run 'batch.py sample 14 --count 20000 --jobs 2' --run 'batch.py sample 16 --count 20000 --jobs 3' --run 'batch.py sample 18 --count 8000 --jobs 3' --run 'batch.py sample 20 --count 4000 --jobs 2'"
```

**How much of the leftover rate is the cap.** Same draws, a much larger walk;
the two ledgers are separate and the difference between them is the slack.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py sample 15 --count 4000 --orbit-limit 100000 --jobs 7' --run 'batch.py sample 13 --count 20000 --orbit-limit 100000 --jobs 7'"
```

**A second seed, for an error bar.**

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py sample 15 --count 8000 --seed 1 --jobs 7' --run 'batch.py sample 15 --count 8000 --seed 2 --jobs 7'"
```

**Structure, for H-017.** Small count, deep search. The first of these continues
the 879 draws already on disk.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py sample 15 --count 6000 --depth 4 --jobs 7' --run 'batch.py sample 17 --count 1500 --depth 4 --jobs 7'"
```

---

## Menu 3 — the older jobs

These have their own names in `overnight.py` and take no parameters:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --only merges10 merges11 classify10"
```

| name | what it is for |
|---|---|
| `merges10`, `merges11` | H-013, which orbits merge at depth |
| `classify10` | an independent route to the `n = 10` class count |
| `sample14`, `sample16` | H-019 at a fixed count |
| `cores15`, `cores13` | the census E-045 ran |

---

## In the morning

Every task answers to `--summary`, which reads the ledger and does no work. The
flags must match the run, because they are what picks the ledger:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python batch.py cores 16 --max-word 4 --summary"
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python batch.py sample 18 --count 8000 --summary"
```

`overnight.py` prints the right line for every job it ran, at the end of its own
output.

Then **write it down**, whichever way it went, in `research/EXPERIMENTS.md` as
an `E-nnn` — a run that found nothing is recorded precisely so it is not
repeated. `research/README.md` has the conventions.

---

## Things that bite

**A run that stops on its budget exits 2 and is not a failure.** Rerun the same
command and it continues; the ledger holds everything already done.

**The flags pick the ledger.** `--max-word`, `--pair-word`, `--max-arrows`,
`--gaps`, `--orbit-limit`, `--join-limit` and, for `sample`, `--seed` and
`--depth` are all in the filename, because each of them changes what an answer
means. `--cores`, `--core-limit`, `--count` and `--jobs` are not, because they
only change how much of the same question gets asked.

**Ledger rows from before 2026-09-20 mean something slightly different.** The
orbit walk now stops at the first certificate instead of enumerating to its cap,
so `orbit` is rows *walked* and not the orbit's size, and `movesTo` /
`certificate` is the first certificate met rather than the smallest one in the
orbit. New rows carry `stopped` (`cores`) or `orbitClosed` (`sample`); rows
without those fields are the old ones. The verdicts are the same either way.

**`outside` is about this move set.** It says the forward orbit closed without
holding an almost separate LNA. It does not say the algebra is outside the
classification, and no run here can say that.
