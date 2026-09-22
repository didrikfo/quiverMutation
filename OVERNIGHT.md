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
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 13 --max-word 4 --jobs 2' --run 'batch.py cores 15 --max-word 4 --jobs 3' --run 'batch.py cores 17 --max-word 2 --pair-word 2 --gaps 5,6 --jobs 3' --run 'batch.py sample 15 --count 20000 --orbit-limit 100000 --jobs 4' --run 'batch.py sample 17 --count 20000 --jobs 2'"
```

Written 2026-09-22, after the nights of E-046 to E-048. It finishes the two
`--max-word 4` censuses that are still partial, so H-020 has every length from
11 to 16 complete; asks the two-cluster question at the gaps where the clusters
are actually separated (E-047: a gap of 1 or 2 never is); walks `n = 15`'s draws
five times further, which is the one thing that can say whether the leftover
rate still rises past 15 (H-019); and adds draws at `n = 17`.

The previous line -- cores 14, 16, 17 (`--core-limit 250`), the `n = 15`
two-cluster shapes, samples at 13 and 17 -- has been run; E-046 to E-048.

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
| `--gaps` | zeros between the two halves | `5,6` for separated clusters, `1,2,3` (default), empty for none |
| `--cores` | run only these words | `45,504` — `245,2045,2245,2555,2556,3344` |
| `--core-limit` | only the first N words of the catalogue | `40` `60` `120` `400` |
| `--orbit-limit` | rows before a placement is undecided | `20000` (default), `60000`, `200000` |
| `--join-limit` | rows per side for the two-ended join | `6000` (default), `40000` |

**A gap is zeros, not free arrows.** The last relation of the first half
reaches over the zeros, so with every relation two arrows or more a gap of 1 or
2 **never** separates the clusters, gap 3 separates one placement in ten, and
gap 6 every one (E-047). A run meant for two clusters with a free arrow between
them — F-040's shape — wants `--gaps 5,6`; the small gaps ask about clusters
that touch, which is a different and also real question.

`--cores` and `--core-limit` are **filters**: they narrow what a run does
without changing what an answer means, so they share the ledger with the full
census and two nights can split one census between them. Everything else is in
the ledger's name, because it changes what a verdict says.

**The same `--core-limit` asks about the same cores at every length** — the
catalogue is built the same way regardless of `n`. That is what makes two
part-finished censuses comparable, and it is what the last run got wrong.

**`--core-limit` counts catalogue words, not cores.** The catalogue lists every
word the digit ranges allow, sorted by length, and most longer words are not
LNAs at all. `--core-limit 250` at `--max-word 4` is **89** cores; `--core-limit
900` with `--pair-word 3` is 135, every one of them a gap-1 pair. `--plan` shows
the placements, which is the number to size by.

### What a census costs, on one core

Measured on the nights of E-046 and E-047, in core-hours — divide by `--jobs`
for the wall clock:

| command | placements | core-hours | of which undecided |
|---|---|---|---|
| `cores 11 --max-word 4` | 859 | 0.3 | 0 |
| `cores 12 --max-word 4` | 1262 | 0.7 | 0 |
| `cores 14 --max-word 4` | 2148 | 9.9 | 0 |
| `cores 16 --max-word 4` | 3034 | **34.6** | 20.3 |
| `cores 17 --max-word 4 --core-limit 250` | 856 | 11.3 | 6.5 |
| `cores 15 --max-word 2 --pair-word 2 --gaps 1,2,3,4` | 1548 | 18.1 | 3.1 |
| `cores 17 --max-word 2 --pair-word 2 --gaps 1,2,3,4` | 2256 | 41.9 | 24.3 |
| `cores 16 --max-word 2 --pair-word 3 --gaps 1,2,3 --core-limit 900` | 859 | 17.9 | 9.7 |

The table this replaces, estimated from a random sample of placements, said
8.7 core-hours for `n = 16`; the whole census took four times that. **The cap
is what it missed**: at `n = 16` the 157 undecided placements took 20 of the 35
hours, 7.8 minutes each, and every one of them is almost certainly outside
(E-046). Inside answers are now nearly free; outside ones cost what closing the
orbit costs; undecided ones cost the whole cap. So the price of a length is set
by how many placements sit near the cap, which grows with `n`.

By interpolation, what is left of the two partial censuses is a few core-hours
at `n = 13` and about a dozen at `n = 15`; the rest of `n = 17` is of the order
of 50 to 60, and `n = 18` more than that. Before `n = 18`, build the orbit cache
F-051 describes: the census walks the same closed orbits dozens of times — one
orbit 72 times at `n = 14` — and caching them would cut the outside rows about
tenfold.

**Always `--plan` first at a new length.** The number that matters is
`left`, and the cost per unit is what the table above is for.

### Nights worth running

**A run of lengths, for H-020.** The head and the tail should not move with `n`.
*Run 2026-09-21 (E-046): 11, 12, 14 and 16 are complete, and for a single
cluster they do not move. 13 and 15 are still partial; the "tonight" line
finishes them.*

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 11 --max-word 4 --jobs 2' --run 'batch.py cores 12 --max-word 4 --jobs 2' --run 'batch.py cores 14 --max-word 4 --jobs 3' --run 'batch.py cores 16 --max-word 4 --jobs 5'"
```

**Two clusters with a free arrow between them, for F-040 and H-018.** The
first attempt (E-047) ran `--gaps 1,2,3,4` at 15 and 17 and `--pair-word 3
--gaps 1,2,3 --core-limit 900` at 16, and most of it asked about clusters that
touch: a gap of 1 or 2 is never a free arrow, and the `--pair-word 3` cut held
nothing but gap-1 pairs. What it found -- an outside pair always has a half
that is outside alone, but a half can be *rescued* by its neighbour -- wants
asking again at the gaps where the clusters really are apart:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 17 --max-word 2 --pair-word 2 --gaps 5,6 --jobs 5' --run 'batch.py cores 18 --max-word 2 --pair-word 2 --gaps 5,6 --jobs 5'"
```

`--pair-word 3 --gaps 5,6` is 4632 placements at `n = 16`, all of them
separated; at the 75 s a placement the last `--pair-word 3` night averaged, that
is several nights. `--plan` it, and cut it with `--core-limit` knowing the cut
is of words.

**The six undecideds, sharpened.** *Run 2026-09-21 (E-046): all six are
outside, closed orbits of 21709 rows, just past the default cap; heads and tails
unchanged.* Not worth repeating for the undecideds at 16 and 17: every slide
holding one is consistent with the other lengths only if it is outside, and
the cap is spent where F-051 says it will be -- on the biggest outside orbit,
next to the inside end.

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
search and leaves the orbit walk, which is still nearly all of the cost: a
depth-0 draw is **14 s** at `n = 13` and **152 s** at `n = 17` (E-048), almost
all of it spent closing leftover orbits. One worker is about 2300 draws a night
at 13 and about 210 at 17. Run depth 0 at a large count for the *rate*, and
depth 4 at a small count for the *structure*.

Neither is worth guessing at: start a new length with a short foreground run and
read the per-unit rate off the progress lines.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python batch.py sample 17 --count 40 --jobs 4"
```

### Nights worth running

**The rate, at every length, cheaply.** This is the series 0.7, 5.4, 16, 28, …
measured by one instrument instead of three. *`n = 13` has 4575 draws (E-048):
39.9% +- 0.7, which is enough; rerunning it only continues it.*

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py sample 12 --count 50000 --jobs 2' --run 'batch.py sample 13 --count 50000 --jobs 2' --run 'batch.py sample 14 --count 20000 --jobs 2' --run 'batch.py sample 16 --count 20000 --jobs 3' --run 'batch.py sample 18 --count 8000 --jobs 3' --run 'batch.py sample 20 --count 4000 --jobs 2'"
```

**How much of the leftover rate is the cap.** Same draws, a much larger walk;
the two ledgers are separate and the difference between them is the slack.
This is now **the** sampling question (H-019): the closed-orbit rate is flat
from 15 to 17 (48.8%, 47.2%) while the capped share goes from 9% to 22%. `n =
13` has no capped leftovers at all, so there is nothing to sharpen there.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py sample 15 --count 20000 --orbit-limit 100000 --jobs 7' --run 'batch.py sample 17 --count 20000 --orbit-limit 100000 --jobs 7'"
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
