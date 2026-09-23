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
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --dry-run --hours 9 --run 'batch.py cores 16 --max-word 4 --jobs 1'"
```

**Sixteen cores.** Keep the `--jobs` across all the night's jobs at **14 or
fewer**, or they fight each other and the machine gets slower, not faster. A
shared census is one job of `--jobs 1` and holds what it has learned in memory:
187 MB at most for all of `n = 14`. Ten of them side by side is a few GB at
worst at the lengths timed so far; watch it the first time `n = 18` runs.

**Size a job before you pick it.** `--plan` prints the ledger it would write,
how many units it is, and how many of them are already done. It runs nothing:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python batch.py cores 17 --max-word 4 --plan"
```

---

## Tonight, if you have no particular question

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 15 --max-word 4 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 16 --max-word 4 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 17 --max-word 4 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 18 --max-word 4 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 17 --max-word 2 --pair-word 2 --gaps 5,6 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 18 --max-word 2 --pair-word 2 --gaps 5,6 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 15 --max-word 5 --min-word 5 --gaps , --orbit-limit 500000 --jobs 1' --run 'batch.py cores 14 --max-word 6 --min-word 6 --gaps , --orbit-limit 500000 --jobs 1' --run 'batch.py cores 17 --max-word 3 --max-arrows 8 --gaps , --orbit-limit 500000 --jobs 1' --run 'batch.py cores 16 --max-word 2 --pair-word 3 --gaps 5,6 --orbit-limit 500000 --jobs 1' --run 'batch.py sample 17 --count 20000 --orbit-limit 500000 --jobs 1' --run 'batch.py sample 15 --count 20000 --walk shared --orbit-limit 500000 --jobs 1'"
```

Written 2026-09-23, after E-051 and E-052. **The walk is thirty to eighty times
faster than on any earlier night** (E-051: the rule lookup was nine tenths of
every walk), so the cap goes from 20000 to 500000 rows -- enough to close every
orbit of `45` up to `n = 18` -- and a census that last night did 410 placements
in nine hours now does its first 1500 in twelve minutes. Twelve jobs, one
worker each:

* **H-020 finished under the reduced walk**, 15 to 18. Last night's
  `-o20000...-shared` ledgers are not resumed (the limit is in the name) but
  **reused**: every decided verdict in them is taken for free, and only their
  undecided and unasked placements are walked.
* **Two clusters apart** (gaps 5,6) at 17 and 18, finishing last night's, and
  the first `--pair-word 3` run whose pairs are all genuinely apart (16). The
  question is whether any half is rescued under the reduced walk (E-051 says
  none so far).
* **New shapes nothing has looked at**: five-letter cores at 15, six-letter
  cores at 14 (`--min-word` asks only the new width -- the narrower words are
  in other ledgers already), and relations of seven and eight arrows at 17.
* **H-019**: `sample 17` at the new cap, which should close last night's 22.9
  percent capped share as `n = 15`'s closed (E-051); and the first sample under
  the reduced walk, at 15, since every number in that series is the plain
  walk's.

Memory is the thing to watch, not cores. WSL has about 6.8 GB of this laptop's
13.7, and a census at 500000 peaked near 280 MB in the trial; every ledger row
now carries `maxRssMB`. If a job is killed for memory, the ledger has
everything it finished -- rerun it alone, or give WSL more with a
`%USERPROFILE%\.wslconfig` holding `[wsl2]` and `memory=10GB`.

The line before this one (E-051's second night) ran ten shared censuses at the
20000 cap and spent most of the night on placements that hit it.

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
| `--max-word` | how many vertices one core spans | `2` `3` `4` `5` `6` |
| `--min-word` | leave out words shorter than this (a filter) | `5` with `--max-word 5`: only the new width |
| `--max-arrows` | the longest relation in a core | `6` (default), `8`, `10` |
| `--pair-word` | each half of a two-cluster core | `2` `3` |
| `--gaps` | zeros between the two halves | `5,6` for separated clusters, `1,2,3` (default), empty for none |
| `--cores` | run only these words | `45,504` — `45,555,556,3344` |
| `--walk` | how the free move is walked | `shared` (default for `cores`), `reduced`, `plain` — see below |
| `--no-mirror` | ask both halves of each mirror pair | off by default; the mirror halves are read off the other |
| `--no-reuse` | walk every placement even if another ledger of the length has it | off by default -- see below |
| `--core-limit` | only the first N words of the catalogue | `40` `60` `120` `400` |
| `--orbit-limit` | rows before a placement is undecided | `500000` from `n = 15` (E-051); `20000` is the default and too small there |
| `--join-limit` | rows per side for the two-ended join | `6000` (default), `40000` |

**Three walks, and the census wants the shared one.** Every move, the free move
and the mirror is an equivalence, so what one placement's walk settles is
settled for its whole class (GLOSSARY, "Equivalences that save work").

| `--walk` | what it does | ledger | at `n = 11` / `12` |
|---|---|---|---|
| `plain` | deletes two-arrow relations, never adds one; each unit alone | the old names | 464 s / 1516 s, and one derived class can get two verdicts (E-049) |
| `reduced` | adds or deletes them, so `245@0` and `45@1` are one state; each unit alone | `-reduced` | 950 s / 5628 s; 60 more placements inside |
| `shared` | plain then reduced, and **shares** every class it settles with every later unit in the same process | `-shared` | **33 s / 170 s**; exactly the reduced verdicts, every placement (E-050) |

(Seconds of one core on the cloud machine these were measured on, which runs
the plain census about 1.7 times faster than this laptop does.)

**Run a shared census with `--jobs 1`.** The sharing is within one process:
a walk stops at the first row of a class some earlier walk put inside, and
never expands a row an earlier closed orbit has already shown to hold nothing.
Spread over four workers, each re-walks the big orbits for itself, and the
census took 147 s of wall clock against 170 s on one (E-050). So give each
census one worker and use the cores for **more censuses**: lengths, widths and
gap sets side by side, one `--run` each.

**A census reuses every other census of its length** (E-051). Before walking a
placement it looks for the row, or its mirror, in every other `cores-n<n>-*`
ledger: an `inside` from any of them is taken (a certificate is a certificate),
an `outside` only from a shared or reduced ledger (the plain walk's is weaker),
an `undecided` only from one at limits at least as large. So raising
`--orbit-limit`, which starts a new ledger, re-walks only what the old one left
undecided or unasked, and a wider catalogue re-walks nothing the narrower one
did. Reused rows say `"by": "reused"` and name the ledger they came from.
`--no-reuse` walks everything, which is only for timing a census.

**A resumed shared census starts with less memory.** The ledger is still the
ledger, and every `inside` in the other ledgers of the length is handed to the
walk before it starts; what is lost is the closed orbits' rows, which are not on
disk. A census that fits in one night is still the cheap way to run it.

**Do not run two catalogues of one length side by side.** A `--max-word 5`
catalogue starts with the whole `--max-word 4` one, so the two walk the same
placements at the same time and neither can reuse the other: on E-051's night
the two `n = 16` jobs did the same 410 placements verdict for verdict. Run the
narrow one, or give the wide one `--min-word` so it asks only what is new.

**Aliases stay in a shared catalogue.** Under `--walk shared` a word with a `2`
in it costs next to nothing once its class is settled, so the catalogue keeps
them and `--summary` prints their slides too. Under `--walk reduced` it drops
them. `--walk plain` is only for finishing or reading the census ledgers
written before 2026-09-22.

**The mirror is asked once.** `45` at `o` and `504` at `n - 7 - o` are one
question (E-049), so a census asks the first of each pair in catalogue order
and `--summary` prints both slides. It is a filter like `--cores`: the ledger is
the same, and a census already on disk simply has fewer units left.
`--no-mirror` asks both.

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

**Since E-051** (the indexed rule lookup), on this laptop, one worker each:

| command | placements | wall clock |
|---|---|---|
| `cores 13 --max-word 4` | 1457 | 64 s (32 s of CPU), 139 MB |
| `cores 14 --max-word 4` | 1850 | 2.5 min (117 s of CPU), 154 MB |
| `cores 16 --max-word 4 --orbit-limit 500000` | 2637 | 1476 done in 12 min, 1018 of them reused, none undecided; 235 MB |
| `cores 18 --max-word 4 --orbit-limit 500000` | 3424 | `45@2` alone closes at 355328 rows in 510 s; 280 MB |

A plain walk is now about 5600 rows a second at `n = 17` and a reduced one about
2500. The cap is still what sets the price of a length, but it is now a price
in minutes: a placement that walks 500000 rows and fails costs three or four of
them, where one at 20000 cost ten to fifteen before.

**Before E-051, under the shared walk**, measured on the cloud machine of E-050
(about 1.7 times faster than this laptop), one worker each:

| command | placements | wall clock |
|---|---|---|
| `cores 11 --max-word 4` | 709 | 33 s |
| `cores 12 --max-word 4` | 1063 | 170 s |
| `cores 14 --max-word 4` | 1850 | 46 min, 187 MB at most |

On this laptop, allow 1.7 times those: `n = 14` is about an hour and a half on
one core, against 9.9 core-hours plain. The shared census grew sixteenfold from
12 to 14, about as the plain one did (fourteenfold), so the saving held at
about 7.6x. Nothing past 14 has been timed shared. If 14 to 16 grows as it did
plain (3.5x), `n = 16` is four to five hours here; 17 and 18 may not finish in
nine, and that is what the budget is for -- they resume, a little slower for
the empty memory. `--plan` says how many units; a short foreground run says how
fast.

**Under the plain walk**, measured on the nights of E-046 and E-047, in
core-hours — divide by `--jobs` for the wall clock:

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
of 50 to 60, and `n = 18` more than that, under the plain walk. The orbit cache
F-051 asked for before `n = 18` is part of the shared walk now.

**Always `--plan` first at a new length.** The number that matters is
`left`, and the cost per unit is what the table above is for. `--plan` does not
know what reuse will take for free; a five-minute foreground run does.

### Nights worth running

**A run of lengths, for H-020.** The head and the tail should not move with `n`.
*Run 2026-09-21 (E-046) under the plain walk: 11, 12, 14 and 16 are complete,
and for a single cluster they do not move. 13 and 15 are still partial, and are
not worth finishing plain: the next night below asks every length again,
shared. The command now runs shared, one worker per length; add `--walk plain`
to reach the old ledgers.*

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 11 --max-word 4 --jobs 1' --run 'batch.py cores 12 --max-word 4 --jobs 1' --run 'batch.py cores 14 --max-word 4 --jobs 1' --run 'batch.py cores 16 --max-word 4 --jobs 1'"
```

**H-020 again, with the free move walked both ways.** The law was measured by
the plain walk, and at `n = 11` and 12 the reduced walk changes slides it rests
on (`4056` is `oii` at 12; F-052). Before building on H-020, ask it again at
the lengths where it is claimed. *Run 2026-09-23 (E-051): 13 and 14 complete
and the law holds; 15 to 18 stalled on the 20000 cap. Tonight's line finishes
them at 500000, reusing what is done.* On its own, one core per length:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 11 --max-word 4 --jobs 1' --run 'batch.py cores 12 --max-word 4 --jobs 1' --run 'batch.py cores 13 --max-word 4 --jobs 1' --run 'batch.py cores 14 --max-word 4 --jobs 1' --run 'batch.py cores 15 --max-word 4 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 16 --max-word 4 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 17 --max-word 4 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 18 --max-word 4 --orbit-limit 500000 --jobs 1'"
```

The morning's comparison is the one E-046 made -- each single-cluster core's
slide at every pair of lengths -- read off the `-shared` ledgers.

**Two clusters with a free arrow between them, for F-040 and H-018.** The
first attempt (E-047) ran `--gaps 1,2,3,4` at 15 and 17 and `--pair-word 3
--gaps 1,2,3 --core-limit 900` at 16, and most of it asked about clusters that
touch: a gap of 1 or 2 is never a free arrow, and the `--pair-word 3` cut held
nothing but gap-1 pairs. What it found -- an outside pair always has a half
that is outside alone, but a half can be *rescued* by its neighbour -- wants
asking again at the gaps where the clusters really are apart. *Run 2026-09-22
plain (17, complete) and 2026-09-23 shared (17 and 18, partial), E-051: every
rescue on record is a `35` or `36` half, which the reduced walk places alone,
and under the shared walk nothing is rescued so far.*

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 17 --max-word 2 --pair-word 2 --gaps 5,6 --orbit-limit 500000 --jobs 1' --run 'batch.py cores 18 --max-word 2 --pair-word 2 --gaps 5,6 --orbit-limit 500000 --jobs 1'"
```

`--pair-word 3 --gaps 5,6` is 4487 placements at `n = 16` and 8502 at 17, all
of them separated. At the old 75 s a placement it was several nights; since
E-051 it should be one. It is in tonight's line at 16.

**The six undecideds, sharpened.** *Run 2026-09-21 (E-046): all six are
outside, closed orbits of 21709 rows, just past the default cap; heads and tails
unchanged.* Not worth repeating for the undecideds at 16 and 17: every slide
holding one is consistent with the other lengths only if it is outside, and
the cap is spent where F-051 says it will be -- on the biggest outside orbit,
next to the inside end.

**Wider cores, one length.** Everything before 2026-09-23 was words of at most
four vertices. *Run 2026-09-23 (E-051): `n = 14` at `--max-word 5` is complete,
3330 placements. The `n = 16` run did nothing but repeat the `--max-word 4`
census beside it.* Ask only the new width with `--min-word`:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 15 --max-word 5 --min-word 5 --gaps , --orbit-limit 500000 --jobs 1' --run 'batch.py cores 16 --max-word 5 --min-word 5 --gaps , --orbit-limit 500000 --jobs 1' --run 'batch.py cores 14 --max-word 6 --min-word 6 --gaps , --orbit-limit 500000 --jobs 1'"
```

(`--gaps ,` leaves the pairs out, so the night is single cores only.)

**Longer relations.** `--max-arrows` has never been above 6, and a relation of
eight arrows needs a long line before it means anything.

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py cores 17 --max-word 3 --max-arrows 8 --gaps , --jobs 1' --run 'batch.py cores 18 --max-word 3 --max-arrows 8 --gaps , --jobs 1'"
```

**A handful of cores, every length.** Cheap enough to run in the foreground, and
the fastest way to see whether a head or a tail moves. One length per `--run`,
so it is the same command shape as everything else:

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 4 --run 'batch.py cores 11 --max-word 4 --cores 45,504,455,605,3344 --jobs 1' --run 'batch.py cores 13 --max-word 4 --cores 45,504,455,605,3344 --jobs 1' --run 'batch.py cores 15 --max-word 4 --cores 45,504,455,605,3344 --jobs 1' --run 'batch.py cores 17 --max-word 4 --cores 45,504,455,605,3344 --jobs 1' --run 'batch.py cores 18 --max-word 4 --cores 45,504,455,605,3344 --jobs 1'"
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
| `--orbit-limit` | rows before a draw is a leftover | `20000` (default), `100000`, `500000` |
| `--walk` | as for `cores` | `plain` (default, the series so far), `shared` |

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
depth 4 at a small count for the *structure*. *All of these are the timings
before E-051, whose rule lookup makes the orbit walk thirty to eighty times
faster; the search at depth 4 runs the mutation engine and is not affected.*

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
*Run at `n = 15` and 100000, 1444 draws (E-051): every leftover closed, 57.4%
+- 1.3, so the capped share there was all outside.* `n = 17` is what is left,
at 500000 (in tonight's line).

```bash
wsl -e bash -lc "cd /mnt/c/Users/didri/kode/quiverMutation && .venv/bin/python overnight.py --hours 9 --run 'batch.py sample 17 --count 20000 --orbit-limit 500000 --jobs 2' --run 'batch.py sample 19 --count 20000 --orbit-limit 500000 --jobs 2'"
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
means. `--cores`, `--core-limit`, `--min-word`, `--no-reuse`, `--count` and
`--jobs` are not, because they only change how much of the same question gets
asked, or how. Tonight's censuses run at `--orbit-limit 500000`, so their
`--summary` needs it too.

**A wider catalogue holds the narrower one.** Catalogues are sorted by word
length, so `--max-word 5` begins with all of `--max-word 4`. Two of them at one
length on one night do the same work twice (E-051); use `--min-word`.

**Memory, not cores, is the limit at `--orbit-limit 500000`.** WSL has about
6.8 GB. Each ledger row carries `maxRssMB`, the process's peak so far: read the
last row of each ledger in the morning, and run fewer jobs, or raise WSL's
share in `.wslconfig`, if they add up to more than about five.

**Ledger rows from before 2026-09-20 mean something slightly different.** The
orbit walk now stops at the first certificate instead of enumerating to its cap,
so `orbit` is rows *walked* and not the orbit's size, and `movesTo` /
`certificate` is the first certificate met rather than the smallest one in the
orbit. New rows carry `stopped` (`cores`) or `orbitClosed` (`sample`); rows
without those fields are the old ones. The verdicts are the same either way.

**`outside` is about this move set.** It says the forward orbit closed without
holding an almost separate LNA. It does not say the algebra is outside the
classification, and no run here can say that.
