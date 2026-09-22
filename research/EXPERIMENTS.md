# Experiments

Runs made, newest first, with parameters and outcome — including runs that found
nothing, which are recorded precisely so they are not repeated. See
[`README.md`](README.md).

---

## E-049 — A core and the same core with a two-arrow relation, walked both ways
*2026-09-22* · **the plain walk gave one derived class two verdicts 19 times at `n = 11` and 12; walked as one state, 60 placements move from outside to inside and none the other way; the census also asks every mirror pair twice** → F-052, H-020, F-051

**The question.** F-051 lists `45` at offset 1 and `245` at offset 0 as two of
"four different cores" in one orbit. They are one LNA up to a relation of two
arrows, so one derived class by the free move (F-028). Does the census treat
them as one, and what does asking both cost?

**How often it happens.** Under `--max-word 4 --gaps 1,2,3` a word holding a
`2` is its stripped word at another offset, and that stripped row is always
another placement of the same catalogue:

| n | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|
| placements | 859 | 1262 | 1705 | 2148 | 2591 | 3034 | 3477 |
| holding a two-arrow relation | 189 | 249 | 309 | 369 | 429 | 489 | 549 |
| stripped row also in the census | 189 | 249 | 309 | 369 | 429 | 489 | 549 |

Five of E-045's six undecideds -- `245`, `2045`, `2245`, `2555`, `2556` -- are
such words: `45` at one or two offsets further on, `555` and `556` at one.

**The two did not get the same verdict.** Every alias and its stripped
placement were run through `_verdictFor` at `n = 11` and 12 with the default
limits. **6 of 189 and 13 of 249 disagree**, always the same way round: the
alias `inside`, the stripped row `outside`. The smallest: `404` at offset 1 of
`n = 11` has **no move at all** (a closed orbit of one row), while `2404` at
offset 0 reaches an almost separate row in nine. The cause is that
`freeMoves.movesFrom(free = True)` only deletes two-arrow relations. The free
move is symmetric and the walk made it one-way, and the added relation is
exactly the spectator a rule needs (F-023).

**Walked as one state.** `freeMoves.REDUCED`: every state is a reduced row, and a
step is a move out of it or out of it with one two-arrow relation added,
stripped again. Every reduced placement at `n = 11` and 12 was run under both
walks and compared with the best verdict any of its aliases had under the plain
walk:

| n | placements, plain | reduced | reduced worse than the best alias | outside → inside | core-seconds, plain | reduced |
|---|---|---|---|---|---|---|
| 11 | 859 | 670 | **0** | 20 | 464 | 950 |
| 12 | 1262 | 1013 | **0** | 40 | 1516 | 5628 |

So the reduced walk loses nothing any alias had and places 60 placements the
plain walk called outside. The `45` and `504` slides are unchanged at both
lengths (`iooii`, `ioooii`; `iiooi`, `iioooi`), so F-042 stands. The ones that
move include `404` (`ioooi` → `iiiii` at 11, `iooooi` → `iiiiii` at 12),
`405`, `5004`, `36`, `6006` (`ooo` → `iii` at 12), and `4056` (`oio` →
`oii` at 12), the one shape H-020 had to amend its statement for.

**It costs more, and that is not a saving.** A reduced step tries every vertex
a two-arrow relation can be added at, so a row costs about twice as much, and
the outside orbits it has to close are the dearest part: 307 outside placements
took 4145 core-seconds at `n = 12`, against 1311 for 389 under the plain walk.
Walking plain first and reduced only where plain did not find a certificate
barely helps (5454 against 5628), because it is the outside orbits and not the
inside answers that cost. The factor went from 2.0 at 11 to 3.7 at 12;
nothing longer has been measured.

**The mirror is a second duplication.** The relation dual sends a row to a row
of the same class and the move set is closed under it (F-026), so the verdicts
of a row and its mirror should agree. Over the same runs: **300 of 300** mirror
pairs agree at `n = 11` under the plain walk, 288 of 288 under the reduced one,
398 of 398 and 384 of 384 at `n = 12`. Asking one of each pair removes 12 to 17
percent of a census (3034 → 2637 placements at `n = 16` under the plain walk,
2545 → 2159 under the reduced one).

**What changed in the code.** `batch.py cores` and `batch.py sample` take
`--walk plain|reduced`; the default is `plain`, so every command in
`OVERNIGHT.md` means what it meant, and `reduced` writes a ledger whose name
ends `-reduced`. Under `reduced` the catalogue drops every word holding a `2`.
Under both, the census asks only the first of each mirror pair in catalogue
order and `--summary` reads the other half of each slide in the mirror;
`--no-mirror` turns that off. The mirror is a filter, not an instrument, so
the ledger is the same one.

Reproduce: the comparison is `_verdictFor(n, word, offset, 20000, 6000, free =
...)` over `batch._placements` with `--walk plain` and `--walk reduced`, which
at `n = 12` is about 2 core-hours. `tests/test_reduced_walk.py` pins the `404`
case, the catalogue and the mirror.

---

## E-048 — The leftover rate at `n = 13` and `n = 17`, by one instrument
*2026-09-22* · **`n = 13` is 39.9% +- 0.7 with no cap in it; `n = 17` is somewhere between 47% and 70%** → H-019

Two depth-0 samples from `OVERNIGHT.md`'s first "tonight" line, run on the night
of 2026-09-20 alongside the censuses of E-046 and E-047. Both stopped on the
9-hour budget, which for a sample is a prefix and still uniform.

| n | drawn | core-h | s/draw | theorem | moves | **leftover, orbit closed** | leftover, orbit capped |
|---|---|---|---|---|---|---|---|
| 13 | 4575 | 18.0 | 14.1 | 13.5% | 46.6% | **39.9% +- 0.7** | **0** |
| 15 *(E-045, depth 4, old rows)* | 879 | 62.7 | 257 | 7.7% | 34.2% | **48.8% +- 1.7** | 9.2% |
| 17 | 214 | 9.0 | 152 | 3.7% | 26.6% | **47.2% +- 3.4** | 22.4% |

**`n = 13` is the first clean point above `n = 11`.** Every one of its 1825
leftovers had a forward orbit that *closed*, so the rate carries no cap slack at
all -- it is exactly "the moves do not reach an almost separate row", the same
statement as F-032's exhaustive 0.6, 5.4 and 16. The series by one definition is
now 0.7, 5.4, 16, 28 (60 draws), **39.9**.

**Above that the cap decides the answer.** The closed-orbit fraction, which is a
lower bound on the leftover rate, reads 48.8% at `n = 15` and 47.2% +- 3.4 at
`n = 17`: flat. The capped fraction, which could go either way, goes from 9% to
22%. So "still rising past 15" and "levelling off near a half" are both
consistent with what is on disk, and only a larger `--orbit-limit` separates
them. The `n = 15` row is E-045's, whose rows predate the closed/capped split;
its 81 capped leftovers were counted there from the orbit sizes.

**What it costs, which the menu had wrong.** `OVERNIGHT.md` said a depth-0 draw
is about a millisecond below `n = 13`. At `n = 13` it is **14 s**, and 17.5 of
the 18 core-hours went on leftovers -- closing an orbit is the whole cost, and
at `n = 13` that is most draws. At `n = 17` a draw is 152 s, and the 48 capped
ones took 4.7 of the 9 core-hours between them. One worker at `n = 17` is 214
draws a night; the `--count 20000` on the menu was two orders of magnitude past
what a night does, which is harmless (the budget cuts it) but should not be read
as a plan.

**By overlap, at `n = 13`**, the leftovers are 254 at overlap 2, 696 at 3, 509
at 4, 251 at 5 and 115 above. Overlap 2 -- F-042's cell `(2, 1)` and its
neighbours -- is one leftover in seven, as it should be if no bound on the overlap
cuts the leftovers off.

**A correction to E-044 found on the way.** Its table gives `n = 12` 208012
LNAs. That is Catalan(12), the count at `n = 13` (this run's summary prints it
for `n = 13`); `n = 12` has Catalan(11) = 58786. The rate there is unaffected,
since the draw is by index and the count is only printed.

Reproduce:

```bash
python batch.py sample 13 --count 50000 --summary
python batch.py sample 17 --count 20000 --summary
```

---

## E-047 — Two heavy clusters in one line, at lengths with room for them
*2026-09-22* · **outsiders with a free arrow between the clusters exist from `n = 14`; every one has a half that is outside on its own** → H-018, F-040

F-040 found no LNA outside a quipu class with two heavy clusters and a free arrow
between them, at every length up to 11, and noted the shape barely fits there.
Three runs from `OVERNIGHT.md` asked the same of the move closure at lengths
with room:

| run | placements | core-h | inside | outside | undecided | free-gap rows |
|---|---|---|---|---|---|---|
| `cores 15 --max-word 2 --pair-word 2 --gaps 1,2,3,4` | 1548 | 18.1 | 776 | 751 | 21 | 130 |
| `cores 17 --max-word 2 --pair-word 2 --gaps 1,2,3,4` | 2256 | 41.9 | 1017 | 1029 | 210 | 210 |
| `cores 16 --max-word 2 --pair-word 3 --gaps 1,2,3 --core-limit 900` | 859 | 17.9 | 391 | 401 | 67 | **0** |

The single-core censuses of E-046 hold pair words too (their default `--gaps
1,2,3`), and contribute 4, 10, 30 and 50 free-gap rows at `n = 11, 12, 14, 16`.

**The count F-040 made, at lengths it could not reach.** Over all 434 rows with
two heavy clusters and a free arrow between them: **388 inside, 32 outside, 14
undecided.** None at `n = 11` or 12 is outside. The first is at `n = 14`:
`330004500000`, a `33` against the source and a `45` five arrows in, whose
forward orbit closes at 1895 rows without an almost separate member.

**Every outsider is explained by one of its halves.** Looking each half up alone,
at the same offset in the same length:

| left half alone | right half alone | the pair | rows |
|---|---|---|---|
| inside | inside | inside | **364** |
| inside | inside | outside | **0** |
| inside | outside | outside | 32 |
| inside | outside | **inside** | **24** |
| inside | outside | undecided | 14 |

So with a free arrow between them, two placeable clusters make a placeable pair
every time (364 of 364), and each of the 32 outside pairs has a half that is
outside where it sits. The `45` in `330004500000` sits in its own outside band
(`45` at `n = 14` is `ioooooii`; it is at offset 5). What looked like the case
that would break the single-cluster reading is a single-cluster outsider with a
harmless neighbour.

**But the halves are not independent, and the 24 say so.** A `33`, `34` or `44`
at or near the source, three or four arrows before a `35` or `36`, carries the
pair inside where the `35` or `36` alone is outside: `3300035` at offsets 0 to 2
at `n = 17`, `33000035` at 0 and 1, and so on through `34` and `44`, at `n = 15,
16, 17`. The rescue reaches less far as the gap widens -- three offsets at gap 3,
two at gap 4 -- and every right half involved is one whose own slide has an
inside head (`35`: `ii…`, `36`: `i…`). One reading is that the left cluster,
pushed off the source, lets the right one reach its own head; nothing here
checks that.

**Without a free arrow the interaction goes the other way too.** Among the rows
where the two clusters touch through a one-arrow overlap, eight have both halves
inside alone and the pair outside -- all of them `35 0^g xx` at offset 1
(`3500033` at `n = 11`, `3500034` and `3500044` at 12, `3500036`, `3500046`,
`3500056`, `3500066` at 14, `3500035` at 15), where the `5` reaches over the gap
into the second cluster.

**The run that answered nothing, and why.** `--pair-word 3 --core-limit 900` at
`n = 16` has no free-gap row at all. The catalogue is sorted by word length, so
its first 900 words are the gap-1 pairs -- and **a gap of one or two never has a
free arrow** when every relation is two arrows or more: the last relation of the
first half reaches over it. Measured on the catalogue at `n = 16`:

| pair-word 2, gap | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| placements | 406 | 496 | 500 | 400 | 300 | 200 |
| with a free arrow | 0 | 0 | 50 | 120 | 180 | **200** |

and for `--pair-word 3`, 0 of 6135 at gap 1, 170 of 7520 at gap 2, and all 1447
at gap 6. The night's 859 placements at `--pair-word 3` are real rows about
touching clusters and are kept, but the free-gap question needs `--gaps 5,6`.

Reproduce:

```bash
python batch.py cores 15 --max-word 2 --pair-word 2 --gaps 1,2,3,4 --summary
python batch.py cores 17 --max-word 2 --pair-word 2 --gaps 1,2,3,4 --summary
python batch.py cores 16 --max-word 2 --pair-word 3 --gaps 1,2,3 --summary
```

The half-by-half table is not something `--summary` prints; it came from joining
each pair row to the single-core rows of the same length on `(core, offset)`.

---

## E-046 — The core census at lengths 11 to 17, and the six undecideds resolved
*2026-09-22* · **H-020 holds for every single cluster from `n = 13` to 17; the undecideds were all outside; the census walks the same orbits over and over** → H-020, F-051

Four `overnight.py` invocations between 2026-09-20 19:15 and 2026-09-21 23:17,
all from `OVERNIGHT.md`: the first "tonight" line, the "run of lengths", the
"two-cluster shapes" (E-047) and the "six undecideds, sharpened". The
`--max-word 4` ledgers now stand at:

| n | placements | cores | inside | outside | undecided | complete | core-h |
|---|---|---|---|---|---|---|---|
| 11 | 859 | 337 | 675 | 184 | 0 | yes | 0.3 |
| 12 | 1262 | 403 | 873 | 389 | 0 | yes | 0.7 |
| 13 | 367 | 62 | 287 | 80 | 0 | no (E-045) | -- |
| 14 | 2148 | 443 | 1252 | 896 | 0 | **yes** | 9.9 |
| 15 | 932 | 125 | 666 | 260 | 6 | no (E-045) | -- |
| 16 | 3034 | 443 | 1599 | 1278 | 157 | **yes** | 34.6 |
| 17 | 856 | 89 | 636 | 196 | 24 | `--core-limit 250` | 11.3 |

**H-020, asked properly.** Every core with one heavy cluster that is decided at
two lengths from 13 to 17 was compared between them: **975 comparisons over 283
cores, no failure** of "the longer slide is the shorter one with the interior
verdict repeated". The `45` core, F-042's example, reads `iooii`, `ioooii`, …,
`iooooooooii` from `n = 11` to 17 -- head 1, tail 2 at every length -- and its
opposite `504` is head 2, tail 1 at every length. At `n = 14` the 343
single-cluster cores split 118 inside everywhere, 97 outside everywhere, 127
inside at the ends and outside between, and **one** with a shape the head/tail
reading does not cover: `4056`, `oio` / `oooio` / `oooooio` at 12, 14, 16, which
is inside one offset from the sink and outside at the sink itself. Its suffix
`io` is as fixed as any tail; it is a word, not a count. No single-cluster slide
is inside in its interior and outside at an end.

**Below `n = 13` the law is not visible yet, as F-042 predicted.** Twelve
single-cluster cores change shape between 11 and 12 or 12 and 14 -- `44066`
`io` → `oooo`, `6600044` `oi` → `oooo`, `550066` `i` → `ooo` -- and every one of
them is a long word whose slide at the shorter length is three offsets or
fewer, where a placement is near both ends at once.

**Two-cluster cores do not obey it, and should not.** 31 incompatibilities
between lengths, all from pair words: `3500035` is `iio` at 14 and `ooooo` at 16,
because its right `35` moves away from the sink as the line grows while its left
one stays at the source. E-047 is what they obey instead.

**The six undecideds of E-045, at three times the orbit limit.** `--orbit-limit
60000 --join-limit 40000 --cores 245,2045,2245,2555,2556,3344` at `n = 15`, 43
placements, 2.3 core-hours: **all six are outside**, each a closed orbit of
**21709 rows** -- just past the 20000 the census had. The heads and tails E-045
read off are unchanged. And across the whole census, the 92 slides that hold an
undecided are every one consistent with the other lengths when `?` is read as
`o`, and two of them also when read as `i`; none only as `i`. The undecideds
are the largest outside orbits, which sit next to the inside end (F-051), and
they are almost certainly outside.

**Where the time went -- and it is the opposite of E-045.** The inside answers
are now nearly free: 1.8 of `n = 16`'s 34.6 core-hours. Outside took 12.5, and
**the 157 undecided took 20.3**, 7.8 minutes each, for placements every other
line of evidence says are outside. The cost table in `OVERNIGHT.md` estimated
8.7 core-hours for `n = 16`; the real figure is four times that, and it is the
cap that the estimate missed.

**The census walks the same orbits again and again.** Placements with the same
closed-orbit size are frequent -- at `n = 14`, the 896 outside placements have
**107** distinct orbit sizes between them, and one size, 7393, was walked 72
times for 3.7 core-hours. Walking a handful of them again and comparing the sets
(F-051): same size is same orbit, or a mirror pair of orbits. If each distinct
orbit were walked once, the outside rows would cost 0.5 core-hours instead of
9.2 at `n = 14`, 1.0 instead of 12.5 at `n = 16`, and 2.1 instead of 13.9 for
the `n = 17` pairs.

**A catalogue word is not a core.** `--core-limit 250` at `n = 17` gave 89
cores, not 250: the catalogue lists every word the digit ranges allow, and most
four-letter words are not LNAs at all (their relations' ends do not increase)
and have no placement at any length. The cut is still the same at every length,
which is what it is for, but it is a cut of the catalogue, not a count of cores.

Reproduce:

```bash
python batch.py cores 14 --max-word 4 --summary
python batch.py cores 16 --max-word 4 --summary
python batch.py cores 15 --max-word 4 --orbit-limit 60000 --join-limit 40000 --summary
```

The cross-length comparison is not in `--summary`: it reads every
`cores-n*-w4p2a6g123-o20000j6000.jsonl` ledger, builds each core's slide per
length, and checks that each pair of decided slides differs by a run of the
interior verdict.

---

## E-045 — The length-15 night: two censuses and a sample, all three cut off
*2026-09-20* · **no census finished; the instrument was the bottleneck, not the machine** → H-018, H-019, H-020

The night `overnight.py --hours 9 --only sample15 cores15 cores13` was left to
run. All three jobs used their whole budget and all three stopped on it, with
the work they had done on disk:

| job | units done | of | core-hours | what it was for |
|---|---|---|---|---|
| `cores15` | 932 | 2591 | 53.8 | H-018, the census at a length with room |
| `cores13` | 367 | 1705 | 18.0 | the control length F-042 already covered |
| `sample15` | 879 | 6000 | 62.7 | H-019, the leftover rate at `n = 15` |

The machine was not the problem: 134 core-hours came back from 15 workers in 9
hours, which is a saturated machine. **The instrument was.** Three things came
out of reading the ledgers, and the first is much the largest.

**1. The walk was paying for the whole orbit to answer a membership question.**
Of `cores15`'s 932 placements, the 666 that came back `inside` cost 49 of the
53.8 core-hours; the 260 `outside` ones cost 3.8. That is backwards -- `inside`
is the *easy* answer -- and the reason was that `_verdictFor` called
`freeMoves.orbitOf`, which enumerates to its 20000-row cap, and only then looked
through the result for an almost separate row. Re-running a random twelve of
those `inside` placements with the walk stopping at the first almost separate
row it meets: every one of them found its certificate within 351 rows, nine of
the twelve within 100, and the twelve together took **5.8 seconds against the
2818 seconds they cost on the night** -- 486x. `sampling.probe` had the same
shape and the same 15.5 core-hours of `moves` rows to show for it.

`freeMoves.orbitReport` is the fix, with a `stopWhen` predicate and, because a
walk that stops early must not be mistaken for one that closed, an explicit
`stoppedBy` of `closed`, `found` or `cap`. The verdicts are unchanged; what
changed is the price. Measured afterwards, single core, over a random sample of
each catalogue:

| census | placements | median | mean | one core |
|---|---|---|---|---|
| `cores 13 --max-word 4` | 1705 | 0.18s | 2.35s | **1.1 h** |
| `cores 15 --max-word 4` | 2591 | 0.28s | 5.27s | **3.8 h** |
| `cores 16 --max-word 4` | 3034 | 0.61s | 10.35s | **8.7 h** |
| `cores 17 --max-word 4` | 3477 | 2.66s | 38.82s | **37.5 h** |

`cores15` was four nights of work and is now half an hour on eight workers. The
catalogue grows slowly and the price of a placement does not: `n = 17` is still
a night, `n = 18` at this width is several.

**2. The two censuses walked their catalogues at different speeds, so only
their overlap could be read.** `cores13` reached 62 cores and `cores15` reached
125, both from the front of the same catalogue. A census read against another
length is the entire point of running one, and the budget, not the question,
decided where each stopped. `--core-limit` and `--cores` now cut the catalogue
deliberately instead, and neither is in the ledger's name, so two nights can
split one census between them.

**3. The shapes the length was chosen for got no units at all.** `--gaps` puts
two-cluster cores in the catalogue after every single core, and at the rate the
night ran neither length reached them: "rows with two heavy clusters and a free
arrow between them: 0" in both summaries. F-040's question -- the one that needs
`n >= 13` to be askable -- was not asked. It needs its own run, which
`--max-word 2 --pair-word 2` gives, and that is now in `OVERNIGHT.md`.

**What the partial data looks like, recorded and not concluded from.** Of the 62
cores done at both lengths, 25 have an `outside` somewhere in their slide. Every
one of those 25 has the same shape at both lengths: some number of inside
offsets at the source end, some number at the sink end, and outside everywhere
between, with **both counts identical at `n = 13` and `n = 15`** and the outside
middle absorbing the two extra offsets. Not one of the 25 has an inside in its
interior. H-020 is that observation written down as something to test; it rests
on two lengths, one of them a third finished, and is not a finding.

Two smaller things worth having on the record:

* **Six placements came back `undecided`, and five of the six sit exactly at the
  boundary** between the outside middle and the inside tail (`ooooo?ii`,
  `oooo?ii`, `ooooo?i`), the sixth at the head boundary (`ii?ooooo`). The hard
  cases are where the answer changes, which is where a sharper limit is worth
  spending on: `--cores 245,2045,2245,2555,2556,3344 --join-limit 40000`.
* **81 of the 510 sampled leftovers, one in six, had walks that hit the cap**
  rather than closing. `settledBy == 'leftover'` was covering both, so a
  leftover *rate* at `n = 15` is part rate and part budget. `probe` now records
  `orbitClosed` and the summary splits them, and `--orbit-limit` has joined the
  sampler's ledger name, which it should always have been in -- the same
  mistake `--pair-word` made in `cores` before the first night.

Reproduce, or continue:

```bash
python batch.py cores 15 --max-word 4 --summary
python batch.py cores 13 --max-word 4 --summary
python batch.py sample 15 --count 6000 --depth 4 --summary
```

The ledgers from the night are kept and the new runs resume from them. Their
rows carry no `stopped` or `orbitClosed` field, which is how a row written
before the walk learned to stop early can be told from one written after. The
sample ledger was renamed from `sample-n15-s0-d4.jsonl` to
`sample-n15-s0-d4-o20000.jsonl`, because `--orbit-limit` was missing from the
name and two runs under different limits disagree about what a leftover is.

`OVERNIGHT.md` is the menu this run produced: what each parameter changes, what
a census of each length costs, and the runs worth making next.

---

## E-044 — The sampler, calibrated against the lengths whose answer is known
*2026-09-19* · **agrees at `n = 9`, 10 and 11; `n = 12` is 28% leftover** → H-019

`batch.py sample` is only worth running at `n = 16` if it gives the right answer
at `n = 10`, where the answer is known exhaustively. This is that check. The
draw is uniform over all Catalan(n - 1) LNAs of the length; each is put through
the cheap pipeline one row at a time (`sampling.probe`), where F-032 measured the
same thing a length at a time.

| n | LNAs | drawn | theorem | moves | **leftover** | exhaustive leftover |
|---|---|---|---|---|---|---|
| 9 | 1430 | 60 | 40% | 58% | **1.7% +- 1.7** | **0.70%** (F-032: 0.6%) |
| 10 | 4862 | 400 | 31.2% | 64.0% | **4.8% +- 1.1** | **5.4%** (F-032) |
| 11 | 16796 | 6 | 33% | 50% | 17% +- 15 | 16% (F-032) |
| 12 | 208012 *(sic: 58786 -- E-048)* | 60 | 27% | 45% | **28% +- 5.8** | not known |

**`n = 10` is the load-bearing row**: 400 draws put the leftover rate at
4.8% +- 1.1, and the exhaustive answer is 5.4%. `n = 11` drew only 6 before the
run's time limit — the orbit walk is much slower there — so its agreement is
suggestive and no more.

**The exhaustive `n = 9` pass was run here too**, over all 1430 rows, as a check
that `probe` asked one row at a time agrees with `overlaps.py` asked a length at
a time: 610 by the theorem, 810 by the moves, **10 leftovers, 0.70%**. It does.

**`n = 12` is the new number.** 28% +- 5.8, against 0.7% at `n = 9` — the rate
is rising steeply and the sampler is the only way to see it above `n = 11`.
Preliminary: 60 draws, one seed. H-019 is the hypothesis this feeds and says
what would settle it.

**A leftover rate from a sample is an upper bound, not an estimate of the
truth.** `freeMoves.orbitOf` walks forwards only and stops at `--orbit-limit`,
so a row it fails to place may still be placeable — by a longer walk, or by a
walk from the other direction, which is the asymmetry E-037 recorded. The rates
above are therefore "not placed by this much walking", and the exhaustive
figures they are checked against carry the same caveat.

**The run that was cut off resumed correctly**, which was not the point of the
experiment but is worth recording: `n = 11` stopped at 6 of 400 draws on a time
limit, and the ledger holds those 6, so the same command continues from draw 7.

Reproduce:

```bash
python batch.py sample 10 --count 400 --jobs 4
python batch.py sample 10 --summary
```

---

## E-043 — The deduplicated walk, checked against the plain walk
*2026-09-19* · **same answers everywhere it was checked; 2.2x to 6.1x faster** → F-049, F-050

The dedup of E-042 is only worth having if it changes no answer. Checked by
running both walks over **every LNA of a length** and comparing three things at
once: the set of lines collected, the set of hereditary forms collected, and the
set of canonical keys the visitor was shown.

| length | depth | LNAs | mismatches | nodes plain | nodes deduped | ratio |
|---|---|---|---|---|---|---|
| 5 | 4 | 14 | **0** | 847 | 674 | 1.3x |
| 6 | 5 | 42 | **0** | 14701 | 7427 | 2.0x |
| 7 | 5 | 132 | **0** | 94498 | 39280 | 2.4x |
| 8 | 4 | 429 | **0** | 143068 | 78295 | 1.8x |

**It failed first, and the failure was real.** Before the sign gauge was
quotiented out, `n = 7`, `30300`, depth 5 lost one node: the walk reached one
algebra under two presentations differing only in the sign of a zero relation,
which are the same ideal, so the dedup was right to identify them — and the
*procedure* then gave two different answers from them. That is F-050, and it was
found by this check rather than by reading.

**Wall clock, on the searches the merge hunt runs** (`n = 9`, one container, one
process):

| start | depth | plain | deduped | speedup |
|---|---|---|---|---|
| `3345000` | 4 | 8.0 s | 3.6 s | 2.2x |
| `3345000` | 5 | 38.4 s | 10.7 s | 3.6x |
| `3345000` | 6 | 184.9 s | 30.5 s | **6.1x** |
| `3033030` | 4 | 4.9 s | 2.4 s | 2.1x |
| `3033030` | 5 | 18.4 s | 5.9 s | 3.1x |
| `3033030` | 6 | 77.1 s | 15.2 s | **5.1x** |

The speedup is below the node ratio because the key costs something at every
node; it grows with depth for the same reason the node ratio does.

**What the dedup does not preserve is the number of *visits*.** At `n = 9`,
`3345000`, depth 6 the plain walk collects 54 line records and the deduped walk
13 — and the **set** of LNAs reached is identical, 4 either way, at depths 4 and
5 as well. Every consumer in the repo builds a set, so this is invisible to all
of them; a caller that counted records would be counting routes, which was never
a meaningful number.

**And end to end, through `merges.py` itself.** The same run at `n = 9`, depth 4,
`--all-groups`, once with the dedup and once with `--no-dedupe`: all 9 members
searched, **every `lnasReached`, `reachedOrbits`, `coveredOrbits` and `alarms`
identical**, 121.5 s against 215.1 s of search time — 1.8x. That is the shallow
end; the depths the overnight job runs at are 5 to 8, where the single-search
measurement above is 3.6x to 6.1x.

`merges.py` now passes a `fingerprint.Visited` by default and records what it
skipped in the checkpoint's `walk` field, so a run's cost can be read back
rather than re-measured. `--no-dedupe` restores the old walk, for measuring what
the dedup changes rather than for producing answers.

Reproduce: `tests/test_fingerprint.py::test_dedup_reaches_exactly_what_the_plain_walk_reaches`
pins the n = 5 to 7 cases at the depths above.

```bash
python merges.py 9 --depths 4 --all-groups --checkpoint logs/a.jsonl
python merges.py 9 --depths 4 --all-groups --no-dedupe --checkpoint logs/b.jsonl
```

---

## E-042 — How much of a mutation search is repeated work
*2026-09-19* · **3.6x at depth 4, 6.3x at 5, 11.4x at 6, and compounding** → F-049

`mutationSearchDepthFirst` walks mutation *sequences* and has never deduplicated
the algebras they reach. This counts what that costs, over the two leftover
orbits at `n = 9` that `merges.py` searches from, by keying every node the
visitor is shown with `search.quiverKey`.

| start | depth | nodes | distinct | ratio | parallel-arrow nodes |
|---|---|---|---|---|---|
| `3345000` | 4 | 779 | 216 | 3.6x | 0 |
| `3345000` | 5 | 3887 | 615 | 6.3x | 0 |
| `3345000` | 6 | 19483 | 1708 | **11.4x** | 40 (0.2%) |
| `3033030` | 4 | 350 | 108 | 3.2x | 4 (1.1%) |
| `3033030` | 5 | 1435 | 262 | 5.5x | 44 (3.1%) |
| `3033030` | 6 | 6088 | 639 | **9.5x** | 284 (4.7%) |

**The ratio roughly doubles per level**, which is the number that matters: the
waste is not a constant overhead to be shrugged at but the dominant term at the
depths the merge hunt wants. Extrapolating the two columns, depth 8 is 30x to
40x, and depth 8 at `n = 10` is exactly what `overnight.py` is running.

**Why there are so many routes.** A mutation is invertible and mutations at
distant vertices commute, so the number of sequences reaching a given algebra
grows with the depth faster than the number of algebras does. Nothing about
this is specific to these two starts.

**The parallel-arrow census, which decides whether a canonical form is hard.**
Over the two depth-5 walks, 44 of 5322 nodes have parallel arrows, and **every
one of them has exactly one bundle of exactly two arrows**. So the arrow-key
ambiguity that NOTES.md has recorded as an open want since F-039 costs a
minimisation over 2 relabelings, not a search. There was no node with two
bundles, and none with a bundle of three.

**Vertex labels do not move under mutation**, so none of this is graph
isomorphism testing: two algebras reached from one start are equal on the nose
or not at all. That is why the key is exact and cheap, and why the probabilistic
fingerprint this measurement was meant to justify turned out not to be needed
for identification at all — only, in `fingerprint.digest`, for storing a very
large visited set in less memory.

Reproduce:

```python
from quivermutation import nakayama as nk, search, fingerprint
alg = nk.LinearNakayamaAlgebra(9, [3, 5, 0, 5, 0, 0, 0])
visited = fingerprint.Visited()
search.mutationSearchDepthFirst(alg, 6, [], 'x', printOutput = False)   # plain
print(visited.summarise())
```

---

## E-041 — The two separators of the sweep, checked
*2026-09-19* · **both hold; one reproduces the quipu boundary exactly, the other certifies 619 rows at `n = 11`** → F-047, F-048

E-040 ended by naming two criteria as its highest-value unverified output, and
warning that a criterion quoted out of a paper should be assumed to be missing a
hypothesis until it has been run against a known answer. This is that run.

### (a) The `Z`-congruence invariant (math/0610685 Cor. 3.15) → F-047

Profile: the Smith normal form of `g(Φ)` for each irreducible factor `g` of the
Coxeter polynomial.

* **Soundness, on every member of every orbit** rather than a sample: 1430 LNAs in
  48 orbits at `n = 9`, 4862 in 113 at `n = 10`. **Zero orbits split.**
* **Sharpness:** splits 1 of the 11 cospectral orbit-groups at `n = 9`, 3 of 25 at
  `n = 10`.
* **The `n = 9` split was then named by the quipu theorem**, which is the check
  that turns a refinement into a separation: all four orbits on one profile are
  `P^(1,4)_(1,0,1)`, both on the other are `P^(1,2)_(1,1,2)`, no orbit on the
  wrong side. The invariant reproduces the true class boundary on the exact pair
  the Coxeter polynomial cannot see.

Cost: about three minutes for `n = 9`, fifteen for `n = 10`, in sympy. The Smith
normal form is the expensive part and it is per irreducible factor, so this scales
with the factorisation rather than with `n`.

### (b) The periodicity criterion (math/0611201 Thm. 3.4) → F-048

`Φ` periodic and the Euler form indefinite certifies not piecewise hereditary.

* Fires on **0** LNAs at `n = 5` to `9`, **3** at `n = 10`, **638** at `n = 11`.
* At `n = 10` the three are `34504030`, `50505000`, `45050400` and
  `piecewiseHereditary` certifies **none** of them — the rows backlog 20 wants
  `lemma:taupathimpliesnotpwh` for. At `n = 11`, 619 of the 638 are new.
* **Falsification test:** an almost separate LNA is piecewise hereditary by the
  quipu theorem, so the criterion must never fire on one. Over 30648 rows at
  `n = 5` to `11` it fires on 641, **not one almost separate**.

This is the criterion E-040's warning was about, in its usable form. The version
that misfired there was de la Peña's, which needs "not of Dynkin module type"
supplied separately; Ladkani's asks for an indefinite Euler form instead, and
indefiniteness rules the Dynkin and Euclidean cases out by itself. **The same
mathematics, and one statement of it is safe to implement while the other is
not** — which is the concrete lesson, rather than the general caution.

### What is still not done

Both are verified as *criteria*; neither is in the code. F-047's profile has been
checked at two lengths and F-048's at seven, so porting them is now a matter of
writing them into `invariants` and `piecewiseHereditary` rather than of deciding
whether they are true. Backlog 32 is updated.

---

## E-040 — A literature sweep aimed at merging classes of LNAs
*2026-09-19* · **21 summaries; one criterion that classifies the tame half outright, and three merges at `n = 11` no move of ours makes** → F-043, F-044, F-045, F-046

The question put to the literature was narrow on purpose: **what merges two
LNAs?** Not what classifies Nakayama algebras, not what invariants exist — what
would let two of our orbits be joined, or be proved distinct. Anything that could
not be tied to that in a sentence was rejected.

### How it was searched

Three passes, and the third is the one that paid.

1. **Backwards**, through the reference lists of the three papers the project
   already uses (arXiv:2112.08129, 2305.06642, 2310.08346). Seventeen distinct
   references between them — a small enough set to read in full. This is what the
   `literature/README.md` candidate list was built from.
2. **Outwards**, from those into the authors' own corpora: Ladkani's 27 arXiv
   papers, the Happel school, the silting-mutation line.
3. **Forwards**, by citation. Semantic Scholar's graph API on the three papers'
   arXiv ids, plus `export.arxiv.org` title and abstract search on "Nakayama" ×
   "derived equivalence". **This found the two best papers in the sweep, and no
   backward reference list could have**: arXiv:2302.02880 (Ueda) and
   arXiv:2203.15735 (Dong–Lin–Ruan) are both later than everything we cite, and
   Brüstle came in as a reference of arXiv:1910.01494, which itself was only found
   forwards. *Do the forward pass first next time.*

Roughly 60 papers screened on abstracts, 21 read closely enough for a file, 22 of
Ladkani's rejected with a recorded one-line reason each so they are not re-screened.

### What came back, sorted by what it does

**Merges.** `research/literature/` now holds four sources that produce derived
equivalences between LNAs: Ueda's Cor. 1.3 (F-043), the Happel–Seidel symmetry and
its extension in Lenzing–Meltzer–Ruan Prop. 4.1, Dong–Lin–Ruan Prop. 4.5, and
Ladkani's `A(mn, m+1) ≃ kA_m ⊗ kA_n` (0911.5137 Cor. 1.2). Every one of them is
about **radical powers** `kA_n/rad^r` or a tensor of two lines — the thinnest
family of rows we have. Nothing found speaks about an LNA with relations of
several different lengths, which is the open case.

**A classification.** Brüstle's Theorem 1.2, which decides the derived class of
any LNA with non-negative Euler form from three numbers: F-045.

**Separators.** Ladkani's Cartan-matrix-up-to-`Z`-congruence (math/0610685
Cor. 3.13), reported to split cospectral quipu groups the Coxeter polynomial
cannot; and two periodicity criteria (math/0611201 Thm. 3.4; de la Peña,
arXiv:1310.1557) certifying non-piecewise-heredity, one of which is reported to
certify `34504030`, `50505000` and `45050400` at `n = 10`, which our own criteria
miss. **Neither has been re-verified here** — they are the obvious next thing to
check.

**Two lines closed.** `HH*(A) = k` for every LNA (arXiv:2312.14699), so idea 22
is dead; and no extension of the Avella-Alaminos–Geiß invariant to string algebras
exists, so R-008's line stays closed — arXiv:1910.01494's skewed-gentle conclusion
does not reach us, because its hypothesis is "no simple projective module" and
every LNA has one (`P_n = S_n`, from the sink).

### What was checked against the code, and what it cost

Everything below was run here rather than taken on trust, which is the only reason
the findings above are findings.

| check | result |
|---|---|
| Ueda Cor. 1.3, 15 instances at `n ≤ 16`, against `movesJoin` | all 15 joined — but by the **double mutation**, not the rule table (F-043) |
| Happel–Seidel Table 1, 11 star types, against `treeCoxeterKey` | 11/11, and an invented 12th row correctly fails (F-044) |
| Happel–Seidel Table 1, 12 sheaf types, against `canonicalWeightType` | 12/12 (F-044) |
| Brüstle Thm. 1.2, reimplemented, at `n = 9, 10, 11` | reproduces F-011's tame partition; 1 new merge per length; 0 contradictions (F-045) |
| LMR Prop. 4.1, 21 instances, against `movesJoin` | 15 joined, 6 not — all six with both orbits **exhausted**, not capped, at `n = 11, 13, 15` (F-046) |
| the three `n = 11` merges, against `derivedOrbits(11)` | four of our orbits merge into two (F-046) |
| de la Peña's periodicity criterion, naive reading, at `n = 9` | **certifies 273 LNAs including `A_9` itself** — the Dynkin exclusion is the whole criterion, caveat recorded in the file |

That last row is the one to remember: a criterion quoted out of a proof, applied
without its exclusions, certified the hereditary line as non-piecewise-hereditary.
Every summary in this sweep that states an implementable criterion should be
assumed to be missing a hypothesis until it has been run against something whose
answer we already know.

### What was not done

* The two separators above are unverified here (Ladkani's congruence invariant and
  the two periodicity criteria). They are the highest-value follow-up, because a
  separator is what F-010 has wanted since R-008.
* `n = 12` and up for Brüstle: the Smith normal form is cheap but `derivedOrbits`
  is not, so there is nothing to compare against past `n = 11`.
* Three papers are summarised **from secondary sources** — Happel–Seidel, Rickard,
  Assem–Happel — because they are journal-only and pre-arXiv. Happel–Seidel's has
  been checked (F-044); the other two have not.
* The cellular-automaton sweep of `literature/README.md` (H-009) was not touched;
  this sweep was about merging, not about the move rules as a rewriting system.

---

## E-039 — Ueda's radical-power equivalence, checked against the move table
*2026-09-19* · **15 instances at `n <= 16`, all 15 joined by the moves alone** → F-043

Prompted by the literature sweep: arXiv:2302.02880 Corollary 1.3 gives a triangle
equivalence `per N(n, l+1) -> per N(n, l)` for `n = p(p+1)q + p(p-1)r`,
`l = (p+1)q + pr`. Enumerating `p` in 2..5, `q` in 1..4, `r` in 0..4 plus the
half-integer case `p = 2`, and keeping `4 <= n <= 16`, gives 15 distinct `(n, l)`.

For each: build both LNAs, compare `coxeterKey`, and ask `freeMoves.movesJoin`
whether the moves join them, with the meeting row then checked reachable from
**both** ends by `orbitOf(..., target = meet)` -- a `movesJoin` result alone is
one walk from each side and worth confirming when the conclusion is that a paper
adds nothing.

* Coxeter keys agree in all 15, as they must.
* All 15 joined with `free = True` (derived equivalence, the same relation Ueda
  proves).
* All 15 joined again with **`free = False`**, so they are joined by *mutation*
  moves alone -- a stronger statement than the paper's, for these instances.
* Both-ends check passed in all 15.

Cost: under a minute for the whole sweep, `limit = 20000` rows per walk, never
approached.

**What was not done.** `n > 16` was not tried: the parameter grid thins out fast
(the next instances are at `n = 18` and `n = 20`) and the orbits grow, and the
point was made. The *other* corollary of the paper -- an equivalence from every
`N(n,l)` to an algebra of global dimension at most 2 -- was not checked, because
the target is not an LNA and the pipeline has nothing to compare it against.

Reproduce: `ueda2.py` in the merge session's scratchpad.

---

## E-038 — The branch's free-move table, re-measured under the guarded search
*2026-09-19* · **both rows reproduce exactly, and the merge changes no number** → E-037, F-041

E-037 was measured on `claude/quiver-nakayama-investigation` before F-038 and
F-039 landed on `main`, so it was measured with a search that took steps which
are not derived equivalences and with a Cartan matrix that counted a parallel
pair of arrows as one path. Merging the branch puts its measurements on top of
three changes at once:

* `search.mutationSearchDepthFirst` now refuses a step whose Coxeter key moves
  (`coxeterGuard`, on by default), so the walk is strictly shorter-reaching;
* `invariants.integerCartanMatrix` counts **arrow** paths and takes the rank of
  the ideal where the cheap count is not provably right, so the guard is reading
  a different number than it would have;
* `search.quiverKey` returns `None` at a quiver with parallel arrows, so those
  nodes are no longer meeting points (below).

Any of the three could have moved the table. Re-run of E-037 section 2, same
depth, `alsoDual = True`:

| | `n = 8` | `n = 9` |
|---|---|---|
| single deletions | 572 | 2002 |
| joined by the known moves | 562 | 1937 |
| joined by meeting at depth 3 | 10 | 52 |
| left | **0** | **13** |

**Identical to E-037 in every cell**, and the 13 left at `n = 9` are the same 13
pairs: `2302230/2300230`, `2302300/2300300`, `2302330/2300330`,
`2302302/2300302`, `2302030/2300030`, `2303020/2303000`, `3302302/3300302`,
`3022302/3020302`, `3002302/3000302`, `0230300/0030300`, `0230302/0030302`,
`0302302/0300302`, `0303020/0303000`. So F-041 stands as written: the guard
removed nothing the branch had counted, which is the same answer F-039 got for
the sweeps it re-ran — at these lengths the unguarded search was not reaching
anything the guarded one cannot.

`n = 10` was **not** re-run: 36s at `n = 8` and 296s at `n = 9` extrapolates past
what this was worth, and the two rows that were re-run are the two the finding's
claim rests on. The `n = 10` row of E-037 is therefore still a pre-guard
measurement and is marked as such there.

**One soundness fix the merge did need.** `quiverKey` keyed a node on
`pathAlg.rels`, and F-039 established that `rels` is a *lossy projection* once a
quiver has parallel arrows -- a path named by its vertices no longer says which
of two arrows it runs along, so two different algebras write down the same
`rels`. Before the merge this could not bite, because `procedure.isMutable`
refused to mutate a quiver with a parallel pair at all and the search never
descended past one. R-013 removed that refusal. Two searches could then have
"met" at a key that is not an algebra, which would be a fabricated proof of
derived equivalence -- the precise failure F-038 is about, reintroduced by a
different door. `arrowRels` cannot serve as the key either: the procedure hands
out arrow keys in the order it builds them, so the same algebra reached by two
routes carries different ones and there is no canonical form to compare across
routes (F-039 says this in as many words). So a parallel-arrow node is now not a
meeting point in either direction; the search still walks through it. No meeting
recorded at `n = 8` or `n = 9` was at such a node, which is why the table above
did not move.

Reproduce: the script is `reverify.py` in the merge session's scratchpad, and
the two rows it prints are `tests/test_other_families.py`'s
`test_every_single_two_arrow_deletion_at_length_eight_is_a_mutation` generalised
to a second length. Full suite after the merge: **3015 passed, 1 xfailed**.

---

## E-037 — What the short quivers can show, and two questions asked properly
*2026-09-17* · **the evidence base is narrower than it looks; the free move settles at `n = 8`; the barricade does not hold and bounded overlap is not the pattern** → F-040, F-041, F-042, H-017

*Renumbered at merge from `E-032`, which was taken on `main` first by an unrelated entry while this branch was open. Session logs and commit messages from the branch use the old identifier.*

Prompted from outside the code, and the prompt was the useful part: everything
known about the classes the quipu theorem misses comes from lengths 9 to 11,
where an LNA has very little room -- at `n = 11` no vertex is more than five from
an end -- so patterns found there may be patterns of the small cases rather than
of the problem.

### 1. How narrow the evidence is (F-040)

Counting heavy clusters -- runs of relations linked by overlaps of two arrows or
more, which is what puts an LNA outside the theorem:

| | `n = 9` | `n = 10` | `n = 11` |
|---|---|---|---|
| LNAs outside a quipu class | 9 | 262 | 2647 |
| of those, with two heavy clusters | 0 | 2 | 42 |
| with two clusters and a **free arrow between them** | **0** | **0** | **0** |
| with the cluster touching an end | 8 | 221 | 2119 |

So every unclassifiable LNA at every length worked on here is **one overlapping
cluster, usually against an end**. Two clusters with a free arrow between them
first fit at `n = 10`, and every LNA that has them is in a quipu class. The
*barricade* -- two clusters walling in a two-arrow relation -- needs 12 arrows and
first fits at `n = 13`.

### 2. The free move, asked one relation at a time (F-041)

*Measured before the Coxeter guard of F-038 existed. Re-run after the merge at
`n = 8` and `n = 9`, both rows identical; the `n = 10` row below has not been
re-measured. E-038.*

H-012 asks whether deleting a two-arrow relation is a mutation equivalence. Two
changes to how it is asked:

* **one deletion at a time.** The whole strip is a composition of single
  deletions, and each of those is a much shorter journey.
* **meeting in the middle.** `search.meetingPoints`: two searches that reach the
  same quiver have joined their algebras, at twice the depth for the same cost.
  Nothing in the pipeline did this -- `resolveMergeCandidates` keeps only the
  lines a search lands on and throws the rest of the tree away.

| | `n = 8` | `n = 9` | `n = 10` |
|---|---|---|---|
| single deletions | 572 | 2002 | 7072 |
| joined by the known moves | 562 | 1937 | 6768 |
| joined by meeting at depth 3 | 10 | 52 | 206 |
| left | **0** | 13 | 98 |

`n = 8` is settled outright. Of the 13 left at `n = 9`, eleven have both sides
almost separate with the same quipu, so the theorem already calls them one
derived class; only `2302330 / 2300330` and `3302302 / 3300302` are outside
everything, and they do not meet at depth 4 either. A depth-4 pass over the 98
left at `n = 10` joins none of them, at 1060s: depth 4 buys nothing over depth 3
here, on either length, which says the remaining pairs are either much further
apart than 4 + 4 mutations or not joined at all.

### 3. The barricade, built on purpose

The shape H-012's doubts are about, built at the lengths where it first fits: a
heavy cluster, a free arrow, a two-arrow relation, a free arrow, another heavy
cluster. `freeMoves.orbitOf` walks the moves out of one row, which is what makes
a length-13-to-16 question affordable at all -- `derivedOrbits` would partition
208012 rows to answer it about one.

**The moves get it out, in all 95 shapes tried.** Every pairing of six heavy
clusters on the left and right, at gaps of one to three arrows on each side and
lengths 13 to 16: `33000200330` at `n = 13` reaches its strip in an orbit of a
thousand rows, and so does every wider version. **The barricade does not trap a
two-arrow relation.**

**Two false starts, and both were the measurement rather than the mathematics.**

* Run first with the rule table left out -- F-032 having found that the double
  mutation subsumes it at `n <= 10` -- the same barricades came out **not**
  joined, with orbits of 53 to 89 rows against 1073 with the table. The table is
  not subsumed at these lengths.
* With the table in, 49 of the 95 still came out not joined -- and every one of
  them had an orbit of *exactly the 20000-row cap*, so what was measured was the
  budget. Walking from both ends instead (`freeMoves.movesJoin`, the move-level
  twin of `search.meetingPoints`) joins **all 49 in 45 seconds**, most of them
  instantly. A one-way walk that runs out of budget says nothing whatever, and
  the orbit size is the tell: if it equals the cap, there is no result.

So the shape H-012's doubts are about does not hold the relation in, at the
lengths where it first exists. What is untested is the user's fuller version --
four clusters at `n = 30` to `50` -- and the shapes here are the minimal ones.

### 4. Bounded overlap is not the pattern (F-042)

The generalisation asked for was a weaker version of "almost separate": overlaps
of at most two, or at most so many overlaps above one. Crossing those coordinates
against membership of a quipu class kills both. The cell `(max overlap 2, exactly
one of them)` -- the smallest possible step past the theorem -- already holds
outsiders at `n = 9` (`3033030`), `n = 10` (12 of them) and `n = 11` (84).

Chasing what those outsiders have in common instead: the ones with the most free
arrows are the same two little cores at every length, `45` and `504`, and sliding
`45` along the quiver gives a clean law -- inside when it touches the source or
comes within one arrow of the sink, outside everywhere in between, with the band
growing by one place per vertex. `0450000` at `n = 9` is inside, `04500000` at
`n = 10` is not. The condition that decides it is **where the cluster sits**, and
no condition on the relations by themselves can see the difference. F-042.

### 5. What the quipu members look like (H-017)

Walking every LNA outside a quipu class at `n = 9` to depth 5 and grouping the
quipu algebras reached by (cords, relations): the pairs that occur are (1,2)
(1,3) (1,4) (1,5) (2,3) (2,4) (2,5) (2,6) (3,4) (3,5), and **never one with
relations at most cords**. Reading `relations - cords` as a defect, a quipu class
is defect `<= 0` and these are all `>= 1`; the minimum over a class is 2 for
`3033030` and 1 for the eight-member class. That is the beginning of a normal
form, and it makes a sharp prediction about the polynomial-only candidates of
F-034 that sit below the diagonal. H-017.

---

## E-036 — What conditional deeper probing costs, and what it reaches at n = 9
*2026-09-18* · **the parallel-arrow region is reached by exactly one of the nine leftover members at n = 9; giving its branches two more mutations costs 2-4% of the run, and nine mutations into the region it still reaches nothing but its own dual** → H-016

A depth-bounded search gives every branch the same budget. `search.DeeperWhen`
gives extra mutations to the branches that reach a quiver meeting a condition,
with a per-branch budget so the walk still terminates; `merges.py --deeper-on`
asks for it from a command line. This measures the two things that decide
whether it is worth using: what it costs, and whether the extra depth it buys
reaches anything.

### 1. Cost and gain at n = 9

Every member of every leftover orbit at `n = 9` (9 members: the eight of
`3345000` and `3033030`), searched from itself and from its relation dual, plain
against `parallel-arrows:2` — two extra mutations for a branch that reaches a
quiver with a parallel pair, once per branch. Four cores; the seconds are the
sum over the nine members, not wall clock.

| depth | firings | grants | LNAs gained | LNAs lost | plain | probed |
|---|---|---|---|---|---|---|
| 4 | 58 | 4 | 0 | 0 | 89s | 91s |
| 5 | 600 | 32 | 0 | 0 | 418s | 435s |

**Every firing is `3033030`.** The other eight members never reach a quiver with
parallel arrows at all at these depths, so the probe is free for them and the
whole cost is one member's: 19s to 36s at depth 5, where raising the depth for
everyone from 5 to 7 would cost about twenty-five times the run. That ratio is
the case for the mechanism, and it is the one thing here that is not about
parallel arrows in particular.

**Nothing is gained, to depth 5.** The same LNAs are reached either way. This is
a weaker negative than E-035's: `3033030` is alone in its orbit *and* alone in
its Coxeter polynomial group, so the only outcome that would show here is it
reaching a **seeded** LNA — a leftover turning out to be in a quipu class after
all — and at depth 5 it reaches one LNA, its own dual.

### 2. One pass against two

Whether granting depth inside one search differs from recording the interesting
quivers and searching from those afterwards. The second contains the first, and
the containment can only be strict for a branch that leaves the region and comes
back, since a firing *below* the one that bought the depth is already inside the
subtree the grant paid for, at exactly the remaining depth a re-search would
give it.

| start | depth | extra | plain | one pass | two passes |
|---|---|---|---|---|---|
| `3030` | 4 | 3 | 50 | 105 | 105 |
| `3030` | 5 | 2 | 118 | 195 | 195 |
| `30330` | 4 | 2 | 67 | 86 | 86 |
| `30330` | 5 | 2 | 158 | 237 | 237 |

Nodes reached, counted by the quiver with its arrow names and its relations.
**Equal in every case**: no branch at these sizes leaves the region and returns
within the depth searched. So the one-go run is not the weaker of the two in
practice, which is the practical answer — the two-pass route's advantage is that
the count of firings is visible before the second round is paid for, not that it
reaches more.

### 3. The one member that reaches the region, pushed to depth 9 inside it

`3033030` is the only member of either leftover orbit whose walk ever reaches a
quiver with parallel arrows, so it is the whole of the `n = 9` test and it is
cheap. From it and its relation dual, plain against `parallel-arrows:2`:

| depth | firings | grants | deepest firing | reaches | seconds |
|---|---|---|---|---|---|
| 6 | — | — | — | itself | 98s |
| 6 + 2 | 4322 | 154 | 8 | itself | 221s |
| 7 | — | — | — | itself | 350s |
| 7 + 2 | 29122 | 826 | 9 | itself | 1200s |

The four ran together on four cores and the last two shared the machine with a
test run, so the seconds are an upper bound and the ratio between them is the
part worth reading.

"Deepest firing" is the length of the longest mutation path at which the
condition still held, so the last row walked **nine** mutations into the region.
It reaches nothing but its own relation dual, which is what depth 5 already
reached.

This is the `n = 9` half of what H-016 asks for, and past the depth it asks for.
It does not settle H-016: `3033030` is alone in its Coxeter polynomial group, so
the only thing it *could* show is a leftover turning out to be in a quipu class,
and one member at one length is not the hypothesis. `n = 10` and `n = 11`, where
H-013's leftover orbits sit, have not been looked at this way. But it is a
negative at a depth and a length E-035 could not reach, and it cost 20 minutes
on one core rather than the run over every member that a uniform depth 9 would
have been.

### 4. Over a whole classification

`classify.py --deeper-on` gives one probe to all three searching steps. At
`n = 8`, the shortest length whose classification needs a search at all:

| | classes | rows | firings | grants | seconds |
|---|---|---|---|---|---|
| plain | 11 | 429 | — | — | 38s |
| `parallel-arrows:2` | 11 | 429 | 44 | 4 | 39s |

**The answer does not move**, which is the check that matters: extra depth may
place a row the depth could not reach, and may never place one differently. The
published table of arXiv:2305.06642 is what both are checked against, as
`tests/test_classify_end_to_end.py` does. The condition does fire here, four
times buying depth, so this is the probe running over a real classification
rather than a no-op.

### Reproducing

Section 4 is `classify.py 8 --quiet` with and without
`--deeper-on parallel-arrows:2`.

Section 1 is `merges.py 9 --depths 4 5 --all-groups` run twice, once with
`--deeper-on parallel-arrows:2` and once without, comparing `lnasReached` and
`deeperFirings` per checkpoint record. The two runs can share one checkpoint
file: it records the condition, and a plain record does not count as covering a
probed search or the other way round.

Section 2 is the computation of
`tests/test_deeper_probing.py::test_recording_and_searching_again_contains_deepening_in_one_pass`
at the four settings in the table.

Section 3 is `merges.searchFrom((9, (3, 0, 3, 3, 0, 3, 0), depth, spec))` for
each row.

---

## E-035 — Lifting the parallel-arrow restriction, and re-measuring E-033
*2026-09-18* · **every wrong-key node at n = 6 and n = 7 to depth 5 was a mis-count; the new region reaches nothing new at these sizes** → F-039, R-013, H-016

E-033 walked the search tree with the Coxeter key in hand and split the nodes
where it had moved into "parallel arrows, harmless" and "clean, the real fault".
This is the same sweep after the relations were moved onto paths that name their
arrows, so a parallel pair can be stated, counted and mutated at.

### What was changed

* `arrowPaths` — an arrow is `(tail, head, key)`, a path is a tuple of arrows.
* `procedure` — steps 1 to 7 per arrow; step 5 divides by the arrow, step 7 reads
  each candidate's first arrow back as its relation and its tail back into the
  old quiver, step 6 and the carried-past relations name the composite they use.
* `procedure.isMutable` — no longer refuses a quiver with a parallel pair.
* `invariants` — the Cartan matrix counts arrow paths, exactly and cheaply.
* `search` — the illegal-relation check is over arrow relations.
* `arrowPaths.homDimensionByClosure` — the cheap count closes the commutativity
  relations to a fixed point, where `paths.numberOfPathsUpToRels` applied each of
  a subset once.
* `invariants.integerCartanMatrix` — and the key is **exact** wherever the cheap
  count is not provably right, which is wherever the ideal is not monomial. No
  closure makes the cheap count see a relation of three or more paths, and step 4
  produces one at every vertex with three arrows out.

### 1. The sweep, before and after

`mutationSearchDepthFirst` from every LNA of the length, `coxeterGuard = False`
so the corrupt region is visible, `coxeterKey` compared at every node against the
start's. One core.

| | | nodes | wrong key | parallel | clean | lines | seconds |
|---|---|---|---|---|---|---|---|
| `n = 6`, depth 5 | before | 14,693 | 4 | 4 | 0 | 1,789 | 28 |
| | after | 14,701 | **0** | 0 | 0 | 1,789 | 27 |
| `n = 7`, depth 5 | before | 94,446 | 79 | 75 | 4 | 7,175 | 254 |
| | after | 94,498 | **0** | 0 | 0 | 7,175 | 268 |
| `n = 7`, depth 6 | before | 336,760 | 339 | 277 | 62 | 20,683 | 866 |
| | after | 337,360 | **10** | 0 | 10 | 20,683 | 949 |

The node count rises by 8, 52 and 600: a parallel-arrow node is no longer
terminal. Lines are unchanged, exactly, at all three. The cost is 10% at depth 6,
which is the exact Cartan matrix against the cheap one, and it is not optional.

Three mechanisms, all three mis-counts:

* the *parallel* nodes -- 4, 75 and 277 of them -- had the key computed over
  vertex sequences, so a parallel pair contributed 1 to the Cartan matrix where
  it contributes 2;
* some *clean* nodes are the incomplete closure. The four at `n = 7` depth 5 are
  `40030` by `[1, 4, 2, 5, 2]`, by `[1, 2, 4, 5, 2]`, by `[1, 1, 5, 4]` and one
  more;
* the rest are relations of three or more paths, which the cheap count has no
  reading of and ignores, e.g. `34400` by `[1, 3, 4, 2, 2, 1]`.

The ten that survive at depth 6 are **one bad step and its descendants**: all ten
are reached from `30330` along the `[4, 1, 3, 1, ...]` family. Replaying `[1, 4, 2, 5, 2]` in both engines gives the **same quiver and
  the same `rels`**, and the key holds in one and moves in the other, which is
  what says it is the measurement.

### 2. What still moves the key

From the relation dual of `33030` at `n = 7` by `[4, 1, 3, 1, 3, 3]`, F-038's own
smallest clean case, the key moves at step 6 under the arrow model too, and there
the cheap and the exact Cartan matrices agree -- so it is the algebra. R-012
stands and the guard stays.

### 3. The classification is unchanged

`classifyLength` on both engines, same machine:

| | classes | rows | before | after |
|---|---|---|---|---|
| `n = 6` | 4 | 42 | 0.3 s | 0.3 s |
| `n = 7` | 6 | 132 | 1.2 s | 1.3 s |
| `n = 8` | 11 | 429 | 29.9 s | 33.4 s |

Class for class, size for size, identical at all three, and about 10% dearer --
the exact Cartan matrix where the ideal is not monomial, against the cheap count
everywhere.

### 4. What the region beyond a parallel pair looks like

From `3030` at `n = 6`, mutated at `[1, 3, 4, 1, 4]`, the quiver has two arrows
`1 -> 6` and one relation between the two parallel paths `5 -> 1 -> 6`. The old
gate refused all six vertices. The paper's criterion admits three, every one
keeps the Coxeter key, and mutating at 3 comes back out to a quiver with **no**
parallel arrows and a genuine commutativity relation
`5 -> 1 -> 3 -> 6 = 5 -> 1 -> 6` -- a quiver the search could not reach at any
depth before.

### 5. Does the new region reach anything? Not at these sizes

The point of walking through a parallel-arrow node is what lies beyond it, so:
for every LNA of the length **and its relation dual**, the set of LNAs the
guarded search reaches, with the gate allowing parallel arrows and with it
refusing them, at the same depth.

| | starts | lines reached, allowing | refusing | starts gaining | losing |
|---|---|---|---|---|---|
| `n = 6`, depth 6 | 74 | 801 | 801 | 0 | 0 |
| `n = 7`, depth 6 | 244 | 3,390 | 3,390 | 0 | 0 |

(The counts are the sum over starts of how many LNAs that start reaches, so a
line reached from two starts counts twice; what matters is that no start gained
or lost one.)

**So the region is reachable and walkable and yields nothing new here.** Two
reasons not to read that as "it never will". The comparison holds the *depth*
fixed, and entering the region and returning from it costs steps, so at depth 6
the part of it that can come back to a line at all is thin -- the one walk
looked at by hand, `3030` by `[1, 3, 4, 1, 4]` then 3, takes six mutations to
get back out to a quiver with no parallel pair, and that quiver is not a line.
And `n <= 7` is fully covered by the move rules with no search at all (F-021),
so there is nothing left for a search to find at these lengths whatever it walks
through. H-016.

Reproduce: `tests/test_parallel_arrows.py`; the sweep is a `visitor` on
`search.mutationSearchDepthFirst` comparing `search._coxeterKeyOrNone` at each
node, and the before column is the same script against `git show 78328e7`. The
reachability comparison patches `procedure.isMutable` to pass
`allowParallelArrows = False` and compares `search.linesReachedFrom` either way;
it is 84 minutes at `n = 7` depth 6 on one core.

---

## E-034 — Can the guarded search be fooled where the polynomial is known to fail?
*2026-09-18* · **not at depth 6, at the smallest collision** → H-015

F-038's guard refuses any step that moves the Coxeter polynomial. That is
necessary for a derived equivalence; H-015 asks whether it is sufficient. The
place to look is where the polynomial is known to be blind, and F-010 says
exactly where that is: cospectral quipus, the smallest collision at order 9,

    P^(1,4)_(1,0,1) = A_{9,(1,3)}^{(3,6)}   class 3060000
    P^(1,2)_(1,1,2) = A_{9,(1,4)}^{(3,4)}   class 3004000

Different trees, so **not derived equivalent**, yet one Coxeter polynomial —
confirmed here, both `(1, 1, -1, -3, -4, -4, -3, -1, 1, 1)`. If a guarded search
out of one reaches anything the other reaches, then a step the guard admits is
not a derived equivalence and H-015 falls.

Searching from each, from the member and from the relation dual:

| depth | guard | LNAs from `3060000` | from `3004000` | **shared** | time |
|---|---|---|---|---|---|
| 5 | on | 2 | 7 | **0** | 85 s |
| 5 | off | 2 | 7 | **0** | 45 s |
| 6 | on | 2 | 8 | **0** | 395 s |
| 6 | off | 2 | 8 | **0** | 190 s |

**Nothing shared, either way.** `3060000` is remarkably rigid — it reaches only
`6000030`, its own dual, at either depth — while `3004000` moves to seven or
eight. The guarded and unguarded searches return *identical* sets here, which is
consistent with F-038: at this depth the corrupt region is entered but does not
come back round to a line.

**What this is worth, and what it is not.** It is the sharpest single test
available and H-015 survives it. But a negative is a depth bound, not a proof,
and this probes **one** collision: order 10 has two collision groups and order 11
has four (F-010), none of them tried. It also cannot detect a guard-passing step
between two algebras that are cospectral for some *other* reason than being
quipus.

Reproduce: the pair is `python classify.py 9 --collisions`; the search is
`search.mutationSearchDepthFirst` from each with `coxeterGuard` both ways.

---

## E-033 — Is the mutation search sound?
*2026-09-18* · **no, and the fix costs 1.85× and changes no answer** → F-038, R-012, F-037

Prompted by the two ALARMs of E-032. The question is narrow and had never been
asked directly: **does `mutationSearchDepthFirst` stay inside one derived
equivalence class?** Every step is meant to be a tilting mutation, so the Coxeter
polynomial must be constant over the whole search tree — at every quiver reached,
not only at the lines it records.

### 1. Walk the tree with the invariant in hand

A `visitor` that computes `coxeterKey` at each node and compares it to the start,
aborting at the first mismatch. **It fails at `n = 6` in 7 seconds**: from
`3030`, the path `[1, 3, 4, 1, 4]` reaches a quiver whose key has moved from
`(1,1,-1,-2,-1,1,1)` to `(1,1,0,-1,0,1,1)`. Replaying it one mutation at a time
shows the last step producing a **parallel arrow**, `1 -> 6` twice.

### 2. Classify every wrong node

| | nodes | wrong key | parallel | cyclic | **clean** | lines | lines wrong |
|---|---|---|---|---|---|---|---|
| `n = 6`, depth 5 | 25,398 | 4 | 4 | 0 | **0** | 3,263 | 0 |
| `n = 7`, depth 6 | 609,474 | 604 | 507 | 0 | **97** | 37,911 | 0 |
| `n = 8`, depth 5 | 1,093,976 | 1,204 | 1,030 | 0 | **174** | 55,175 | 0 |

The parallel-arrow ones are a limitation of the model and are **harmless**:
`procedure.isMutable` refuses every vertex of a quiver that has parallel arrows
anywhere, so the node is terminal, and such a quiver can never be mistaken for a
line. Not one oriented cycle appeared at all.

**The clean ones are the fault** — acyclic, no parallel arrows, key moved,
nothing to stop the search descending. Smallest: `n = 7`, from the relation dual
of `33030`, path `[4, 1, 3, 1, 3, 3]`, where step 6 loses a commutativity
relation outright. Written up as F-038, retracted as R-012.

**No answer was wrong at these sizes.** All 96,349 lines collected carried the
starting key. The corrupt region exists but had not reached an answer.

### 3. The guard

`search.mutationSearchDepthFirst(..., coxeterGuard = True)`, now the default,
refuses a step whose key differs from the start's — R-005's third requirement,
applied per step rather than per rule.

| | without the guard | with it |
|---|---|---|
| wrong-key nodes, `n = 6` depth 5 | 4 | **0** |
| lines reached, `n = 6` depth 5 | 3,263 | 3,263 |
| core-seconds, `n = 6` depth 5 | 75 | 138 (**1.85×**) |
| core-seconds, `n = 7` depth 5 | 759 | 1,403 (**1.85×**) |
| LNAs losing a reached line, `n = 6`, `n = 7` | — | **0** |
| LNAs gaining one | — | **0** |

Comparing the *sets* of lines reached, from every LNA and from its dual, at
`n = 6` and `n = 7` to depth 5: not one line lost, not one gained. **Existing
results at these sizes stand unchanged.**

A cyclic quiver has no unimodular Cartan matrix, so `coxeterKey` raises there;
`_coxeterKeyOrNone` returns None and such a step is let through, because the
search does not descend from a cycle anyway. Without that, starting a search at a
cyclic algebra crashed — `tests/test_cycles.py` caught it.

### 4. The links of E-032 re-derived and replayed

A link is a positive claim, so each was found again and checked move by move:
every step admissible, no illegal relation, **no parallel arrow or oriented
cycle**, and the key held. All five pass.

| link | paths into the target | checked | steps |
|---|---|---|---|
| `n = 10` `34504030 -> 50505000` | 19 of 20 lines collected | 5 | 7 |
| `n = 10` `05040330 -> 33460000` | 9 | 5 | 7 |
| `n = 11` `030233030 -> 300330400` | 3 | 3 | 3–4 |
| `n = 11` `346004030 -> 060040400` | 4 | 4 | 3–4 |
| `n = 11` `302340030 -> 300403030` | 3 | 3 | 3–4 |

The first is F-037.

### 5. The ALARM itself, reproduced

The depth-8 search from `03033030` that raised it takes 8257 s on one core, so it
was split by first mutation and the branches run in parallel. It reproduces
exactly, and it is the clean mechanism of part 2:

* From the **member**, ten branches, **not one line reached at all** — everything
  the overnight run recorded from this side is the start itself.
* From the **relation dual**, branch `9` reaches one line correctly, and branch
  `4` reaches **four lines, all four wrong**, all of them `30233330`, which is a
  member of orbit `00330400`. No parallel arrows, no oriented cycle.

Bisecting `[4, 6, 4, 6, 9, 4, 4, 6]`, every step admissible and acyclic
throughout:

    steps 1-6   key (1, 1, -2, -3, 1, 4, 1, -3, -2, 1, 1)   held
    step 7 at vertex 4   -> (1, 1, -1, -2, -1, 0, -1, -2, -1, 1, 1)   moved
    step 8 at vertex 6   lands on the line 30233330, outside the class

So the ALARM was neither a broken orbit (H-013's guess) nor the alarm test's own
`polyOf` defect (E-032, part 4): it is one admissible-but-not-derived-equivalent
mutation at step 7 of eight, and the guard refuses it.

Re-running that one branch both ways settles it:

    coxeterGuard = False   4 lines, 4 WRONG, 534s   ['30233330']
    coxeterGuard = True    0 lines, 0 WRONG, 941s   []

1.76× here, and the false answer is gone.

**This is why the small-`n` sweeps found nothing.** The corruption needs a deep
enough tree to come back round to a line: seven clean steps, one bad one, then
one more. At `n ≤ 8` and depth ≤ 6 the corrupt region is reached but never
returns to a line, which is exactly what part 2's "lines wrong: 0" column says.

### 6. Which existing results the fix disturbs: none found

The guard only ever *refuses* a step, so it can only shrink what a search
reaches. Every **negative** in the record — "reaches nothing seeded", "stayed
apart to depth 8" — is therefore untouched. Every **positive** needed checking,
and the shallow ones are where most of them live: F-034 and H-014 rest on
depth-3 walks, F-036 and E-031 on depth-4 searches.

Comparing the sets of lines reached with the guard and without, from each LNA and
from its dual:

| | sampled | LNAs where the guard changes what is reached |
|---|---|---|
| `n = 10`, depth 3 | 300 of 4862 | **0** |
| `n = 10`, depth 4 | 150 of 4862 | **0** |
| `n = 9`, depth 4 | 300 of 1430 | **0** |
| `n = 6`, depth 5 | all 42 | **0** |
| `n = 7`, depth 5 | all 132 | **0** |

Together with parts 2 and 4 this says the damage was confined to depth 7 and
beyond, and the only wrong answer the record ever contained is the ALARM itself,
which was already excluded from E-032's unions by the alarm test. **E-032's
conclusions stand unchanged: 43–46 classes at `n = 10`, 84–115 at `n = 11`.**

### 7. The `break` that meant `continue`

`search.py` abandoned every remaining vertex at a node as soon as one mutation
there produced an illegal relation, rather than just that vertex — and the loop
runs over `reversed(vertices)`, so it lost every lower-numbered one. Fixed.
**It never fired in E-032**: `isIllegalRelation` prints when it triggers and all
three overnight logs contain zero such lines over 121 core-hours.

---

## E-032 — The overnight run: H-013 at `n = 10` and `n = 11`
*2026-09-18* · **`n = 10` settled to 43–46 classes; one prediction wrong; and an ALARM that indicts the search rather than the orbits** → E-033

`python overnight.py`, started 2026-09-17 15:46, budget 9 h, on 16 cores under
WSL. Three jobs, and all three reached a definite state:

| job | command | outcome |
|---|---|---|
| `merges10` | `merges.py 10 --depths 5 6 7 8 --jobs 7` | **finished**, exit 0, 00:38 |
| `merges11` | `merges.py 11 --depths 4 5 6 --jobs 7` | stopped on its budget, exit 2, 00:51 |
| `classify10` | `classify.py 10 --resume` | terminated at the deadline, 01:07 |

279 searches at `n = 10` (57.7 core-hours) and 3133 at `n = 11` (63.3), by depth
`{5: 122, 6: 122, 7: 26, 8: 9}` and `{4: 2415, 5: 718}`. Depth 7 and 8 are thin
because `settled(poly)` stops searching a group once its orbits have merged. The
slowest single search was 9691 s, at depth 8.

### 1. `n = 10` is finished, and H-013 was right twice and wrong once

| polynomial group | orbits | H-013 predicted | what happened |
|---|---|---|---|
| `(λ-1)²(λ+1)²(λ²+1)(λ⁴+λ³+λ²+λ+1)`, `C(2,4,5)`'s | 2 (69, 42) | merge, depth ≤ 7 | **merged at depth 6**, 9 searches found it |
| `(λ-1)²(λ+1)²(λ²+λ+1)(λ⁴-λ²+1)` | 4 (4, 2, 2, 1) | at least two stay apart | **all four stayed apart to depth 8** |
| `(λ+1)²(λ²-λ+1)(λ⁶-λ³+1)` = `T¹⁰+T⁹+T+1` | 2 (1, 1) | **no link to depth 8** | **merged at depth 7**, from both sides |

So 12 orbits fall to **at most 10** non-quipu classes, and the derived classes at
`n = 10` number between `36 + 7 = 43` and `36 + 10 = 46`, where H-013 said 43–48.

The wrong prediction is the one worth keeping. `34504030` and `50505000` are the
two the Coxeter polynomial can never separate — the pair `remark:Coxeter` of
arXiv:2310.08346 makes its point with, and the ones H-013 called singletons under
every move known. A depth-7 search links them in both directions, in 553 s from
`34504030` and 1951 s from `50505000`. **A link is a positive claim and the
search is now known to be unsound in some places, so this one is checked
separately in E-033 rather than believed here.**

### 2. `n = 11` got through depth 4 and most of depth 5

54 orbits in 20 polynomial groups, down to **at most 51** non-quipu classes, so
between `64 + 20 = 84` and `64 + 51 = 115`. Three merges, all first seen at depth
4 and all confirmed from both sides:

    030033030 <-> 300330400     28 searches
    040344030 <-> 060040400     12
    300340030 <-> 300403030     12

No group beyond those merged at depth 5, and the two largest groups — 9 orbits
each — stayed entirely apart. Rerunning the same command resumes from the
checkpoint.

### 3. `classify10` got nowhere worth keeping

1010 of 4862 rows attempted in 9 hours, 31 resolved and 5 still unresolved at
depth 6. It is an independent route to the `n = 10` count and it is far too slow
to be one; `merges.py` answered the same question in a fraction of the time
because it only searches what the moves leave over. **Do not schedule
`classify.py 10` again as a whole-range run.**

### 4. The ALARM

Two depth-8 searches, from `03033030` and `30330300`, reported reaching orbit
`00330400`, which carries a different Coxeter polynomial. H-013 said an alarm
"would refute F-032's orbits". It does not: both orbits are internally
consistent, each carrying exactly one polynomial over all its members (4 and 142
of them). What the alarm indicts is the **search**. That is E-033.

Two things about the alarm test itself, found while reading it:

* `polyOf` is built from the leftover orbits only, so `polyOf.get(other)` is
  `None` for any orbit a seed reaches. A leftover linking to a **quipu** orbit —
  which would be the most interesting result the run could produce, an LNA the
  theorem misses turning out to be in a theorem class after all — is therefore
  reported as an ALARM and **excluded from the union**, not recorded as a merge.
  It did not happen here, but the test would have hidden it.
* An alarm does not stop the run or mark the group, so the two alarms sat in the
  log for five hours.

---

## E-031 — Is reorientation a mutation, and does it help the search?
*2026-09-17* · **yes, and it is the merge step rather than the search that needed it** → F-036

Suggested from outside the code: all orientations of a relation-free tree should
be mutation equivalent, by mutating at a source -- which only flips the outgoing
arrows -- and then at each newly created source, with left mutation at sinks for
the other direction. The classification has been merging classes on a shared
hereditary form all along, which is a statement about mutation classes resting on
a fact about derived ones, so this is the step in between.

### 1. It is true, and the sequence can be written down (F-036)

Right mutation at a source of a relation-free tree quiver reverses exactly the
arrows there, creates no relations and does not renumber -- 2339 cases over every
tree and orientation at orders 3 to 7, no failures, and the same for left
mutation at a sink. `reflections.reflectionSequence` turns one orientation into
another by flipping every vertex on one side of a differing edge, in topological
order; verified against the engine at orders 4 to 10, every step admissible, 5568
reorientations, no failures.

Right mutations **alone** already connect every orientation of every tree up to
order 8 -- which matters because the search walks only those -- but the worst
distance grows: 4, 6, 9, 12, 16 at orders 4 to 8, against 3, 3, 6, 6, 10 when
sinks are allowed too.

### 2. It does *not* widen a bounded search

Allowing a search to jump to any orientation whenever it reaches a relation-free
quiver (`reflections.linesReachedThroughReflections`), against the plain search
from the same start:

| start | depth | plain | with reorientation |
|---|---|---|---|
| `kA_7` | 2 | 5 LNAs | 5 |
| `kA_7` | 3 | 8 | 8 |
| `kA_8` | 2 | 5 | 5 |
| `kA_8` | 3 | 9 | 9 |

Not one new LNA, and the reason is in part 1: a right mutation at a source *is* a
reflection, so the plain search already performs them where they are cheap. What
the lemma adds is the reflections that are **not** cheap -- and those do not lead
anywhere new within the depth either. Seeding backwards, from every one of the
128 orientations of each quipu of order 8 at depth 2, reaches 3 LNAs for the line
and 0 or 1 for every other quipu: the hereditary side is simply a long way from
any LNA.

### 3. It is the merge step that needed it

`mergeReport` merges two classes when both reach the same tree. Searching every
LNA to depth 4 and pairing the ones that reach the same tree:

| | `n = 7` | `n = 8` |
|---|---|---|
| LNAs reaching a relation-free quiver | 23 of 132 | 26 of 429 |
| pairs reaching the same tree | 79 | 80 |
| of those, pairs reaching **isomorphic quivers** | 61 | 65 |
| pairs whose orientations must be joined | **18** | **15** |

Those merges are derived equivalences until the orientations are joined, and
joining them takes 4 to 11 mutations, which is beyond any depth the pipeline runs
at. `reflections.mutationBridge` produces the whole path and the engine confirms
where it lands.

**A false start worth recording.** The first version of
`relationFreeQuiversReached` recorded what the search of the *opposite* algebra
found without carrying it back through the opposite, so half the orientations in
each list were reversed. It made the answer to part 3 come out as 4 of 79 at
`n = 7` and 0 of 80 at `n = 8`, where the true counts are 18 and 15 -- the
measurement looked like a nearly-empty result when it is a fifth of the pairs.

### 4. What it does not do: H-012

The natural hope was that the bridge would close H-012's gaps -- the 8 LNAs at
`n = 8` and 44 at `n = 9` that the known moves do not join to their stripped
form. It cannot: **none of those 52 rows reaches a relation-free quiver at all at
depth 4**, on either side of the pair, so there is no tree to bridge through.
Recorded against H-012.

---

## E-030 — Two other families, measured by Coxeter polynomial
*2026-09-17* · **no tree outside the quipu shape; quipus *with* relations carry every class the theorem misses** → F-033, F-034, F-035, H-014

The question: the quipu theorem names a class by a tree with no relations, and
every LNA it does not cover has to be classified some other way. Is there a
second family that does for those what quipus do for the almost separate ones?
Two candidates were measured, both through the Coxeter polynomial, which is a
derived invariant and forces the number of simples — so a candidate can only
match an LNA of its own length, and a difference settles it.

### The instrument

`invariants.coxeterCoefficients`. For a quiver with no oriented cycles the Cartan
matrix is unimodular, so

    det(lambda I - Phi) = det(lambda C^T + C),

which is a determinant of integers and `lambda` with no inversion in it. Taking
it at `n + 1` points and interpolating gives the polynomial as an integer
coefficient tuple: exact, hashable, and a hundred times faster than the symbolic
route — 16796 LNAs at `n = 11` in 19 s.

`coxeterTables` turns that into the three tables a search matches against: every
LNA's polynomial, every quipu's, and every LNA's **status** — in a quipu class by
F-032's moves (QUIPU), in none because no quipu of the order carries its
polynomial (NOT_QUIPU), or neither (UNPLACED). It reproduces F-032 exactly:
1421/9/0 at `n = 9`, 4600/260/2 at `n = 10`, 14149/2631/16 at `n = 11`.

### 1. Every tree, not only the quipus (F-033)

`python families.py trees 9 10 11` — and 12.

| order | trees | not quipus | cospectral with a quipu | **leads** |
|---|---|---|---|---|
| 9 | 47 | 29 | 3 | **0** |
| 10 | 106 | 70 | 0 | **0** |
| 11 | 235 | 171 | 7 | **0** |
| 12 | 551 | 424 | 15 | **0** |

F-031 asked this of the trees of maximum degree three; these are all of them,
694 non-quipu trees over the four orders. 25 of them do share a polynomial with
some LNA — but in every case every LNA under that polynomial is one the moves
place in a quipu class, and the tree is cospectral with that quipu rather than
isomorphic to it. Two tree algebras are derived equivalent exactly when the trees
are isomorphic, so those are refuted outright, with no search.

### 2. Quipus that carry relations (F-034)

`python families.py quipus 9 --min-arrows 2`, and the same at 10 and 11.

Enumerated: each quipu of the order, each orientation of its edges up to the
tree's automorphisms, each admissible monomial ideal — an antichain of directed
paths under "is a contiguous subpath of". The relation-free ideal and the
linearly oriented line are left out, being the quipu theorem's own case and the
LNAs themselves.

| order | shortest relation | ideals walked | matching | polynomials covered |
|---|---|---|---|---|
| 9 | 2 arrows | 370 483 | 3677 (2820 up to isomorphism) | **2 of 2** |
| 10 | 2 arrows | 3 411 263 | 306 624 | **7 of 7** |
| 11 | **3 arrows** | 5 465 194 | 1 246 011 | **20 of 20** |

Every Coxeter polynomial of an LNA outside a quipu class is carried by quipu
algebras with relations, in quantity — 1746 of them on `3033030` alone at
`n = 9`, over 16 different quipu shapes.

**Confirmed from the other side.** A polynomial match is necessary and not
sufficient, so the classes were also walked: every quiver a mutation search out
of an LNA reaches is in its class by construction, and `reachedQuipuAlgebras`
keeps the ones that are quipus with monomial relations. Walked from **every** LNA
outside a quipu class, and its dual, to depth 3:

| | `n = 9` | `n = 10` | `n = 11` (sample of 200) |
|---|---|---|---|
| LNAs outside a quipu class | 9 | 262 | 2647 |
| reaching a quipu with relations | **9** | **262** | **200 of 200** |
| reaching none | 0 | 0 | 0 |
| confirmed per LNA: min / median / max | 8 / 16 / 18 | 3 / 20 / 84 | 5 / 27 / 117 |
| distinct algebras confirmed | 178 (depth 4) | 3510 | -- |

Every algebra reached this way is in the enumeration — the only things reached
and not enumerated were the linearly oriented lines, which are the LNAs. The
`n = 10` walk takes 7 minutes and `n = 11` would take an hour and a half, so
`n = 11` is a random sample of 200 of its 2647 rows, seeded, and every one of
them reaches something too.

**And what they reach is not arbitrary.** `python families.py members 9` prints
the confirmed members per class, simplest first, and seven of the nine LNAs at
`n = 9` reach the *same quipu* — `P^(6)_(1,1)`, the line on eight vertices with a
pendant at the second — carrying their own relations shifted by one vertex with
one absorbed into the branch. That is the raw material for H-014's second part,
which is the part that would make this a theorem rather than a census.

**The smallest instance of the phenomenon is at order 4.** `D_4` with a single
two-arrow relation has the Coxeter polynomial of `kA_4`, and **one** mutation
takes it there. So a quipu with relations being derived equivalent to a line is
not exotic; what is new is that it happens for the lines the theorem cannot name.

### 3. A relation of two arrows is not free on a quipu (F-035)

`python families.py free 4 5 6 7 8`. Deleting every two-arrow relation and
comparing the polynomial:

| order | ideals with one | polynomial kept | **changed** |
|---|---|---|---|
| 4 | 11 | 7 | **4** |
| 5 | 72 | 48 | **24** |
| 6 | 543 | 300 | **243** |
| 7 | 4160 | 2138 | **2022** |
| 8 | 34938 | 15337 | **19601** |

The control passes: on the **linearly oriented line** the polynomial is kept
every time, which is `corollary:lengthtworelations` of arXiv:2310.08346 and
F-028. Off the line it fails immediately — the order-4 counterexample above is
the smallest. So `--min-arrows 3`, which is what makes order 11 affordable, is a
real restriction of the family and not a normalisation, and the order-11 run is a
statement about ideals whose relations all have three arrows or more.

### 4. Relation-free sightings, as instrumentation

`search.relationFreeSightings` records every quiver a search reaches with no
relations left, with whether its underlying graph is a tree, whether it is a
quipu, whether the quiver has an oriented cycle, and the path that got there;
`classify.py --sightings FILE` writes them as JSON lines. Nothing is recorded
unless a sink is open, so it costs nothing when it is not asked for. It answers a
question nobody had asked of the searches: the hereditary algebras they pass
through are collected and all but the first thrown away, and whether any of them
is *not* a tree has never been looked at.

**First measurement: `python classify.py 9 --sightings` records none at all.**
That is the expected answer and worth having written down. The quipu classes are
named by the theorem without a search, and the only rows the classification does
search at `n = 9` are the nine in no quipu class -- which contain no hereditary
algebra, so there is nothing for a search to find. The instrument will only have
something to say where a search runs into a class that *does* have one, which
means `--form-depth` at a length with an unnamed class, or the deeper resolving
runs of `merges.py`.

---

## E-029 — Reading `proposition:doubleMutation`, and what it does to coverage
*2026-09-17* · **a mechanism, found by reading; n = 9 needs no search, n = 10 is 16 orbits** → F-032

The first session run locally, with the `.tex` sources of all three papers to
hand. NOTES item 5 said to look for mechanisms before rules; this is the third,
and the largest.

1. **Read** the proposition and its proof (`main.tex` lines 331–497 of
   arXiv:2310.08346v1). Replaced the second-hand lead in the literature summary.
2. **Implemented** it as an interval rewrite, dual through F-026, `s = 1`
   extension behind `allowSource`. Reproduced `example:mutationToA11_5` steps
   `L_8`, `L_9` and `example:A13tworelations` steps `R_1`, `R_1^2` exactly.
3. **Verified** against the engine at `n = 5..10`, 14 processes: 17556
   confirmations, 0 failures, 1 min 22 s for `n = 9, 10`.
4. **Measured** coverage with and without the table, the free move, and the
   extension. The extension changes nothing (`doubles, paper only` gives the same
   numbers at `n = 7, 8, 9`). Interior-only (`s > 1` and `t < n`) gives 273 / 429
   and 770 / 1430 at `n = 8, 9`, with or without the free move.
5. **Named what is left** by Coxeter polynomial against the quipu polynomials of
   the order (34 at `n = 10`, 60 at `n = 11`, matching F-010's collision counts),
   and by the two implemented non-piecewise-hereditary certificates.
6. **H-012, by known moves only**: with the double mutation, edge moves and the
   rule table (all mutations), 8 of 429 LNAs at `n = 8` and 44 of 1430 at `n = 9`
   are not joined to their stripped form. That is a statement about the known
   moves, not a counterexample; nothing settled.
7. **H-011's sharp question**: of the 262 rows left at `n = 10`, 190 have a
   relation at both source and sink (73 %); at `n = 11`, 1743 of 2647. Not a
   characterisation — but the question has changed, since what is left is no
   longer a gap in the rules (F-032).
8. **What an interior application does** (`s > 1`, `t < n`), `n = 6..10`: at
   `n = 10`, 7164 applications; maximum overlap down in 1416, up in 1416, same in
   4320; relation count −1 in 1430, +1 in 1430. So interior moves *do* lower
   overlap when bystanders cross `r`, while on an isolated pair they only slide
   it. A claim to the contrary was written into H-010 and corrected before
   commit.
9. **Relation dual.** The first `merges.py` smoke test at `n = 10`, depth 3,
   found 18 links and every one was a relation dual, reached at depth 0 because
   the search starts from each member's dual. Closing the orbits under it is free
   and takes `n = 10` from 16 leftover orbits to **12 in 7 polynomial groups**,
   `n = 11` from 86 to **54 in 20**. So the derived classes number 43 to 48 at
   `n = 10` and 84 to 118 at `n = 11` (H-013).
10. **Reading, not yet used**: `lemma:taupathimpliesnotpwh` certifies both members
   of `example:A10double`, which A9, A13 and vertex deletion all miss.

11. **Sweep for further mechanisms, negative.** arXiv:2310.08346 §2 has exactly
    two derived equivalences (the free move and this one). arXiv:2305.06642's
    `algorithm:CRswap` needs almost separate relations and is what the seeding
    already encodes; `corollary:2Rels` is the free move restricted to that case.
    arXiv:2112.08129 is the procedure itself. No further *merge* mechanism in the
    three papers; the one unused tool is a *certificate*,
    `lemma:taupathimpliesnotpwh`.

Reproduce: `python overlaps.py 9 10 11 --free --doubles --no-rules`, and
`pytest tests/test_double_mutation.py` (the `n = 9, 10` engine checks are
`slow`).

**What this makes obsolete in the handover plan.** The raised-bound rule
discovery (`--max-arrows 7 --max-width 8` interior, `--extend` again): its purpose
was coverage, and the table adds nothing at `n = 10` on top of this. Not run.

---

## E-028 — Checkpointing the classification, and what a resumed n = 8 gives
*2026-09-17* · **the published table, out of two interrupted halves**

E-008 records that a long classification does not survive the night, and F-014
records the cost: the n = 10 run named 61 classes and kept none of them, because
only the search step wrote anything as it went. The naming and resolving steps
now write the table after every class and record what they finished in a sidecar
JSON file beside the CSV.

Checked by interruption rather than by argument. `classifyLength(8, ...)` with a
budget of zero seconds stops inside the search with 10 of the 429 rows still
unplaced and exits 2; resumed, it finishes and gives

    133, 65, 64, 64, 40, 26, 13, 10, 9, 4, 1

which is arXiv:2305.06642's n = 8 table exactly, with nothing left as a
candidate. Ten tests in `tests/test_checkpointing.py`, ~54 s.

**One thing the first version got wrong.** `stoppedEarly` was computed as "the
deadline has passed by the time the run returns", which is not the same as "a
step broke out". A zero budget at n = 7 expires before the first check and still
leaves a complete classification, because seeding places that whole length
outright -- and the run reported itself as stopped and unfinished, which would
send someone to resume a run with nothing in it. It now means a step actually
broke out of its loop, and is pinned as a test.

**What is still not checkpointed, and deliberately.** `probe.py` holds its
search in memory and prints at the end, so the deep run H-010 asks for -- seven
mutations at clearance 9, where six already took 43 minutes -- either completes
or is lost. `overnight.sh` runs it under a hard timeout for that reason rather
than pretending otherwise. Making it resumable is the obvious next piece of work
if that probe is going to be run repeatedly.

---

## E-027 — The free move, the square walked instead of searched, and the trees
*2026-09-16* · **the free move beats the whole table; two new families; no non-quipu tree reached** → F-028, F-029, F-030, F-031

Four threads, all suggested from outside the search, and the cheapest of them is
the largest result this project has had.

### 1. Relations of two arrows are free (F-028)

`corollary:lengthtworelations` of arXiv:2310.08346 has been sitting in
`research/literature/` unused since that paper was read. Deleting every
two-arrow relation keeps the Coxeter polynomial in all 4861 cases at `n = 3..10`
and keeps the quipu name in all 44320 cases at `n = 3..13`. Added to the orbit
computation it merges 20052 pairs at `n = 12` that the 1794 verified rules of the
table at the time do not, and cuts what a search must still place by 76% at
`n = 9` and 90% at `n = 8`.

The reduced space is exactly the LNAs of one fewer vertex, by shortening every
relation by one arrow — checked as a bijection, not just a count, for
`n = 4..13`.

### 2. Walking the square instead of searching for it

NOTES backlog 26. Open a relation with one mutation — F-027 says that always
gives a 2-by-k square — then mutate at consecutive vertices along the side it
opens. The cost is linear in the relation's length, where a search is exponential
in the depth, so a family's later members cost no more than its first. The
instrument was validated against F-020's lone slide, which it reproduced out to
`d = 9` in 13 seconds; discovery had reached `d = 7` at far greater cost.

**The march must be allowed to repeat its first vertex.** A first version
advanced one vertex per mutation and found nothing at all, because the two known
end families are `[1, 1]` and `[2, 2]`. That is worth recording: a walk that
cannot stand still cannot see either of them.

| planted | where | outcome |
|---|---|---|
| a lone relation, `l = 3..6` | interior | nothing — the square closes only by undoing itself, at every length, not just F-027's `l = 5` |
| a long relation and a two-arrow one, all gaps | interior | nothing but the two-arrow relation sliding past; the long one is a pure spectator |
| two relations of ≥ 3 arrows, `L, m = 3..7`, all gaps | interior | nothing, in 305 configurations |
| the same, all gaps ≥ 2 | either end | nothing |
| an overlapping pair, gap 1 or 2 | interior | **F-030**, the pair-to-triple family |
| a relation at an end | either end | **F-029**, the doubling |

So the square is a real mechanism and a cheap one, and what it finds is
concentrated exactly where F-022 and F-024 said the action was: on relations
that overlap, or against an end. A relation with room around it does nothing,
whatever its length and whatever it is next to.

### 3. The two families, and what they cost the story

F-029 is not expressible as a table rule at all — the honest conclusion is that
the encoding needs widening, which is NOTES backlog 27. F-030 is expressible, and
the table already held its first two members and none of the rest, which is
H-008's shape for the third time.

Between them and the free move, **`A_8` is fully covered with no search**: 21
orbits, nothing left. `A_9` falls from 380 orbits and 222 rows needing a search
to 77 and 37.

**A false start worth recording.** The collapse directions of F-029 were first
written by inverting the doubling's condition and sequence by hand. Both were
wrong: the sequence, because the procedure relabels and a left mutation at a
vertex is not undone by a right mutation there — the collapse is two right
mutations at the *source*; and the condition, which fired on 97 cases where 282
were available. Defining each collapse as "the LNA whose doubling is this one"
fixed both at once, and all four moves then fire 907 times apiece over `n = 7..10`
with no failures.

### 4. What is left at n = 9, now that it is small enough to read

37 rows in 77 orbits, and they have a property in common: **every one has a
relation at the source and a relation at the sink**, and every one is already
reduced. Only 1 of the 37 is certified non-piecewise-hereditary. At `n = 10` the
same two counts are 670 and 660 of 887, so the characterisation is strong there
but not complete. Recorded against H-011, whose mechanism it is the natural limit
of: an LNA with both ends occupied has no free end to walk a run to.

### 5. Trees that are not quipus (F-031)

Every tree of maximum degree three up to order 12, tested for quipu-ness and then
compared by Coxeter polynomial against every LNA of the same length. The first
non-quipu appears at order 10 and is unique — the centre with three neighbours,
each carrying two leaves, which is the tree the question was asked about. Eleven
non-quipu trees over orders 10, 11 and 12; 80444 LNAs compared; not one shared
Coxeter polynomial.

---

## E-026 — The dual as a mirror, and what a rule walks through
*2026-09-16* · **410 duals, all holding; every intermediate a square with a side of two; no shortcut survives** → F-026, F-027, R-011

Four things suggested from outside the search, all cheap, and the first of them
corrects a finding made the same day.

### 1. The proper mirror of a rule

F-025 read an asymmetry between the two ends of the quiver off a comparison
between a rule at the sink and *the same pattern* at the source. That is not the
mirror. The mirror is the relation dual: reverse every arrow **and** exchange
right mutation for left.

| | |
|---|---|
| at the sink | `(0:3) (1:6) -> (0:2) (1:6)` via `[2, 2]` -- 4 confirmations, 0 failures |
| the same pattern at the source | 1 confirmation, **3 failures** |
| its **dual** at the source | `(0:6) (4:3) -> (0:6) (5:2)` via `[-7, -7]` -- 4 confirmations, **0 failures** |

Checked at `(l, m)` = (3,6), (4,7), (3,4), (5,6). The pattern's dual is a pair
sharing an *end*, not a pair sharing a *start*, which is why comparing a pattern
with itself at the other end says nothing. R-011.

### 2. Closing the table under it

Whether the sequence's **order** reverses under the dual was the one thing not
obvious. Over 50 rules sampled from both halves of the table, order **kept**
works for all 50 (8 of them exclusively; the other 42 have sequences symmetric
enough that both work) and order reversed works for none exclusively. So the
order is kept.

Then, over the whole table: 1384 rules, **none self-dual**, **410 duals
missing**, verified at up to four lengths each --

> **410 hold, 0 fail, 0 never apply**, in 106 s on four processes.

The table is generated closed now (364 floating, 1430 anchored). Coverage barely
moves -- 7 rows at n = 10, none at n = 9 -- because those orbits were already
joined another way. The value is that it is free, that it halves what a search
has to look for, and that it is what corrected F-025. F-026.

### 3. What a rule walks through

Every multi-mutation rule in the table run one step at a time, each intermediate
quiver classified:

| | |
|---|---|
| commutative squares | **1364**, sides `2 x k` for k = 2..8 |
| squares with a short side other than 2 | **0** |
| a line again, mid-sequence | 336 |
| anything else | 5 |

So a mutation leaves the line only into a square with a side of exactly two, and
a rule is: open a relation into such a square, do something along its long side,
close it back. Exactly the structure the search has been finding by brute force.

**Tried and it does not immediately give a construction.** Opening the lone
relation of `00500000` in A_10 into its 2-by-4 square and searching four
mutations over the square's vertices finds one way back to a line: `[-3]`, the
undo. A square with nothing to interact with closes only onto itself, which is
F-023's lesson again -- the companion relation is the whole point. F-027.

### 4. Whether mixing directions shortens the rules already known

**Mixed sequences are not an unexplored region.** `localMutationSequences` tries
both signs at every vertex at every step, so every run so far has been free to
find them, and 274 of the table's 1794 rules do mix a left mutation with a right
one -- 82 floating and 192 anchored, almost all of them two or three mutations
long.

**And no rule in the table got shorter.** 150 rules of three mutations or more
were sampled; for each, one LNA it matches, and a search at one mutation fewer
over the vertices of its window. Three came back with a shorter sequence at that
LNA:

| rule | listed | shorter, at one LNA |
|---|---|---|
| `(4:2) -> (0:2)` | `[5, 4, 3, 2]` | `[-2, -7, -8]` |
| `(5:2) -> (0:2)` | `[6, 5, 4, 3, 2]` | `[-2, -8, -9]` |
| `(0:2) -> (6:2)` | `[-3, …, -8]` | `[1, 9, 8]` |

All three are the lone short-relation slide of F-020, which would have been the
interesting place to find a shortcut -- and **none of the three is a rule**.
Stated as floating rewrites they give 1 confirmation and 21 failures apiece;
anchored to the left end, 1 confirmation and 8 failures. They work at the single
LNA the search tried, where the window is flush against both ends of the
shortest quiver its width fits in, and nowhere else. Checked for every d from 1
to 7 with the same answer.

**So F-020's one mutation per arrow travelled stands**, and this is the evidence
against the obvious objection to it. A run that finds a shorter sequence for one
LNA has found nothing until `verifyMove` says otherwise.

**The limits of this, stated so it can be redone properly.** One matching LNA per
rule, and mutations only within one vertex of the window. A shortcut that needs
to reach further out, or that only applies at some positions, would not show up
here.

---

## E-025 — Discovery with the bounds raised, and probes deep and wide
*2026-09-16* · **3045 rules, n = 8 to 98%, and H-010 tested two deeper** → F-024, F-025

E-024's residue said what to do: of the LNAs at n = 9 for which no rule had the
pattern at all, 32 of 33 contained a relation of five arrows or more, and every
run so far had stopped at five arrows in a six-arrow window, where such a
relation cannot sit beside another. So: raise the bounds, and probe the
configurations the findings actually turn on.

### The discovery run

```bash
python discover.py --anchor both --max-arrows 7 --max-width 8 \
    --anchor-lengths 13,14 --jobs 3 --verify-cap 12
```

| stage | |
|---|---|
| patterns x ends x lengths | 412 x 2 x 2 = 1648 searches, 3 mutations each |
| rewrites described | 5807, in 6146 s |
| recurring at both lengths | 3826 |
| verified with no failures | **3045**, in 579 s |
| neither listed nor a floating rule restricted to an end | 2929 |
| changing the orbit partition at lengths 6 to **11** | **728** |

Two hours of search, ten minutes of verification. 1469 of the fresh rules are at
the source and 1467 at the sink, which is the consistency check the relation
dual demands.

**The criterion had to change with the bounds, and this is the trap.** A window
of nine arrows does not fit in A_9 at all, so judging by what changes the
partition at n <= 9 -- which is what E-023 and E-024 did -- *cannot* select a
wide rule however useful it is. Judged that way this run yields 209 rules, all
of window 7 or 8. Judged at lengths 6 to 11 it yields **728**, of which 208 have
a window of nine arrows or more. The earlier curations should be read with that
in mind: they were not wrong at the lengths they measured, but they could not
see past them.

**Coverage, with no mutation search at all:**

| n | LNAs | before this run | after |
|---|---|---|---|
| 7 | 132 | 100% | 100% |
| 8 | 429 | 95% | **98%** -- 10 rows left |
| 9 | 1430 | 73% | **84%** |
| 10 | 4862 | 55% | **63%** |
| 11 | 16796 | 43% | **47%** |

The measurement now reaches n = 10 and n = 11, which it had not before; at 16 s
for n = 10 and about four minutes for n = 11 there was never a reason not to.

### The probes

`probe.py`, written for this run, plants one named pattern and enumerates what
the mutations near it reach, reporting the arrows a run can actually rewrite.

**The frozen pair, deeper and honestly interior (F-024).** The earlier probes
allowed mutations at every vertex of A_13; re-run in A_21 where the ends are out
of reach, `(1:3) (2:3)` reaches 2 LNAs at three mutations, 4 at four, 4 at five
and 6 at six -- against 8, 14, 22 and 36 with an end in reach. None of them
lowers the overlap. Six mutations *with* an end in reach does lower it, to zero,
by walking the relation down to arrow 1; that sequence is not translation
invariant and fails at every shift tried.

**Long pairs, in the interior.** `(1:3) (2:6)`, `(1:5) (2:6)`, `(1:5) (2:7)`,
`(1:3) (2:7)` at four mutations, and the equal pairs `(1:6) (2:6)` and
`(1:7) (2:7)`: every one frozen, and the overlap goes *up* in a third to a half
of what they reach.

**Long pairs at the ends, which is where the new family came from (F-025).** The
same four against each end at three mutations: nothing at all moves at the
source, and at the sink every one loses an arrow off the shorter relation --
`(1:3) (2:6)` and `(1:3) (2:7)` going to overlap 0 outright.
`sinkShortRelationShrinkRules` generates that family and it is verified for
every 3 <= l < m <= 9, 21 members, 4 confirmations apiece, no failures. The
mirror at the source gives 1 confirmation and 3 failures.

**A run of three no longer dissolves when the relations are long.** F-022 had it
that a run of three heavily overlapping relations comes apart where a run of two
does not, on the evidence of `(1:3) (2:3) (3:3)` and two others. At three
mutations `(1:3) (2:6) (3:7)` and `(1:5) (2:6) (3:7)` reach two LNAs each and
neither lowers the overlap. So the dissolution of a run of three is not a
property of the run; it is a property of the *short* runs that were tested, and
what it probably costs is mutations -- F-020's one-per-arrow-travelled again.
Do not quote F-022's run-of-three line without that qualification.

### What to do next

The bounds can go up again -- `--max-arrows 9 --max-width 10`, and the interior
run at the same bounds, which this session did not get to. And the whole
judgement should now be made at lengths 6 to 11 as a matter of course, since it
is affordable and the alternative silently discards every wide rule.

---

## E-024 — Widening the rules to tolerate a bystander
*2026-09-16* · **625 verified, and A_7 needs no search at all** → F-023, H-011

F-023 said what stops a known rule from firing is a relation in its window that
it does not touch. This is that read as a construction rather than a diagnosis:
take each rule in the table, put one untouched relation -- a *spectator* --
somewhere in its window, growing the window by up to three arrows to make room,
and let `verifyMove` decide.

```bash
python discover.py --extend --jobs 4
```

| stage | |
|---|---|
| widenings generated from the 368 rules then in the table | 10609, in 13 s |
| of those, firing on an LNA no search had placed at n <= 9 | **1008** |
| verified with no failures | **625**, in 474 s |
| changing the orbit partition at n <= 9 | **270** -- 125 floating, 145 anchored |

Eight minutes, no search, no mutation run speculatively. The filter is what makes
it affordable and is worth keeping: generate freely, throw away everything that
would not fire on a row still needing a search, then verify.

**Coverage with no search at all.**

| n | theorem | before this batch | after |
|---|---|---|---|
| 6 | 81% | 100% | 100% |
| 7 | 67% | 96% | **100%** |
| 8 | 54% | 81% | **95%** -- 23 rows left of 429 |
| 9 | 43% | 60% | **73%** -- 392 of 1430 |

A classification of A_7 is now a table lookup. At n = 8 twenty-three rows need a
mutation search and at n = 9, 392.

**The number H-011 said to watch stayed small.** Re-running the blocked-rule
diagnostic -- for each unplaced LNA, the rule whose left-hand pattern is present
with the fewest extra relations in its window:

| | n = 8 before | n = 8 after | n = 9 after |
|---|---|---|---|
| blocked by a bystander | 126 | 18 | 359 |
| **no rule has this pattern at all** | 29 | **5** | **33** |

So the mechanical half is still the whole story: what is left is overwhelmingly
more of the same, and another widening pass is the obvious next run.

**And the 33 turn out to say the same thing.** Looked at by hand, they are
almost all a pair of relations one of which is *long*: `(1:3) (2:6)`,
`(1:5) (2:6)`, `(1:5) (2:7)` and the like. 32 of the 33 contain a relation of
five arrows or more and 23 contain one of six or more -- and discovery has never
been given a pattern like that. Both runs used `--max-arrows 5 --max-width 6`,
so a six-arrow relation could not appear beside another at all. This is F-013's
lesson again: *absence of a pattern from the table is evidence about the search,
not about the mathematics*. Raise the bounds before concluding anything about
these.

**What this does not say.** Every one of these rules was verified, but 355 of
the 625 are not listed, and the 270 that are were chosen for changing the
partition at n <= 9. That is a curation against the lengths measured, not a
claim about n >= 10. Re-run the command.

---

## E-023 — Discovery against the ends of the quiver
*2026-09-16* · **630 anchored rules, and coverage at n = 6 becomes complete** → F-023

The first run of `discoverAnchoredMoves`, on the framework F-022 added.

```bash
python discover.py --anchor both --max-arrows 5 --max-width 6 --jobs 4 --verify-cap 12
```

74 patterns x 2 ends x 2 lengths (A_11 and A_12) = 296 searches at three
mutations, margin 3.

| stage | |
|---|---|
| rewrites described | 1344, in 286 s |
| recurring at both lengths | 892 |
| verified with no failures | 724, in 70 s |
| a floating rule restricted to an end | 94 |
| genuinely anchored | **630** -- 315 at each end |
| changing the orbit partition at n <= 9 | **229**, and those are what is listed |

**A first attempt at the same run had to be abandoned**, and why is worth
recording. `verifyMove` enumerated the LNAs by building a path algebra for each,
which at length 12 is 58786 of them and six seconds -- per rule, and there were
886 to check. The enumeration is now cached as relation-length rows
(`nakayama.allRelationLengths`), the algebras built only where a rule actually
matches: 0.3 s instead of 6.5, and the verification of all 886 fell from hours
to 70 seconds. Anything that verifies many rules over the same lengths should go
through that function.

**What the run cost and bought.** Ten minutes end to end. Coverage with no
search at all: 100% at n = 6 (from 83%), 96% at n = 7 (72%), 81% at n = 8 (57%),
60% at n = 9 (45%).

**What is still not reached, and the next question.** 567 rows at n = 9, all of
them heavily overlapping, 306 at overlap 2. The diagnostic that pointed at the
spectators -- for each unplaced LNA, the rule whose left-hand pattern is present
with the fewest extra relations in the window -- says at n = 8 that 126 of 155
are blocked by a bystander and 29 by having no rule at all. Run it again after
the next batch: the count of "no rule has this pattern" is the one to watch,
because it is the part more discovery cannot fix.

---

## E-022 — Whether the relation dual widens the move orbits
*2026-09-16* · **it halves the orbit count and adds no coverage**

The relation dual -- reverse every arrow, renumber -- is one of the three
class-preserving operations of arXiv:2305.06642 and holds for *any* LNA, not
only an almost separate one (`nakayama.relationDual`). It is free, it is not in
the move orbit, and the obvious thought is that adding it would carry rows
across the overlap line for nothing. It does not.

| n | floating | + anchored | + anchored + dual |
|---|---|---|---|
| 7 | 95 covered, 84 orbits | 107, 71 | 107, **45** |
| 8 | 246, 310 | 274, 277 | 274, **156** |
| 9 | 644, 1106 | 726, 1019 | 726, **542** |

The orbit count roughly halves at every length and the covered count does not
move by one row. The reason is structural rather than accidental: the almost
separate condition is itself dual-symmetric, so the dual maps seeded to seeded,
and the rule table already contains the mirror of every rule it contains, so it
maps orbit to orbit. The dual therefore identifies orbits pairwise and never
joins a covered one to an uncovered one.

**Worth knowing, and worth not repeating.** Halving the orbit count is real and
would be worth having if orbits were the expensive object; they are not, the
uncovered rows are. Do not reach for the dual again expecting coverage.

---

## E-021 — What can be done to a heavily overlapping run, in the interior
*2026-09-16* · **the pair is frozen, a run of three is not** → F-022, R-010

The experiment H-003 asked for, aimed where F-021 says to aim it. Each pattern
planted in the middle of A_13 at offset 4 -- four arrows of empty quiver on the
left, six on the right -- with `lnaMoves.localMutationSequences` enumerating
every admissible sequence at vertices within the margin, and the reached LNAs
reported by maximum overlap.

**The isolated pair, at four settings.**

| pattern | start overlap | mutations | margin | reached | any lower |
|---|---|---|---|---|---|
| `(1:3) (2:3)` | 2 | 3 | 3 | 8 | no |
| `(1:3) (2:3)` | 2 | 4 | 3 | 14 | no |
| `(1:3) (2:3)` | 2 | 5 | 3 | 22 | no |
| `(1:3) (2:3)` | 2 | 4 | 6 | 34 | no |
| `(1:4) (2:4)` | 3 | 3 | 3 | 17 | no |
| `(1:5) (2:5)` | 4 | 3 | 3 | 16 | no |
| `(1:3) (2:4)` | 2 | 3 | 3 | 16 | no (3 of them go **up** to 3) |
| `(1:4) (3:3)` | 2 | 3 | 3 | 16 | no (3 go up to 3) |

About a quarter of an hour in total, the depth-5 probe a third of it. Neither depth nor margin
is the dial: doubling the margin at four mutations reaches 34 LNAs instead of
14 and not one of them has a smaller overlap.

**The margin-6 row is stronger than it was written as, and the description was
wrong.** A margin of 6 around a pattern at arrows 5 to 8 of A_13 admits the
vertices 1 to 13 -- *every vertex of the quiver*, both ends included. So that
row is not a probe of the interior at all: it says that from `00003300000`,
**four mutations anywhere in A_13** reach 34 LNAs and none of them has a smaller
overlap. That is a claim about the LNA rather than about locality, and it is the
stronger one. It was recorded here as an interior probe with a wide margin,
which it was not. `probe.py --allow-ends` is how to ask that question on
purpose; without the flag the quiver is lengthened to keep the ends out of
reach, so an interior probe stays one.

**A third relation, and which third relations count.**

| pattern | overlapping run | three mutations |
|---|---|---|
| `(1:3) (2:3) (3:3)` | 3 | down to **0**, via `[6, 5, 6]` |
| `(1:3) (2:4) (3:4)` | 3 | down to **0** |
| `(1:4) (2:4) (4:3)` | 3 | down to **0** |
| `(1:4) (2:4) (3:4)` | 3 | down to 2, from 3 |
| `(1:2) (2:3) (3:3)` | 2 | 31 reached, none lower |
| `(1:3) (2:3) (4:2)` | 2 | 31 reached, none lower |
| `(1:3) (2:3) (5:2)` | 2 | 35 reached, none lower |

The parameter is the length of the run of relations linked by an overlap of two
or more, not the number of relations present: `(1:2) (2:3) (3:3)` has three
relations and is as frozen as the bare pair, because its first shares one arrow
and not two. The rewrites that dissolve a run of three were already in the
table, found at length 8 and in E-011 -- so nothing here is a new rule, and that
is the result. **Do not run a deeper interior search for a rule that pulls an
isolated pair apart**; four probes at three settings of depth and two of margin
say there is none to find, and F-022 says where the pair does come apart.

Reproduce with `lnaMoves.localMutationSequences(13, relLengths, lo, hi, steps,
margin)` on `lnaMoves.embedPattern(13, pattern, 4)`.

---

## E-020 — What the theorem and the move orbits reach, by relation overlap
*2026-09-16* · **the gap is exactly overlap two and above** → F-021

The measurement H-003 has been asking for since 2026-09-13, now that there is a
coordinate to make it in. Every LNA of a length partitioned into orbits under
the verified rules -- applied as rewrites on the relation lengths, with no
mutation computed, which F-017 licenses -- and an orbit called covered when it
contains one the quipu theorem names.

With the 123 floating rules:

| n | LNAs | overlap 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|---|
| 6 | 42 | 16/16 | 18/18 | 1/7 | 0/1 | | | |
| 7 | 132 | 32/32 | 57/57 | 4/33 | 2/9 | 0/1 | | |
| 8 | 429 | 64/64 | 169/169 | 9/132 | 4/52 | 0/11 | 0/1 | |
| 9 | 1430 | 128/128 | 482/482 | 24/484 | 10/247 | 0/75 | 0/13 | 0/1 |

covered over total at each maximum overlap. Two readings, and both matter.
Everything at overlap 0 or 1 is covered, at every length -- which is the almost
separate set exactly, so the theorem's reach is not merely *mostly* the low
overlap rows, it is precisely them. And above the line the table reaches 34 rows
out of 820 at n = 9, none at all past overlap 3.

The leftovers' heavily overlapping runs, commonest first at n = 9: `(1:3) (2:3)`
391 times, `(1:4) (2:4)` 198, `(1:3) (2:4)` and `(1:4) (3:3)` 144 each. By the
longest run in the LNA, 434 of the 786 have nothing longer than a pair.

Seconds per length. `python overlaps.py 6 7 8 9 --cores`, and the same numbers
are pinned in `tests/test_overlap.py`.

---

## E-019 — Two more families, and whether a rule's inverse is free
*2026-09-15* · **two families confirmed, the inverse shortcut refuted** → F-020

Acting on F-020's own lesson rather than raising the search bound.

**The families.** Two candidates read straight off consecutive window widths in
the enlarged table, then verified at three lengths each with `verifyMove`:

| family | d | result |
|---|---|---|
| `(0:2) (2:2)` → `(1:2) (d+2:2)` via `[-3, -5, …]` | 1–6 | 8 confirmations each, no failures |
| `(1:2) (4:2)` → `(0:2) (d+4:2)` via `[2, -7, …]` | 1–5 | 8 confirmations each, no failures |

Discovery had found d = 1 and 2 of each; d = 3 needs four mutations and d = 6
needs seven, so the rest were out of reach of any search run so far. Minutes to
check, against the hours a four-mutation run costs. Both are generated now.

**The inverse shortcut, and it does not hold.** Family A's listed left slide is
exactly its right slide with the sequence reversed, each vertex negated, and each
then moved one step toward zero — which looked like it might be a property of the
window's numbering and so give every rule's inverse for nothing. Applied to all
96 listed rules and verified:

| sequence | inverts | fails | inverse leaves the window |
|---|---|---|---|
| all one direction | 12 | 20 | 28 |
| mixed directions | 0 | 5 | 31 |

So **12 of 96**. The transform works for the slide families because a slide's
sequence is a single uniform run; it is not a general fact about the table, and
the spreading pair -- whose sequence mixes a right mutation with left ones -- is
the counterexample closest to hand. Do not try it again.

---

## E-018 — How far the gentle condition reaches
*2026-09-15* · **one class, the hereditary one** → F-019, R-008

Idea 22's premise, checked before implementing anything.

Every LNA of lengths 5 to 9 tested for `all(arrows in (0, 2))`, then grouped by
the class the classification puts it in:

| n | gentle LNAs | 2^(n-2) | classes containing one |
|---|---|---|---|
| 5 | 8 | 8 | `P^(0)_(0,4)` only |
| 6 | 16 | 16 | `P^(0)_(0,5)` only |
| 7 | 32 | 32 | `P^(0)_(0,6)` only |
| 8 | 64 | 64 | `P^(0)_(0,7)` only |
| 9 | 128 | 128 | `P^(0)_(0,8)` only |

So every gentle LNA is in the class of the path algebra of A_n, and the other 21
classes at n = 9 contain none -- including both members of the cospectral pair,
which have 18 members each and not one gentle among them. Seconds to run, against
the days an AAG implementation would have taken.

The reason is a theorem, not a coincidence: all relations of length 2 implies
almost separate relations, and operation 2 of `cor:EquivNakayamaAlgebras` drops
such a relation without changing the class, so dropping them all leaves the
hereditary algebra. Pinned for n = 4 to 10 in
`test_every_gentle_lna_is_the_hereditary_one`.

Also established while looking: `WebSearch` reaches the literature from the
session sandbox even though `curl` and `WebFetch` to arxiv.org are refused by the
egress proxy. Enough to find and identify a paper, not to read one.

---

## E-017 — n = 9 on the corrected engine
*2026-09-15* · **22 classes, then 20** → F-018

`classify.py 9` on the engine of F-015 with the gate of F-016, default depths
(`--depth 6 --resolve-depth 6`). About 90 minutes.

The run placed all 1430 rows and left nothing a candidate, but at **22** classes:
18 quipus plus `C(2,3,5)` (46 LNAs), `C(2,2,6)` (13), `C(2,4,4)` (8) and one not
piecewise hereditary. It reported three groups "proved distinct despite sharing a
Coxeter polynomial", two of which were `C(2,3,5)` against `P^(5)_(1,2)` and
`C(2,2,6)` against `P^(1,1)_(1,3,1)`.

**Diagnosis, in this order.**

1. The `C(...)` names come from `canonicalWeightType`, which reads them off the
   class' Coxeter polynomial — so they cannot separate two classes that share
   one. Circular.
2. A direct probe: iterative deepening from every member of each of the two
   classes, and from every member's relation dual, reporting every foreign class
   reached. Both merged **at depth 2**, from the first member tried and from its
   dual as well. Seconds, against the 90 minutes of the run.
3. (2,3,5) and (2,2,6) are domestic weight types, and their extended Dynkin trees
   are `P^(5)_(1,2)` and `P^(1,1)_(1,3,1)` — the two classes they were separated
   from. Computed, not asserted.

**Re-run of the post-search half only** (the rows were sound; only the merge step
was wrong), on the same table: both merged at the first depth tried, then
`3033030` certified not piecewise hereditary and `3345000` named tubular
`C(2,4,4)`. **1430 LNAs, 20 classes, 0 candidates, 1 separated** — the separated
group being the cospectral pair of F-010, exactly F-011. Under a minute.

**Then re-run whole, from scratch, on the fixed pipeline: the same 20**, with the
same sizes class for class. The naming order is visible in the log -- the theorem
names 18 classes, `resolveMergeCandidates` then merges 11 away (`2233030` into
`P^(1,1)_(1,3,1)` and `2334400` into `P^(5)_(1,2)` among them, both at depth 6),
and only then do the fallbacks name what is left: `3033030` not piecewise
hereditary by Proposition A9, `3345000` the tubular `C(2,4,4)`. One separated
group, the cospectral pair of F-010. About 35 minutes, 10 of them the search.

---

## E-016 — Are the move rules local?
*2026-09-15* · **yes, both halves** → F-017

H-009's own caveat, checked before anything else was built on it.

1. **Applicability.** `matchesAt` against a predicate reading only the window's
   cells plus one bit (a relation covering the window's first arrow having
   started earlier), over every rule in `VERIFIED_MOVES` x every admissible LNA
   x every window position:

   | n | comparisons | matches | disagreements |
   |---|---|---|---|
   | 5, 6, 7 | 67,712 | 150 | 0 |
   | 8, 9 | 924,352 | 1084 | 0 |

2. **Legality.** The whole table re-verified where each rule fits: **1218
   confirmations, zero failures** -- every match is an admissible sequence
   landing on the predicted LNA with the Coxeter polynomial kept.

Cheap: seconds for lengths 5 to 7, a couple of minutes for 8 and 9, and about
four minutes for the legality half. Rows 1 and 2 are tests now
(`test_whether_a_move_applies_is_a_local_condition`,
`test_each_rule_holds_wherever_it_applies`), so **do not repeat them by hand**.

The result that was not the question: the state has to be the **arrow** row, not
the vertex row -- see F-017. Anyone starting the CA literature sweep should start
there rather than from `relLengths`.

---

## E-015 — Every mutation the loosened gate newly allows
*2026-09-15* · **280 of them, all Coxeter-preserving** → F-016

Before switching the search's gate from the strict reading to the paper's
criterion, every mutation the switch would newly allow was enumerated and
checked. Walking out of every LNA of the length, at every vertex of every quiver
reached, comparing the old gate against the new one and computing the Coxeter
polynomial wherever they disagreed:

| n | depth | allowed by both | newly allowed | Coxeter moved | new gate narrower |
|---|---|---|---|---|---|
| 5 | 3 | 304 | 22 | 0 | 0 |
| 6 | 3 | 1450 | 138 | 0 | 0 |
| 7 | 2 | 1938 | 120 | 0 | 0 |

The last column matters as much as the others: a criterion that was *narrower*
anywhere would have meant the switch loses a mutation the published runs used,
and it never is.

Then the classifications, which are the acceptance test: n = 6, 7 and 8 all give
the same classes with the same sizes as before, and n = 7 dropped from 38
seconds to 20.

The old criterion is kept in `tests/test_procedure.py` as `strictlyMutable`,
which is what makes the comparison re-runnable; the first two rows are a test
now. **Do not repeat the n = 7 row** — about four minutes, and it says the same
thing as the other two.

---

## E-014 — The procedure on coefficients, against the one it replaced
*2026-09-14* · **agreement everywhere but two cases, which are R-007** → F-015

Five runs, all gated on `mutationIsPossibleAtVertex` so both implementations walk
the same mutations:

1. **One mutation, every admissible vertex, every LNA of n = 4..8.** 45 + 126 +
   462 + 1716 = 2349 comparisons, **zero** differences. Minutes.
2. **Depth-3 walks, n = 5 and 6.** 1446 and 7496 step comparisons, **zero**
   differences.
3. **Depth-3 walks, n = 7.** 37470 step comparisons, **2** differences, both
   after three mutations, both a relation the old implementation did not
   produce. These are the whole of R-007.
4. **The exact cleanup on the old steps' output**, n = 5 and 6 at depth 3 and
   n = 7 at depth 2: 1446 + 7496 + 4710 = 13652 comparisons, **zero**
   differences. Worth having separately, because it says the disagreement is in
   step 7 and not in the cleanup.
5. **Coefficients against the guess**, over every quiver within depth 3 of every
   LNA of n = 5 and 6: 1239 Cartan matrices, **zero** differences.

**A mixed engine is not an option, and this is how that was learned.** Running
the old steps 1-7 with the exact cleanup passed run 4 above and then reached
*two* different hereditary forms from `A_6` `3030` at depth 7 — a degree-4 tree
alongside `P^(1,1)_(1,0,1)`, which cannot both be one class. The exact cleanup
expects the relations step 7 produces; with step 7's output missing a relation it
cuts the wrong generators. Use one engine or the other, whole.

**Timings**, n = 7 over the 462 admissible single mutations: procedure 0.26 s
against 0.94 s, admissibility 0.14 s against 1.74 s. The exact versions are
3.6x and 12x *faster*.

**Do not repeat runs 1, 2, 4 and 5** — they are `tests/test_procedure.py` now.
Run 3 at n = 7 depth 3 takes about eight minutes and is worth re-running only if
step 7 changes.

---

## E-013 — Audit of the quipu symmetry, after R-006 was challenged
*2026-09-14* · **no defect found** → F-014

Four runs, in increasing cost:

1. **Canonicalisation against `networkx.is_isomorphic`**, over every quipu
   parameter pair of orders 3–11 (12 names at order 3 up to 28656 at order 11).
   Same canonical parameters iff isomorphic graphs, both directions, zero
   exceptions. Seconds.
2. **The paper's class-preserving operations against the quipu fibres**, over
   every LNA of lengths 4–10 with almost separate relations and no length-2
   relation. Orbits equal fibres exactly at every length; largest orbit 8, the
   paper's bound. Seconds. Now `tests/test_quipu_symmetry.py`.
3. **Tree enumeration against parameter enumeration.** The old
   `generateAllQuipus` (enumerate non-isomorphic trees, test the degrees) and
   `quipuForms.allQuipusOfOrder` (enumerate the P^(m)_(k) parameters,
   canonicalise) agree on the counts for orders 4–12, and no tree the first
   accepts is rejected by `quipuForms.isQuipu`. The old function only tested the
   "degree-3 vertices lie on one path" condition when there were more than three
   of them, which looked like a hole — but three branch vertices in a tree of
   maximum degree 3 always do lie on one path, since a path through two of them
   passes through the third, so four is the smallest number that can fail.

   Kept, since it is a genuinely independent route: it is now
   `quipuForms.quipusByTreeEnumeration`, with the degree test written out as
   `isQuipuByDegrees`, and the agreement is a test rather than a note here.
4. **Hereditary form by mutation search**, from all four long-relation members of
   `P^(1,4)_(1,0,1)` (`0003030`, `3030000`, `3060000`, `6000030`) and both of
   `P^(1,2)_(1,1,2)` (`0400030`, `3004000`), at depth 6, plus iterative
   deepening 2–6 from `3060000` and `3004000`. **Nothing reached** — no
   relation-free quiver from any of them. Tens of minutes.

Run 4 is the one that would have been independent of `thm:QuipuToAn`, and it is
simply out of range here, the same way `A_{7,(2,4)}^{(3,3)}` is (see the test
`test_the_theorem_answers_where_the_search_gives_up`). **Do not repeat it at
depth 6 or less.** Depth 7+ at n = 9 was not attempted and is expected to be
hours; the cheaper route to an independent check is a derived invariant computed
from the algebra, not a deeper search.

---

## E-012 — Pair slide at relation lengths 2 to 7
*2026-09-14* · **confirmed a family**

`lnaMoves.verifyMove` on the pair-slide rewrite for each `l`, both directions,
over lengths `l+3 .. l+6`. 22 confirmations per direction per length, zero
failures throughout. → F-013, confirming H-001.

Seconds to run. Should have been the first thing tried after finding the rule at
`l = 3`.

---

## E-011 — Interior discovery, three mutations
*2026-09-14, concluded 2026-09-15* · **44 rules, and 30 false ones caught** → F-020, R-009

`lnaMoves.discoverLocalMoves`, 26 patterns of up to 3 relations spanning ≤ 5
arrows, planted at offset 4 in A_13 and offset 5 in A_14, `maxSteps=3`,
`margin=3`. Re-run as

    python discover.py --jobs 2

after the original `interior.py` turned out never to have been committed.

**Discovery.** 52 searches, 315 s on two cores. 336 rewrites described, **166
recurring across both embeddings**, 134 of them not already in the table.

**Verification, first attempt — wrong, and instructively so.** All 134 checked at
the fixed lengths 7 to 10, which E-010 had used: 74 passed. But the lengths have
to follow the window, and a window of 9 arrows fits in A_10 at exactly one
position, flush against both ends. Each window-9 rule therefore got one
confirmation from one length.

**Verification, redone per rule at `width + 1 .. width + 4`.** 44 survive; **all
30 window-9 rules fail at length 11**, where the window can sit clear of the
ends — wrong rules, not thin ones. R-009.

| window | rules | lengths checked | confirmations |
|---|---|---|---|
| 5 | 2 | 6, 7, 8, 9 | 22 |
| 6 | 16 | 7, 8, 9, 10 | 22 |
| 7 | 18 | 8, 9 (+11, 12) | 3 (+14, 42) |
| 8 | 8 | 9, 10 (+11, 12) | 3 (+2, 8) |
| 9 | 0 | 10, 11 | **all 30 failed at 11** |

The window-7 and window-8 survivors were then checked at lengths 11 and 12 as
well, since two lengths is the minimum that rules out an end effect and those had
only two: **52 checks, no failures** (14 and 42 confirmations at 11 and 12 for a
window of 7; 2 and 8 for a window of 8). All 44 are in `VERIFIED_MOVES` now,
taking the listed table from 52 rules to 96 and the table with families from 64
to 116.

**And the rule that mattered was not one of the 44.** Among them,
`(0:2) -> (3:2)` via three left mutations, next to E-010's one- and two-mutation
versions, is the third member of a family whose `d`-th member needs `d`
mutations — so discovery at any bounded depth sees only an initial segment of it.
Generating the family instead gives every member: F-020, and H-008 confirmed.
That is the return on this run, more than the 44.

Cost: roughly 30 s per (pattern, embedding) at `maxSteps=3`; the re-verification
is the expensive half, since a window of 8 wants length 12 and its 58786 LNAs.
Tests H-007, and H-007 bit back.

---

## E-010 — Whole-quiver discovery at length 8
*2026-09-13* · **34 rules**

`discoverMoves([8], maxSteps=2)`: 101 candidates, 85 not already known, 34
verified. All of window 6 — the width that first has room to sit clear of both
ends at that length.

Each was confirmed only 3 times at lengths 7–8, which is thin, so all 34 were
**re-verified over lengths 7 to 10**: all survived, 22 confirmations each. Keep
doing this for wide rules found at short lengths.

---

## E-009 — Whole-quiver discovery at length 7, three mutations
*2026-09-13* · **2 rules** — poor yield

`discoverMoves([7], maxSteps=3)`: 61 candidates, 45 new, **2** verified.

The yield is low because `describeLink` only admits a *local* rewrite — one whose
window contains every relation it touches — and at length 7 a three-mutation
sequence usually disturbs the whole quiver, so nothing recurs across positions.
This is the experiment that motivated interior embedding (H-007). **Do not repeat
at this length.**

---

## E-008 — Classification of n = 10
*2026-09-13, updated 2026-09-14* · **unfinished — resume it**

`classifyLength(10)`, several attempts, none yet complete. Furthest reached:
about 1900 of 4862 rows.

**Long runs do not survive.** Three separate causes, all worth knowing:

1. two attempts were killed by over-broad `pkill -f` patterns issued by the
   session itself — a pattern that also matches the shell issuing it kills the
   shell, and anything sharing its process group;
2. one died under `setsid` when the machine went away between sittings;
3. n = 10 takes hours, so any of the above is likely to happen at least once.

**So: resume rather than restart.** The table is written after every class
searched, and `classify.py --resume` continues from the existing CSV:

    python classify.py 10 --resume

Before assuming a long job is still running, check `ps` — a stalled row count
looks the same as a dead process.

---

## E-007 — Certificate propagation by vertex deletion
*2026-09-13* · **0 / 0 / 1 / 24 / 308**

`notPiecewiseHereditaryByDeletion` over every LNA of lengths 4 to 11 → F-012.
The zero below length 9 is the correctness check, not an absence of result.

---

## E-006 — Cospectral quipu enumeration to order 13
*2026-09-13* · **the collision map**

`quipuForms.cospectralQuipuGroups(n)` for n = 4..13, cross-checked against equal
Coxeter polynomials computed through each algebra's Cartan matrix for n = 4..11.
The two agree exactly. → F-010.

Seconds to run, no mutation search involved. `python classify.py <n> --collisions`.

---

## E-005 — Orbit verification under the move table
*2026-09-13* · **clean**

Every LNA of lengths 5 to 9: compute its orbit under the verified moves, apply
each recorded mutation sequence to check it reaches the class it claims, and check
the Coxeter polynomial is constant on the orbit. 1764 orbit members, zero
failures.

Run **before** trusting any change to the rule table — it is what caught R-005.

---

## E-004 — Reduction preserves the Cartan matrix
*2026-09-13* · **clean, after R-003**

Every legal mutation of depth ≤ 3 out of every LNA of lengths 5 to 8: 38095
reductions, zero changes. Exact and heuristic Cartan matrices agree throughout.
→ F-008.

First run gave 9 apparent failures; all were R-003, not the reduction.

---

## E-003 — Exact against heuristic Cartan matrix
*2026-09-13* · **agree everywhere tested**

All 624 LNAs of length ≤ 8, and all 8101 quivers reached by walking every legal
mutation of depth ≤ 3 out of all 188 LNAs of lengths 5 to 7. No disagreement, so
no published Coxeter polynomial moves. The shapes where the two models differ
(F-004) have not turned up in an LNA search.

---

## E-002 — Classifications of n = 5 to 9
*2026-09-12 – 2026-09-13* · **match the published table**

| n | LNAs | classes | time |
|---|---|---|---|
| 6 | 42 | 4 | ~13 s |
| 7 | 132 | 6 | ~37 s |
| 8 | 429 | 11 | ~4 min |
| 9 | 1430 | 20 | ~56 min |

n = 6, 7, 8 match arXiv:2305.06642 exactly. n = 9 → F-011. Lengths 6–8 are pinned
as `slow` tests.

---

## E-001 — Reproducing the published n ≤ 8 classification
*2026-09-12* · **the baseline**

The first cross-check of the restored code against the papers: relation-set counts
against the Catalan numbers, the worked example of arXiv:2112.08129 step by step,
Coxeter polynomials of A_n and D_n, and the class membership of the n ≤ 8 table.
Everything agreed once F-001 was fixed.
