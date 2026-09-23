# Hypotheses

Things suspected but not established, newest first, each with what would settle
it. Status is one of `OPEN`, `SUPPORTED`, `CONFIRMED → F-nnn`, `REFUTED → R-nnn`,
`PARKED`. See [`README.md`](README.md).

---

## H-021 — A core's slide is a palindrome when its class is self-dual, and the head/tail difference is the reflection's shortfall
*2026-09-23* · **OPEN** *(eight cores pair their offsets by a reflection at `n = 13` and 14, `3346` does not -- E-052, F-053)*

**The conjecture.** For a single-cluster core `c` there is a centre `s(c)` with
`s(c) = n - k(c)` such that, under the reduced walk, `c@o` and `c@(s(c) - o)`
are always in one orbit, **exactly when** the orbit of some placement of `c`
holds the mirror of a placement of `c`. Then the slide of `c` read on offsets
`0 .. s(c)` is a palindrome, and the offsets beyond `s(c)` -- one for `45`,
where `s = n - 8` and the last offset is `n - 7` -- are what make the head and
the tail of H-020 differ.

**Why it is worth asking.** H-020 says the head and the tail are two numbers per
core; this would make them one number and a shift, and would say which cores
can have a head that differs from their tail by more than the shift. It would
also halve every census of such a core at no risk: half its offsets are the
other half's reflections, which the shared walk already exploits when an orbit
closes and now, since E-051, when one caps.

**What would settle it.** For every single-cluster core of `--max-word 4` at
`n = 14`: walk each offset under the reduced walk to closure (a few seconds
each since E-051), record which offsets and which mirrors each orbit holds, and
check the two sides of "exactly when". `3346` is the known case without a
pairing and should have no orbit holding a mirror of its own placements; a core
with a mirror in its orbit and no reflection pairing refutes it.

**What it would not explain.** Cores with a different centre for their inside
and outside offsets, and the two-cluster words, where each cluster has its own
distance to its own end (H-020's third amendment).

---

## H-020 — Where a core may sit is fixed by its distance to the two ends, not by the length
*2026-09-20* · **SUPPORTED** *(for a single heavy cluster, 283 cores at `n = 13` to 17 without an exception under the plain walk -- E-046; under the reduced walk 1186 of 1192 comparisons from 13 to 18, the six failures all words whose slide at 13 is one or two offsets -- E-051; `45` closed at every offset to `n = 17` -- E-052; two clusters do not obey it -- E-047)*

**The conjecture.** For an overlapping core placed alone in a line, whether the
moves carry it to an almost separate LNA depends only on how far it sits from
the source and how far from the sink. There are two numbers `head(c)` and
`tail(c)`, belonging to the core and not to the length, such that the core at
offset `k` in a line of `n` is placeable exactly when `k < head(c)` or the core
ends within `tail(c)` of the sink, and is outside everywhere between.

**Where it comes from.** F-042 found one instance of this by hand: the `45` core
is placed against the source or within one arrow of the sink and nowhere else,
at every length from 9 to 13, and the law is invisible at `n = 9` where all
three offsets are inside. E-045's census is the first run to ask it of many
cores at once. Of the 62 cores it completed at both `n = 13` and `n = 15`, the
25 that have an `outside` anywhere all read as `i…i o…o i…i` with the head and
the tail **the same at both lengths** and the outside middle taking up both
extra offsets. Not one has an inside in its interior.

**Why that is not yet a finding.** It is two lengths, and the shorter of the two
is the one F-042 already covered, so `n = 15` is carrying the claim alone. Both
censuses stopped a third of the way through their catalogue and both stopped at
the same place in it, so the 62 cores are the *front* of the catalogue — short
words, small overlaps — and not a sample of it. Six placements at `n = 15` are
undecided, and five of them sit exactly on a head or tail boundary, which is
precisely where a wrong `head` or `tail` would come from. And `outside` is a
statement about this move set, not about derived equivalence.

**What would settle it.** The same census at a run of lengths, complete rather
than truncated, which E-045's 486x now makes a night's work rather than a
month's:

```
python overnight.py --hours 9   --run "batch.py cores 12 --max-word 4 --jobs 4"   --run "batch.py cores 14 --max-word 4 --jobs 5"   --run "batch.py cores 16 --max-word 4 --jobs 5"
python batch.py cores 14 --max-word 4 --summary
```

The summary prints `head` and `tail` per core. The hypothesis says those two
columns are identical at every length; a single core whose head or tail moves
with `n` refutes it as stated, and a single core with an inside in its interior
refutes the shape. Resolving the six undecideds is worth doing first, since they
sit where the answer changes:

```
python batch.py cores 15 --max-word 4 --join-limit 40000   --cores 245,2045,2245,2555,2556,3344 --jobs 6
```

**What it would mean if it held.** Placeability would be a boundary condition
and nothing else — the interior of a long line would be uniform, and the whole
question would reduce to two finite numbers per core. It would also say what is
*not* happening: no interaction between the core and the length, and nothing
that appears only at some particular `n`. H-018 asks the same question the other
way round, in terms of overlap rather than placement.

**2026-09-22: asked at a run of lengths, and it holds for one cluster** (E-046).
The censuses at `n = 11, 12, 14, 16` are complete and `13, 15, 17` are partial.
Every single-cluster core decided at two lengths from 13 to 17 was compared
between them: **975 comparisons over 283 cores, no failure**. `45` is head 1,
tail 2 at every length from 11 to 17 and `504` head 2, tail 1. The six
undecideds this entry said to resolve first are all outside at three times the
limit, with the heads and tails unchanged. Three amendments to the statement:

* **The ends are words, not counts.** `4056` is `oio`, `oooio`, `oooooio` at
  12, 14, 16: inside one offset from the sink and outside at the sink itself.
  Its suffix `io` does not move with `n`, which is the conjecture's substance,
  but "inside exactly when `k < head(c)` or within `tail(c)` of the sink" is too
  narrow a way to say it. The statement that survives: *there are two fixed words
  `P(c)` and `S(c)` and a verdict `v(c)` such that the slide at every long
  enough `n` is `P(c) v(c)^m S(c)`.* Every mixed single-cluster slide seen has
  `v(c) = o`; none is inside in the interior and outside at an end.
* **"Long enough" is about 13.** Twelve single-cluster cores change shape
  between 11 and 12 or 12 and 14, all of them words long enough that their
  slide at the shorter length is three offsets or fewer and every placement
  is near both ends.
* **One cluster only.** Pair words break it 31 times, as they must: each of
  two clusters has its own distance to its own end, and the one far from the
  sink drifts further away as the line grows (`3500035`: `iio` at 14, `ooooo`
  at 16). What two clusters obey is H-018's question, answered in part by E-047.

**2026-09-22: the census above was measured with a one-way free move** (F-052,
E-049). Walked with the relation added as well as deleted, 60 placements at
`n = 11` and 12 move from outside to inside. `45` and `504` do not move, but
`4056` becomes `oii` at 12 (outside at the source, inside in the interior),
which is the shape the amended statement says is never seen, and `404`, `405`,
`5004` lose their outside bands altogether at those lengths. Twelve is below
the "long enough" of the second amendment, so this refutes nothing yet, but the
975 comparisons are comparisons of the plain walk and the law has to be asked
again of the reduced one, at 13 and above, before it is read as a property of
the moves:

```
python batch.py cores 13 --max-word 4 --walk reduced --jobs 4
python batch.py cores 14 --max-word 4 --walk reduced --jobs 4
```

F-051 is the mechanism, seen from the orbits: the outside interior of `45` at
`n = 14` is one closed orbit, containing offsets 1 and 5, so it can only have
one verdict.

**What would settle it now.** A proof, not more lengths: the moves carry a
cluster through the interior of the line (F-051), so the claim is that the rule
table acts the same way on every interior position, and only the anchored and
edge moves see the ends. The two partial lengths are worth finishing for the
record -- `cores 13` and `cores 15` at `--max-word 4` are a few core-hours each
-- and `n = 18` is the length at which a counterexample would first have room
beyond anything checked.

**2026-09-23: asked again under the reduced walk** (E-051, E-052). The shared
censuses at `n = 13` and 14 are complete, and 15 to 18 are partial. No
single-cluster slide at 13 or 14 is inside in its interior. Every pair of
complete slides from 13 to 18 was compared: 1186 of 1192 hold, and the six that
do not are at `n = 13`, for words of six or seven letters whose slide there has
one or two offsets -- the "long enough" of the second amendment is better read
as *three offsets or more* than as a length. `45` walked to closure at every
offset from 12 to 17 is `i o^(n-9) i i`, head 1, tail 2, as under the plain walk.
F-053 gives a mechanism for part of it: the outside band of `45` is one closed
orbit per **reflected pair** of offsets, so the slide is a palindrome up to one
offset at the sink; H-021 asks whether that is the general shape.

**What would settle it now.** The census at 15 to 18 under the reduced walk,
finished, which E-051's faster walk and larger cap make a night's work.

---

## H-019 — The leftovers do not thin out with the length, and past `n = 12` they are most of it
*2026-09-19* · **SUPPORTED**

**The conjecture.** The fraction of LNAs that neither the quipu theorem names nor
the moves carry to one it names does not tend to 0 with the length. It rises,
and the classification's present machinery covers a share of `A_n` that goes to
nothing.

**Why it is worth stating as a hypothesis rather than assuming either way.**
Every complete answer the repo has is at a length where **every vertex is within
three arrows of an end** — `n = 8` is the last length classified without a
search — and the anchored rules are exactly the ones an end is in reach of. So
the short lengths are not a small version of the long ones, they are the case
where the mechanism that does the work is always available. Nothing measured so
far distinguishes "the leftovers are a boundary effect that thins out" from "the
leftovers are the generic case and the short lengths are the exception".

**What is already known, exhaustively** (F-032, by coverage rather than by
sampling): the moves with no rule table leave 0.6% of `n = 9`, 5.4% of `n = 10`
and 16% of `n = 11` unplaced. So the series so far is 0.6, 5.4, 16, 28 — rising,
and rising faster than linearly in the first three steps.

**What a first sample says** (E-044). `batch.py sample`, uniform over all the
LNAs of the length: **28% +- 5.8 leftovers at `n = 12`**, from 60 draws. The
sampler is calibrated where the answer is known — 400 draws at `n = 10` give
4.8% +- 1.1 against the exhaustive 5.4%, and an exhaustive pass at `n = 9` gives
0.70% — so the `n = 12` figure is a measurement and not a guess. Preliminary
nonetheless: 60 draws, one seed. The orbit walk is forward only and capped, so a
row it fails to place may still be placeable and every rate here is an **upper**
bound on the true leftover fraction.

**What would settle it.** Samples of a few thousand at `n = 12, 14, 16, 18, 20`,
which is what the `sample` task exists to run and what a machine can be left
doing:

```
python batch.py sample 16 --count 4000 --jobs 7 --budget-hours 9
python batch.py sample 16 --summary
```

The rate is the first thing to read off. The second, and the reason each
leftover is recorded with its overlap profile and relation count, is **whether
long leftovers take shapes the short lengths have no room for** — a
configuration needing, say, six clear arrows on both sides cannot exist below
`n = 14` at all, and if the shapes at `n = 16` are the `n = 11` ones scaled up
then the mechanism is understood and only the coverage is missing. If they are
not, there is something new to find, and H-018's "pushed to an end" is the
reading to test it against.

**What would refute it.** A rate that levels off or falls between `n = 14` and
`n = 20`. That is a real possibility and not a formality: the moves get more
room to act as the quiver lengthens, and the double mutation of F-032 in
particular carries every relation crossing the one it slides.

**2026-09-20: `n = 15` is set up and not yet run**, with a mutation search out
of every leftover, which no sampling run has done before -- `--depth` was added
with the task and E-044 used depth 0:

```
python batch.py sample 15 --count 6000 --depth 4 --jobs 7
python batch.py sample 15 --count 6000 --depth 4 --summary
```

A pilot of 12 draws timed while setting it up came out 1 by the theorem, 6 by
the moves, **5 leftover**, which is the right order for the series 0.7, 5.4, 16,
28 to be continuing but is 12 draws and means nothing on its own. The costs, for
whoever sizes the next one: a draw is about 52 s of one core, ranging 0 to 150,
and a depth-4 deduplicated search out of a leftover is 27 to 65 s.

The search records `reached`, `hereditary` and the walk's dedup ratio per
leftover. That is H-017's census -- relations against cords over what each
leftover reaches -- asked at a length where F-040's caveat about single clusters
near an end may not hold.

**2026-09-20: 879 of the 6000 draws are in, and the number they give is not yet
comparable with the series above** (E-045). The raw split is 7.7% theorem, 34.2%
moves, **58.0% +- 1.7 leftover**. Before that is read against 0.7, 5.4, 16, 28,
two things have to be settled, and neither is a matter of drawing more:

* **One leftover in six was the cap and not the moves.** 81 of the 510 had
  walks that hit the 20000-row limit rather than closing, so the figure is an
  upper bound with a known and sizeable slack in it, and the slack grows with
  the length -- which is exactly the direction that would manufacture a rising
  series out of a flat one. `probe` now records `orbitClosed` and the summary
  splits the two. The number to compare with the exhaustive figures is the
  closed one; the way to shrink the gap is a larger `--orbit-limit`, which is
  now in the ledger's name so the two runs cannot be mixed.
* **The earlier points in the series were measured a different way.** 0.6, 5.4
  and 16 are F-032's exhaustive coverage counts, and 28 is 60 sampled draws at
  `n = 12`. A clean series wants the same instrument at every length, which is
  cheap now that the walk stops early: `--depth 0` is about a millisecond a draw
  below `n = 13`, so `n = 12, 13, 14` can be re-measured by sampling in minutes
  and compared with the exhaustive answers they already have.

The cost note above is superseded. A draw at `n = 15` with `--depth 4` cost 257
s of one core over the 879, which is five times the depth-0 pilot estimate --
the depth-4 search runs on every leftover, and leftovers are most of the draws
at this length. Depth 0 and depth 4 are worth running as two different jobs at
two different counts rather than one.

**2026-09-22: `n = 13` and `n = 17` by the same instrument** (E-048). `n = 13`,
4575 draws at depth 0: **39.9% +- 0.7 leftover, every one of them a closed
orbit**, so no cap in the number at all. The series by one definition is 0.7,
5.4, 16, 28, 39.9 -- still rising, by about twelve points a length. `n = 17`, 214
draws: 47.2% +- 3.4 leftover with a closed orbit and 22.4% more whose walk hit
the cap; E-045's `n = 15` rows split the same way give 48.8% closed and 9.2%
capped.

So the rise is established through `n = 13` and **not above it**. The closed
fraction, a lower bound, is flat from 15 to 17; the capped fraction, which could
land on either side, more than doubles. "Most of it past `n = 12`" is true of the
upper bound from `n = 15` and not yet of the lower one anywhere. The one run that
separates the readings is the same draws at a larger walk:

```
python batch.py sample 15 --count 20000 --orbit-limit 100000 --jobs 4
python batch.py sample 15 --count 20000 --orbit-limit 100000 --summary
```

and the same at `n = 17` if `n = 15`'s capped rows turn out to close. A depth-0
draw is 14 s at `n = 13` and 152 s at `n = 17`, not the millisecond the entry
above guessed -- the closing walk is the whole cost.

**2026-09-23: the capped share at `n = 15` was all outside** (E-051). The same
draws at `--orbit-limit 100000`, 1444 of them: **57.4% +- 1.3** leftover, and
every leftover orbit closed -- the largest at 60852 rows. E-045's 48.8% closed
plus 9.2% capped was 58.0%, so the cap was hiding no placements at all. The
series by one instrument, closed orbits only, is now 0.7, 5.4, 16, 28, 39.9, --,
**57.4** at `n = 11` to 15: still rising, and "most of it past `n = 12`" holds
from 15. `n = 17` (1131 draws at 20000) is 49.6% closed and 22.9% capped; if
its capped draws close as `n = 15`'s did, it is 72.5%. Read as the fraction the
moves *place*, 84, 72, 60, 43 and perhaps 27 percent from 11 to 17: a fall of
about a sixth per unit of length, steady enough to be worth fitting once `n = 17`
is clean.

**What would settle it now.** `sample 17` at `--orbit-limit 500000`, which
E-051's faster walk makes a few seconds a draw; and the same at 15 under
`--walk shared`, since every number in this series is the plain walk's and the
reduced walk places more (E-049).

---

## H-018 — What escapes the quipu theorem is a placement, not an overlap
*2026-09-18* · **OPEN**

*Renumbered at merge from `H-016`, which was taken on `main` first by an unrelated entry while this branch was open. Session logs and commit messages from the branch use the old identifier.*

F-042 measured that no bound on the overlap separates the LNAs the quipu theorem
and the move table place from the ones they do not: the smallest step past
"almost separate" -- one pair of relations sharing two arrows -- already contains
outsiders at every length from 9 up. What the same measurement suggests instead:

**The conjecture.** An LNA is carried to an almost separate one by the moves iff
every one of its heavy clusters can be **pushed to an end**. The evidence is that
the cores with room -- `0^a 4 5 0^b` and its opposite `0^a 5 0 4 0^b` -- are
outside exactly when `a >= 1` and `b >= 2`, that is, exactly when neither end is
within reach, while `33`, `44`, `34`, `43`, `54`, `55` and `333` are inside at
every placement because each of those *can* be pushed out.

**What is odd about it, and is the part to explain.** `55` overlaps in four arrows
and is always inside; `45` overlaps in three and is not. Whatever "can be pushed
to an end" means precisely, it is not monotone in the overlap, and the only
structural handle so far is that the two escaping cores are opposite algebras of
each other.

**What would settle it.** A census of every core at `n = 14` and `n = 15` -- if a
third family appears, this phrasing is already too narrow. A move-by-move account
of why `55` clears and `45` does not; whichever move does the clearing should say
what the obstruction is. And the sharper version: whether the escapees are outside
the *classification* or only outside this move set, which the quiver-level search
can answer where `orbitOf` cannot.

**Why it matters more than the pattern asked for.** If placement is the
coordinate, then a generalisation of the quipu theorem cannot be a condition on
the relation profile alone, and the short quivers could never have shown this --
at `n = 9` every placement of the `45` core is inside. F-042, E-037.

**2026-09-20: the census this asks for now has a command, and has not been run.**
`batch.py cores` slides every heavily overlapping core word along a line at every
offset, which is F-042's hand slide as a resumable job. Written against this
entry and against F-040's count of separated clusters:

```
python overnight.py --hours 9 --only sample15 cores15 cores13
python batch.py cores 15 --max-word 4 --summary
python batch.py cores 13 --max-word 4 --summary
```

Three things about it are worth stating before the answer is in, so that they
are not read into the answer afterwards.

* **`n = 13` is run as a control, and is the more important half.** A census at
  one length is a list of rows; the law in F-042 is what *changes* between
  lengths, and the table there is five rows for one core. The two runs give the
  same table for every core at two lengths.
* **A verdict of "outside" here means the forward move orbit closed**, which is
  a statement about this move set and not about derived equivalence -- the same
  caveat `orbitOf` has always carried. A placement whose orbit hit its cap comes
  back `undecided` and not `outside`, which is the E-037 lesson built into the
  instrument: 49 barricades were once recorded as failures by a walk that had
  run out of rows, and `movesJoin` joined them in 45 seconds.
* **What would already be new.** F-040 found no LNA outside a quipu class with
  two heavy clusters separated by a free arrow, at any length it could count,
  and said the shape barely fits below `n = 12`. The summary counts exactly
  those rows. Either they are all placed, which extends F-040's count to a
  length with room, or one is not, and then the single-cluster reading that
  every finding since F-040 rests on is too narrow.

The instrument is checked against F-042's published slide at lengths 9 to 12 and
against the mirrored `504` slide, in `tests/test_cores_task.py`.

**2026-09-20: it was run, it was cut off a third of the way through both
lengths, and none of the three things above got an answer** (E-045). What came
back is 62 cores at both lengths, from the front of the catalogue -- short words
and small overlaps -- which is the part of it F-042 had already looked at. In
particular the count this entry says would already be new got **no units at
all**: the two-cluster words are appended to the catalogue after every single
core, and neither length reached them. They need a run of their own, and
`--max-word 2 --pair-word 2` is it.

The partial table does say one thing, and H-020 is that thing written down
separately rather than folded in here, since it is a different shape of claim:
the head and the tail of every slide that has an outside in it are the same at
`n = 13` as at `n = 15`. This entry asks *which* configurations escape; H-020
asks where an escaping one may sit.

**2026-09-22: the two-cluster shapes, run at `n = 14` to 17, and the count
this entry said would already be new** (E-047). Of 434 rows with two heavy
clusters and a free arrow between them, **32 are outside** -- so the shape is
not all placed, and the first outsider is at `n = 14`: `330004500000`. But
looking each half up alone at the same offset of the same length:

* **no outside pair has two placeable halves** -- 364 pairs whose halves are
  each inside alone are all inside, and all 32 outsiders have a half that is
  outside alone (in `330004500000` it is the `45`, sitting in its own outside
  band);
* **the halves are not independent**: 24 pairs are inside although one half is
  outside alone -- a `33`, `34` or `44` at the source end, three or four arrows
  before a `35` or `36`, carries the pair in;
* **without a free arrow it goes the other way too**: eight pairs touching
  through a one-arrow overlap (`35 0^g xx` at offset 1) are outside with both
  halves inside alone.

For the conjecture as stated, "every heavy cluster can be pushed to an end":
the direction "some cluster cannot, so the LNA is outside" fails if "cannot" is
read as "cannot on its own" -- the 24 are counterexamples to that reading. The
other direction survives everything run. The single-cluster reading that the
findings since F-040 rest on survives in the weaker form: at these lengths a
separated two-cluster LNA is outside only where one of its clusters would be.

What E-047 could not ask: whether the rescue survives a wider gap (it reaches
three offsets at gap 3 and two at gap 4, at `n = 17`), and anything about
`--pair-word 3`, whose run spent its catalogue cut on gap-1 pairs that can never
have a free arrow. Both want `--gaps 5,6`.

**2026-09-23: the rescues were the plain walk's** (E-051). With a free arrow
between the clusters (gaps 5 and 6 at `n = 17`, complete, plain walk) two
placeable halves made a placeable pair 330 times in 330. The 62 pairs rescued
from an outside half all had a `35` or a `36` as that half, and the reduced walk
places `35` and `36` at every offset from 13 to 18. Under the shared walk at 17
and 18 (partial) no half is rescued: every pair that is not inside has a half
that is outside or undecided alone. So far the two-cluster verdict is the worse
of its halves' and nothing else, which is the answer this hypothesis predicts
for clusters that are genuinely apart.

---

## H-017 — A quipu in the class always carries more relations than it has cords
*2026-09-17* · **OPEN**

*Renumbered at merge from `H-015`, which was taken on `main` first by an unrelated entry while this branch was open. Session logs and commit messages from the branch use the old identifier.*

F-034 says every LNA outside a quipu class reaches a quipu quiver *with*
relations. The question this asks is what those quipus look like, because a
theorem in the shape of `thm:QuipuToAn` needs a normal form and not a census.

**The observation.** Walking every LNA outside a quipu class at `n = 9` to depth
5 and grouping what it reaches by (cords, relations):

| LNA | excess overlap | (cords, relations) reached |
|---|---|---|
| `3033030` | 1 | (1,3) (1,4) (1,5) (2,4) (2,5) (2,6) (3,5) |
| `3345000` | 4 | (1,2) (1,3) (1,4) (1,5) (2,3) (2,4) (2,5) (2,6) (3,5) |
| `3505000` | 3 | (1,2) (1,3) (1,4) (1,5) (2,3) (2,4) (2,5) (2,6) (3,4) (3,5) |
| `4444400` | 8 | (1,3) (1,4) (1,5) (2,4) (2,5) (2,6) |

and so on for the rest. Over all nine LNAs and every member reached, **the number
of relations exceeds the number of cords, every time**: the pairs seen are
(1,2) (1,3) (1,4) (1,5) (2,3) (2,4) (2,5) (2,6) (3,4) (3,5), and never (1,1),
(2,2) or anything below the diagonal.

**Why that would be the right shape.** A quipu class is exactly the case where
the relations can all be turned into cords, so its members include one with
relations = 0. Reading `relations - cords` as a **defect**, the theorem's case is
defect `<= 0` and everything this project cannot classify has defect `>= 1`. The
minimum over a class is a class invariant by construction; at `n = 9` it is 2 for
`3033030` and 1 for the eight-member class, so it is not merely "not zero" -- it
separates the two classes the theorem misses.

**A sharp prediction, and the reason to care.** The polynomial enumeration of
F-034 offers candidates *below* the diagonal -- `3033030`'s polynomial is carried
by `P^(1,1,1)_(1,0,1,1)` with three cords and a **single** relation. If the
observation is a law, that candidate is not in the class and the polynomial match
is a coincidence. That is a falsifiable statement about a specific algebra, and
settling it either finds the normal form or kills the pattern.

**What would settle it.** Deeper walks -- depth 6 and 7 at `n = 9`, where the
depth-5 counts have stopped moving for most rows -- and the same census at
`n = 10` and `n = 11`. A member below the diagonal refutes it outright. A proof
would want an invariant that counts relations against cords; the Euler form is
the obvious place to look, since for a tree quiver of global dimension 2 the
number of relations is read off it, and what breaks that here is exactly the
higher `Ext` the overlapping relations create.

**The caveat that applies to all of this.** F-040: at these lengths every LNA
outside a quipu class is a *single* overlapping cluster, usually against an end.
A normal form fitted to those may say nothing about a quiver long enough to hold
two clusters far from both ends.

---

## H-016 — Walking through a parallel-arrow quiver reaches a merge nothing else does
*2026-09-18* · **OPEN** *(no gain at n ≤ 7 to depth 6, which is where a gain could not show anyway — E-035; at n = 9 one member of nine reaches the region at all, and it reaches nothing new nine mutations into it — E-036)*

The procedure produces quivers with parallel arrows and, since F-039, the engine
can state them, the gate admits them and the Coxeter key over them is right. So
there is a region of the mutation graph that no search has ever entered. **Is
anything in it?**

The reason to think so is the shape of the one walk looked at by hand. From
`3030` at `n = 6` by `[1, 3, 4, 1, 4]` the quiver has two arrows `1 -> 6` and a
commutativity relation between the two parallel paths `5 -> 1 -> 6`; mutating at
3 comes back out to a quiver with **no** parallel pair carrying
`5 -> 1 -> 3 -> 6 = 5 -> 1 -> 6`. That quiver was unreachable at any depth
before. Whether such an exit ever lands on a *line* — which is what a merge
needs — is the question.

**What is measured.** Nothing gained at `n = 6` or `n = 7` to depth 6: the same
LNAs are reached with the gate allowing parallel arrows and with it refusing
them, from every LNA and every relation dual, not one start gaining or losing a
line (E-035).

**That is a weak negative and should not be read as an answer**, for two reasons
that are both about where it was measured.

* **Depth.** Entering the region and returning costs mutations. The one exit
  known takes six to reach a parallel-free quiver and that quiver is not a line,
  so a depth-6 comparison cannot see a return to a line at all. The test needs
  depth 8 or more, which is why it was not run here.
* **Length.** `n <= 7` is covered 100% by the move rules with no search at all
  (F-021). There is nothing at those lengths for a search to find, whatever it
  walks through. The lengths where a search still has to place rows are `n >= 9`,
  and `n = 10`, `n = 11` are where H-013's leftover orbits sit.

**What would settle it.** Run the comparison at `n = 9` or `n = 10` and depth 8,
against the orbits `merges.py` leaves as singletons — the same targets as H-013.
A single LNA pair joined only through a parallel-arrow node settles it yes; a
clean negative at that depth and length is worth having either way, because it
would say the region is a detour rather than a shortcut and the classification
need never enter it.

**Do not run it at `n <= 8` again.** E-035 is that run and it found nothing, for
reasons that are about the sizes and not about the region.

**2026-09-18, amended: the test is much cheaper than it looked, and the first
part of it is done.** `search.DeeperWhen` gives the extra depth only to the
branches that reach the region, so the comparison does not cost a whole extra
level of search. At `n = 9` it costs 2% of the run at depth 4 and 4% at depth 5,
because **exactly one of the nine leftover members reaches a parallel-arrow
quiver at all** -- `3033030`, whose walk fires 600 times at depth 5 where the
other eight fire never (E-036). Nothing is gained at either depth.

That narrows the hypothesis rather than answering it. `3033030` is alone in its
orbit and alone in its Coxeter polynomial group, so the guard forbids it reaching
any other leftover: the only outcome visible from it is reaching a **seeded** LNA,
which would say a leftover is in a quipu class after all. And because it is the
only member that enters the region, the depth-8 run this hypothesis asks for is a
run of *one* member, not of nine -- which is what makes it affordable. `n = 10`
and `n = 11` have not been looked at this way and are where the leftover orbits
are that H-013 cares about.

**That run has since been made, and it is a negative.** `3033030` from itself and
its dual at depth 7 with two extra mutations for the region: 29122 firings, 826
grants, the condition still holding **nine** mutations in, and the only LNA
reached is its own relation dual -- which depth 5 already reached. Twenty minutes
on one core. Depth 6 + 2 is the same answer (E-036).

So the `n = 9` half is done, past the depth this hypothesis asked for, and the
region is a detour there. What it does **not** do is settle the hypothesis:
`3033030` is alone in its Coxeter polynomial group, so the only outcome it could
ever have shown is a leftover turning out to be in a quipu class, and one member
at one length is not the claim. **What is left is `n = 10` and `n = 11`**, where
H-013's leftover orbits sit several to a polynomial group and a merge between two
of them is a thing the search can actually find. Find which of their members
reach the region first -- if it is again a handful, the run is again cheap.

---

## H-015 — The Coxeter guard is sufficient, not merely necessary
*2026-09-18* · **SUPPORTED** *(survives the sharpest test available, at one collision, to depth 6 — E-034)*

**2026-09-18, first test.** The two cospectral quipus of order 9 — `3060000` and
`3004000`, different trees and one polynomial, so the guard is blind between them
— are **not linked**, at depth 5 or 6, with the guard on or off. `3060000`
reaches only its own dual; `3004000` reaches eight LNAs; the two sets are
disjoint. E-034. This is the place a guard-passing non-equivalence would show
first, and it does not show.

F-038 fixed the search by requiring that every step keep the Coxeter polynomial.
That is a *necessary* condition for a derived equivalence, and the whole repo now
leans on it as though it were sufficient: a path the guarded search finds is
taken as proof that two algebras are in one class (F-037 is exactly such a
proof, and every merge of E-032 rests on it).

**The suspicion is that within this setting it is sufficient** — that an
admissible mutation which holds the Coxeter polynomial fixed really is a tilting
mutation here, so the guarded search proves what it claims.

**Why it is not obvious, and might well be false.** The polynomial is not a
complete invariant, and F-010 says precisely where it fails: cospectral quipus,
the first pair at order 9. If two algebras can be cospectral without being
derived equivalent, then in principle a single mutation could step between them,
hold the key, and pass the guard. Nothing rules that out. The guard removes the
600-odd failures per `n = 7` search tree that F-038 measured; it does not prove
there are none left.

**Evidence for, such as it is.** The condition is the same one R-005 imposed on
rule verification, and under it 16 of 67 candidate rules survived and none has
since been found false. Every one of the five links of E-032 was replayed
step by step under it, and independently each pair is already known to share a
polynomial by a different route. And the guard changes no answer at `n = 6` or
`n = 7` — the corrupt region it removes never contributed a line at those sizes —
so it is not doing violence to results that were right.

**What would settle it.** A second invariant applied along a guarded path:
`τ`-periodicity data, Hochschild cohomology, or the Avella-Alaminos–Geiß
invariant where it applies (R-008 says it does not apply directly here). Cheaper
and worth doing first: take every step the guard admits at `n = 6` and `n = 7`,
and check whether the two algebras have the same *hereditary form* where both
reach one — F-036 makes that a mutation invariant, and it is independent of the
polynomial. A step passing the guard and changing the hereditary form would
refute this outright.

**Why it matters.** If it is false, then the guard is an improvement and not a
fix, and every positive claim the search makes — F-037, the merges of E-032,
F-034's walks — needs a second invariant before it can be believed. If it is
true, the search is sound and can be trusted at depth, which is where the
remaining classification questions live.

---

## H-014 — Every class outside the quipu theorem has a quipu-with-relations member, and one of them is canonical
*2026-09-17* · **SUPPORTED**

F-034 establishes two things at `n = 9`, `n = 10` and `n = 11`: every Coxeter
polynomial carried by an LNA in no quipu class is also carried by quipu algebras
with relations, and at `n = 9` and `n = 10` a mutation walk proves that every
such LNA really does reach one, within three mutations. The hypothesis is the
general statement, in two parts.

**Part 1, the existence.** Every derived equivalence class of LNAs contains an
algebra whose quiver is a quipu — with relations where the theorem's quipus have
none. The quipu theorem is then the case where the relations can be cleared away
entirely.

**Part 2, the normal form, which is what would make it a theorem.**
`thm:QuipuToAn` is useful because it is a *bijection*: one quipu per class, named
by parameters read off the LNA. What F-034 has is the opposite — 1746 quipu
algebras on one class at `n = 9`, over 16 of the order's 18 quipu shapes. A
theorem in the shape of the quipu theorem needs a canonical one, and the
hypothesis is that some rule (fewest relations? shortest total relation length? a
particular orientation?) picks it out, and that its parameters are a function of
the LNA the way `k` and `m` are.

**Part 2 has a shape already, at `n = 9`.** `python families.py members 9
--depth 3` prints what the walks actually reach, and what they reach first is one
quipu over and over: **`P^(6)_(1,1)`**, the line on eight vertices with a single
pendant vertex at the second — one mutation's worth of quipu. `3033030`, whose
relations are `1-2-3-4`, `3-4-5-6`, `4-5-6-7`, `6-7-8-9`, appears there as

    1 -> 2 -> ... -> 8  with  2 -> 9,  relations  2-3-4-5, 3-4-5-6, 5-6-7-8

— the relations carried along one vertex, one of them absorbed into the branch,
and the count down by one. Seven of the nine LNAs at `n = 9` reach that same
shape, the other two reach `P^(5)_(1,2)`. Whether the correspondence is a
function of the LNA, and what it does to the relations in general, is exactly
what part 2 is asking; the material to read it off is what `members` prints.

**Part 1 is now measured, not just suspected.** Every one of the 9 LNAs outside a
quipu class at `n = 9`, every one of the 262 at `n = 10`, and every one of a
random 200 of the 2647 at `n = 11` reaches a quipu with relations **within three
mutations** — none reaches none, and the median LNA reaches 16, 20 and 27 of them
(F-034). What is open in part 1 is the rest of `n = 11` and beyond, and whether
it is a theorem rather than a run of small cases.

**Evidence for.** The counts of F-034 and the walks behind them; and the fact
that the phenomenon is old and small — `D_4` with one two-arrow relation is
`kA_4` after a single mutation (F-035), so quipu-with-relations to line is a
mechanism the procedure performs routinely rather than a coincidence of the
polynomial.

**What is in the way.** A polynomial match is necessary and not sufficient, and
at `n = 10` one of the seven polynomials -- `T^10 + T^9 + T + 1`, which
`34504030` and `50505000` carry -- is shared with a quipu class, so a match
against it says nothing at all about those two (the paper's `remark:Coxeter`
makes the same point). Those two are the `UNPLACED` rows, and they are exactly
the ones a polynomial can never settle.

**What would settle it.** For part 1: `quipuRelations.reachedQuipuAlgebras` from
every LNA outside a quipu class at `n = 11` and `n = 12` — `n = 10` is done and
takes 7 minutes, `n = 11` about an hour and a half. An LNA that reaches none at a depth where its neighbours reach
twenty is the interesting outcome and is where a counterexample would show. For part 2: take the confirmed members of one
class and look at what they have in common — the `--verify` path of
`families.py` produces them, and `classpage`-style drawing would make a family
visible faster than a table will.

**Why it matters.** It would be the first statement about the classes the quipu
theorem misses that is about *shape* rather than about search. Every tool the
project has for those classes at the moment is negative -- a certificate that an
algebra is in no quipu class -- and F-033 closes the hereditary route, so a
family with relations is what is left.

---

## H-013 — The leftover orbits sharing a polynomial are few classes, and the search can say which
*2026-09-17, written before the overnight run* · **SUPPORTED** *(2026-09-18: `n = 10` answered, two predictions right and one wrong; the search's own soundness is now the open question — E-032, E-033)*

**2026-09-18, what the run returned.** `n = 10` finished all four depths. Twelve
orbits fall to **at most 10** classes, so the derived classes at `n = 10` number
between **43 and 46**, not the 43–48 written below. `n = 11` got through depth 4
and most of depth 5: 54 orbits to at most 51, so between 84 and 115. The full
table of predictions against outcomes is in E-032.

The prediction that failed is the one to keep. `34504030` and `50505000` —
the pair the Coxeter polynomial can never separate, which this hypothesis called
singletons under every move known — **are linked by a depth-7 mutation path, from
both sides**. If that path is real, the invariant that could not separate them
did not need to: they are one class. The "what would settle it" line below said a
real separation needs an invariant beyond the polynomial; what this shows is that
the pair needing one may simply not exist.

**What the alarm clause got wrong.** It said an ALARM "would refute F-032's
orbits". Two fired, and they refute nothing about the orbits: each of the two
orbits involved carries exactly one Coxeter polynomial across all of its members.
The alarm is about the **search**, which is a possibility this hypothesis did not
consider at all. E-033.

**Still open.** A negative remains a depth bound rather than a separation, and
nothing here changes that.

The `break`-for-`continue` bug in `mutationSearchDepthFirst` (E-033) would have
made the depths weaker than they read, by abandoning the remaining vertices at a
node whenever one mutation there yields an illegal relation — but it **did not
fire once** in this run: `isIllegalRelation` prints when it triggers, and all
three logs contain zero such lines over 121 core-hours. The bug is latent here,
not active, so the depths below stand as searched.

**The original statement and predictions, unchanged:**

F-032 leaves, at `n = 10`, 12 orbits in 7 Coxeter-polynomial groups once the
orbits are closed under the relation dual (which is derived equivalence and
which neither the free move nor the double mutation performs — E-029 caught a
first search rediscovering only that). So the derived classes at `n = 10`
number between **36 + 7 = 43** and **36 + 12 = 48**. Three groups have more than
one orbit:

| polynomial | orbits (members) | certified not p.h. | prediction, before the run |
|---|---|---|---|
| `(λ-1)²(λ+1)²(λ²+1)(λ⁴+λ³+λ²+λ+1)`, `C(2,4,5)`'s | 2 (69, 42) | none | **merge**, at depth ≤ 7. Three weights have no parameters, so if both are piecewise hereditary they are the one canonical algebra |
| `(λ-1)²(λ+1)²(λ²+λ+1)(λ⁴-λ²+1)` | 4 (all tiny) | every member | **at least two stay apart** to the depth reached — their members barely move |
| `(λ+1)²(λ²-λ+1)(λ⁶-λ³+1)` = `T¹⁰+T⁹+T+1` | 2 (1, 1): `34504030`, `50505000` | neither, by our criteria; both, by the paper's τ-path | **no link** to depth 8 — both are singletons under every move known. Whether they are derived equivalent is open and the paper does not say |

At `n = 11`: 54 orbits in 20 groups, so between 84 and 118 derived classes.

**What would settle it.** `python merges.py 10 --depths 5 6 7 8` and
`python merges.py 11 --depths 4 5 6`. A link merges; its absence is only a depth
bound, and a real separation needs an invariant — `τ`-periodicity data or
Hochschild cohomology, not the Coxeter polynomial.

**An ALARM in either run** — a link between orbits of different polynomials —
would refute F-032's orbits, and would matter more than any merge.

---

## H-012 — Every free move is also a mutation equivalence
*2026-09-16* · **SUPPORTED** -- settled at `n = 8`, F-041

F-028 establishes that deleting a relation of two arrows keeps the *derived*
equivalence class. It says nothing about the mutation class, and the two are not
the same question. The suspicion is that they coincide here: that an LNA is
always mutation equivalent to its stripped form, so the free move is a shortcut
through the mutation graph rather than a step outside it.

**2026-09-17.** By the known mutation moves only — double mutation, edge moves
and the whole rule table — 8 of 429 LNAs at `n = 8` and 44 of 1430 at `n = 9`
are not joined to their stripped form (E-029). That is a gap in the known moves,
not a counterexample; the test this hypothesis asks for still needs the
mutation classes, i.e. a finished classification with its search paths.

**Evidence for.** Against an end it is a single mutation: E-027 found
`(l, …, 2) → (l)` by one left mutation at the sink, for every `l` it tried, and
F-029's collapse is the same thing with a spectator. And a two-arrow relation
alone in its window slides any distance (F-020's lone slide), so one that can
reach an end can be deleted there.

**Evidence against, or at least in the way.** Sliding needs room. The table
already shows a two-arrow relation that can be slid to the sink and still not
deleted there — `(3,0,0,0,2,0)` in `A_8` reaches `(3,0,0,0,0,2)` and never
`(3,0,0,0,0,0)` — because the collapse rule's window is not clean. F-029 removes
that particular obstruction, but only for a companion starting on the window's
last arrow. Whether every configuration can be cleared is open.

**2026-09-17, one route tried and closed.** If both an LNA and its strip reached
a relation-free quiver with the same underlying tree, F-036 would join them by an
explicit sequence and settle the pair. They do not: of the 8 gap rows at `n = 8`
and the 44 at `n = 9`, **not one reaches a relation-free quiver at depth 4**, on
either side. The bridge has nothing to work with here. E-031.

**2026-09-17, and then a route that works.** Asked one relation at a time -- the
whole strip is a composition of single deletions -- and settled by meeting in the
middle rather than by one search reaching the other: **every single two-arrow
deletion at `n = 8` is a mutation equivalence**, 562 by the known moves and the
last 10 by a meeting at 3 + 3 mutations. So the hypothesis holds outright at
`n = 8`. At `n = 9` it holds for 1989 of 2002 deletions and at `n = 10` for 6974
of 7072; of the 13 left at `n = 9`, eleven have both sides almost separate with
the same quipu, so only two are outside the theorem's reach as well. F-041.

**2026-09-18, the barricade.** The shape the hypothesis is doubted for -- two
heavy clusters walling a two-arrow relation in, with free arrows between -- built
at the lengths where it first fits. All 95 shapes at `n = 13` to `16`, six cluster
types on each side and gaps of one to three arrows: **the moves strip the
two-arrow relation out of every one of them**. 46 by a one-way walk and the other
49 by meeting in the middle (`freeMoves.movesJoin`), those 49 having been reported
as failures by a walk that had merely run out of budget. So the barricade does not
trap the relation at the lengths where it first exists. It remains untested at the
lengths the doubt was raised for -- four clusters at 30 to 50 vertices -- where
nothing can be enumerated. E-037.

**What would settle it.** For each `n` where the mutation classes are known,
check whether every LNA and its strip share one. A single pair that does not,
with the mutation classes verified, would be the more interesting outcome: it
would be a derived equivalence that is not a mutation equivalence, which the
classification has never yet had to handle.

**Why it matters either way.** If it holds, the free move can be used anywhere
the mutation class is what is wanted, and the rule table is simply missing rows.
If it fails, then `classification` has a real distinction to maintain and
`freeMoves`' warning is not a formality.

---

## H-011 — Every class is reached by walking a heavily overlapping run to an end
*2026-09-16* · **OPEN**

F-022 gives the one mechanism known to reduce the overlap of an isolated pair:
walk it to an end of the quiver with the pair slide and collapse it there. The
suspicion is that this is not one mechanism among several but **the** mechanism
-- that the derived equivalence class of any LNA is reached from an almost
separate one by a sequence of interior moves that carry its heavily overlapping
runs to an end, and collapses there.

**2026-09-17 — the mechanism is in the literature, and the residue is no longer
a gap (F-032).** `proposition:doubleMutation` of arXiv:2310.08346 *is* this
hypothesis's walk-and-collapse, with any number of bystanders allowed: an
interior double mutation slides a relation while carrying everything crossing
it, and the same move at `t = n` is the collapse. With the free move it places
every LNA in a quipu class at `n = 9, 10, 11`; what is left is provably not a
quipu class. The end still matters — the interior half alone reaches 770 of
1430 at `n = 9` against 1421 — so the hypothesis's *shape* is confirmed. Its
literal statement ("every class is reached from an almost separate one") was
always false for the non-quipu classes, which contain no almost separate LNA;
read it as quipu classes only. The sharp question below is overtaken: the rows
left at `n = 10` are not waiting for rules (190 of 262 have both ends occupied,
and they are non-quipu classes whatever their ends look like).

**Why it is worth stating.** If it holds, the classification needs no search at
all: the seeding, the slide families and the anchored collapses generate
everything, and n = 10 and beyond become a matter of counting rather than of
mutation. If it fails, the LNA it fails on is the first evidence of a genuinely
different obstruction, which is worth more than another rule.

**Evidence for, and it is now substantial.** Coverage at n = 9 has gone from 45%
to **60%** on nothing but rules anchored to an end -- 100% at n = 6, where a
classification needs no search at all any more (F-022, F-023). Every rule that
crosses the almost separate line does so at an end; not one floating rule
reduces the overlap of an isolated pair.

**What the first end-discovery run settled.** Of the two readings this
hypothesis offered for why the runs of three were stuck, the second is right.
At n = 8, of the 155 LNAs then unplaced, **126 had a rule whose left-hand
pattern was present and which did not fire**, because further relations were
sitting in its window; only 29 had no rule with their pattern at all. The rules
were not too narrow, they were too **clean** -- and 188 of the 229 anchored
rules now listed carry a bystander they step around (F-023).

**The residue was the search bound, as predicted, and raising it paid.**
E-024's 33 unreachable LNAs at n = 9 were almost all pairs with a long relation,
which no run had ever planted. With `--max-arrows 7 --max-width 8`, discovery
against the ends verified 3045 more rules and took n = 8 to 98% and n = 9 to
84%; measured for the first time, n = 10 is 63% and n = 11 is 47% (E-025). So
the mechanical reading keeps being the right one, and the residue keeps being
about the bounds rather than about a new obstruction.

**And the first widening batch behaved exactly as the hypothesis predicts.**
`discover.py --extend` verified 625 widened rules in eight minutes: coverage went
to **100% at n = 7**, 95% at n = 8 and 73% at n = 9, and the number this
hypothesis said to watch -- LNAs for which *no* rule has the pattern at all, the
ones more widening cannot fix -- fell from 29 to **5** at n = 8 and stands at
**33** of 392 at n = 9 (E-024).

**Two moves found from outside the search have now taken n = 8 to the end of it
(E-027).** With F-028's free move and F-029's end doubling added, `A_8` needs
**no search at all**: 21 orbits and nothing left over, where the rule table alone
left 10 rows. `A_9` falls from 380 orbits and 222 rows to **77 and 37**. That is
the strongest evidence this hypothesis has had, and it came from two mechanisms
the discovery runs could not express rather than from more searching. It also
sharpens what is left, and the residue at `n = 9` turns out to have one property
in common, which is the sharpest form this hypothesis has yet taken:

**All 37 have a relation at the source *and* a relation at the sink.** Every one
of them is also already reduced, so F-028 has nothing left to give them. If the
mechanism this hypothesis names is walking a heavily overlapping run to an end
and collapsing it there, then an LNA with both ends already occupied is exactly
the case where there is no end to walk to — and that is precisely, and only,
what is left. Just 1 of the 37 is certified non-piecewise-hereditary, so the
other 36 are a gap in the rules rather than algebras outside any quipu class.

The property is an `n = 9` fact and not yet a general one: at `n = 10`, 670 of
the 887 rows left have both ends occupied and 660 are reduced, so the
characterisation is strong but not complete there. Whether it becomes complete
once `n = 10` has the rules `n = 9` has is the question to ask next, and it is
the concrete form of this hypothesis to try to break.

**Evidence against, and it is weaker than it looks.** 392 rows at n = 9 are
still not placed. The 33 with no rule at all were the candidate for a fourth
configuration, and they are not one: 32 of them contain a relation of five
arrows or more, and discovery has only ever been run at `--max-arrows 5
--max-width 6`, where such a relation cannot appear beside another. They are a
search bound, not a phenomenon.

**What would settle it.** Keep re-running the blocked-rule diagnostic after each
batch and watch that count as a *share*. While it stays small the hypothesis is
holding and the remaining work is mechanical -- another widening pass, a wider
spectator margin, two spectators instead of one. If it grows, there is a
configuration the whole approach does not reach, and the LNA it first appears on
is worth more than another hundred rules. Raise `--max-arrows` and
`--max-width` before reading anything into the current residue -- looking at it
by hand is cheap and, this time, it was the search bound.

---

## H-010 — Overlap is reducible only at an end, and that is a theorem about the procedure
*2026-09-16, strengthened the same day* · **SUPPORTED**

F-022 is an empirical statement: no interior sequence found so far reduces the
overlap of an isolated pair. The suspicion is that it is exact -- that **no**
interior mutation sequence does, whatever its length -- and that the reason is
visible in the procedure rather than in the search.

**2026-09-17 — not contradicted by the paper, but "isolated" now carries all
the weight.** The lead said `proposition:doubleMutation` might be an interior
overlap-reducing statement. On an *isolated* pair it is not: the companion
lengthens onto `r`, drops, and `s+1 → t+1` replaces it — the pair slide, nothing
more. But with bystanders the interior move (`s > 1`, `t < n`) **does** lower
the maximum overlap of the whole LNA: 1416 of 7164 interior applications at
`n = 10`, with the relation count going down by one in 1430 (E-029). So an
interior mechanism reduces overlap when something else crosses `r`, and this
hypothesis survives only as a statement about a pair with nothing crossing it.
Coverage still depends on the ends — interior applications alone place 770 of
1430 at `n = 9` against 1421 with them — so the ends are not incidental. The deep
probe (`probe.py 1:3,2:3 --steps 7 --clearance 9`) was **not** run on 2026-09-17;
it tests the isolated case only, which the proposition does not reach.

**It now rests on better evidence than it did.** The probes F-022 quoted allowed
mutations at every vertex of A_13, so they were not testing interiority at all.
Re-run where the ends are genuinely out of reach, `(1:3) (2:3)` reaches 2, 4, 4
and 6 LNAs at three, four, five and six mutations, and not one of them has a
smaller overlap (F-024). Six is two deeper than before, and the interior orbit
turns out to be five times smaller than the earlier numbers suggested -- most of
what those probes reached, they reached with an end's help.

**And the one sequence that ever lowered it is the boundary in disguise.** With
an end in reach, six mutations take `00003300000` to `30000020000` in A_13. The
middle of that sequence walks the relation down the quiver one vertex at a time
until it is at arrow 1, and the sequence is not translation invariant: shifted
by one to four, it does not even produce an LNA. So it is not a counterexample
to this hypothesis; it is H-011's mechanism, arrived at in one sequence rather
than as a composition of rules.

**Why it should be provable rather than searched for.** The collapse at the end
uses the one thing an end has: the source of the line has no arrow into it. The
procedure's steps at a vertex are about the paths through it, and at the source
there are none coming in, so the relation ending there has nothing to be
re-formed against and is dropped. In the interior the incoming arrow puts it
back. If that argument can be made properly it is a **conserved quantity**
statement -- the overlap of an isolated run is invariant under interior
mutation -- and it says the move table is complete in a direction, which no
amount of searching can say.

**What would settle it.** Either a proof from the procedure's step 7, or a
counterexample: an interior sequence, of any length, that lowers the overlap of
an isolated pair. `python probe.py 1:3,2:3 --steps 7 --clearance 9` is the next
search, and it is expensive -- six mutations took 43 minutes. The probes already
run are in E-021 and E-025 and should not be repeated; in particular **do not
re-run them without a clearance**, which is what made the earlier ones measure
the wrong thing.

**A cheaper line than searching deeper.** F-025 found that the two ends of the
quiver behave differently for an unequal pair -- the sink shortens it, the source
does nothing. Whatever asymmetry in the procedure explains *that* is likely the
same one that explains this, and it is a question about one mutation rather than
about seven.

**Caveat.** "Isolated" is doing work. A pair with a third heavily overlapping
relation was thought not to be invariant -- `(1:3) (2:3) (3:3)` dissolves in
three mutations -- but that turns out to hold only for short runs: at three
mutations `(1:3) (2:6) (3:7)` is as frozen as a bare pair (E-025). So the
statement is not simply "runs of two are invariant and longer runs are not", and
the right form of it is still open.

---

## H-009 — The move rules are a one-dimensional cellular automaton, and its theory applies
*2026-09-14, first check 2026-09-15, parked 2026-09-16* · **PARKED**

**Parked, deliberately, in favour of H-010 and H-011.** Not because anything
below is wrong -- F-017 stands, and the caveat it answers was a real one -- but
because the reading has to be fitted to the part of the problem that is still
open, and the part that is still open has just been identified precisely
(F-021): the heavily overlapping LNAs, and among them the isolated overlapping
pair. Its own last paragraph already says the translation to an arrow row is
finite-state only under *almost separate* relations, which is the part that is
already easy. A theory of the easy part, arrived at by fitting this problem to
another domain, is the likely yield, and the cost of getting it is a literature
sweep.

The two findings since are also the wrong shape for classical CA. F-022's rule
is anchored to an end of the quiver, so the boundary is not a perturbation of
the interior rule but where the interesting behaviour lives; and the empirical
invariant of H-010 -- the overlap of an isolated run, unchanged by any interior
sequence -- is a conserved quantity that wants proving from the mutation
procedure, not recognising in a rule table.

**When to unpark.** If H-010 is proved and the proof reads as a conservation law
rather than a computation, the CA literature on additive invariants is then
looking at something known to be there, which is a different proposition from
looking for one. Start from the arrow row, as F-017 says.


An LNA of length `n` is a row of `n - 2` cells, cell `i` holding the number of
arrows in the relation starting at vertex `i + 1`. Every move rule found so far
is then literally a **local rewrite on that row**: `lnaMoves.describeLink` only
admits a rewrite whose *window* contains every relation it touches, `matchesAt`
slides that window along the row, and the rules themselves read as
neighbourhood-to-neighbourhood maps — the pair slide (F-013) translates a
two-cell pattern one place along; the other verified rules lengthen or shorten a
relation depending on what overlaps it on either side.

That is the setting of **one-dimensional cellular automata**: a finite alphabet
(relation lengths, bounded by `n`), a finite neighbourhood (the window width), a
local transition rule, and boundary behaviour at the two ends that differs from
the interior — which is exactly the phenomenon H-007 is about.

**Why it might pay.** Questions we are currently answering by brute-force search
are standard questions there, with machinery behind them:

* *which rows are reachable from which* — the orbit/reachability problem for a
  rewriting system, and the injectivity/surjectivity theory of CA maps;
* *do the rules generate everything, or are there invariant classes* — additive
  invariants and conserved quantities of a local rule, which is the CA way of
  saying "a derived invariant the moves preserve";
* *when does a family of rules parameterised by window width collapse to one
  statement* — H-008's question, and the block/rescaling constructions are built
  for it;
* *how much does the boundary matter* — the difference between a CA on `Z` and on
  a finite interval, which is well studied and is H-007's question.

The nearest formal fit is probably not classical CA (synchronous, everywhere at
once) but **asynchronous CA** or a **one-dimensional rewriting / subshift**
presentation, since a mutation applies at one place at a time. Sand-pile and
chip-firing models are the closest-looking relatives: local, order-independent
in the right circumstances, with a well-developed theory of reachability and
invariants.

**What would settle whether it is worth pursuing.** A literature sweep first —
asynchronous CA and local rewriting on finite words, reachability under a finite
set of local rewrites, conserved quantities of local rules — then one concrete
attempt: state the verified move table of `lnaMoves.VERIFIED_MOVES` as a rule set
in that language and ask whether any standard result gives the orbit structure we
are currently getting by search.

**Caveat that would sink it.** A move is only valid when every mutation in its
sequence is admissible (R-005), and admissibility is a condition on the *algebra*,
not on the row of numbers. If the admissibility side conditions cannot be written
as part of the local neighbourhood, the CA picture describes something strictly
larger than the moves and its conclusions do not transfer.

**That check is done, and it passes → F-017.** Whether a rule applies comes out
of the window's cells plus one bit -- whether a relation covers the window's
first arrow having started earlier -- over 991,064 comparisons at lengths 5 to 9
with no disagreement; and wherever a rule matches, its mutations are legal, 1218
confirmations and no failures. So the side conditions *are* local and the CA
picture is about the right object.

**What the check also settled, which was not the question asked.** The state has
to be indexed by **arrows**, not vertices. A per-vertex cell holds a relation
length, which is unbounded in `n`, so the alphabet is unbounded and a relation
reaches arbitrarily far right; a per-arrow row carrying "covered / starts /
ends" has a fixed alphabet and makes that one bit a property of the cell at the
window's edge. Any attempt at this should start from the arrow row.

**Where it will strain.** Translating a per-vertex row into a per-arrow one needs
the number of relations open at each arrow. That is bounded by two under *almost
separate* relations and unbounded otherwise -- and the heavily overlapping LNAs
are exactly the ones the classification still has to search for (H-003). So the
CA reading may be exactly a theory of the part that is already easy.

---

## H-008 — Rule families are parameterised by relation length and overlap, and some need more than three mutations
*2026-09-14, settled 2026-09-15* · **CONFIRMED → F-020**

**Confirmed in both halves.** The lone short-relation slide is one statement for
every `d`, and it needs `d` mutations -- so its members from d = 4 on are exactly
the rules a three-mutation search cannot see, however simple they are.  Verified
for d = 1 to 7, both directions, 63 confirmations each with no failures, and
generated by `lnaMoves.shortRelationSlideRules` rather than waiting for discovery
to reach them.  The pair slide (F-013) is the other half: a family whose mutation
count does *not* grow.  The note below is what was suspected.

The rules found so far are individually verified but individually unilluminating.
The suspicion is that they are members of a handful of families parameterised by
the lengths of the relations involved and the size of their overlap — the shape
the quipu result has — and that **a family can be simple to state while needing
more mutations for larger parameters**. A rule requiring four or five mutations
would be invisible to a search bounded at three, however simple its statement.

H-001 is the proof of concept: the pair slide is one statement for every relation
length, and the search only found two members of it.

**Deliberately not being pursued yet.** Generalising from a three-mutation search
risks fitting families to an artefact of the search bound rather than to the
mathematics. Deepen the search first — see H-007.

**What would settle it.** Search deeper (four or five mutations, interior
embeddings), then look for families across the enlarged table; check whether the
mutation count of a family grows with its parameters.

---

## H-007 — The rules that matter live in the interior of long quivers
*2026-09-14* · **SUPPORTED**

On a quiver of length 6 or 7 every vertex is within a step or two of an end, so
the special cases that apply near a boundary apply almost everywhere, and a rule
that is really about the interior cannot be told apart from one that depends on
an end being nearby. The interesting interaction effects — between the number of
relations, their lengths, and how much they overlap — need room to appear.

**Evidence so far.** Discovery at length 8, where a six-arrow window first has
room to sit clear of both ends, found 34 rules that lengths 6 and 7 did not.

**What would settle it.** `lnaMoves.discoverLocalMoves` plants a pattern in the
middle of a length-13 or 14 quiver and mutates only near it; comparing its yield
against whole-quiver discovery at the same sequence length measures the effect
directly.

---

## H-006 — Every class of LNAs is a quipu class, a canonical-type class, or non-piecewise-hereditary
*2026-09-13* · **SUPPORTED**

By Happel's classification a hereditary abelian category is either the module
category of a hereditary algebra or derived equivalent to a canonical algebra.
So a piecewise hereditary LNA should be derived equivalent either to a tree
algebra — a quipu, in this setting — or to a canonical algebra; and what is left
is the non-piecewise-hereditary case.

**Evidence.** Exactly true at n = 9: 18 + 1 + 1 = 20 (F-011).

**Caveat that would break it.** "Tree algebra" is not the same as "quipu algebra".
A hereditary algebra of a tree that is *not* a quipu — maximum degree 4, or
degree-3 vertices not on one path — would be a fourth case. Whether such a tree
can be derived equivalent to an LNA is not known here.

**What would settle it.** Check at n = 10 and 11, where there is room for a
non-quipu tree; `quipuForms.canonicalUndirectedForm` already reports the tree for
any relation-free quiver a search reaches, so a non-quipu tree would be visible
as a canonical form with no quipu notation.

---

## H-005 — The class count is (quipus of order n) + (canonical classes) + (non-piecewise-hereditary classes)
*2026-09-13* · **OPEN**

A counting form of H-006. If it holds, the number of classes is computable
without any mutation search once the three counts are known — and the quipu count
already is (F-010's table).

**Evidence.** Holds at n = 9. Untested above.

**Caveat.** Needs every quipu of order n to actually occur as the class of some
LNA of length n, which is true up to 9 but not proved in general.

---

## H-004 — Non-piecewise-hereditary classes become the common case as n grows
*2026-09-13* · **SUPPORTED**

The certified counts are 0, 0, 1, 24, 308 for n ≤ 8, 9, 10, 11 (F-012). The paper
argues the scarcity at small n is an artefact of the sizes that computer search
can reach.

**What would settle it.** These are *certified* counts, a lower bound on the
truth — the criteria are sufficient, not a characterisation. The gap between
certified and actual is unknown and worth measuring at n = 10 by another route.

---

## H-003 — Wider rules are what unlock the heavily overlapping LNAs
*2026-09-13, settled 2026-09-16* · **REFUTED → R-010**

**Its question is answered and its diagnosis is wrong.** The measurement it asked
for is F-021: the rows a search still has to place are exactly the LNAs whose
consecutive relations share two or more arrows -- nothing below that line is ever
left over, and of the 820 above it at n = 9 the whole rule table reaches 34. So
the diagnosis below is right that the rules keep an LNA inside the almost
separate set.

What is wrong is the remedy. The commonest blocking configuration, by a factor
of four, is an *isolated pair* of relations sharing two or more arrows, and its
overlap cannot be reduced by any interior sequence at any width tried -- three,
four and five mutations, and a margin of six (F-022, E-021). It comes apart at an
**end** of the quiver instead, under a rule the framework could not state until
it grew anchored descriptions. Aiming discovery at bigger interior patterns,
which is what this hypothesis asked for, would not have found it. R-010.


Seeding by the quipu theorem and expanding along move orbits places 72% of the
n = 7 table, 57% of n = 8, 45% of n = 9 — and the 34 rules added at length 8
barely moved those numbers. The diagnosis is that the rules found so far mostly
keep an LNA *inside* the almost-separate set the theorem already covers, while
the rows still needing a search are the heavily overlapping ones.

**What would settle it.** Measure, for the rows a search still has to place, what
relation patterns they have; then aim discovery at exactly those patterns rather
than at small ones.

---

## H-002 — The Coxeter polynomial plus the hereditary form settles every merge
*2026-09-12* · **SUPPORTED**

The merge step needs to decide, for classes sharing a Coxeter polynomial, whether
they are one class. The hereditary form does this where it exists, and F-010 says
where it is needed.

**Evidence.** At n = 6, 7, 8, 9 the combination leaves nothing as a candidate.

**Caveat.** It relies on every class having a form. A class with neither a quipu
nor a canonical type nor a certificate would be left as a candidate, and none has
been seen — which is H-006 again.

---

## H-001 — The pair slide holds for every relation length
*2026-09-14* · **CONFIRMED → F-013**

Discovery found the pair slide at relation lengths 3 and 4 only. The suspicion was
that this is an artefact of which window widths fit at the lengths searched, and
that it holds for every length with the same two mutations.

Confirmed the same day: `l = 2..7`, 22 confirmations per direction per length, no
failures. The family is now generated by `lnaMoves.pairSlideRules` rather than
waiting for discovery to reach each member.

**The lesson, which generalises:** a rule family shows up in discovery only at the
widths the search reaches. Absence of a family member from the table is evidence
about the search, not about the mathematics.
