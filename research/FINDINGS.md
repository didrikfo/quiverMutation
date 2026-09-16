# Findings

Established results, newest first. Each carries the evidence it rests on.
See [`README.md`](README.md) for conventions.

---

## F-027 — Leaving the line, a mutation of an LNA goes to a square with a side of two
*2026-09-16*

A rule of two or more mutations passes through quivers that are not lines. They
are not arbitrary. Walking every multi-mutation rule in the table one step at a
time and classifying each intermediate:

| intermediate | count |
|---|---|
| commutative square, sides 2 and 2 | 450 |
| sides 2 and 3 | 351 |
| sides 2 and 4 | 297 |
| sides 2 and 5 | 216 |
| sides 2 and 6 | 31 |
| sides 2 and 7 | 9 |
| sides 2 and 8 | 10 |
| a line again | 336 |
| anything else | **5** |

1296 rules, 1705 intermediates. **Every square has a short side of exactly two**
-- 1364 of them, and not one with a short side of three. The five exceptions
have a branching vertex and no matching join, so they are not squares at all.

**What the square is.** Mutating at the source vertex of a relation of k arrows
takes the relation `v -> v+1 -> ... -> v+k` and replaces it with a commutative
square: a new two-arrow path from a vertex to the relation's end, commuting with
the k-1 arrows still on the line. Right mutation at the start does it one way
round and left mutation at the end the other:

```
A_10, one relation on the arrows 3..7, right mutation at 3
    relations  (2,4,3)  and  (4,3,8) = (4,5,6,7,8)
                              ^^^^^ two arrows          ^^^^^^^ four arrows
```

so the zero relation has become a *commutativity* relation between a side of two
and a side of k - 1.

**Why it matters, and it is a construction rather than an observation.** The
rules were all found by search -- enumerate mutation sequences, describe what
recurs. The square says what the search is walking through: a rule is *open the
relation into a square, do something along its long side, close it back onto a
line*. That is the same shape as F-020's families, where the mutation count grows
one per arrow travelled, and it suggests building rules directly instead of
finding them: open at a chosen relation, walk the long side, and read off where
it closes.

**The caveat, from trying it.** A lone relation opened into a square closes only
by undoing itself -- searching four mutations from the 2-by-4 square of
`00500000` in A_10 finds nothing but `[-3]` back to where it started. The square
has to have something to interact with, which is the same lesson as F-023's
spectators: the interesting rules are the ones with a second relation in the
window. Constructing rules this way therefore means opening a square *and*
choosing the companion, which is a smaller search than the one being run now but
not a formula.

E-026.

---

## F-026 — A rule's dual is a rule, and the table was missing 410 of them
*2026-09-16*

The relation dual -- reverse every arrow of the line and renumber -- preserves
the derived equivalence class of any LNA, and left mutation at a vertex is right
mutation at that vertex of the dual. So a rule must carry over to the dual
picture, and the transform is mechanical.

**The transform.** For a rewrite on a window of `width` arrows:

* a relation covering the arrows `s .. s+a-1` covers `width-s-a .. width-s-1`;
* the vertex at offset `o` becomes the one at `width-o+2`, and a **right**
  mutation there becomes a **left** one, and the other way about;
* the sequence keeps its order, the dual being applied step by step;
* an anchor to one end becomes an anchor to the other.

`lnaMoves.dualRule`, and it is an involution.

**It holds.** Of the 1384 rules the table then held, **none** was its own dual and
**410** had a dual that was not in the table. Verified where each fits, at up to
four lengths apiece: **410 hold, 0 fail, 0 never apply.** The table is generated
closed under the dual now -- 364 floating rules and 1430 anchored, 1794 in all.

**It is not the transform E-019 refuted**, and the difference is worth keeping
straight. That one tried to read a rule's **inverse** off its window by reversing
the sequence and negating: it worked for 12 of 96 rules and was abandoned. This
is a symmetry of the problem rather than a shortcut, and it works for all of
them.

**What it is worth, honestly: very little coverage.** Closing the table adds 410
rules and moves the count by **7 rows at n = 10 and none at n = 9**. The orbits
those duals join were already joined another way. Its value is elsewhere:

* it is free, and a table that is not closed under a symmetry of the problem is
  wrong to leave that way;
* it halves what a search has to look for -- discovery could plant patterns at
  one end only and dual the results, at half the cost of E-025's two hours;
* and it is the correction to F-025, which claimed an asymmetry between the two
  ends of the quiver on the strength of comparing a pattern with itself rather
  than with its dual (R-011).

E-026. Tests: `test_the_dual_of_a_rule_reverses_the_window_and_turns_the_mutations_round`,
`test_the_table_is_closed_under_the_dual`.

---

## F-025 — The two ends of the quiver are not the same end
*2026-09-16* · **RETRACTED 2026-09-16 → R-011**

**The asymmetry is the pattern's, not the quiver's.** The mirror of a rule is its
relation dual -- reverse the arrows *and* exchange right mutation for left -- and
under that transform the sink rule below holds perfectly well at the source, on
the dual pattern. What was compared with it was the same pattern at the other
end, which is a different configuration. The rules and the family below stand and
are in the table; the conclusion drawn from them does not. R-011, F-026.


`endPairCollapseRules` (F-022) collapses a pair of relations of **equal** length
against either end, and the two directions are mirror images, as the relation
dual says they must be. For a pair of **unequal** lengths that symmetry breaks,
and only one end works.

**The probe.** Each of the four long unequal pairs the residue of E-024 is made
of, planted flush against each end of A_13 or A_14 and mutated within three
vertices, at three mutations:

| pattern | at the source | at the sink |
|---|---|---|
| `(1:3) (2:6)` | 2 LNAs, **nothing lower** | overlap 2 → **0**, via `[7, 6, 7]` |
| `(1:3) (2:7)` | 2 LNAs, **nothing lower** | overlap 2 → **0**, via `[7, 6, 7]` |
| `(1:5) (2:6)` | 2 LNAs, **nothing lower** | overlap 4 → 3, via `[7, 7]` |
| `(1:5) (2:7)` | 2 LNAs, **nothing lower** | overlap 4 → 3, via `[7, 7]` |

At the source nothing moves at all. At the sink the **shorter** relation, which
is the one that starts first, loses an arrow — and where it had three arrows to
begin with, losing one takes the overlap to zero and the LNA into the quipu
theorem's reach outright.

**Two of them as rules, and the anchor is not decoration.**

```
window 7 arrows at the right end:  (0:3) (1:6)  ->  (0:2) (1:6)     via [2, 2]
window 8 arrows at the right end:  (1:3) (2:6)  ->  (0:2) (2:6)     via [3, 2, 3]
```

4 confirmations each over the three lengths their windows fit in, no failures.
The identical rewrites stated as floating rules give **4 confirmations and 4
failures** apiece: they are true against the sink and false elsewhere, checked
rather than assumed.

**Why the asymmetry is not surprising once stated.** The pair is `(n, l)` and
`(n+1, m)` with `l < m`, so the two relations start one vertex apart but end
`m - l + 1` arrows apart. Flush against the sink it is the *ends* that are
pinned, and the two relations end at different places, so the configuration
there is genuinely different from the one at the source, where it is the
*starts* that are pinned and they are one apart either way. The equal-length
pair is exactly the case where the two descriptions coincide, which is why
F-022's family is symmetric and this is not.

**What it opens.** The rule is stated here for the two pairs probed, not for a
family in `(l, m)`. Sweeping `l` and `m` is what E-025's discovery run is for,
and the family is the thing to look for in its output. E-025.

---

## F-024 — The interior is emptier than it looked, and the boundary was doing the work
*2026-09-16*

E-021 measured what a heavily overlapping pair can be turned into and found it
frozen. The measurement was right and the label on it was wrong: those probes
allowed mutations at every vertex of A_13, ends included, because a margin of 6
around arrows 5 to 8 of a 13-vertex quiver reaches both of them and even a margin
of 3 reaches vertex 2, which rewrites arrow 1. Re-run with the quiver lengthened
so that the ends are genuinely out of reach, the picture changes in size but not
in conclusion — and the difference is the finding.

**The same pattern, the same margin, with and without an end in reach.**

`(1:3) (2:3)`, margin 3, in A_13 at offset 4 (arrows 1 to 12 rewritable, so both
ends in reach) against A_21 at offset 8 (arrows 5 to 14 only):

| mutations | LNAs reached, ends in reach | LNAs reached, genuine interior |
|---|---|---|
| 3 | 8 | **2** |
| 4 | 14 | **4** |
| 5 | 22 | **4** |
| 6 | 36 | **6** |

In a genuine interior the pair is not so much frozen as nearly immobile: six
mutations reach six LNAs, and the fifth mutation buys nothing at all. Five of
every six LNAs the earlier probes reported were reached with the help of an end.

**And the overlap never comes down in the interior, now checked to six.** Every
one of those 2, 4, 4 and 6 still has two relations sharing two arrows, the
depth-6 run taking 2583 seconds to say so. That is H-010's claim tested two mutations
deeper than before, and it survives.

**With an end in reach, six mutations do pull the pair apart** — the one thing
that has ever done so:

```
A_13:  00003300000  ->  30000020000   via [-8, 5, 4, 3, 2, -6]
```

overlap 2 down to 0. The middle of that sequence is `5, 4, 3, 2`: four mutations
walking down the quiver, one vertex at a time, until the relation is at arrow 1
and there is no further to go. It is **not translation invariant** — shifted by
1, 2, 3 or 4 vertices in a quiver lengthened to match, it does not even produce
an LNA, let alone the shifted answer. So it is not a rule that happens to need
six mutations; it is the boundary, reached the long way round.

**What this settles and what it costs.** It settles that the escape route is the
one H-011 names -- walk the run to an end -- and that it is a single mutation
sequence, not only a composition of table rules. It costs E-021 its framing: the
rows there labelled interior were whole-quiver rows, which made them stronger
claims about the *LNA* and weaker ones about *locality*, and the distinction
matters because locality is what a move rule is. `probe.py` reports the arrows a
run can rewrite for that reason, and lengthens the quiver unless told not to.

E-025.

---

## F-023 — What a rule needs to fire is a spectator, and the ends are where they are
*2026-09-16*

F-022 put sixteen rules against the ends of the quiver and coverage rose more
than the whole floating table had bought. Running discovery there properly says
how much more there is, and why the rules found until now so rarely fire.

**Discovery against an end.** `discover.py --anchor both --max-arrows 5
--max-width 6` plants each of the 74 patterns flush against each end of A_11 and
A_12, mutates within three vertices of it, keeps the rewrites described at both
lengths, and verifies each at the four lengths its window fits in:

| | |
|---|---|
| rewrites described | 1344, in 286 s |
| recurring at both lengths | 892 |
| verified, no failures | **724** |
| of those, a floating rule restricted to an end | 94 |
| genuinely anchored | **630**, 315 at each end |

315 at each end is a consistency check worth noticing: the relation dual
exchanges the two ends, so a rule at one has a mirror at the other, and the
counts had to come out equal.

**Most of them carry a spectator, and that is the finding.** A rule's window has
until now held nothing but the relations it rewrites -- `matchesAt` refuses a
position where any other relation reaches in. That is what makes a rule true and
it is why so few of them match anything: of the 155 LNAs left unplaced at n = 8,
**126 have a rule whose left-hand pattern is present and which does not fire**,
because one to four further relations are sitting in the window doing nothing.
Only 29 have no rule with their pattern at all. Among the 229 anchored rules
that change the orbit partition, **188 carry at least one relation that appears
unchanged on both sides** -- a bystander the rewrite steps around.

So the blockage was never that the patterns were too small. It was that they
were too clean.

**They do not compress into families.** Setting the spectators aside leaves 190
distinct rewrites among the 229, so unlike the slide families (F-013, F-020)
there is no statement covering many at once, and `quivermutation/endMoves.py`
lists them. 38 need one mutation, 99 two, 92 three.

**Two of them are worth reading on their own.**

```
window 2 arrows at the left end:   (0:2)  ->  -    via [1]
window 2 arrows at the right end:  (0:2)  ->  -    via [-3]
```

A lone relation of two arrows at an end of the quiver is simply deleted, by a
single mutation at the end vertex. That is operation 2 of
`cor:EquivNakayamaAlgebras` -- "a relation of two arrows does not change the
class" -- appearing as a mutation rather than as a theorem, and it is the rule
NOTES warned about: stated as a *floating* two-arrow window it holds 63 times
and fails 130, which is exactly right, because away from an end it is false.

**And widening the rules to admit a spectator is what the diagnosis was for.**
`lnaMoves.spectatorExtensions` puts one untouched relation into a rule's window,
growing the window by up to three arrows to make room, and `verifyMove` decides;
10609 such widenings produced **625** verified rules in eight minutes of no
searching at all, of which `spectatorMoves.SPECTATOR_MOVES` lists the 270 that
change the orbit partition -- 125 that float and 145 that need an end (E-024).

**What the two batches are worth.** LNAs placed with no mutation search:

| n | LNAs | theorem | + floating | + anchored | + widened |
|---|---|---|---|---|---|
| 6 | 42 | 34 (81%) | 35 (83%) | 42 (100%) | 42 (100%) |
| 7 | 132 | 89 (67%) | 95 (72%) | 127 (96%) | **132 (100%)** |
| 8 | 429 | 233 (54%) | 246 (57%) | 347 (81%) | **406 (95%)** |
| 9 | 1430 | 610 (43%) | 644 (45%) | 863 (60%) | **1038 (73%)** |

(E-025 has since raised the last column again, to 98% at n = 8 and 84% at n = 9,
and measured n = 10 and n = 11 for the first time: 63% and 47%.)

A classification of A_6 or A_7 is now a table lookup; A_8 needs a search for 23
rows and A_9 for 392. The rows still left are still exactly the heavily
overlapping ones, so F-021's reading is unchanged; there is simply much less of
it. And the diagnostic says the same thing about them as before: at n = 9, 359
of the 392 have a rule whose pattern is present and blocked by a bystander, and
only **33** have no rule with their pattern at all. Those 33 are what to look at
by hand.

**The curation, stated plainly.** `endMoves.DISCOVERED_END_MOVES` lists the 229
rules that change the orbit partition at n <= 9, not all 630 verified ones, and
`spectatorMoves.SPECTATOR_MOVES` the 270 of 625 the same way. A rule left out
reaches nothing the listed ones do not *at the lengths measured*, and could in
principle be the one that matters at n >= 10; re-running the commands gets them
all back. E-023, E-024.

---

## F-022 — An overlapping pair is frozen in the interior and comes apart at an end
*2026-09-16*

F-021 says the whole gap is the LNAs whose relations overlap in two or more
arrows, and that the commonest blocking configuration by a wide margin is a
**pair** of relations sharing two or more arrows. This is why the move rules
cannot place them, and what does.

**In the interior the overlap of an isolated pair does not move.** Plant
`(1:3) (2:3)` -- two relations of three arrows at consecutive vertices, sharing
two -- in the middle of A_13, with four arrows of empty quiver on the left and
six on the right, and enumerate every admissible mutation sequence near it:

| mutations | margin | LNAs reached | with a smaller maximum overlap |
|---|---|---|---|
| 3 | 3 | 8 | **0** |
| 4 | 3 | 14 | **0** |
| 5 | 3 | 22 | **0** |
| 4 | 6 | 34 | **0** |

The last row is stronger than the others and was first recorded as though it
were not: a margin of 6 around arrows 5 to 8 of A_13 admits *every vertex of the
quiver*, so it says that four mutations anywhere in A_13 -- ends included --
leave the pair intact.

Every LNA reachable still has two relations sharing two arrows. The same holds
for `(1:4) (2:4)` (overlap 3, 17 reached at three mutations, none below 3),
`(1:5) (2:5)` (overlap 4, 16 reached, none below 4) and for the unequal pairs
`(1:3) (2:4)` and `(1:4) (3:3)`, where the overlap goes *up* to 3 in three of
the sixteen but never down. Neither depth nor width is the obstacle: depth 5 and
a margin of 6 reach further into the quiver and find nothing new.

**A third heavily overlapping relation unlocks it; anything else does not.**

| pattern | run of | three mutations reach |
|---|---|---|
| `(1:3) (2:3)` | 2 | overlap 2, all 8 |
| `(1:3) (2:3) (3:3)` | 3 | **overlap 0**, via `[6, 5, 6]` |
| `(1:3) (2:4) (3:4)` | 3 | **overlap 0** |
| `(1:4) (2:4) (4:3)` | 3 | **overlap 0** |
| `(1:4) (2:4) (3:4)` | 3 | overlap 2, down from 3 |
| `(1:2) (2:3) (3:3)` | 2 | overlap 2, all 31 |
| `(1:3) (2:3) (4:2)` | 2 | overlap 2, all 31 |
| `(1:3) (2:3) (5:2)` | 2 | overlap 2, all 35 |

A third relation only helps when it *also* shares two or more arrows with the
pair: `(1:2) (2:3) (3:3)` has three relations and is as stuck as the bare pair,
because its first relation shares only one arrow. So the parameter is the length
of the **overlapping run** -- maximal relations linked by an overlap of two or
more, `overlap.overlapRuns` -- and a run of two is frozen where a run of three is
not.

**QUALIFIED 2026-09-16 by E-025.** "A run of three is not frozen" is true of the
short runs in the table above and false in general: `(1:3) (2:6) (3:7)` and
`(1:5) (2:6) (3:7)` reach two LNAs each at three mutations and neither lowers
the overlap. What the length of the run buys is probably mutations rather than
freedom -- F-020's one per arrow travelled -- so a long run may well dissolve
deeper down. Do not quote the sentence above without this. The rules that dissolve a run of three were already in the table; nothing in
it dissolves a run of two, and E-021 says why nothing was ever going to be found.

**At an end of the quiver the pair collapses in two mutations.** The source of
the line has no arrow into it, so a mutation there is not the mutation the same
rewrite would be in the interior. Where the pair starts at vertex 1, two *right*
mutations at vertex 1 delete the second relation outright; at the sink, two left
mutations at vertex n delete the first:

```
window l + 1 arrows at the left end:   (0:l) (1:l)  ->  (0:l)    via [1, 1]
window l + 1 arrows at the right end:  (0:l) (1:l)  ->  (1:l)    via [-(l+2), -(l+2)]
```

Verified for `l = 2` to 7, both ends, at the four lengths `l + 2 .. l + 5` each:
**9 confirmations apiece, no failures**, a confirmation being an admissible
sequence landing on the predicted LNA with the Coxeter polynomial kept. The
window is exactly the pair's span, so no other relation may touch it.
`lnaMoves.endPairCollapseRules` generates the family.

**Why this needed the framework to grow a notion it did not have.** Every rule
until now was a rewrite holding at *every* window position. The collapse holds at
one position and is false at all the others -- checked, not assumed:
`(0:3) (1:3) -> (0:3)` via `[1, 1]` stated as a floating rule fails, which is
`test_the_end_pair_collapse_is_false_in_the_interior`. So a description now
carries an optional anchor, `'left'` or `'right'`; `lnaMoves.windowStartsFor` is
the single gate every caller slides a rule through, and an anchored rule is
offered only its own position. Stating such a rule as though it floated is
exactly how R-009's false rules arose.

**What it buys, and it is more than sixteen rules should.** With the anchored
family and no other change:

| n | LNAs | theorem | + floating orbits | + anchored |
|---|---|---|---|---|
| 6 | 42 | 34 | 35 (83%) | **38 (90%)** |
| 7 | 132 | 89 | 95 (72%) | **107 (81%)** |
| 8 | 429 | 233 | 246 (57%) | **274 (64%)** |
| 9 | 1430 | 610 | 644 (45%) | **726 (51%)** |

Sixteen anchored rules place 82 rows at n = 9 where all 123 floating rules
placed 34. The mechanism is the pair slide (F-013) walking a pair to an end and
the collapse taking it from there -- `closureUnderMoves(10, [3,3,0,0,0,0,0,0])`
is the six positions of the pair *and* the two LNAs where it has lost a
relation, which is
`test_the_pair_slide_walks_a_pair_along_the_quiver_and_off_each_end`.

E-021, E-023. Tests: `tests/test_overlap.py`.

---

## F-021 — What a classification search still has to find is exactly the heavily overlapping LNAs
*2026-09-16*

H-003 asked what relation patterns the rows still needing a search actually
have. They have one, and it is sharp: **every LNA the quipu theorem and the move
orbits fail to place has two consecutive relations sharing two or more arrows,
and almost every LNA that has two such relations is one of them.**

**Overlap is the right coordinate because the theorem's condition is a bound on
it.** Consecutive relations `(n_i, l_i)`, `(n_{i+1}, l_{i+1})` share
`max(0, n_i + l_i - n_{i+1})` arrows, and *almost separate* -- the hypothesis of
`thm:QuipuToAn` -- is exactly that this never exceeds one. So "what the theorem
misses" and "what overlaps by two or more" are the same set, not merely
correlated ones; `test_almost_separate_is_exactly_overlap_at_most_one` checks the
two predicates against each other over every LNA of lengths 4 to 8.

**The measurement.** For each length, partition the LNAs into orbits under the
verified move rules and call an LNA *covered* when its orbit contains one the
theorem names -- which is what `seedTableFromQuipuTheorem` fills in before a
single mutation is computed. With the 123 floating rules:

| n | LNAs | overlap ≤ 1 | overlap ≥ 2 | of those, covered | left |
|---|---|---|---|---|---|
| 6 | 42 | 34 | 8 | 1 | 7 |
| 7 | 132 | 89 | 43 | 6 | 37 |
| 8 | 429 | 233 | 196 | 13 | 183 |
| 9 | 1430 | 610 | 820 | 34 | 786 |

Two things at once. Nothing at overlap ≤ 1 is ever left over -- so there is no
second phenomenon hiding among the rows the theorem does cover, and the search's
whole remaining job is the overlapping ones. And of the 820 heavily overlapping
LNAs at n = 9 the entire rule table reaches **34**: the rules found so far are
almost exactly the rules that keep an LNA where it already was.

**What the leftovers look like.** Their heavily overlapping runs -- maximal
relations linked by an overlap of two or more -- are overwhelmingly *pairs*, and
overwhelmingly the shortest pair there is:

| run | n = 7 | n = 8 | n = 9 |
|---|---|---|---|
| `(1:3) (2:3)` | 21 | 96 | 391 |
| `(1:4) (2:4)` | 6 | 40 | 198 |
| `(1:3) (2:4)` | 5 | 31 | 144 |
| `(1:4) (3:3)` | 5 | 31 | 144 |

and by the length of the longest run in the LNA, 434 of the 786 left at n = 9
have no run longer than two. That pair is what F-022 is about.

**The instrument.** `quivermutation/overlap.py` -- `overlapProfile`,
`maxOverlap`, `overlapRuns`, `coverage`, `blockingCores` -- and `overlaps.py`
over it:

```bash
python overlaps.py 6 7 8 9              # the table above
python overlaps.py 9 --cores            # what is left, by overlapping run
python overlaps.py 9 --floating         # without the rules anchored to an end
```

The whole thing is seconds per length because the orbits are computed as
rewrites on the relation lengths with no mutation run at all -- legitimate
exactly because every rule has been checked against the engine wherever it
applies (F-017).

H-003 → its question answered, its diagnosis refuted (R-010). E-020.

---

## F-020 — Three rule families whose mutation count grows with their parameter
*2026-09-15*

H-008 suspected that the move rules are members of families parameterised by
relation length and overlap, and that **a family can be simple to state while
needing more mutations for larger parameters** — which would make its later
members invisible to a search bounded at three. That is now a family, not a
suspicion.

**The lone short-relation slide.** A relation of two arrows with nothing else in
its window travels `d` arrows right under the `d` left mutations at the window's
vertices 3, 4, …, d + 2, and back under the `d` right mutations at d + 1, d, …, 2.
The window is `d + 2` arrows wide.

| d | window | sequence | found by |
|---|---|---|---|
| 1 | 3 | `[-3]` | E-010, two mutations |
| 2 | 4 | `[-3, -4]` | E-010 |
| 3 | 5 | `[-3, -4, -5]` | E-011, three mutations |
| 4 | 6 | `[-3, -4, -5, -6]` | **nothing** — needs four |
| 5 | 7 | `[-3, …, -7]` | needs five |
| 6 | 8 | `[-3, …, -8]` | needs six |
| 7 | 9 | `[-3, …, -9]` | needs seven |

Verified for d = 1 to 7, both directions, at the four lengths `d+3 .. d+6` each —
so up to A_14 and its 742900 LNAs: **63 confirmations apiece, no failures**, where
a confirmation is an admissible sequence landing on the predicted LNA with the
Coxeter polynomial kept. The count
being 63 at every d is itself a consistency check — the pattern is one relation
alone in its window, so the number of matching LNAs does not depend on how far it
travels.

**Two more of the same shape, read off consecutive widths in the enlarged table.**
Once the first family was recognised, looking for others cost minutes rather than
the hours a four-mutation search would have:

| family | rewrite on the window | sequence | mutations | verified |
|---|---|---|---|---|
| lone slide | `(0:2)` → `(d:2)` | `[-3, …, -(d+2)]` | d | d = 1–7, both directions |
| trailing walk | `(0:2) (2:2)` → `(1:2) (d+2:2)` | `[-3, -5, …, -(d+4)]` | d + 1 | d = 1–6 |
| spreading pair | `(1:2) (4:2)` → `(0:2) (d+4:2)` | `[2, -7, …, -(d+6)]` | d + 1 | d = 1–5 |

8 confirmations per member for the latter two, at three lengths each, no
failures; the spreading pair reaches A_14. Discovery had found d = 1 and 2 of each
and could not have found more — d = 3 of either needs four mutations. **One
mutation per arrow travelled** is the shape all three share: the relation that
moves furthest pays for each arrow, and any companion relation costs one more.

**Contrast with the pair slide (F-013).** That family is *two* mutations for every
relation length: only the window grows. Here the mutation count grows with the
parameter, so the two kinds are the two halves of H-008's statement, and this kind
explains why discovery keeps finding "new" rules that are the same rule.

**The practical consequence, and it is the point.** Discovery can only ever find
an **initial segment** of such a family — `d <= maxSteps` — so raising the search
bound by one buys one more member at multiplying cost, while recognising the
family gives every member at once. `shortRelationSlideRules`,
`trailingRelationWalkRules` and `spreadingPairRules` generate the three, as
`pairSlideRules` does for the other kind, taking the table from 96 listed rules
to 123. Before spending a four-mutation search, look at what the three-mutation
one found for a family whose members would be out of reach.

**What does not generalise.** The three slide families' inverses come for free:
reverse the sequence, negate each vertex, move each one step toward zero. That is
not a property of the table, though — it holds for only 12 of the 96 listed
rules, and fails or leaves the window for the rest (E-019). It works here because
a slide's sequence is one uniform run of mutations in one direction.

H-008 → CONFIRMED. E-011, E-019.

---

## F-019 — The gentle LNAs are one class, so gentle invariants separate nothing here
*2026-09-15*

Idea 22 (and NOTES idea 14) proposed the Avella-Alaminos–Geiss invariant as the
independent route to separating classes the quipu theorem does not reach, on the
stated ground that "LNAs are gentle". They are not, and even where they are the
route is empty.

**An LNA is gentle exactly when every relation has two arrows.** A gentle algebra
is a monomial algebra whose ideal is generated by paths of **length two**, plus
degree conditions on each vertex and two conditions on pairs of arrows. On a line
quiver every vertex has at most one arrow in and one out, so all of that is
automatic and only the length of the relations is at issue — and a relation of
three or more arrows is not generated by paths of length two. Such an LNA is
still a *string* algebra, but the AAG invariant is a gentle-algebra invariant.

**And every gentle LNA is the hereditary one.** Two facts collide:

1. all relations of length 2 implies almost separate relations — the condition
   `n_{i+1} >= n_i + l_i - 1` reads `n_{i+1} >= n_i + 1` when every `l_i` is 2,
   which distinct starts already give; and
2. operation 2 of `cor:EquivNakayamaAlgebras` drops a relation of two arrows
   without changing the derived equivalence class.

Drop them all and nothing is left. So the `2^(n-2)` gentle LNAs of length n are
all in the single class `P^(0)_(0,n-1)`, the class of the path algebra of A_n.
Checked against the classifications: 8, 16, 32, 64, 128 gentle LNAs at n = 5 to 9,
every one of them in that class and no other class containing any, and the
theorem-level statement pinned for n = 4 to 10 in
`test_every_gentle_lna_is_the_hereditary_one`.

**Consequence.** No gentle-algebra invariant can separate two LNA classes, because
no class but the hereditary one contains a gentle algebra at all. In particular
neither member of the cospectral pair `P^(1,4)_(1,0,1)` / `P^(1,2)_(1,1,2)` has a
single gentle member among its 18 — checked directly. The AAG route is closed;
what is left on idea 22's list is Hochschild cohomology, which is a derived
invariant of any algebra and does not care about the relation lengths. R-008.

---

## F-018 — A name read off the Coxeter polynomial must never separate two classes
*2026-09-15*

The n = 9 classification, re-run on the corrected engine (F-015) and the loosened
gate (F-016), came out at **22** classes rather than F-011's 20. The engine was
not at fault. The two extra classes were a **circular use of the Coxeter
polynomial** in the merge step, and once that is removed the partition is F-011's
20 again, on the very same table of rows.

**The mechanism.** A class the quipu theorem does not name gets a fallback name
from `canonicalWeightType`, which searches the weight types of the right order
for one whose canonical algebra has *exactly this class' Coxeter polynomial*. So
`C(2,3,5)` is a restatement of the polynomial and nothing more. But it is a
different *string* from `P^(5)_(1,2)`, and `mergeReport` separated two classes
whenever their form strings differed. The class named `C(2,3,5)` shared its
polynomial with `P^(5)_(1,2)` — which is how it got that name — and was declared
distinct from it on the strength of the very polynomial the two have in common.
The same happened to `C(2,2,6)` against `P^(1,1)_(1,3,1)`. Worse, the label also
stopped the class ever being searched again: `resolveMergeCandidates` only looked
at classes with *no* form.

**Both pairs really are one class, two independent ways.**

1. **Mutation.** `1;2;3;4;5|2;3;4;5;6;7`, of the class the run called `2334400`,
   reaches `1;2;3;4;5|2;3;4;5;6;7|6;7;8` of `P^(5)_(1,2)` in **two** mutations;
   `2;3;4;5|3;4;5;6|6;7;8;9`, of `2233030`, reaches `1;2;3;4|2;3;4;5|6;7;8;9` of
   `P^(1,1)_(1,3,1)` in **two**. Both also merge from the relation dual. These
   were available to the old resolution step at any depth it ran; it never looked.
2. **Theory.** A canonical algebra is derived equivalent to a hereditary algebra
   exactly when its weight type is *domestic* — `(p,q)`, `(2,2,n)`, `(2,3,3)`,
   `(2,3,4)`, `(2,3,5)` — with the partner of extended Dynkin type Ã, D̃, Ẽ6,
   Ẽ7, Ẽ8. Both (2,3,5) and (2,2,6) are domestic, and the trees come out at
   exactly the two quipus in question:

   | weight type | affine type | vertices | tree |
   |---|---|---|---|
   | (2,3,3) | Ẽ6 | 7 | `P^(2)_(2,2)` |
   | (2,3,4) | Ẽ7 | 8 | `P^(3)_(1,3)` |
   | (2,3,5) | Ẽ8 | 9 | `P^(5)_(1,2)` |
   | (2,2,6) | D̃8 | 9 | `P^(1,1)_(1,3,1)` |

   So a **domestic** `C(...)` can never be a class of its own: the class is a
   quipu class, and the quipu theorem has already named it under another name in
   the same table. A domestic weight type in that column is always a merge nobody
   found. `C(2,4,4)`, the n = 9 class that survives, is **tubular** — the boundary
   past which a canonical algebra is derived equivalent to no hereditary algebra
   — so that one is genuine.

**What the column means now.** Three kinds of value, and only two of them decide
anything:

* a quipu name or a tree encoding — **proved**, by a mutation path to a
  relation-free quiver or by the theorem. Two classes with different proved forms
  are different classes; two with the same one are the same class.
* `not piecewise hereditary` — **proved**, negatively. It separates such a class
  from every quipu class and merges nothing.
* `C(...)` — **not proved**. `isCoxeterDerivedForm` marks it; it may not separate
  and it may not merge, and a class carrying one stays a merge candidate.

**Why the pipeline also had to be reordered.** The cheap proof (a bounded search
for a path to an already-named class) now runs *before* the fallback that reads
the polynomial. Run the weak name first and it looks like an answer: a class one
mutation away from a quipu class gets labelled and is never searched again, which
is exactly what happened. `nameClassesFromTheorem` → `resolveMergeCandidates` →
`nameRemainingClasses`.

**One more thing this settles.** Every quipu of order n is realised by an LNA with
almost separate relations, so the quipu theorem names *every* quipu class in the
table. A class it leaves unnamed is therefore not a quipu class, and if it turns
out to be one after all, it must share a Coxeter polynomial with the named copy —
so it is always in a merge-candidate group, and the expensive search for a
relation-free quiver of its own (`--form-depth`) can never be what finds a quipu.
It is off by default for that reason.

E-017. Tests: `tests/test_merge_decisions.py`.

---

## F-017 — A move rule is a local rewrite, and the right encoding is per arrow
*2026-09-15*

H-009 suspects the move table is a one-dimensional cellular automaton, and names
its own caveat: a rule whose applicability depended on the row far away would not
be local at all, whatever it looked like. It does not. Both halves of that check
pass.

**Applicability is local.** `lnaMoves.matchesAt` scans the whole relation-length
row, but it does not need to. The same answer comes out of

  (a) the window's own cells -- a relation starting at an offset inside the
      window and staying inside it has at most `width` arrows, so what a cell
      inside can say is bounded by the width; and
  (b) **one bit**: whether any relation covers the window's first arrow having
      started strictly before it.

Checked exhaustively over every rule in `VERIFIED_MOVES`, every admissible LNA
and every window position at lengths 5 to 9: **991,064 comparisons, 1234 of them
an actual match, and no disagreement anywhere.**

**And the mutations a local match licenses are legal.** Re-verifying the whole
table where each rule fits gives **1218 confirmations and zero failures** --
every position where a rule matches, the mutation sequence is admissible at every
step, lands on the predicted relation lengths, and keeps the Coxeter polynomial.
So the admissibility side condition is not an extra non-local hypothesis riding
along; wherever the local pattern holds, the rewrite is legitimate.

**What the one bit says about the encoding, and this is the useful part.** Bit
(b) is exactly what the per-vertex encoding *cannot* supply locally: cell `i`
holds the number of arrows in the relation starting at vertex `i + 1`, and that
number is unbounded in `n`, so a relation can reach arbitrarily far to the right
and no fixed neighbourhood of cells sees it coming. Two consequences:

* the per-vertex row is **not** a good CA state -- unbounded alphabet, unbounded
  reach;
* a row indexed by **arrows** rather than vertices, each carrying whether it is
  covered and whether a relation starts or ends there, has a **fixed alphabet**
  and makes bit (b) a property of the cell at the window's edge. Over that
  encoding the move table is local with a margin of one cell.

So the CA reading is about the right object, provided the state is the arrow row.
That also says where the analogy will strain: converting a per-vertex row to a
per-arrow one needs to know how many relations are open at each arrow, which is
a counter, and the count is bounded only under *almost separate* relations
(overlap at most one arrow, so at most two). For the heavily overlapping LNAs --
which is exactly where the classification still needs a search (H-003) -- the
translation is not finite-state.

**Evidence.** `tests/test_lna_moves.py`,
`test_whether_a_move_applies_is_a_local_condition` at lengths 5 to 7 with 8 and
9 marked slow, and `test_each_rule_holds_wherever_it_applies` for the legality
half. E-016.

---

## F-016 — The search's gate is the paper's criterion now, and it costs nothing
*2026-09-15*

The admissibility gate was a **stricter** reading of the paper's theorem than the
theorem states, and it was applied twice:

* `mutation.mutationIsPossibleAtVertex` walked the relations and refused a vertex
  on any minimal zero relation whose last arrow left it and whose truncation was
  not itself written as a relation;
* `search.mutationSearchDepthFirst` then counted, for every predecessor `v` and
  every arrow `i -> w`, the paths `v -> i` against the paths `v -> w`, and
  refused the vertex if any one arrow lost a path.

Both are the same idea and both are "**every** arrow out of the vertex must keep
every nonzero path nonzero". The paper's theorem rules mutation out only when a
nonzero path ending at the vertex dies against **every** arrow out of it — so
where the two differ, the old gate was refusing mutations the paper allows. Both
also decided "nonzero" syntactically: a zero relation written inside the path, or
a path count taken up to the commutativity relations.

The gate is now `procedure.isMutable`, the theorem read exactly with "nonzero"
decided over the ideal, and it is the only gate — the inline path count is gone.

**What it changed.**

| | |
|---|---|
| mutations newly allowed, n = 5 and 6 to depth 3, n = 7 to depth 2 | **280** |
| of those, Coxeter polynomial preserved | **280** — all of them |
| mutations the new gate refuses that the old one allowed | **0** |
| n = 6 classification | same 4 classes, same sizes |
| n = 7 classification | same 6 classes, same sizes |
| n = 8 classification | same 11 classes, same sizes |
| n = 7 classification time | **38 s → 20 s** |

Faster, which is the opposite of what a more permissive gate suggests: the old
gate cost two syntactic sweeps per vertex — `allRelsInPathAlgebra` at every
search node, then a path count per (predecessor, successor) pair — where the new
one answers from the ideal directly. Removing the sweeps more than pays for the
extra branches.

**Why every newly allowed mutation had to be checked.** The criterion is
*necessary and not sufficient*: the theorem's hypothesis is
`Hom(P_i*[1], Lambda) = 0`, and the paper says explicitly that this is in
general not equivalent to a condition on the quiver. So a mutation the gate
newly allows could in principle fail to be a derived equivalence, and the
Coxeter polynomial would move across it. None does. This is the check R-005
exists to insist on, and being stricter than the paper was the old gate's way of
avoiding having to make it.

**What did not change.** `A_{7,(2,4)}^{(3,3)}` still reaches no relation-free
quiver within depth 8 — the loosened gate does not rescue it, so the point of
`test_the_theorem_answers_where_the_search_gives_up` stands: there are classes
the theorem names outright that no affordable search reaches.

**Evidence.** `tests/test_procedure.py` holds the old criterion, as
`strictlyMutable`, and checks the two against each other: they agree on every
vertex of every LNA of lengths 4 to 6 (a line's vertices have one arrow out,
which is where they coincide), and after one mutation every disagreement is in
the permissive direction, at a branching vertex, with the Coxeter polynomial
preserved. The wider sweep is E-015.

---

## F-015 — The procedure on coefficients, and the two relations the old one missed
*2026-09-14*

`procedure.py` is steps 1-7 of arXiv:2112.08129 and the cleanup after step 7,
with relations as `relationAlgebra` combinations instead of sets of paths. The
coefficients come out of the steps rather than being guessed back, and step 7 is
computed as the kernel the paper says it is. It is what
`mutation.quiverMutationAtVertex` and `reduction.reducePathAlgebra` now run.

**Agreement with the implementation it replaced**, gated on
`mutationIsPossibleAtVertex` throughout so both walk the same mutations:

| check | cases | differences |
|---|---|---|
| one mutation at every admissible vertex of every LNA, n = 4..8 | 2349 | 0 |
| depth-3 walks, n = 5 | 1446 steps | 0 |
| depth-3 walks, n = 6 | 7496 steps | 0 |
| depth-3 walks, n = 7 | 37470 steps | **2** |
| the exact cleanup on the old steps' output, n = 5,6 depth 3 and n = 7 depth 2 | 13652 | 0 |

**The two differences are relations the old implementation missed**, and the
paper's step 7 says so — see R-007. Both are at n = 7 after three mutations:

    A_{7,(1,4)}^{(4,3)} = 40030, mutated at 1, 4, 2
      old: 1;2;5;6;7 | 3;2;5;4 | 5;4;7 = 5;6;7
      new:   2;5;6;7 | 3;2;5;4 | 5;4;7 = 5;6;7     (and 2;5;6;7 implies the old one)

    A_{7,(1,4,4)}^{(4,3)}-ish = 44030, mutated at 1, 4, 2
      new has 2;5;6;7 = 0 in addition to everything the old one has

In the first, `2;5;6;7` is **not** in the ideal the old answer generates, so the
two are different algebras: total dimension 23 against 22, where the algebra
being mutated has 22.

**Neither is caught by the checks that were in place.** Both preserve the
Coxeter polynomial, both survive a there-and-back left mutation, and both reach
`P^(1,2)_(1,0,1)` — the quipu the theorem names for that class — by a depth-7
search. Which is worth recording on its own: *the Coxeter polynomial, the round
trip and the hereditary form can all three be satisfied by an algebra that is
not the mutation*, because all three see only the derived equivalence class, and
a missing relation can leave the class unchanged.

**The coefficients never change a Cartan matrix.** Computing it from the true
combinations and from `fromPathSet`'s guess over the same algebra agrees on every
quiver reached within depth 3 of every LNA of lengths 5 and 6 — 1239 of them,
zero differences. So no published Coxeter polynomial was ever wrong because of
the guess; what the guess cost was the reasoning, as R-003 records.

**It is faster, not slower.** At n = 7 over the 462 admissible single mutations:
the procedure 0.26 s against 0.94 s, the admissibility condition 0.14 s against
1.74 s. Exact linear algebra over the ideal beats the hand-rolled list surgery it
replaces, by 3.6x and 12x. That was the opposite of what was expected, and is
why the switch was affordable.

**The classification is unchanged.** `python classify.py 8` on the new procedure
gives the same 11 classes with the same sizes --
133 + 65 + 64 + 64 + 40 + 26 + 13 + 10 + 9 + 4 + 1 = 429 -- as the partition
recorded in NOTES.md for the old one, and n = 7 gives the same
54 + 32 + 29 + 7 + 6 + 4 = 132. So the two missing relations did not change a
published class: they were lost on quivers the search passes through, not on the
LNAs it records. n = 6 now takes 5 seconds against 13, and n = 7 40 seconds
against 37.

**Evidence.** `tests/test_procedure.py` (agreement at n = 4..6, with n = 7 and 8
marked slow; step 5's minus sign; step 7 as a kernel; the cleanup; the two
recovered relations). The wider runs are E-014.

---

## F-014 — The quipu's end exchanges are already exactly right, and they do not merge the n = 9 pair
*2026-09-14*

The quipu notation is not unique, and one of the re-readings is the one that
matters most here: at the **outermost** foot the main-string end segment and the
cord are the only two branches, so exchanging them — `k_0` with `m_0`, or
`k_{r+1}` with `m_r` — is an isomorphism of the tree and the two LNAs it names
are derived equivalent. **This is not a symmetry the code was missing.** It was
already subsumed by canonicalisation, and it is now checked directly, from both
sides:

1. **Canonicalisation is exactly tree isomorphism.** Over *every* parameter pair
   of each order 3–11 — 4180 names at order 9, 28656 at order 11 —
   `quipuForms.quipuParameters` gives two names the same parameters **if and only
   if** `networkx` says their graphs are isomorphic. Both directions, no
   exceptions. So no end exchange (nor the backwards reading, nor any other
   relabelling) can split a class, and no pair of genuinely different trees is
   being run together.
2. **The LNA-side operations agree with it exactly.** `cor:EquivNakayamaAlgebras`
   is now implemented directly on the relations, independently of `quipuForms`,
   as `LinearNakayamaAlgebra.swapFirstRelation`, `swapLastRelation`,
   `relationDual` and `withoutShortRelations`, with `classPreservingOrbit` for
   their closure. For every LNA of length 4–10 with almost separate relations and
   no length-2 relation, **the orbit of the paper's operations equals the set of
   LNAs that `quipu` names with the same quipu** — 64 algebras in 18 orbits at
   n = 9, 128 in 36 at n = 10, and the largest orbit is 8, which is the paper's
   own bound.

The exchange is at the ends only. Applied to an interior gap `k_i`,
`0 < i < r + 1`, it changes the tree, because the foot there has a third branch
running on along the main string: `P^(1,2,1)_(1,1,2,1)` exchanged at `k_1` is a
different tree of the same order.

**Consequence for the n = 9 pair.** The two classes R-006 says the workbook
over-merged are not related by any of this:

    P^(1,4)_(1,0,1)   (3060000)   diameter 6   -- exchanges to P^(1,1)_(1,0,4)
    P^(1,2)_(1,1,2)   (3004000)   diameter 5   -- both exchanges fix it

The end exchange at the last foot of the first one does fire, and lands on
`P^(1,1)_(1,0,4)` = `3030000`, which canonicalises straight back to
`P^(1,4)_(1,0,1)` — so it adds `3030000` and `0003030` to that class, where they
already are. It does not reach the other quipu, and nothing can: the trees have
different diameters, so they are not isomorphic, and for hereditary algebras of
tree type the underlying tree *is* the derived equivalence class.

The same split is **already in the published table** one vertex down:
`P_(1,0,3)^(1,1)` and `P_(1,1,2)^(1,1)` are two separate rows of the n <= 8
classification in arXiv:2305.06642, and adding a vertex to the last cord of each
gives exactly this pair. n = 9 is not a new claim — it is the published one at
the first order where the Coxeter polynomial can no longer see it (F-010).

**Evidence.** `tests/test_quipu_symmetry.py`, 44 tests at orders/lengths up to 9
plus three marked slow at 10 and 11. Reproduce with

    .venv/bin/python -m pytest tests/test_quipu_symmetry.py -q

**What this does not settle.** Both routes above descend from `thm:QuipuToAn`
being correctly inverted, and the third route — reaching a relation-free quiver
by mutation — is out of range for this pair: searches from all four long-relation
members of `P^(1,4)_(1,0,1)` and both of `P^(1,2)_(1,1,2)` reach no hereditary
quiver at depth 6 (E-013). A genuinely theorem-free separation would need a
derived invariant computed from the algebra itself; the Avella-Alaminos–Geiss
invariant for gentle algebras is the candidate (NOTES idea 22).

---

## F-013 — The pair slide holds for every relation length, with two mutations
*2026-09-14*

Two relations of equal length `l` starting at consecutive vertices (maximally
overlapping) slide one arrow along the quiver, provided no other relation shares
an arrow with their span:

- two **right** mutations at the first relation's source move the pair one arrow
  **left**;
- two **left** mutations at the second relation's target move it one arrow
  **right**.

The mutation count is **two whatever `l` is**; only the window widens, to
`l + 2` arrows.

**Evidence.** `lnaMoves.verifyMove` at every window position of every LNA, for
`l = 2..7`, over four lengths each (`l+3` to `l+6`): 22 confirmations per
direction per length, no failures. Generated rather than listed, by
`lnaMoves.pairSlideRules`.

*Amended 2026-09-14.* The evidence above stands, but the **regression test** for
it did not: `test_each_rule_holds_wherever_it_applies` asked for lengths 5 to 8
whatever the rule, which is `l+3 .. l+6` only for `l = 2`. The window of the
pair slide is `l + 2` arrows wide and needs `l + 3` vertices to sit in, so the
rules for `l = 6..9` got **zero confirmations** at those lengths and the test
failed on all eight of them — for asking where they cannot occur, not for
anything wrong with the rules. It had been failing since the family was
generated up to `l = 9`; it is in the slow set, which is why it went unseen. The
test now derives its lengths from each rule's own width, and a second test
asserts that it always can.

Note the direction: right mutations slide the pair **left**. See R-004.

Discovery found only `l = 3` and `l = 4` on its own, because those are the
window widths that fit at the lengths being searched — see H-001, which this
confirms.

---

## F-012 — Certificates of non-piecewise-heredity propagate by deleting vertices
*2026-09-13*

Corollary `removevertex` of [arXiv:2310.08346](literature/2310.08346-non-piecewise-hereditary-nakayama.md):
if a Nakayama algebra is piecewise hereditary, so is the algebra obtained by
deleting any one vertex. Contrapositively, if some one-vertex deletion is not
piecewise hereditary, neither is the algebra.

Reach, against the paper's two direct criteria alone:

| n | LNAs | direct | with deletion |
|---|---|---|---|
| ≤ 8 | | 0 | **0** |
| 9 | 1430 | 1 | 1 |
| 10 | 4862 | 22 | 24 |
| 11 | 16796 | 265 | **308** |

**Evidence.** `piecewiseHereditary.notPiecewiseHereditaryByDeletion` over every
LNA at each length. The zero below length 9 is the load-bearing check: the paper
states every LNA of length ≤ 8 is piecewise hereditary, so any certificate there
would mean the deletion construction drops or extends the wrong relation.

**What this does not give.** Deleting a vertex preserves piecewise heredity
*only* — not the derived equivalence class, the Coxeter polynomial, or the quipu.
The certificate says exactly "this class is not a quipu class", so it can
separate a class from every quipu class but can never merge two classes.

---

## F-011 — Every n = 9 class is named, and there are exactly 20
*2026-09-13*

*Amended 2026-09-15: confirmed on the corrected engine, after F-018. The first
re-run gave 22, from a circular use of the Coxeter polynomial in the merge step
and not from the engine; the two extra classes were `C(2,3,5)` and `C(2,2,6)`,
which are the domestic weight types of `P^(5)_(1,2)` and `P^(1,1)_(1,3,1)`. The
count and the partition below are unchanged.*

The 1430 LNAs of length 9 fall into **20** derived equivalence classes:

- **18** quipu classes, one per quipu of order 9;
- **1** class (`3345000`, 8 members) of canonical type, tubular weight `(2,4,4)`;
- **1** class (`3033030`, 1 member) not piecewise hereditary at all.

**Evidence, three independent routes agreeing.**

1. `classifyLength(9)` assigns all 1430 rows and leaves nothing as a candidate.
2. The count: 18 quipus of order 9 exist; arXiv:2310.08346 states exactly one LNA
   of length 9 is not piecewise hereditary; the remaining class is of canonical
   type. 18 + 1 + 1 = 20.
3. **20 is exact, not a lower bound.** Two classes can only merge if they share a
   Coxeter polynomial. Across the 20 the polynomials are distinct except for the
   cospectral pair of F-010, which is provably two classes.

`3033030` is precisely the quiver `(**)` of arXiv:2310.08346 — relations
1→4, 3→6, 4→7, 6→9 — checked relation by relation.

---

## F-010 — The Coxeter polynomial fails exactly at cospectral quipus
*2026-09-13*

The Coxeter polynomial of the path algebra of a tree is determined by the tree's
adjacency spectrum. So two **cospectral** non-isomorphic quipus give algebras that
are not derived equivalent yet share a Coxeter polynomial — and among these
classes that is the *only* way it can fail.

| order | quipus | collision groups | quipus involved |
|---|---|---|---|
| ≤ 8 | 2..11 | **0** | 0 |
| 9 | 18 | 1 | 2 |
| 10 | 36 | 2 | 4 |
| 11 | 64 | 4 | 8 |
| 12 | 127 | 13 | 27 |
| 13 | 241 | 30 | 61 |

**Evidence.** `quipuForms.cospectralQuipuGroups(n)`. For orders 4 to 11 the
groups it finds are *exactly* the groups with equal Coxeter polynomials computed
the expensive way, through each algebra's Cartan matrix.

**This is why the published n ≤ 8 table is clean**: below order 9 there are no
cospectral quipus, so grouping by Coxeter polynomial is right there and nowhere
else. The smallest collision, at order 9:

    P^(1,4)_(1,0,1)  =  A_{9,(1,3)}^{(3,6)}   (class 3060000)
    P^(1,2)_(1,1,2)  =  A_{9,(1,4)}^{(3,4)}   (class 3004000)

Reproduce with `python classify.py 9 --collisions`.

---

## F-009 — The length-12 crash was cycles, not length
*2026-09-13*

`allRelsBetweenVertices` and `extendRel` recurse along the arrows out of a vertex
and tracked no visited set, so any cycle in the quiver made them descend forever.
`mutationSearchDepthFirst` then made that fatal rather than local by computing all
relations of a quiver at the top of every node *before* testing it for cycles, so
the first mutation producing a cyclic quiver killed the whole search on the next
node.

This explains the shape of the reported failure exactly: the crash appeared only
when running length 12, yet the offending quiver "didn't necessarily consist of
12 vertices". What mattered was that a mutation at that length finally produced
a cycle.

**Evidence.** `allRelsBetweenVertices` raises `RecursionError` on the three-cycle
`1→2→3→1` with `[1,2,3] = 0`. Fixed by bounding both recursions to simple paths
and testing for cycles before enumerating relations; the full n = 7 search is
byte-identical before and after, so nothing else moved. Pinned in
`tests/test_cycles.py`.

---

## F-008 — `reducePathAlgebra` preserves the algebra
*2026-09-13*

Reduction changes the quiver — it cancels arrows against inadmissible relations —
but must present the same algebra, and the Cartan matrix is the sharpest cheap
witness since it is indexed by vertices that do not move.

**Evidence.** Every legal mutation of depth ≤ 3 out of every LNA of lengths 5 to
8: **38095 reductions, zero changes** to the Cartan matrix, and the exact and
heuristic Cartan matrices agree on all of them.

---

## F-007 — Mutation reachability at bounded depth is directional
*2026-09-13*

`mutationSearchDepthFirst` walks only **right** mutations, so A can reach B at
depth `d` while B reaches nothing at that depth. Since
`rightMutate(dual(P)) = dual(leftMutate(P))` and the relation dual of an LNA is
derived equivalent to it, a right-mutation search out of `dual(X)` covers the
left-mutation directions out of `X`, and everything it reaches is in `X`'s class.

**Evidence.** At n = 8, seeding leaves class `340030` = A_{8,(1,2,5)}^{(3,4,3)}
reaching nothing classified at depth 7; searching from its relation dual finds
the link at depth 6 and closes the classification.

---

## F-006 — Two-path relations are differences, not sums
*2026-09-13*

A relation written as two paths means commutativity throughout this codebase —
`applyRelSetToPath` substitutes one for the other — so it means `p - q`, not
`p + q`.

Not cosmetic: three commutativity relations among three parallel paths read as
sums give `p = -q`, `r = -q` and `p + r = -2q`, forcing `q = 0` and collapsing a
Hom space that should be one-dimensional.

**Evidence.** Found via F-008: nine of 38095 reductions appeared to change the
Cartan matrix, and every one was this sign reading rather than a fault in the
reduction. See R-002.

---

## F-005 — Signs alone do not close; coefficients must be integers
*2026-09-13*

From `p + q + r = 0` and `p - q = 0` follows `2p + r = 0`, which has no expression
with coefficients in `{-1, 0, +1}`. The first step that combines two relations
leaves the sign-only world.

Integer coefficients are the smallest choice that closes and cost nothing over
signs — the same dict with a wider value type. The row reduction in
`relationAlgebra` is over the rationals already, so moving to a field would change
nothing.

---

## F-004 — The set-of-paths model misses zero relations through chained squares
*2026-09-13*

In the 2×2 commutative grid, the two commutativity relations make all three paths
from 1 to 6 equal, so adding `[1,2,3,6] = 0` kills all three. `pathHasZeroRel`
recognises only `[1,2,3,6]`, because it looks for a zero relation sitting
contiguously inside the path. `relationAlgebra.isInIdeal` gets all three.

A second symptom: `reducePathAlgebra` turns `{p,q,r}` plus `{p,q}` into `{p,q}`
and `{r}`, valid for `p+q+r=0, p+q=0` but not for `p+q+r=0, p-q=0` where the
answer is `2p+r=0`; `numberOfPathsUpToRels` reports 2 for that same algebra, so
the two halves of the code disagree about one object.

**Nothing published is affected**: see F-008.

---

## F-003 — The quipu theorem inverts, naming a class in O(1)
*2026-09-12*

Theorem `thm:QuipuToAn` of arXiv:2305.06642 read backwards: an LNA
`A_{n,(n_0..n_r)}^{(l_0..l_r)}` with almost separate relations comes from the
quipu with

    m_i = l_i - 2,  k_0 = n_0,  k_i = n_i - n_{i-1} - m_{i-1} - 1,
    k_{r+1} = n - n_r - m_r - 1

so it names its own class with no searching. Relations of length 2 are dropped
first, since they do not change the class.

**Evidence.** Round-trips on every quipu the paper names; constant on every
published class; agrees with what a mutation search reaches wherever the search
can reach an answer. Coverage falls with length (100% of rows at n = 4, 54% at
n = 8, 19% at n = 12) but the set of quipus it names does not: it finds every
class at every length checked.

---

## F-002 — Step 3's cyclic case is not implemented
*2026-09-12*

For a minimal relation `r: i ⇢ i` the procedure calls for one arrow
`(α r̄): i* → t(α)` per arrow `α` out of `i`; the code adds a single arrow from
`r`'s source to its target, which on a cycle is a loop on `i*`. The second worked
example of arXiv:2112.08129 does not reproduce.

Pinned as a strict xfail in `tests/test_mutation_procedure.py`. The LNA search
never reaches it, since it stops descending at the first cycle, so no published
result depends on it.

---

## F-001 — A NetworkX change had silently broken every Coxeter polynomial
*2026-09-12*

Since NetworkX 3.1, `all_simple_paths(G, v, v)` yields the trivial path `[v]`
where it used to yield nothing. `numberOfPathsUpToRels` counted it, putting 2
instead of 1 on the Cartan diagonal, so `det(C) = 2^n` instead of 1 and every
Coxeter polynomial came out with fractional coefficients.

**Evidence.** After the fix, A_5 gives `1 + λ + … + λ^5` and D_5 gives
`(λ+1)(λ^4+1)`, as they must.
