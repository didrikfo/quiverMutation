# Hypotheses

Things suspected but not established, newest first, each with what would settle
it. Status is one of `OPEN`, `SUPPORTED`, `CONFIRMED → F-nnn`, `REFUTED → R-nnn`,
`PARKED`. See [`README.md`](README.md).

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
quipu class at `n = 9` and every one of the 262 at `n = 10` reaches a quipu with
relations **within three mutations** — none reaches none, and the median LNA
reaches 16 and 20 of them respectively (F-034). What is open in part 1 is `n = 11`
and beyond, and whether it is a theorem rather than a run of small cases.

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
every LNA outside a quipu class at `n = 11` and `n = 12` — `n = 10` is done, and
takes 7 minutes. An LNA that reaches none at a depth where its neighbours reach
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
*2026-09-17, written before the overnight run* · **OPEN**

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
*2026-09-16* · **OPEN**

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
