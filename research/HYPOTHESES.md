# Hypotheses

Things suspected but not established, newest first, each with what would settle
it. Status is one of `OPEN`, `SUPPORTED`, `CONFIRMED → F-nnn`, `REFUTED → R-nnn`,
`PARKED`. See [`README.md`](README.md).

---

## H-011 — Every class is reached by walking a heavily overlapping run to an end
*2026-09-16* · **OPEN**

F-022 gives the one mechanism known to reduce the overlap of an isolated pair:
walk it to an end of the quiver with the pair slide and collapse it there. The
suspicion is that this is not one mechanism among several but **the** mechanism
-- that the derived equivalence class of any LNA is reached from an almost
separate one by a sequence of interior moves that carry its heavily overlapping
runs to an end, and collapses there.

**Why it is worth stating.** If it holds, the classification needs no search at
all: the seeding, the slide families and the anchored collapses generate
everything, and n = 10 and beyond become a matter of counting rather than of
mutation. If it fails, the LNA it fails on is the first evidence of a genuinely
different obstruction, which is worth more than another rule.

**Evidence for.** Coverage rises from 45% to 51% at n = 9 on the strength of
sixteen anchored rules, where 123 floating ones had bought 34 rows (F-022).

**Evidence against, and it is the thing to explain.** 704 rows at n = 9 are still
not placed, and 352 of them have a heavily overlapping run of three or more --
configurations the interior rules *do* dissolve in isolation (F-022's table).
So the blockage is not the run itself but what surrounds it: a companion
relation inside the window a rule needs. Two readings, and they want different
work:

1. the runs cannot reach an end because other relations are in the way, and what
   is missing is rules for moving a run *past* a neighbour; or
2. the rules that dissolve a run in isolation have windows too clean to match
   anything real, and what is missing is the same rules with a companion
   relation admitted into the window.

**What would settle it.** For each uncovered LNA, ask which is the case: whether
some sequence of known moves brings a run to an end (then the gap is search
order, not rules), and if not, which relation blocks it. `overlaps.py --cores`
already reports the runs; what it does not yet report is the blocker.

---

## H-010 — Overlap is reducible only at an end, and that is a theorem about the procedure
*2026-09-16* · **SUPPORTED**

F-022 is an empirical statement: no interior sequence found so far reduces the
overlap of an isolated pair, over four probes reaching depth 5 and a margin of 6.
The suspicion is that it is exact -- that **no** interior mutation sequence does,
whatever its length -- and that the reason is visible in the procedure rather
than in the search.

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
an isolated pair. `lnaMoves.localMutationSequences` at a larger `maxSteps` and
`margin` is the search; note that the four probes already run are recorded in
E-021 and should not be repeated.

**Caveat.** "Isolated" is doing work. A pair with a third heavily overlapping
relation is *not* invariant -- it dissolves in three mutations -- so any
statement has to be about the whole overlapping run, and the right form of it is
probably about the run's ends rather than the pair.

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
