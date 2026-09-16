# Hypotheses

Things suspected but not established, newest first, each with what would settle
it. Status is one of `OPEN`, `SUPPORTED`, `CONFIRMED → F-nnn`, `REFUTED → R-nnn`,
`PARKED`. See [`README.md`](README.md).

---

## H-009 — A short mutation sequence breaks the maximally overlapping pair
*2026-09-16* · **OPEN** — this is what the overnight discovery run is for

F-015 measures the obstruction: the LNAs that seeding and the move orbits cannot
place all carry an overlapping run, and the single commonest is the maximally
overlapping pair of length-3 relations, `(1:3) (2:3)` — 26% of the unplaced rows
at n = 8, 20% at n = 9. No rule in the table reduces the overlap of a bare pair.

The hypothesis is that a sequence of **four or five** mutations does, and that it
is invisible to every search run so far because all of them were bounded at
three. H-008 is the general form of this; the pair is the instance worth
settling first, because it is the one that pays.

**Evidence, such as it is.** A four-step search at quiver length 11 does find a
link out of the pair to a less overlapping LNA, `(0:3)(1:3) -> (0:2)(3:3)` via
`[2, -5, -6, 2]`, which three steps do not. That particular rewrite is **false**
as a rule — it fails at lengths 10, 11 and 12 (R-008) — so this is evidence only
that four-step links out of the pair exist at all, not that a valid one does.

**What would settle it.**

    python discover.py 9 --steps 4 --jobs 8 --resume
    python discover.py 9 --steps 5 --jobs 8 --resume

Plant the commonest blocking runs in the interior of A_13 and A_14, walk
sequences of four and then five mutations near them, keep only what lands on a
strictly smaller overlap, and verify each over lengths 7 to 10 with a minimum
confirmation count. A negative result is worth as much as a positive one: if
five mutations do not break the pair either, the obstruction is not a matter of
search depth and the move table is the wrong lever for getting past n = 11.

**What it would be worth.** A single valid rule for `(1:3) (2:3)`, with the
family in relation length that F-013 suggests such a rule would have, would place
a fifth of the rows that currently need a depth-6 search — at n = 9 and, being
local, at every length above it.

---

## H-008 — Rule families are parameterised by relation length and overlap, and some need more than three mutations
*2026-09-14* · **OPEN**

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
*2026-09-13, measured 2026-09-16* · **SUPPORTED, and sharpened — see F-015**

Seeding by the quipu theorem and expanding along move orbits places 72% of the
n = 7 table, 57% of n = 8, 45% of n = 9 — and the 34 rules added at length 8
barely moved those numbers. The diagnosis is that the rules found so far mostly
keep an LNA *inside* the almost-separate set the theorem already covers, while
the rows still needing a search are the heavily overlapping ones.

**Measured** (E-013, F-015, `python unplaced.py 9`). The diagnosis is right and
the measurement makes it precise: every unplaced row carries an overlapping run,
and one shape — the maximally overlapping pair `(1:3) (2:3)` — blocks a fifth to
a quarter of them.

**But "wider" is the wrong word for what is missing.** Eleven of the 64 rules do
reduce overlap; they are simply too specialised to fire, needing a third relation
in the window or relations of length 2, so only 49 of the 820 LNAs of length 9
with an overlap of 2 or more ever reach a smaller one. What is missing is not a
wider rule but one that applies to a **bare overlapping pair**. That is H-009,
which is the successor to this hypothesis.

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
