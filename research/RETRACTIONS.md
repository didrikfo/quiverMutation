# Retractions

Things believed and then found false, newest first. These are kept because a
recorded wrong belief stops the same reasoning being repeated; a tidied-away one
does not. See [`README.md`](README.md).

---

## R-006 — "The workbook's n = 9 classification has 19 classes"
*2026-09-13* · superseded by **F-011**

**Challenged 2026-09-14 and upheld → F-014.** The objection was that the two
quipus are isomorphic after all, being related by the exchange of an end segment
of the main string with the cord at the outermost foot. That exchange is real,
and it is already implemented; it is checked exhaustively in F-014, and it does
not relate these two — their trees have diameters 6 and 5. The analogous pair is
two separate rows of the paper's own n <= 8 table.

The hand-made classification merged two classes of 18 into one of 36. Every other
class matches the computed partition exactly, and the computed partition refines
the workbook's 65 unmerged classes with no contradiction anywhere — so this is one
merge too many, not a disagreement about the rest.

The over-merged pair is the cospectral one of F-010:
`P^(1,4)_(1,0,1)` (containing `3060000` = A_{9,(1,3)}^{(3,6)}) and
`P^(1,2)_(1,1,2)` (containing `3004000` = A_{9,(1,4)}^{(3,4)}). Each contains an
LNA with almost separate relations naming its quipu, the two quipus are
non-isomorphic trees, so the classes are distinct.

**Why it happened:** merging on equal Coxeter polynomials. That is right for every
length up to 8 and first goes wrong at 9 — exactly where F-010 says it must.

---

## R-005 — "A discovered rewrite that produces the predicted LNA is a valid move"
*2026-09-13* · corrected in `lnaMoves.verifyMove`

The first version of move verification checked only that applying a rewrite's
mutation sequence produced the predicted relation lengths. **38 rules passed on
that basis, and the orbits they generated had the wrong Coxeter polynomial 6561
times out of 8388.**

A mutation applied outside the procedure's admissibility condition still returns
a quiver — just not a derived equivalent one. So a rewrite can land on exactly the
predicted LNA and be entirely false.

Verification now requires three things, and all three are needed: the predicted
result, **every step admissible**, and the Coxeter polynomial unchanged. Under
that, 16 of 67 candidates survived. `movesFrom` also only walks admissible
mutations now.

**The general lesson:** in this setting, landing on the right object is not
evidence of having got there legitimately.

---

## R-004 — "Right mutations slide a maximally overlapping pair to the right"
*2026-09-13* · corrected by **F-013**

As originally described, and as I first implemented it. Empirically they slide it
**left**; the rightward slide is two *left* mutations, and at the **target of the
second relation** (`v + l + 1`), not at `v + l` as I first guessed.

Cheap to get wrong and cheap to check, which is the argument for verifying a rule
against the engine rather than reasoning it out.

---

## R-003 — "A two-path relation can be read as a sum"
*2026-09-13* · corrected by **F-006**

`relationAlgebra.fromPathSet` read every relation in the set-of-paths model as a
sum with all coefficients +1. For a two-path relation that is wrong: two paths
written as a relation means commutativity throughout this codebase, so `p - q`.

Cost: nine apparent failures in the reduction check of F-008, every one of them
this reading rather than a fault in the code being checked. Nearly led to
"reducePathAlgebra is broken", which it is not.

---

## R-002 — "An intermediate quiver already seen need not be explored again"
*2026-09-13* · corrected in `lnaMoves.localMutationSequences`

Pruning on mere membership of a `seen` set. A state first reached deep in one
branch then blocked a later branch that reached it with more steps remaining, so
the search **found fewer LNAs at `maxSteps=3` than at `maxSteps=2`** — which is
how it was noticed. The set now records how many steps were still available.

Worth remembering for any future bounded-depth search here: monotonicity in the
depth bound is a cheap self-check, and it caught this immediately.

---

## R-001 — "Verification can enumerate all tuples of relation lengths and filter"
*2026-09-13* · corrected in `lnaMoves.verifyMove`

That is 10^8 tuples at length 10, against the 4862 LNAs that actually exist there.
Enumerate the admissible LNAs directly. Verification at length 10 went from
impractical to about a second, which is what made re-verifying the length-8 rules
at greater length possible at all.
