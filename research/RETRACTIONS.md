# Retractions

Things believed and then found false, newest first. These are kept because a
recorded wrong belief stops the same reasoning being repeated; a tidied-away one
does not. See [`README.md`](README.md).

---

## R-007 — "Two classes with different names are different classes"
*2026-09-14* · corrected in `mutationClassTable.formsAreCompatible`

`mergeReport` treated any two differing identifying forms as proof that the
classes are distinct. That is false when one is a **quipu that is tame
hereditary** and the other is that quipu's **canonical type**: a tame hereditary
algebra is also derived equivalent to a canonical algebra, so `P^(1,1)_(1,4,1)`
and `C(2,2,7)` are two names for one thing.

At n = 10 this reported a false separation — three classes named `C(2,2,7)`
against the quipu `P^(1,1)_(1,4,1)`, whose own canonical weight type is `(2,2,7)`.
They are merge candidates, not separated classes.

**How it was caught.** The cospectral analysis (F-010) predicts *exactly* two
collision groups at order 10. The run reported two separated groups, one of which
was not either of them. A prediction that did not match is what exposed it; the
count alone would have looked right.

The earlier note that "the two identifications agree rather than compete" was
correct and is in `tests/test_nakayama_classes.py` — but the comparison logic did
not honour it. **A property recorded in a test is not enforced anywhere else.**

Does not affect n ≤ 9: the separated pair at n = 9 is two quipus, and no class
below 10 carries a canonical name against a quipu sharing its polynomial.

---

## R-006 — "The workbook's n = 9 classification has 19 classes"
*2026-09-13* · superseded by **F-011**

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
