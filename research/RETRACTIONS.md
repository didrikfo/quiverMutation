# Retractions

Things believed and then found false, newest first. These are kept because a
recorded wrong belief stops the same reasoning being repeated; a tidied-away one
does not. See [`README.md`](README.md).

---

## R-007 — "Step 7 of the mutation procedure is about the new arrows themselves"
*2026-09-14* · corrected against the paper's own words, → **F-015**

The summary in `literature/` had step 7 as

> `Σ_r ε_r r̄ = 0` is a relation iff `Σ_r ε_r (r/α) = 0` is one in `Q`, for every
> `α` out of `i`

which reads as a statement about the arrows `r̄`, with scalar coefficients. It is
not. The paper says the `ε_r` are **linear combinations of paths `t(r) → l`**, for
**any** vertex `l`, and says it as an **if and only if**. So step 7 is about every
path out of `i*` — an `r̄` followed by a tail — and it determines them completely:
the relations out of `i*` to `l` are exactly the kernel of the map that
precomposes with each `α*`.

**What the abridged reading cost.** Half a day, and nearly the wrong decision.
Writing the procedure on coefficients, I implemented step 7 as that kernel, found
it disagreed with the trusted implementation in 2 of 37470 walk steps, and could
not tell which was right: both preserved the Coxeter polynomial, both survived a
there-and-back mutation, and both reached the quipu the theorem names. I was
about to revert the switch and file the disagreement as an open question. The
paper's actual sentence settles it in one reading: the kernel is right, and the
old implementation was **missing relations**.

**The lesson, and it is the point of `literature/` existing.** A summary that
abridges a statement can be worse than no summary, because it reads as
authoritative. The rule going in: **quote the statements the code implements,
verbatim.** Step 7 is now quoted in full, and the two things the abridgement
dropped — that the coefficients are paths, and that it is an iff — are called out
under the quote, because both are what made it misleading.

**Also corrected:** the same summary said the admissibility criterion was a
condition for mutation being *allowed*. The paper gives it as two cases where
mutation is *impossible*, and says explicitly that the homological condition is
in general **not** equivalent to a condition on the quiver. Ruling out is not the
same as ruling in, which is why the stricter of our two implementations stays the
search's gate.

---

## R-006 — "The workbook's n = 9 classification has 19 classes"
*2026-09-13* · superseded by **F-011**

**Challenged 2026-09-14, upheld, challenge withdrawn 2026-09-15 → F-014.**
The objection was that the two quipus are isomorphic after all, being related by
the exchange of an end segment of the main string with the cord at the outermost
foot. That exchange is real and was already implemented; it does not relate
these two, whose trees have diameters 6 and 5, and the analogous pair is two
separate rows of the paper's own n <= 8 table. The objection was withdrawn by
its author the next day as a misread — the quivers had been manipulated in the
head rather than on paper.

Kept, with the outcome, for two reasons. The exchange **is** a real symmetry of
the notation and mistaking its reach is an easy error to repeat — F-014 now pins
where it does and does not apply, including that it fails at an interior gap.
And it is the record of a doubt that was answered rather than left hanging, which
is the more useful half of "nothing is deleted".

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
