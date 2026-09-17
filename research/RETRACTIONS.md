# Retractions

Things believed and then found false, newest first. These are kept because a
recorded wrong belief stops the same reasoning being repeated; a tidied-away one
does not. See [`README.md`](README.md).

---

## R-011 — "The two ends of the quiver are not the same end" (F-025)
*retracted 2026-09-16*

F-025 claimed an asymmetry between the source and the sink of the line: a pair
of relations of unequal length comes apart against the sink and not against the
source. The evidence was that `(0:l) (1:m)` with l < m loses an arrow at the
sink under two right mutations, and that *the same rule written at the left end*
gives 1 confirmation and 3 failures.

**The mirror was the wrong mirror.** Mirroring a rule is not reflecting its
pattern and keeping the mutations; it is the **relation dual** -- reverse every
arrow *and* exchange right mutation for left. Under that transform

```
    right end:  (0:l) (1:m)  ->  (0:l-1) (1:m)         via [2, 2]
    left  end:  (0:m) (m+1-l:l)  ->  (0:m) (m+2-l:l-1) via [-(m+1), -(m+1)]
```

and the left-hand one holds: 4 confirmations, no failures, at every `(l, m)`
tried. What was compared against it instead was the *same pattern* at the other
end, and that is a different configuration -- a pair sharing a start rather than
a pair sharing an end -- so of course it behaves differently. The ends are
mirror images; the pattern was not.

**What survives.** The rule itself, and the family: `(0:l) (1:m)` shortens at the
sink for every 3 <= l < m <= 9, 21 members, 4 confirmations apiece, and its dual
does the same at the source. The probe results behind it are also untouched --
`(1:3) (2:6)` planted at the source really does reach only two LNAs. What is
withdrawn is the *interpretation*: that asymmetry belongs to the pattern, not to
the quiver's ends, and no reading of the procedure is needed to explain it.

**What it cost and what it bought.** The mistake was worth making, because
looking for its cause turned up the transform itself, which the table had never
been closed under: 410 rules were missing their duals, all 410 verify, and they
are generated now (F-026). A wrong mirror is how the right one got written down.
E-026.

---

## R-010 — "Wider rules are what unlock the heavily overlapping LNAs" (H-003)
*retracted 2026-09-16*

H-003's measurement half is right and is now F-021: the rows a classification
search still has to place are exactly the heavily overlapping ones, all 786 of
them at n = 9. Its **diagnosis** -- that the rules found so far are too narrow,
and that wider ones would reach them -- is wrong for the configuration that
dominates the leftovers.

**What was believed.** That discovery kept finding rules which keep an LNA inside
the almost separate set because the windows searched were too small, and that
aiming discovery at bigger patterns would produce rules crossing the line.

**What is true.** An *isolated pair* of relations sharing two or more arrows --
434 of the 786 rows left at n = 9 have no heavily overlapping run longer than
that -- cannot have its overlap reduced by any interior sequence, at any window
width tried. Planted in the middle of A_13, it reaches 8 LNAs at three mutations,
14 at four, 22 at five and 34 at four with the mutations allowed twice as far
out, and **every one of them still has the pair** (F-022, E-021). Widening the
window is not a dial that turns here; it reaches further along the quiver and
finds the same thing.

**Where the diagnosis does hold, and why that misled.** A heavily overlapping run
of *three* relations does dissolve under an interior rule, and the rules that do
it are wide -- the window-5 and window-6 entries found at length 8 and in E-011.
So the belief was confirmed every time it was tested on a triple, and the pair,
which is the commonest core by a factor of four, was never the thing being
tested.

**What corrects it.** The overlap of an isolated pair is reduced at an **end** of
the quiver, not by a wider window: two mutations at the source or sink delete one
of the two relations (F-022). That is a rule the framework could not even state
until it grew anchored descriptions, because it is false at every other position
-- so no amount of searching for *floating* rules, at any width, was ever going
to find it. E-021, E-023.

---

## R-009 — "A rewrite that verifyMove confirms with no failures is a rule"
*retracted 2026-09-15*

`verifyMove` enumerates every admissible LNA of each length it is given, matches
the rule at every window position and checks all three conditions, so
`confirmed > 0 and failures == []` reads like verification. It is not, unless the
lengths are chosen to fit the rule. The lengths have to follow the **window**.

**What it cost.** The three-mutation interior run (E-011) verified everything at
the fixed lengths 7 to 10, the range E-010 had used. A window of 9 arrows fits in
A_10 at exactly one position -- flush against both ends -- so each of its 30 rules
got **one** confirmation from **one** length and was reported as verified.
Re-checked at length 11, where the window can sit clear of the ends, **all 30
failed**: not thin evidence, wrong rules. Example: window 9,
`(1:2) (4:2) -> (0:2) (5:2) (7:2)` via `[2, -7, 9]`, applied to `020020000` at
length 11, does not even land on a line quiver.

This is H-007 again, on the other side. Discovery was moved into the interior of
A_13 and A_14 precisely so that an end effect could not pass for a rule -- and
then the verification put the window back flush against the ends, where every
special case applies at once, and let the end effects through.

**The correction.** Verify at `width + 1 .. width + 4`, so the window has room to
move; require confirmations at **two or more lengths**; and drop a rule that
cannot get them within the affordable range rather than keeping it on one.
`discover.py` derives the lengths from each rule's width (`--verify-span`,
`--verify-cap`), and `tests/test_lna_moves.py::lengthsToCheck` already did this
for the table -- which is why nothing false ever reached `VERIFIED_MOVES`.

**The general lesson.** A count of confirmations is not evidence until you know
how many *positions* produced it. One position is one case, and one case at the
only place a window fits is the worst case there is.

---

## R-008 — "LNAs are gentle, so the Avella-Alaminos-Geiss invariant applies directly"
*retracted 2026-09-15*

Written into `research/literature/README.md`'s candidate list and into NOTES ideas
14 and 22, as the plan for an independent separation of the cospectral pair.

**Wrong twice.** A gentle algebra's ideal is generated by paths of **length two**;
an LNA with a relation of three or more arrows is a string algebra but not a
gentle one. And an LNA that *is* gentle is derived equivalent to the path algebra
of A_n, by operation 2 of `cor:EquivNakayamaAlgebras` applied to every one of its
relations — so the gentle LNAs are a single class and a gentle invariant has
nothing to separate. Neither member of the cospectral pair has a gentle member
among its 18.

**The correction.** F-019. The candidate list and both ideas now say so. What
survives of idea 22 is the fallback it already named: Hochschild cohomology, which
is a derived invariant of any finite-dimensional algebra.

**How it slipped in.** The LNAs are special biserial and monomial on a quiver of
maximum degree 2, which is most of the definition of gentle, and the one
remaining condition is the one that fails. Worth the general lesson: a definition
that is "obviously satisfied except for one clause" is where to look, not where to
stop.

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
