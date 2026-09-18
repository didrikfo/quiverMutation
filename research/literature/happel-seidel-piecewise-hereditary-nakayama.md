# Piecewise hereditary Nakayama algebras

**D. Happel, U. Seidel**, *Piecewise hereditary Nakayama algebras*,
Algebr. Represent. Theory **13** (2010), no. 6, 693–704.
**Not on arXiv, no preprint found.** DOI `10.1007/s10468-009-9169-y`.

> ⚠ **Written from secondary sources, 2026-09-18.** The original was not
> reachable. Everything below is quoted from papers that restate it and reprove
> it — principally **Lenzing–Meltzer–Ruan, arXiv:2112.15587** (whose §4 is an
> independent second proof and whose Theorem 4.4 *is* the classification),
> **Ueda, arXiv:2302.02880**, **Lenzing–de la Peña, arXiv:0805.1018 §5.4**, and
> **arXiv:2310.08346**. Theorem numbers below are those of the citing paper, not
> of Happel–Seidel, except where a Happel–Seidel number is given explicitly.
> Where the sources disagree, or where only one source says something, it is
> flagged.

This is the paper the README's candidate list describes as "a complete answer in
that setting". **That description is correct**, and this file records what the
answer is. The setting is the ideal being a power of the radical: `Λ(n,r) =
kA_n / rad^r`, over an **algebraically closed** field. In our digit strings that
is `"r" * (n − r) + "0" * (r − 2)`, a string of `n − 2` digits.

## The answer

Restated from Lenzing–Meltzer–Ruan Theorem 4.4, which is explicitly "[11] The
complete classification of piecewise hereditary Nakayama categories". Standing
hypotheses `n ≥ 5`, `r ≥ 3` (the two excluded families are stated separately
below). `[a,b,c]` = the hereditary star with branch lengths `a,b,c`; `(a,b,c)` =
the weighted projective line `X(a,b,c)`.

**Piecewise hereditary of module type** — these are the ones in a **quipu class**,
and every one of them is a *star*:

| `Λ(n,r)` | type |
|---|---|
| `N_{r+2}(r)`, `r ≥ 3` | `[2,3,r−1]` |
| `N_{r+3}(r)`, `r ≥ 3` | `[2,3,r]` |
| `N_7(3)` | `[2,3,4]` |
| `N_8(3)` | `[2,3,5]` |
| `N_8(4)` | `[2,4,4]` |
| `N_9(3)`, `N_9(5)` | `[2,3,6]` |
| `N_10(5)` | `[2,3,7]` |

**Piecewise hereditary of sheaf type** — piecewise hereditary but derived
equivalent to a *canonical* algebra, so **not** in any quipu class:

| `Λ(n,r)` | type |
|---|---|
| `N_{r+4}(r)`, `r ≥ 4` | `(2,3,r)` |
| `N_9(3)`, `N_9(6)`, `N_9(7)` | `(2,3,5)` |
| `N_10(3)` | `(2,3,6)` |
| `N_11(3)`, `N_11(6)` | `(2,3,7)` |
| `N_9(4)` | `(2,4,4)` |
| `N_10(4)` | `(2,4,5)` |

The two lists overlap exactly in the tame cases, where the weighted projective
line is domestic and the canonical algebra is derived equivalent to an affine
quiver:

- `(2,3,4) ≃ N_8(4) ≃ [2,4,4]`, affine `Ẽ7`;
- `(2,3,5) ≃ N_9(3) ≃ N_9(5) ≃ N_9(6) ≃ N_9(7) ≃ [2,3,6]`, affine `Ẽ8`.

**Excluded from the table, true separately:** `Λ(n,2) = kA_n/rad^2` is piecewise
hereditary of type `A_n`; `Λ(n, n−1)` is piecewise hereditary of type `D_n` for
`n ≥ 4`. (The first agrees with our own `corollary:lengthtworelations` — all
relations of length two, hence derived equivalent to `kA_n`; the second is the
"one short relation at the far end" case.)

**Everything else is not piecewise hereditary.** That is the content of the word
"complete", and it is what makes this usable as a certificate.

## Why the list is finite — Lemma 4.2

The mechanism, as Lenzing–Meltzer–Ruan reconstruct it:

> If `Λ(n,r)` is piecewise hereditary then so is `Λ(n−1, r)`.

`P_n` is exceptional in `D^b(Λ(n,r))` and its right perpendicular category is
`D^b(Λ(n−1,r))`; heredity passes to perpendicular categories of exceptional
objects. In digit strings this drops one leading copy of `r` and **keeps `r`** —
sharper and simpler than the `removevertex` corollary of arXiv:2310.08346, but
only available in the radical-power family.

Contrapositive: **for each fixed `r`, once some `n₀` fails, every `n ≥ n₀`
fails.** So for each `r` there is a threshold, and the table is the region to the
left of it. Lenzing–Meltzer–Ruan's "red wall" is that threshold: the first
non-piecewise-hereditary member at each `r` is
`N_12(3), N_11(4), N_11(5), N_12(6), N_12(7)`, and `N_{r+5}(r)` for `r ≥ 7`.

Those five are exactly the list arXiv:2310.08346's `liste:baseforHS` recovers by
its own means, **without algebraic closedness**.

## The Happel–Seidel symmetry

The other thing the paper is cited for, and the only *merge* in it:

> For `a, b ≥ 2` and `n = (a − 1)(b − 1)`, `D^b(Λ(n,a)) ≃ D^b(Λ(n,b))`.

Lenzing–Meltzer–Ruan extend it to `n ± 1` (their Proposition 4.1(2),(3)) and Ueda
reproves the `n+1` case without algebraic closedness (his Corollary 5.17). See
[`2112.15587-nakayama-fuchsian-singularities.md`](2112.15587-nakayama-fuchsian-singularities.md)
for the statement we actually use.

Two further Happel–Seidel results cited elsewhere and *not* recovered here in full:

- **Happel–Seidel Prop. 2.3**, cited by Ueda: `per Λ(s+6, s+4) ≃ per Λ(s+6, s+3)`
  for `s ≥ 0`, both being derived equivalent to the star `[2, 3, s+3]`. This is
  the `q = 1` case of Ueda's Corollary 1.3, and it is the `N_{r+2}(r) / N_{r+3}(r)`
  pair of the table read as a merge: `N_{r+2}(r) ≃ [2,3,r−1]` and
  `N_{(r−1)+3}(r−1) ≃ [2,3,r−1]`, so `Λ(n, r) ≃ Λ(n, r−1)` whenever
  `n = r + 2`. Checked cospectral here for `n ≤ 12`.
- **Happel–Seidel Prop. 2.6**, cited by a 2024 paper as giving algebras that are
  *not* piecewise hereditary. Not seen; presumably the threshold argument.

## Caveats and limitations

- **Algebraically closed field.** Happel's trichotomy is the backbone. What
  survives without it is only the part arXiv:2310.08346 reproves by hand.
- **Radical powers only.** The paper says nothing about an LNA with two different
  relation lengths. Every criterion we have for the general case is separate work.
- **"Table 1" is not reproduced here.** The tables above are Lenzing–Meltzer–Ruan's
  restatement. They assert it is Happel–Seidel's classification, and
  Lenzing–de la Peña's independent `r = 3` column (`0805.1018` Prop. 5.8: the types
  of `Λ(n,3)` for `n ≤ 12`) agrees with it entry for entry, which is the best
  cross-check available without the original. **If someone gets the original, check
  the sporadic entries** — those are where a transcription error would hide.
- The two type systems are easy to confuse. **Sheaf type is still piecewise
  hereditary**, so a sheaf-type `Λ(n,r)` is *not* a counterexample to anything
  about piecewise heredity — but it is also **not in a quipu class**, since a quipu
  class means derived equivalent to a hereditary *algebra*. A criterion that only
  asks "piecewise hereditary?" therefore under-certifies for our purposes; the
  useful question is "of module type?".
- The bare statement "`Λ(n,r)` is piecewise hereditary" is the *weaker* half of
  what the table says. The type is the strong half, and it is the half that merges
  classes: every pair in the same row of either table is derived equivalent.

## What it does not give us

- No invariant. The classification is a list, not a test.
- Nothing about relations of unequal length, i.e. nothing about `3033030` or its
  neighbourhood.
- No derived functors; the proofs (as reconstructed) are perpendicular-category
  arguments inside `D^b(coh X)`.

## How it lands on our problem

Two uses, both immediate:

1. **A ground-truth table.** For every `n` and `r ≥ 3`, the table above says
   whether `"r"*(n−r) + "0"*(r−2)` is in a quipu class, which quipu (a star
   `[a,b,c]`), or neither. Our classification must agree with it wherever both are
   defined, and it extends far past where we can search.
2. **Merges within a row.** `N_9(3) ≃ N_9(5) ≃ N_9(6) ≃ N_9(7)`, i.e.
   `3333330 ≃ 5555000 ≃ 6660000 ≃ 7700000` — four LNAs of length 9 in one class,
   named `[2,3,6]`, which is the star with branches 2, 3, 6, a quipu of order 9.
   That single row is the cleanest test case we have for a literature merge that
   no local move produces.
