# On derived equivalences of categories of sheaves over finite posets

**arXiv:** [math/0610685](https://arxiv.org/abs/math/0610685) · S. Ladkani ·
J. Pure Appl. Algebra **212** (2008) 435–451 · *read 2026-09-19 from the arXiv PDF*

Two halves, and we want both. §3 is a list of **derived invariants**, one of which
is strictly finer than the Coxeter polynomial and which we have now checked
separates cospectral quipus. §4 is a **construction**: for any poset `X` and any
closed subset `Y ⊆ X` it builds an algebra `A_Y` with `D^b(X) ≃ D^b(A_Y)`, and
`A_Y` is an algebra **with zero relations** — sometimes an LNA. That is a machine
for producing derived equivalences between LNAs, and it is the machine
Dong–Lin–Ruan (arXiv:2203.15735) run to get their merges.

Throughout, `1_X` is the incidence matrix (`(1_X)_{xy} = 1` iff `x ≤ y`), which is
the **Cartan matrix** of the incidence algebra `kX`; `D^b(X)` is `D^b(mod kX)`.

## Part 1 — the invariants (§3)

### Proposition 3.11 and Corollary 3.13 — the Cartan matrix up to Z-congruence

> `⟨[S_x], [S_y]⟩_X = μ_X(x, y)`, so the matrix of the Euler form on the basis of
> simples is `1_X^{-1}`. Hence if `X ~ Y` then `1_X` and `1_Y` are **congruent
> over Z**: there is `P ∈ GL_n(Z)` with `P^t 1_Y^{-1} P = 1_X^{-1}`.

The proof is one line and uses nothing about posets: a triangulated equivalence
induces an isometry of Grothendieck groups for the Euler form. So it holds for
**any** finite-dimensional algebra of finite global dimension, LNAs included,
with `1_X` replaced by the Cartan matrix. (Ladkani states the same fact again,
with an explicit congruence matrix, as Lemma 5.2 of arXiv:math/0702060:
`P = [[0, I],[-B^{-1}B^t, -B^{-1}C]]` congruates `[[A,0],[C,B]]` to `[[B,0],[C^t,A]]`.)

### Lemma 3.14 and Corollary 3.15 — hence the Coxeter matrix up to Z-conjugacy

> If `M_1, M_2 ∈ GL_n(R)` are congruent then `M_1 M_1^{-t}` and `M_2 M_2^{-t}` are
> **conjugate** in `GL_n(R)`. So `X ~ Y` implies `1_X 1_X^{-t}` and `1_Y 1_Y^{-t}`
> are similar over `Z`, in particular over `Q` and modulo every prime `p`.

`1_X 1_X^{-t}` is, up to sign, the Coxeter matrix. **This is the invariant we
want**: the Coxeter *polynomial* is only the characteristic polynomial of `Φ`,
whereas the `Z`-conjugacy class of `Φ` is the isomorphism class of `Z^n` as a
`Z[Φ]`-module, which is strictly finer. Ladkani's own witness is Example 4.20:
two posets whose Euler forms are equivalent **over Q** but not over `Z`, separated
by Corollary 3.15 at `p = 11`.

**Checked here, 2026-09-19.** Take the profile

    for each irreducible factor g of the Coxeter polynomial of Λ,
        the Smith normal form over Z of g(Φ_Λ)

(basis-independent, since `SNF(P^{-1} g(Φ) P) = SNF(g(Φ))` for `P ∈ GL_n(Z)`).

- *Sound.* Constant on every one of the 73 mutation classes of
  `A_10_mutation_classes.csv` (two members sampled per class): 0 disagreements.
- *Sharp where it must be.* It separates **every** cospectral quipu group at
  order 9 (1 of 1) and order 10 (2 of 2), 1 of 4 at order 11, 8 of 13 at order 12.
  These are the groups F-010 identifies as the exact failure mode of the Coxeter
  polynomial, and they are **provably distinct classes**, so every split is a
  genuine separation. The smallest, at `n = 9`, is `3060000` against `3004000`,
  where the factor `x^2+x+1` gives SNF `(1,1,1,1,1,1,1,0,0)` on one side and
  `(1,1,1,1,1,2,2,0,0)` on the other.
- *Conservative elsewhere.* Among the 73 classes at `n = 10` there are 13 groups
  of classes sharing a Coxeter polynomial; the refinement fully splits exactly
  **one** of them — and that one is `{P^(1,4)_(1,0,2), P^(3,3)_(1,0,1)}`, a
  cospectral quipu pair. It left the other 12 intact, which is the right answer
  if (as H-003 suspects) those groups are search orbits that ought to merge.

So this is the theorem-free separation F-010 asked for and that the
Avella-Alaminos–Geiss route could not supply (R-008). Enriching the profile with
`h(Φ)` for every *product* `h` of irreducible factors adds nothing — checked at
orders 11 and 12, identical results.

### Corollaries 3.21, 3.22 — Betti numbers and the Euler characteristic

> If `X ~ Y` then `β_i(X) = β_i(Y)` for all `i`, and `χ(X) = χ(Y)`, where
> `χ(X) = Σ(-1)^i β_i(X)` is also the **sum of all entries of `1_X^{-1}`**.

**Do not carry this over to LNAs.** The invariance comes from
`HH^i(kX) = H^i(X)` (Theorem 3.19, Cibils / Gerstenhaber–Schack), which is a fact
about *incidence algebras*. `HH^i` is a derived invariant of any algebra, but
"sum of the entries of `C^{-1}`" equals `⟨[k_X],[k_X]⟩` for the *constant sheaf*,
an object a derived equivalence need not preserve. An LNA is not an incidence
algebra, so neither the identification nor the entry-sum shortcut applies.

### Propositions 3.26, 3.29 — operations that preserve the class

> `X ~ Y ⟹ X^op ~ Y^op`. And `X_1 ~ X_2`, `Y_1 ~ Y_2 ⟹ X_1 × Y_1 ~ X_2 × Y_2`,
> because `k(X × Y) = kX ⊗_k kY` (Lemma 3.28).

The opposite statement is the poset form of our relation dual. The product
statement is what makes the rectangle `kA_m ⊗ kA_n` behave.

## Part 2 — the construction (§4)

### Proposition 4.5, Corollary 4.6, Proposition 4.7 — the algebra `A_Y`

Let `X` be a poset, `Y ⊆ X` **closed** (a down-set), `U = X \ Y` its complement.
Put `~P_y = i_*i^{-1}P_y` and `~I_u = j_!j^{-1}I_u` (the projectives truncated to
`Y`, the injectives truncated to `U`). Then:

> **Proposition 4.5.** `E_Y = { ~P_y }_{y∈Y} ∪ { ~I_u[1] }_{u∈U}` is a strongly
> exceptional collection generating `D^b(X)` (order it `U` first, then `Y`).
>
> **Corollary 4.6.** `D^b(X) ≃ D^b(A_Y)` where `A_Y = End_{D^b(X)}(T_Y)`,
> `T_Y = (⊕_y ~P_y) ⊕ (⊕_u ~I_u)[1]`.
>
> **Proposition 4.7.** `A_Y` has `k`-basis `{e_{yy'} : y ≤ y'} ∪ {e_{u'u} : u' ≤ u}
> ∪ {e_{uy} : y < u}` with
> `e_{yy'}e_{y'y''} = e_{yy''}`, `e_{u''u'}e_{u'u} = e_{u''u}`,
> `e_{uy}e_{yy'} = e_{uy'}` if `y' < u` and **0 otherwise**,
> `e_{u'u}e_{uy} = e_{u'y}` if `y < u'` and **0 otherwise**.

The two "0 otherwise" clauses are the whole point: `A_Y` is an algebra **with
zero relations**, and it is exactly an incidence algebra again iff

> **Lemma 4.9 (⋆).** whenever `y ≤ y' ∈ Y`, `u' ≤ u ∈ U` and `y < u`, also
> `y' < u'`.

**Example 4.8** is the smallest instance and shows an LNA coming out: `X` the
poset `1 < 3`, `2 < 3`; `Y = {1}`. Then `A_Y` is `kA_3` on the line `2 → 3 → 1`
modulo the zero relation `(2→3)(3→1) = 0` — our LNA `2` on three vertices.

**How to get an LNA merge out of it.** Run Corollary 4.6 on **one** poset `X`
with **two** different closed subsets `Y, Z`; then `A_Y ≃ A_Z`, and if both are
lines with zero relations that is a merge of two LNAs. That is precisely
Dong–Lin–Ruan's method: their Proposition 4.1 takes `Y = {1..v}` and
`Z = {1..v+u}` in one poset to get `N(2u+v, u+v+1)` on one side and a one-branch
extension of a rectangle on the other, and their Proposition 4.5 takes
`Y = {1..r-1}` in the poset behind `1A(r-1)` to get `N(2r-1, r) ≃ N(2r-1, r+1)`.
See `2203.15735-one-branch-extensions-rectangles.md`.

### Theorem 4.14, Corollaries 4.15, 4.18 — lexicographic sums and BGP

> **Theorem 4.14.** If `S` is a bipartite poset and `X = {X_s}` a collection of
> posets, then `⊕_S X ~ ⊕_{S^op} X`.
>
> **Corollary 4.15.** `X ⊕ Y ~ Y ⊕ X` for any two posets (ordinal sum).
>
> **Corollary 4.18.** `S ~ S^op` for bipartite `S` — i.e. BGP reflection.

Example 4.16 is an APR tilt, Example 4.17 is a source-to-sink BGP reflection; the
theorem is the common generalisation.

### Example 4.20 — and why three summands is different

> With `X` a 3-element antichain, `Y` the poset `• → •`, `•` (three elements,
> one relation) and `Z = X ⊕ Y`, the posets `X ⊕ Y ⊕ Z` and `Y ⊕ X ⊕ Z` are
> **not** derived equivalent — their Euler forms are equivalent over `Q` but not
> over `Z`, shown by Corollary 3.15 at `p = 11`.

Proposition 4.19 shows that permuting `n` summands reduces to the three-summand
case, so this single example kills the general commutativity. Two things to take
from it: the ordinal-sum flip is genuinely a *two*-summand statement, and mod-`p`
conjugacy of the Coxeter matrix is not a formality — it is the only thing that
separated this pair.

## Caveats and limitations

- **Everything in §4 starts from a poset.** `A_Y` is built *from* an incidence
  algebra; you cannot feed an LNA in. So the construction can produce LNA merges
  only when two closed subsets of one poset both give lines — you have to guess
  the poset. It is not a move on digit strings.
- **Corollary 4.6 has no converse.** Not knowing which `A_Y` are lines is the
  whole difficulty; the paper gives no classification of that.
- **`Y` must be closed**, i.e. a down-set. For a chain that forces `Y` to be an
  initial segment, and then `(⋆)` holds and `A_Y` is a chain again — so applying
  §4 to a *line* gives nothing. The interesting input posets are branched.
- **§3's poset-specific invariants do not transfer** (see above on Betti numbers).
  The two that do transfer verbatim — Cor. 3.13 and Cor. 3.15 — are the two we use.
- **Testing Z-congruence is undecidable in practice.** Ladkani says so outright.
  Corollary 3.15 is the usable shadow, and even it is only a necessary condition:
  two algebras with conjugate Coxeter matrices need not be derived equivalent.
- The paper works over a field and with finite-dimensional vector spaces; the
  field-independent ("universal") versions of Theorem 4.14 and Corollary 4.18 are
  in arXiv:0705.0946, which we do not otherwise need.

## What it does not give us

- No decision procedure — the paper says explicitly that no algorithm is known to
  decide derived equivalence of two posets.
- No contact with quipus, piecewise heredity, or relation lengths.
- The `A_Y` construction never produces a *cyclic* quiver and never produces
  relations longer than the zero relations above, so it cannot reach LNAs with
  long relations except through the accidents Dong–Lin–Ruan exploit.

## How it lands on our problem

Implement the §3 refinement as an invariant beside `coxeterKey`: for an LNA `Λ`,
compute `Φ = -C^{-t}C`, factor its characteristic polynomial over `Z`, and record
the Smith normal form of `g(Φ)` for each irreducible factor `g`. It costs one
`n × n` SNF per factor, is sound on every known class at `n = 10`, and is the
first thing in this repo that separates the cospectral quipu pairs of F-010
without appealing to `thm:QuipuToAn`. Where it fails to split a group, that is
weak evidence for a merge.
