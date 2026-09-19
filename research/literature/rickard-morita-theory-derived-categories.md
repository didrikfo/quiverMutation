# Morita theory for derived categories

**J. Rickard**, J. London Math. Soc. (2) **39** (1989) 436–456 · and the sequel
*Derived equivalences as derived functors*, J. London Math. Soc. (2) **43**
(1991) 37–48.

**Written from secondary sources.** Neither paper is on arXiv, and neither was
read here. Everything below is either (a) the modern restatement in
Aihara–Iyama, arXiv:1009.3370 Definition 2.1(b), Example 2.2(a) and Proposition
2.3, which *was* read; (b) an elementary computation done here from those
definitions; or (c) a citation taken at face value from Bobiński–Ciborski,
arXiv:2409.05158, whose reference [21] is the 1991 paper. The **numbering**
"Theorem 6.4" is the one the literature uses for the 1989 result; it was not
verified against the printed paper.

This is the theorem the whole repo stands on: it says what a *proof* that two
LNAs are derived equivalent has to produce. arXiv:2112.08129's procedure and
arXiv:1009.3370's mutation are both machines for producing that object.

## The criterion

A complex `T ∈ K^b(proj A)` is a **tilting complex** when

1. `Hom_{K^b(proj A)}(T, T[i]) = 0` for all `i ≠ 0`; and
2. `thick(add T) = K^b(proj A)` — `T` generates, i.e. the smallest full
   triangulated subcategory closed under summands containing `T` is everything.

> **Theorem** (Rickard 1989, Thm 6.4). For rings `A`, `B` the following are
> equivalent:
> (a) `D^b(Mod A) ≃ D^b(Mod B)` as triangulated categories;
> (b) `K^b(proj A) ≃ K^b(proj B)` as triangulated categories;
> (c) `B ≅ End_{K^b(proj A)}(T)` for some tilting complex `T` over `A`.

In Aihara–Iyama's language: condition 1 is exactly "`add T` is a **tilting**
subcategory", and (c) is their Proposition 2.3 (attributed to Keller): an
algebraic triangulated category with a tilting object `M` is equivalent to
`K^b(proj End(M))`. Note the asymmetry — **silting** only asks `Hom(T, T[i]) = 0`
for `i > 0`, and a silting object's endomorphism ring is *not* derived
equivalent to `A`.

> **Theorem** (Rickard 1991, Thm 3.3). If `A` and `B` are derived equivalent,
> there is a complex `X` of `B`-`A`-bimodules with `X ⊗^L_A −` a derived
> equivalence — a **two-sided tilting complex**.

An equivalence isomorphic to some `X ⊗^L_A −` is called **standard**. Rickard's
own open question — is every derived equivalence standard? — is still open in
general; it is known for triangular algebras, which includes all LNAs.

## The two-term case, computed out

Every mutation we perform produces a **two-term** complex
`T = (P^{-1} --f--> P^0)` concentrated in degrees −1 and 0, with `P^{-1}, P^0`
projective. Condition 1 then collapses. Writing out chain maps and homotopies:

- **`Hom_{K}(T, T[i]) = 0` automatically for every `|i| ≥ 2`.** `T[i]` lives in
  degrees `−1−i, −i`, which for `|i| ≥ 2` is disjoint from `{−1, 0}`, so there is
  no nonzero component to give.
- **`Hom_K(T, T[1]) = coker( Hom(P^0, P^0) --(−∘f)--> Hom(P^{-1}, P^0) )`.** A
  chain map is a single `g : P^{-1} → P^0` (no commutation condition), and the
  null-homotopic ones are exactly `h ∘ f`. So the condition is: **every
  homomorphism `P^{-1} → P^0` factors through `f`.**
- **`Hom_K(T, T[-1]) = { g : P^0 → P^{-1} : g ∘ f = 0 and f ∘ g = 0 }`** — there
  are no nonzero homotopies in this degree, so the chain maps *are* the
  morphisms. The condition is: **no nonzero `g : P^0 → P^{-1}` with
  `gf = 0 = fg`.**

So for a two-term complex, tilting = *these two* vanishing statements + the
generation condition 2. Both are finite linear algebra over `Hom_A(P_a, P_b) =
e_a A e_b`, which for an LNA is 0- or 1-dimensional and is exactly "the path
`a → b` is nonzero in `A`".

The first of the two is the one silting mutation gives for free
(arXiv:1009.3370 Theorem 2.31: any mutation of a silting object is silting). The
second is the one that can fail, and Aihara–Iyama Theorem 2.32(b) is its
approximation-theoretic form. **These are two views of the same condition** —
implement whichever is cheaper.

**Generation (condition 2)** is the one nobody checks and nobody should skip. For
`T = μ⁻_{P_i}(A)`, i.e. `A` with `P_i` replaced by the cone of an approximation
and the other `P_j` untouched, generation holds because the triangle recovers
`P_i` from the rest. It is *not* automatic for an arbitrary two-term complex —
`T = P^{-1} ⊕ P^0` with the wrong summands generates a proper thick subcategory.

## Caveats and limitations

- **The theorem is an existence statement, not an algorithm.** It says a derived
  equivalence is *witnessed* by a tilting complex; it gives no way to find one,
  and no way to prove none exists. It cannot separate two classes. Everything we
  do to prove two LNAs are *not* derived equivalent has to come from invariants,
  never from here.
- **`End(T)` is taken in `K^b(proj A)`, not in the module category.** It is
  homotopy classes of chain maps. Computing it as module endomorphisms of the
  cohomology is wrong and will silently give a different algebra.
- **Condition 1 is `i ≠ 0`, not `i > 0`.** The single most common misreading, and
  it is the exact place where a mutation that "worked" turns out not to be a
  derived equivalence — R-005's 38 admitted rules are that mistake wearing
  different clothes.
- Derived equivalence is **not** Morita equivalence: `B` is determined only up to
  Morita equivalence by `T`, and different tilting complexes over the same `A`
  give genuinely non-isomorphic `B`. Our classes are classes of algebras up to
  derived equivalence, so this is harmless — but a "new" quiver produced by
  mutation may be the same algebra presented differently, and only an
  isomorphism test says so.
- Two derived equivalent algebras need **not** have equivalent module
  categories, equal representation type in the naive sense, or the same quiver.
  They do share: the Grothendieck group with the Euler form (hence the Coxeter
  polynomial), finiteness of global dimension, the centre, and `HH*`.
- The 1991 paper's two-sided statement is what makes "derived equivalence" behave
  like a Morita theory (composition, inverses, a derived Picard group). We do not
  use it, but it is why the relation is an equivalence relation on algebras and
  why our orbits are well defined.

## What it does not give us

No construction, no invariant, no merge. It is the definition of the finish line.

## For LNAs specifically

The usable content is the two-term computation above. `procedure.isMutable`
implements a criterion that arXiv:2112.08129 itself says is **necessary only**;
the two displayed vanishing conditions are **necessary and sufficient**, and both
are computable from the path spaces `e_a A e_b` that `pathAlgebra` already
carries. Checking them at a mutation step would replace F-016's after-the-fact
Coxeter check — which can only detect a *wrong* answer, never certify a right
one — with a certificate.

Worth noting what this means for the repo's standard of proof: a merge of two
classes is established when a tilting complex is exhibited, and a chain of
certified mutations is exactly that, composed. A chain of *uncertified*
mutations, even one that preserves the Coxeter polynomial at every step, is not.
