# On the periodicity of Coxeter transformations and the non-negativity of their Euler forms

**arXiv:** [math/0611201](https://arxiv.org/abs/math/0611201) · S. Ladkani ·
Linear Algebra Appl. **428** (2008) 742–753 · *read 2026-09-19 from the arXiv PDF*

One theorem we can use as a **certificate that an LNA is not piecewise
hereditary**, computable from the Cartan matrix alone in `O(n^3)`, and
independent of the criteria of arXiv:2310.08346. On a first test it certifies
three LNAs at `n = 10` that Propositions A9 and A13 and the deletion recursion
all miss, two of which are exactly the pair arXiv:2310.08346 had to reach with
its Auslander–Reiten-translate lemma.

## Conventions

`⟨x,y⟩_C = x^t C y` for the Euler form with matrix `C`, and the Coxeter matrix is
`Φ_C = -C^{-1}C^t`. For an algebra `Λ` of finite global dimension, `C` is the
matrix of the Euler form on `K_0` over the **simples**, i.e. `C_Cartan^{-1}` up to
transpose, and `Φ_Λ` is the image of the Auslander–Reiten translation on `K_0`.
The form is *positive* if `⟨v,v⟩ > 0` for all `v ≠ 0`, *non-negative* if
`⟨v,v⟩ ≥ 0` for all `v`, *indefinite* otherwise. `Φ` is **periodic** if `Φ^m = I`
for some `m ≥ 1`, **weakly periodic** if `Φ^m - I` is nilpotent for some `m`.

Note the first paragraph of §3: derived equivalence makes the forms `Z`-equivalent,
so **both positivity properties of the Euler form and periodicity properties of
the Coxeter transformation are derived invariants**. That is worth having on its
own: "the Euler form is indefinite" is a derived invariant the Coxeter polynomial
cannot see, since congruence does not preserve eigenvalues (see Example 4.1).

## The result we use

### Theorem 3.4 — periodic Coxeter + piecewise hereditary ⟹ non-negative form

> Let `k` be algebraically closed and `Λ` a finite-dimensional **piecewise
> hereditary** `k`-algebra. If `Φ_Λ` is periodic then `⟨·,·⟩_Λ` is non-negative.

Contrapositive, which is how we use it:

> **If `Φ_Λ^m = I` for some `m` and the Euler form of `Λ` is indefinite, then `Λ`
> is not piecewise hereditary** — so it is derived equivalent to no hereditary
> algebra, hence lies in **no quipu class**, and is not of canonical type either.

Both hypotheses are pure linear algebra over `Z`:

- `Φ^m = I`: compute `Φ = -C^{-t}C` from the Cartan matrix and take powers. (An
  `m` exists only if every eigenvalue is a root of unity *and* `Φ` is
  semisimple; the search can be bounded by the lcm of the orders of the
  cyclotomic factors of the Coxeter polynomial.)
- indefinite: `Q = C^{-1} + C^{-t}` is a symmetric **integer** matrix (`C` is
  unitriangular), so its eigenvalues are real and the number of negative ones is
  the number of sign variations in the coefficients of `det(tI + Q)` — exact, no
  floating point.

The proof is short and is worth knowing because it is where Happel's
classification enters: `Λ` piecewise hereditary ⟹ `D^b(Λ) ≃ D^b(H)` with `H`
hereditary with a tilting object ⟹ `H` is `mod H` for a hereditary algebra or
`mod` of a canonical algebra ⟹ apply Proposition 3.1 or Proposition 3.2. So it is
the same trichotomy H-006 runs on, spent on a different question.

**What it does at `n ≤ 10`** (exact scan, 2026-09-19):

| `n` | LNAs with `Φ` periodic | certified by Theorem 3.4 |
|---|---|---|
| 4–7 | 5, 14, 41, 115 | 0 |
| 8 | 262 | 0 |
| 9 | 282 | 0 |
| 10 | 709 | **3** |

(The same scan at `n = 11` produces **many** certificates — the run was cut off
before it could count them, but every one in the retained tail is an LNA the
repo's `piecewiseHereditary.certificate` returns `None` for. Worth re-running
with a bound on the period search; the `Φ^m = I` loop is what makes it slow.)

and the three are

| LNA | order of `Φ` | Coxeter polynomial | repo `piecewiseHereditary.certificate` |
|---|---|---|---|
| `34504030` | 18 | `λ^10+λ^9+λ+1` | **None** |
| `50505000` | 18 | `λ^10+λ^9+λ+1` | **None** |
| `45050400` | 84 | `λ^10+λ^9+λ^6+λ^5+λ^4+λ+1` | **None** |

The first two are `example:A10double` of arXiv:2310.08346, which that paper
certifies with `lemma:taupathimpliesnotpwh` — the lemma we have not implemented
(E-029) and which was named there as "the natural next certificate". Theorem 3.4
gets them with three matrix operations and no `τ`. The third, `45050400`, is a
singleton class with a Coxeter polynomial shared by no other LNA at `n = 10`, and
is (as far as this sweep found) not in the literature.

All three are classes with an **empty hereditary form** in
`A_10_mutation_classes.csv`, so this closes three of the rows that F-011's
trichotomy left unnamed.

### Proposition 3.1 — the hereditary case, for calibration

> For a connected acyclic quiver `Q`: `Φ_Q` is periodic **iff** `⟨·,·⟩_Q` is
> positive **iff** the underlying graph is Dynkin; `Φ_Q` is weakly periodic
> **iff** the form is non-negative **iff** the graph is Dynkin or extended Dynkin.

So on the hereditary side "periodic and indefinite" is impossible, which is why
Theorem 3.4 has teeth exactly against algebras with relations.

### Proposition 3.2 — the canonical case

> If `Λ` is canonical of type `(p, λ)` and `Φ_Λ` is periodic, then `p` is one of
> `(2,3,6)`, `(2,4,4)`, `(3,3,3)`, `(2,2,2,2)`; in all four the form is
> non-negative.

These are the four **tubular** weight types. Useful next to
`piecewiseHereditary.TUBULAR_TYPES`: if `Φ` is periodic and the algebra *is*
piecewise hereditary of canonical type, its weights must be one of those four.

## Lemmas worth knowing

- **Lemma 3.5.** For a poset `X`, `C_X = 1_X^{-1}` — the Euler-form matrix is the
  inverse incidence matrix, i.e. the Möbius function.
- **Lemma 3.6.** `(Φ_X)_{xy} = -Σ_{z ≥ x} μ_X(y,z)`, and when any two vertices are
  joined by at most one directed path the Möbius function is just
  `1 / -1 / 0` on `=` / Hasse edge / otherwise. That is Boldt's formula, and it is
  a cheap way to write down `Φ` for a tree or line without inverting anything.
- **Lemma 3.7 and Corollary 3.8.** `C_{X×Y} = C_X ⊗ C_Y` and
  `Φ_{X×Y} = -Φ_X ⊗ Φ_Y`; so a product of two posets with periodic Coxeter
  matrices has one. This is how to compute anything about the **rectangle**
  `kA_m ⊗ kA_n`, and hence — via Corollary 1.2 of arXiv:0911.5137 — about the
  radical-power line `A(mn, m+1)`. `Φ_{A_m}` has order `2(m+1)/gcd(…)`, so the
  periodicity of these lines is immediate.
- **§2.2, Theorem 2.9 and Corollary 2.13.** For *any* square integer matrix `A`
  with `2` on the diagonal, the product of the reflections it defines, in any
  order `π`, equals `-A_{π,+}^{-1} A_{π,-}^t`. No generalised-Cartan, bipartite or
  symmetry hypothesis is needed. So the Coxeter matrix of an LNA can be written
  as a product of `n` reflections of the symmetrised form `C + C^t`, in the order
  given by any linear extension — an alternative factorisation of `Φ` if one is
  ever wanted.

## Caveats and limitations

- **Theorem 3.4 needs `k` algebraically closed**, because Happel's classification
  does. arXiv:2310.08346 makes a point of avoiding that hypothesis; this
  certificate does not.
- **It is a sufficient condition, not a characterisation.** It says nothing when
  `Φ` is not periodic, and nothing when the form is non-negative. In particular
  it does **not** catch `3033030`, the unique non-piecewise-hereditary LNA at
  `n = 9`: there `Φ` is periodic of order 8 but the form has no negative
  eigenvalue (it is degenerate, with a two-dimensional radical). So it is
  **orthogonal to** Propositions A9 and A13, not stronger — use all three.
- **Periodicity is restrictive.** At `n = 10` only 709 of 4862 LNAs have `Φ^m = I`
  at all; the other 4153 the criterion cannot look at.
- **"Indefinite" must be decided exactly.** The negative eigenvalues found at
  `n = 10` are around `-0.12` to `-0.46`; a floating-point test is asking for a
  wrong answer. Use sign variations of the characteristic polynomial of
  `C^{-1}+C^{-t}`, or an exact `LDL^t`.
- **Example 4.1 is a trap worth naming.** Four derived equivalent posets (all of
  type `D_5`, so with *the same* Coxeter polynomial `x^5+x^4+x+1`) have
  **different characteristic polynomials of the symmetrised form** `C + C^t`
  (`x^5-10x^4+36x^3-56x^2+34x-4` and three others). Congruence does not preserve
  eigenvalues. So `spec(C + C^t)` is **not** a derived invariant, however natural
  it looks — only its `Z`-congruence class (Smith normal form, signature) is.
- Example 4.2 kills another plausible guess: there is a poset with
  `spec(Φ_X) ⊄ S^1 ∪ R`. For posets on `≤ 7` elements it does not happen.
- The paper's main point is negative — it refutes [7, Prop. 1.2]'s claim that
  triangularity suffices in place of piecewise heredity. Examples 4.3 and 4.4 are
  the counterexamples; Example 4.4 is `A_3 × D_4`, where `Φ` is periodic of order
  12 and even `τ^e ≃ [d]` holds, yet `v^t C v = -1`.

## What it does not give us

- No merge of any kind — there is not a single derived equivalence constructed in
  this paper.
- No way to decide piecewise heredity positively.
- Nothing about relation lengths, overlap or quipus.

## How it lands on our problem

Add a third criterion to `piecewiseHereditary.CRITERIA`:

    periodic(Φ) and indefinite(C^{-1} + C^{-t})  ⟹  not piecewise hereditary

It needs only `invariants.cartanMatrix`, it is exact over `Z`, and on its first
run it names three previously unnamed classes at `n = 10`. Because all our counts
of non-piecewise-hereditary LNAs are lower bounds (H-004), any independent
criterion moves them up; this one is cheap and provably not subsumed by the two
we have, since `3033030` shows A9 catches what it misses and `34504030` shows it
catches what A9 misses.
