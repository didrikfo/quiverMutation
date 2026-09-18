# Literature

One summary per paper: the results we use, the lemmas worth remembering, and the
caveats that bite. **Not** the papers themselves — those live on arXiv and are
linked from each summary.

The point of summarising rather than storing is that a summary can hold a single
useful result from a paper that is otherwise irrelevant to us. That makes a
literature sweep worth doing: anything with one usable lemma earns a file here,
however far from our subject the rest of it is.

| paper | what we take from it |
|---|---|
| [2112.08129](2112.08129-combinatorial-tilting-mutation.md) — A combinatorial procedure for tilting mutation | The mutation rule the whole codebase implements; the admissibility condition |
| [2305.06642](2305.06642-quipu-quivers-nakayama.md) — Quipu quivers and Nakayama algebras with almost separate relations | The quipu ↔ LNA correspondence; the class-preserving operations; the n ≤ 8 table |
| [2310.08346](2310.08346-non-piecewise-hereditary-nakayama.md) — Non-piecewise hereditary Nakayama algebras | Criteria certifying an LNA is in no quipu class; the vertex-deletion corollary |

**Added by the sweep of 2026-09-18/19** (E-040), whose question was: what does the
literature give us for **merging** classes of LNAs?

| paper | what we take from it |
|---|---|
| [2112.15587](2112.15587-nakayama-fuchsian-singularities.md) — Nakayama algebras and Fuchsian singularities | The complete derived classification of the radical-power LNAs, the Happel–Seidel symmetry that merges them, and a periodicity obstruction we do not have |
| [happel-seidel](happel-seidel-piecewise-hereditary-nakayama.md) — Piecewise hereditary Nakayama algebras | Table 1: which `kA_n/rad^r` are piecewise hereditary, and of star or of sheaf type. Checked against our own naming, F-044. *Secondary sources* |
| [2302.02880](2302.02880-ueda-derived-equivalences-nakayama.md) — Ueda, On derived equivalences of Nakayama algebras | `N(n,l+1) ≃ N(n,l)` on an infinite parameter family — the only merge between different relation lengths any paper hands us. Already in our move table at `n ≤ 16`, F-043 |
| [1009.3370](1009.3370-silting-mutation.md) — Aihara–Iyama, Silting mutation in triangulated categories | Theorem 2.32: an **iff** for "this mutation is a derived equivalence", where our gate is necessary only; and transitivity for piecewise hereditary algebras |
| [1504.02617](1504.02617-quivers-for-silting-mutation.md) — Oppermann, Quivers for silting mutation | The same seven-step rewrite on a dg quiver, with no admissibility hypothesis, and the cyclic case F-002 leaves out |
| [rickard](rickard-morita-theory-derived-categories.md) — Rickard, Morita theory for derived categories | What a *proof* of derived equivalence has to produce. *Secondary sources* |
| [2608.08222](2608.08222-iterated-tilted-type-a.md) — Enumerating iterated tilted algebras in type A | Assem–Happel restated: a local four-condition certificate for membership in the `kA_n` class |
| [2312.14699](2312.14699-hochschild-monomial-bardzell.md) — The Hochschild cohomology ring of monomial algebras | Bardzell's resolution in full — and with it, that `HH*(A) = k` for **every** LNA, which closes idea 22 |
| [bruestle](bruestle-derived-tame-tree-algebras.md) — Brüstle, Derived-tame tree algebras | **Theorem 1.2: the derived class of any LNA with non-negative Euler form, from three numbers.** The only two-sided criterion we have, F-045. *Journal PDF* |
| [2203.15735](2203.15735-one-branch-extensions-rectangles.md) — Dong–Lin–Ruan, One-branch extensions of "rectangles" | Prop. 4.5: `N(2r-1,r) ≃ N(2r-1,r+1)`, a merge no move of ours makes from `n = 11` (F-046); and closed-form Coxeter polynomials |
| [0911.5137](0911.5137-lines-rectangles-triangles.md) — Ladkani, Lines, rectangles and triangles | Cor. 1.2: `A(mn, m+1) ≃ kA_m ⊗ kA_n`, so `A(mn,m+1) ≃ A(mn,n+1)` — the result arXiv:2112.08129 recreates |
| [math/0610685](math-0610685-sheaves-over-finite-posets.md) — Ladkani, Sheaves over finite posets | Cor. 3.13: the Cartan matrix up to `Z`-congruence, a derived invariant **finer than the Coxeter polynomial**; and the `A_Y` construction |
| [1001.4765](1001.4765-perverse-equivalences-bb-tilting-mutations.md) — Ladkani, Perverse equivalences, BB-tilting, mutations | Prop. 2.3(c): an **iff** for our own mutation gate; Prop. 3.6: the mutated Cartan matrix in closed form, verified against the repo |
| [math/0611201](math-0611201-coxeter-periodicity-euler-form.md) — Ladkani, Periodicity of Coxeter transformations | Thm. 3.4: periodic `Φ` plus indefinite Euler form certifies non-piecewise-heredity, from the Cartan matrix alone |
| [2509.12983](2509.12983-chz-criterion-derived-equivalences.md) — Pavon, Detecting derived equivalences with the CHZ criterion | Cor. 3.6: an iff for an HRS-tilt to be a derived equivalence, stated for `kA_n/I`; a set-mutation our engine lacks |
| [0805.1018](0805.1018-spectral-analysis-and-singularities.md) — Lenzing–de la Peña, Spectral analysis of finite dimensional algebras | Where the Coxeter polynomial is and is not complete — the two-sided bound on F-010 |
| [1310.1557](1310.1557-algebras-of-cyclotomic-type.md) — de la Peña, Algebras of cyclotomic type | The periodicity obstruction, with its exclusions. **Read the caveat before implementing** |
| [1606.08279](1606.08279-hereditary-triangulated-categories.md) — Chen–Ringel, Hereditary triangulated categories | Cor. 5.4(3): the only *positive* piecewise-heredity certificate in the sweep |
| [1305.5213](1305.5213-strong-global-dimension.md) — Alvares–Le Meur–Marcos | Happel–Zacharia via `s.gl.dim`: finite iff piecewise hereditary, and a height function on the mutation graph. Its **value** is not a derived invariant |
| [1910.01494](1910.01494-derived-tame-nakayama.md) — Bekkert–Giraldo–Vélez-Marulanda, Derived tame Nakayama algebras | Why the skewed-gentle route does **not** reopen R-008: every LNA has a simple projective. The pointer to Brüstle |
| [2509.02375](2509.02375-coxeter-coefficients-trees.md) — Harel–Ladkani, Coefficients of Coxeter polynomials of trees | Thm. 1.1: the second coefficient counts a quipu's cords — the invariant H-017 asks for |

## Writing a summary

Cover, in this order: what the paper is for; the results we actually use, stated
precisely enough to implement; lemmas worth knowing; **caveats and limitations**,
especially anything that would mislead someone who only read the main theorem;
and what it does *not* give us. Record where a result is used in the code.

A summary that only restates the abstract is not worth having.

## Candidates not yet summarised

Worth a file when someone gets to them:

- ~~**Happel–Seidel**, *Piecewise hereditary Nakayama algebras*.~~ **Summarised**
  (from secondary sources — the paper is journal-only and was not reachable), and
  its table checked against our own naming, 12 weight types and 11 tree types, all
  agreeing: F-044.
- **Happel**, the classification of hereditary abelian categories — the result
  behind the trichotomy we lean on (module category of a hereditary algebra, or
  derived equivalent to a canonical algebra).
- ~~**Chen–Ringel**, on hereditary triangulated categories.~~ **Summarised** as
  [1606.08279](1606.08279-hereditary-triangulated-categories.md).
- ~~**Ladkani**, *On derived equivalences of lines, rectangles and triangles*.~~
  **Summarised** as [0911.5137](0911.5137-lines-rectangles-triangles.md), along
  with five more of his papers; the other 22 were screened and rejected, with the
  reason recorded in E-040 so they are not re-screened.
- ~~**Avella-Alaminos–Geiss**, derived invariants of gentle algebras.~~
  **Dropped — R-008.** LNAs are *not* gentle: a gentle algebra's ideal is
  generated by paths of length two, and an LNA with a longer relation is a string
  algebra but not a gentle one. Worse, an LNA that *is* gentle is derived
  equivalent to the path algebra of A_n, so the gentle LNAs are a single class and
  the invariant has nothing to separate (F-019). Read it only if someone finds an
  extension of the invariant to string algebras. **Searched for in E-040 and there
  is none** — the invariant stops at gentle and skew-gentle — and the skewed-gentle
  route does not reach LNAs either, because every LNA has a simple projective
  module and arXiv:1910.01494's theorem excludes those. Consider this line closed
  rather than merely deprioritised.
- ~~**Hochschild cohomology of monomial algebras**, via Bardzell's minimal
  resolution (idea 22).~~ **Dropped.** Summarised as
  [2312.14699](2312.14699-hochschild-monomial-bardzell.md): on a linearly ordered
  quiver there is at most one path between two vertices, every `n`-ambiguity with
  `n ≥ 1` lies in `I`, and the complex collapses. `HH*(A) = k` for **every** LNA,
  with no dependence on `I` at all, so it separates nothing. Read the file before
  reaching for it again.
- **Woo–Neumaier**, where open quipus were introduced, as graphs of small spectral
  radius. Possibly relevant to F-010, since cospectrality is what makes the
  Coxeter polynomial fail.

**Reaching the papers from here.** *Depends where "here" is.* From the sandboxed
sessions the earlier work ran in, outbound HTTPS to `arxiv.org` and
`export.arxiv.org` was refused by the egress proxy, so `curl` and `WebFetch` both
failed and only `WebSearch` worked — enough to identify a paper and read what
secondary sources say about it, not to read the paper. **From a session on the
author's own machine arXiv is reachable in full** (2026-09-19), which is how the
sweep of E-040 was done: `https://arxiv.org/html/<id>v1` renders the paper, the
reference list is `document.querySelector('.ltx_bibliography').innerText`, and
`curl -sL https://arxiv.org/pdf/<id>` plus `pypdf` reads anything with no HTML
build. `https://export.arxiv.org/api/query` answers search queries, and Semantic
Scholar's graph API gives forward citations — which found two of this sweep's best
papers, and no backward reference list would have.

What stays out of reach either way is anything **journal-only and pre-arXiv**:
Happel–Seidel, Rickard, Assem–Happel. Those have files written from secondary
sources and marked as such at the top; check one against the repo's own
computations before relying on it, as F-044 did.

For H-009 — the suspicion that the move rules are a one-dimensional cellular
automaton — the sweep to do, none of it yet read:

- **Asynchronous cellular automata**, where one cell updates at a time, which is
  what a single mutation is. Look for reachability and confluence results.
- **One-dimensional rewriting systems on finite words**, and the reachability
  problem under a finite set of local rewrites. This is the exact shape of
  `lnaMoves.matchesAt` / `applyAt`.
- **Conserved quantities and additive invariants of local rules**, the CA analogue
  of "a derived invariant every move preserves".
- **Chip-firing and sand-pile models**, the closest-looking relatives: local,
  order-independent under the right conditions, with a developed theory of orbits
  and of what the boundary does.
- **Block transformations and rescaling of CA rules**, for H-008's question of
  when a family of rules parameterised by window width is one statement.
