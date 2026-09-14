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

## Writing a summary

Cover, in this order: what the paper is for; the results we actually use, stated
precisely enough to implement; lemmas worth knowing; **caveats and limitations**,
especially anything that would mislead someone who only read the main theorem;
and what it does *not* give us. Record where a result is used in the code.

A summary that only restates the abstract is not worth having.

## Candidates not yet summarised

Worth a file when someone gets to them:

- **Happel–Seidel**, *Piecewise hereditary Nakayama algebras* — solves the case
  where the ideal is a power of the radical, over an algebraically closed field.
  Table 1 is a complete answer in that setting.
- **Happel**, the classification of hereditary abelian categories — the result
  behind the trichotomy we lean on (module category of a hereditary algebra, or
  derived equivalent to a canonical algebra).
- **Chen–Ringel**, on hereditary triangulated categories — cited for the
  characterisation used in the vertex-deletion corollary.
- **Ladkani**, on derived equivalences of Nakayama algebras — the line-to-rectangle
  result recreated in arXiv:2112.08129.
- **Avella-Alaminos–Geiss**, derived invariants of gentle algebras. LNAs are
  gentle, so their invariant applies directly and is a candidate for separating
  classes the quipu theorem does not reach (NOTES idea 14).
- **Woo–Neumaier**, where open quipus were introduced, as graphs of small spectral
  radius. Possibly relevant to F-010, since cospectrality is what makes the
  Coxeter polynomial fail.

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
