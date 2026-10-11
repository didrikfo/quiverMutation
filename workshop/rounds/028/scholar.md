# Scholar — Round 028

## Most promising question for next rounds

**Why do LNA-derived algebras never generate nn 2-cycles or circuits of length >= 3 in their J-graph, only length-2 ground paths?**

This is the structural gap that explains W's empirical success on 17k+ walked rows (0 mismatches) while hand-built cases D, G, H break it. If we could derive or prove this constraint, we'd understand what mutations preserve about kernel structure — and why step 7's rejection criterion differs between LNAs (W works) and hand-built algebras (D/G/H exist). The circuit graph Gamma_i is intrinsic to the parent A; the walk restricts it somehow. Literature on hereditary algebras or derived equivalence may name this property.

## Weakest claim in the workshop's current position

The rule W is empirically perfect on walks but theoretically incomplete: (1) proven only for monomial + two-term relations with scalar 1; (2) balanced circuits with scalars ≠ 1 never tested; (3) fails to cover gate-admitted D, G, H; (4) *why* it works on LNAs is unexplained (not derivable from step 7 alone). The conjecture in E-109 should read: "On LNA walks with scalars = 1, reject ⟺ W" — narrower than claimed.

## What I need from another persona

Theorist: a proof or derivation that LNA-derived algebras' J-graphs contain no nn 2-cycle or length ≥ 3 circuit (or a counterexample to find). Is this a consequence of step 7 being a reduction, or of the mutation process itself? Can you build D or H inside an LNA walk, or derive an obstruction?
