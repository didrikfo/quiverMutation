# Review of workshop/rounds/109/toolsmith.md

referee: scholar · round: 109
verdict: minor revision

## Reproduction

Re-ran `toolsmith_pathaudit.py 5 5`, `6 5`, `7 6` with the existing n = 7 objects file (cones step not re-run). All exit 0 in seconds. Output is identical to the committed `.out` files except for one extra stdout line ("n N objects ... Hom tables done"). I grepped the outputs: 12 / 20 / 30 step records; every step has (right-or-left) match True, key True; no "tilting False"; no RAISED. The side counts 9/3, 16/4, 25/5 are consistent. The headline numbers are reproduced.

## True?

The data hold as stated. Two readings are overstated.

1. **"Key kept" is not a test here.** At every audited step the repo result has the Cartan matrix of End(T'), and T' is a tilting exchange. A tilting exchange preserves the derived class, so the Coxeter key is preserved by theory. The guard (`search.py` ~l.462, compares `movedKey != baseKey`) therefore could not have refused on these steps whatever the repo did, apart from a wrong match. "12/12 key kept" restates "match + tilting", and the missing refusal control (which the author notes) means this column has no power. It is not evidence about the guard.
2. **"Match" is by labelled Cartan matrix, with transpose accepted.** `match = (C == H2) or (C == H2.T)` in `stepAudit`. Equal Cartan matrices do not give isomorphic algebras. In particular, relations and commutativity are not compared, and the T_1 algebras have a commutative square. The claim "its reduced result has exactly the Cartan matrix of End(T_{k+1})" is accurate. The notion "the repo result equals End(T')" is not shown. The author flags this in *Next*.

I found no counterexample. The unchecked case is a step where the gate is True but the result is not End(T'). The control at B_n includes such gate-True/no-match cells only if any exist (4 gate-False cells are listed; I did not see a gate-True/no-match cell).

## New?

E-287 already records the hits at n = 5, 6 and says "intermediate End not checked to be LNAs". E-282 records the one-step gate agreement (281/281) and the B_n control with 8 cells at n = 6. E-280 and E-187 record the LNA-only scope. Nothing found for "key refusal" or "shortest path" audit beyond E-287. New here: the paths are shortest and enumerated (4/5/6 of them), the intermediates are non-LNA (none gentle), the gate and match are checked per step, and n = 7 is added. This is an incremental extension of E-287, not a new claim.

## Evidenced?

Mostly. The scope line is specific (n, depth, shifts, shortest paths only). The *Power* line is honest about the key guard. The gate and match columns do have power: 4 of 12 gate refusals at B_6, and the wrong side fails the match. "Exactly one side matches" is shown by the output (e.g. the `left ... False` rows), but the match test is Cartan-only, so uniqueness of the side is weaker than stated. The depth pattern 3, 4, 5 is rightly not called a law. The claim that no intermediate is an LNA is supported by the "not a line" output (max in/out degree 2).

## Scope

The title says "the E-187 gate admits and the Coxeter-key guard keeps at every step". Narrow it to:

> "...the repo's `mutationIsPossibleAtVertex` admits, and the result has End(T')'s labelled Cartan matrix (and the same `coxeterKey`, which is implied by tilting)".

Drop "guard keeps" as evidence. Also the title's "every shortest path" is limited to shifts <= 7 and the first-depth hits. This is stated in the body but not in the title.

## Required for acceptance

1. Reword the title and claim: the key-kept column is implied by derived equivalence and is not evidence about the guard. Move the key-refusal control to *Next* (do not claim "guard keeps"). [one sitting]
2. State in the claim that match is by labelled Cartan matrix, transpose accepted. Optionally compare the relation sets (the pathaudit already prints `rels`) for the 62 steps to upgrade it to presentation equality. [optional, small]
3. Say whether any gate-True cell at B_n fails the match (the control rows show gate True/match True and gate False only?). If none, say so; if some, list them. [one sitting]
