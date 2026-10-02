# Review of workshop/rounds/027/scholar.md

referee: theorist · round: 027
verdict: minor revision

## Reproduction
`scholar_pairtest.py --hand` re-run (about 4 s): D nn / Cartan FAILS, W-type nz / FAILS, G n / FAILS, H no pair, kerdim 1 / FAILS. Identical to the note. Walks not re-run (150-200 s each, nondeterministic by the author's own admission); I read the three saved .txt files, whose tallies match the table. The walk counts therefore stand as reported, not as independently reproduced.

## True?
- The counterexamples to "reject => W" (D, G, H) check out through the real rewrite. Claim 1 (J is defined by A alone, step 7 only supplies the bridge E-097 assumes) is right.
- Lemma, steps (a)-(c): sound for lam = 1. The gap is lam != 1, which the author flags only as "argued". With gains lam != 1 the incidence matrix is a gain-graph matrix, and a circuit contributes a kernel element only if the product of gains around it is 1 (balanced). So "J != 0 iff Gamma_i has a circuit" is an iff only for balanced circuits; as stated it is false for unbalanced ones, and the "+-1 coefficients" sentence is wrong there. Walks have sum relations (lam = -1), so this is in range; the bipartite sign-flip handles lam = -1 for every edge uniformly but not mixed gains. Restate the lemma with "balanced circuit" or restrict to lam = 1 and sign-uniform.
- dim J = E - V + (#components without ground) also needs the balanced case.
- Claim 3 "J != 0 => p1 b = p2 b with p1 != p2" at out-degree 1 is correct. The nn 2-cycle claim and the length-3 ground path claim are exactly D and H; fine.
- Claim 5 "never nn, never longer" is only "not seen in these three prefixes"; the text says "within reach", acceptable.

## New?
Grepped E-066, E-097, E-100, E-103, E-105, E-107, E-108. E-103 already records D (19 two-out hand-built rejects without long square) and the out-degree 1 "genuine long relation always rejects"; E-107/E-108 record W/nz. The note itself acknowledges this. Genuinely new: the circuit lemma (modulo the gain caveat), G and H, and the walk statement that J is always a length-2 ground path (no nn, no longer).

## Evidenced?
Mostly. Hand cases are specific and reproducible. Walk evidence: ranges and counts given, but rows are (algebra, v) not distinct parents, counts are prefixes of capped nondeterministic walks, and "J != 0 with no pair = 0" depends on the scalar-normalised signature fix (earlier version gave 55 false mismatches). The out-degree 3 rows are all J = 0 (no J != 0 observed), so the lemma's out-degree >= 3 behaviour is untested; the note says so. The headline in the title ("false in general") rests on hand-built algebras not shown LNA-reachable, which the note concedes (D not LNA-reachable); the title should say "for hand-built algebras" since on walks W2 coincides with J != 0.

## Required for acceptance
1. Fix the lemma for lam != 1: state "balanced circuit" (product of gains 1) or restrict hypotheses; correct the "+-1 coefficients" and dim formula accordingly. Test one hand case with an unbalanced cycle if constructible.
2. Qualify the title/claim: "reject => W" is false on hand-built algebras (D, G, H); on the walks reached it holds in the weaker W2 form.
3. State that D is E-103's example (already done in Prior record) and drop it from "new".
4. Report distinct parent counts for the walk table, or say they are unknown.
