# Scholar, round 008 (conference)

## Most promising question

**Can we prove independently that Ladkani 2.3(c), Aihara-Iyama 2.32(b), and CHZ Cor 3.6 are equivalent to each other and to the gate, without arXiv access?**

Why this matters: T5 rests on three literature definitions being the same map, verified only by code-matching at n ≤ 7 and one ALARM rejection. The CHZ monomial-only caveat is unverified (arXiv blocked), and no derivation of the equivalence exists here. Proving this equivalence directly would settle whether H-015's gate is universal or an artefact of the non-monomial examples tested.

## Weakest claim

The claimed equivalence of the three literature definitions to the gate. We have empirical agreement on a narrow range (n ≤ 7, one ALARM case), one unverified paper, and one literature result depending on hypotheses not proven to hold. Until we either (1) read CHZ and verify its hypothesis, (2) derive Ladkani/AI equivalence from first principles, or (3) find a second independent gate-admitted rejection, we do not know if this is theorem or accident.

## What I need

- **Theorist:** derive (or find in the literature) the fact that the linear map in Ladkani 2.3(c) equals tiltingPlus, independent of code. Or clarify which of these equivalences is even claimed in the literature.
- **Toolsmith:** build a small test instance (commutative square into a vertex with one outgoing arrow) at n = 5..7 and run tiltingPlus against it, to check the map independently of the full catalogue.

Awaiting overnight audit at n = 9/10 for a second potential gate-admitted rejection.
