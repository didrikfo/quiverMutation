# Theorist notebook (rewritten round 034)

## What I believe now
- J_i = Hom(S_v, e_iA) = socle copies of S_v in e_iA (E-122). Gate (isMutable) = "no single nonzero path in J" (E-114). Round 034: L1 gives dim J_i <= d_i - 1 for gate-admitted v (d_i = dim e_iAe_v), so J != 0 needs d_i >= 2 and d_i = 2 gives dim J_i = 1. Layered m = 3 hand example has (d,J) = (3,2): dim J = 1 is not forced in general.
- "d_i <= 2 whenever J_i != 0 on walks" is the unproved part (E-116 bound too weak). Walk prefixes n = 7, 8: only (2,1) and (2,0) pairs with d >= 2 appear.
- L2: Hom(N,N[-1]) = {h: P_v -> D' with hg = 0, gh = 0} lies in paths from out-neighbours back to v, so is 0 on acyclic algebras; Hom(D,N[-1]) = 0; Hom(T,T[-1]) = sum J_i. Silting (AI 2.31) cited, not proved; AI text not compared.
- Earlier (031): circuits/nn shapes are not excluded by the kernel structure alone; layered families have no LNA key but E-123 says that is uninformative (0 of 2 704).

## What I tried
- 034: `rounds/034/theorist_dimji.py` (hand: E-078 and layered m=3; walk: (d,J) histogram at all gate-admitted (v,i)).
- 031: `theorist_t3.py`, `theorist_layers.py`. 029: `theorist_circuit.py`, `theorist_keys.py`. 026: W exact match on 17 802 rows.

## Next
- Prove or find a counterexample to d_i <= 2 at J_i != 0 on walks: track how e_iAe_v gains its second path under one mutation (the new arrow composite), try an induction on mutation steps.
- Cyclic quiver case for L2 (does a cycle through v exist in any derived class here? walks skip cyclic algebras).
- Invariant (Coxeter polynomial / Euler form) separating nn / long circuits from LNA keys is still open.
- Blind spots: a miss in a thin family is not "cannot occur"; L1 is a restatement of the gate, do not oversell it as the explanation of d = 2.
