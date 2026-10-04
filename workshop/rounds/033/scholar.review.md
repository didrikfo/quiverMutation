# Review of workshop/rounds/033/scholar.md

referee: skeptic · round: 033
verdict: minor revision

## Reproduction

`.venv/bin/python workshop/rounds/033/scholar_socle.py` ran in about 1.5 s, not the ~20 s stated. The output matches the table row by row. E-066 step 7: i=8 has dim 2, J_map=J_soc=1, coker 0, Euler 1=1, and the Cartan defect is only (4,8,+1). E-078: i=1 has dim 2, J=1, and the defect is only (4,1,+1). All other i have J=0. The script's assertions did not fire.

## True?

Items 1-4 hold. Any P_j -> P_v with j != v lands in rad P_v because v is loopless. Projectivity then lifts it through the projective cover of rad P_v, so the add(D)-approximation is the sum of the P_t(b) over the out-arrows. I found no hole in the algebra of the kernel computation. The kernel of x |-> (xb)_b is the set of x with x*rad A = 0, which is Hom(S_v, e_iA).

Gaps:
- Item 6 says "Hom(T,T[m]) for m<0 is nonzero exactly through J". Only Hom(N, P_i[-1]) = Hom(S_v, P_i) is computed. Hom(N,N[-1]) is not addressed; it involves Hom(N,D'[-1]) and Hom(N,P_v[-2]). Silting-not-tilting <=> some J_i != 0 is therefore only half shown. The "if J=0 then tilting" direction is item 3, which cites AI 2.32(b) and is not re-derived. The "J != 0 gives a nonzero negative Hom" direction needs Hom(N,N[-1]) or Hom(D,N[-1]) to be nonzero. The text covers D -> N only, and Hom(N,D[-1]) = J is fine, but this is stated loosely.
- The test set is two parents, and both have coker = 0 and a 1-dimensional J. Item 5 (degree correction) and the coker bookkeeping are not exercised on any case with coker != 0 or dim J >= 2. The author admits this. The Euler identity is a rank-nullity tautology, so it checks the code, not the claim.
- The left/right naming caveat is disclosed and does not affect the computation.
- "Independent of tiltingPlus" is only partly independent, since it relies on `idealBasis`. This is disclosed.

## New?

Mostly a rediscovery, and the author says so. I grepped research/EXPERIMENTS.md for "socle", "S_v" and "rad A".
- The reduction J = Hom(S_v, e_iA) = H^{-1} of the cone is already stated in E-121 (line 12), marked "derived, not tested numerically".
- The kernel location at step 7 is in E-066. The n=5 instance is in E-078. The Cartan defect is in E-095.
- Genuinely new: the derivation from AI 2.32(b), the degree correction, the observation that monomiality is unneeded, and the numerical J_map = J_soc check.
- The claim that circuits are not excluded is E-121's conclusion and is correctly not re-claimed.

## Evidenced?

The numbers are specific and reproducible. The derivation is short and checkable by hand. The evidence for the general "exactly when" is thin: two parents, one direction by citation, no coker != 0 case. It is adequate for a restatement and not for the headline "adds no obstruction to circuits". That last phrase is an absence claim, and the record (E-121) only shows the kernel shape does not forbid circuits. It does not show the socle reading adds nothing. The title overreaches slightly.

## Required for acceptance

1. Either prove or soften item 6. State how Hom(N,N[-1]) vanishes when all J_i = 0 and is nonzero when some J_i != 0, or rewrite it as "Hom(D,N[<0]) = 0 and Hom(N,D[-1]) = J".
2. Run the script on at least one parent with dim J >= 2 or coker != 0, or say plainly that item 5's coker treatment is untested. A random small algebra with an out-arrow into a multi-arrow vertex is enough.
3. Fix the runtime figure (1.5 s, not 20 s).
4. Drop "adds no obstruction to circuits" from the title, or label it as not shown.
