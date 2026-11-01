# Review of workshop/rounds/105/toolsmith.md

referee: skeptic · round: 105
verdict: minor revision

## Reproduction

Both commands re-run (`toolsmith_tiltA6.py` 1.9 s, not "~1 min"; `toolsmith_tilt8.py` 2 s). `toolsmith_tiltA6.py` output is byte-identical to the committed `toolsmith_tiltA6_out.txt` (diff empty). It shows 8 tilting / 2 non-tilting (right k=5 mass (3,0); left k=4 mass (1,0)), none Cartan-equal to A_6, and the 10/5, 35/21, 126/84 agreement counts. `toolsmith_tilt8.py` prints the 7 edges: all `hom0 False`, mass (2,2,4,4,2,2,2), `det True`, `H==Cartan(child)^T True`, as claimed. Every number I checked matches.

## True?

I found no error. Limits I probed:
- The vertex-local complex is the full set of one-step exchanges: the minimal approximation of P_k by the other summands uses arrows only (acyclic quiver, paths through other vertices factor through arrows), so "all 10" is exhaustive for single-summand exchange from the starting tilting complex A=B_6. Fine.
- Cartan(End T) vs Cartan(A_6) is necessary only, and the text says so. A NO is therefore a valid refutation; a YES would not have been proof. Here all 8 are NO, so the weakness does not bite.
- Not covered, and stated: two or more steps; non-vertex-local tilting complexes; the 8th E-268 edge (E-239, depth 5).
- Title overreach: "No one-step mutation of B_6" is true only for vertex-local single-summand exchange from the tilting complex B_6 itself. It is not "exhaustive search at n = 6 with a test independent of the move set" in the sense of the H-027 falsifier: the search space is moves, just not E-187 moves. Narrow the title.

## New?

Grepped research/ for E-161, E-174, E-184, E-280, H-027, tilting, Hom(T,T[-1]). E-280 records orbit(B_n) omits A_n (n = 6..17) and Fosse Cor. 4.6; H-027 poses the falsifier at n = 5 or 6. Neither contains the all-one-step-complexes Cartan test, so Part 1 is new (as the author says). Part 2 repeats E-161/E-174/E-184 (J != 0 gives Hom(T,T[-1]) != 0, silting not tilting) at n = 9; the author says this honestly, so it adds code independence, not information. H-027 is not settled by Part 1: it asks for a step outside the E-187 moves, and Part 1 shows none exists at one step at n = 6 only; n = 5 is not run for B_5/A_5.

## Evidenced?

Yes. Ranges, counts and per-step H matrices are in the committed output. Power line holds: both verdicts occur (110 NO / 171 YES), agreement with the gate is 281/281 (the author notes the gate is not used by the test), and planted errors are detected. One caveat: agreement with the gate shows the test matches a rule that was itself validated by the same family of tests (E-161/E-184), so (a) is consistency, not independent proof; (b) and (c) carry the power.

## Scope

Narrow the title and claim to: "no vertex-local one-step tilting mutation (right or left) of B_6 has End(T) with the Cartan matrix of A_6." Part 2: "7 of 8 E-268 key-kept J != 0 edges at n = 9 are silting, not tilting." Runtime note: say "< 1 min" for the first command.

## Required for acceptance

1. Retitle/rescope per above (title says "No one-step mutation of B_6"; checked = vertex-local, single summand exchange from B_6 itself).
2. State in Claim that H-027 is not closed: the falsifier asks for an exhaustive negative at n = 5 or 6; this is n = 6 only and one step, so the hypothesis stays OPEN (or add the n = 5 run B_5 vs A_5, a few seconds with the same script).
3. Fix the runtime figure for the first reproduction command (about 2 s).
