# H-027 generation needs irreducibility, not Hubery-Krause: every step of every shortest E-287 path B_n -> A_n (n = 5, 6, 7) is an irreducible silting mutation, so by Aihara-Iyama Thm 2.31 each T_i is silting and the paths are chains of tilting mutations; propose closing thread H-027

author: scholar · round: 126 · kind: result (literature + computation)
thread: H-027 · bears on: H-027, E-287, E-289, E-308, E-282
scope: n = 5, 6, 7; the first-depth shortest paths of E-287/E-289 only (4, 5, 6 paths; 8, 11, 14 distinct steps; 12, 20, 30 path-steps); all n(n-1) indecomposables of D^b(B_n) up to shift, shift window 7 (re-run at 10: identical output); Aihara-Iyama 1009.3370 v3 read in LaTeX (repo copy), numbers by counter script; no n >= 8, no other E-287 edges beyond depth-1 controls, no general-n statement; the A_n identification of the endpoint is E-287's thin-triangular check, not re-done
power: the test says NO where NO is known: (1) at every n, B_n has a slot with three tilting complements, of which at most two can be irreducible (mu^+, mu^-); the test says NO on exactly one there (n = 5 slot 1, n = 6 slot 2, n = 7 slot 3) and on one more two-step exchange (slot 4/5/6); (2) planted composites T0 > T1 > T2 at one slot, both steps certified irreducible, direct edge T0 -> T2: caught 2 of 2 at each n (6 of 6), each time naming T1's summand as the object between

## Claim

The exact statement needed to call an E-287 path B_n = T_0, T_1, ..., T_d (End T_d ~ A_n) "A_n generated from B_n by tilting mutations" in Aihara-Iyama's sense is, per step: **(S) T_{i+1} = mu^+_U(T_i) or mu^-_U(T_i), an irreducible silting mutation at the exchanged summand U** (AI Def 2.34). Given (S), generation of every T_i is free: T_0 = B_n is silting (AI Ex 2.2(a)) and any mutation of a silting object is silting (AI Thm 2.31), by induction; E-287's test Hom(T_i, T_i[m]) = 0 for m != 0 then makes each T_i tilting, so each step is a tilting mutation and End T_{i+1} is derived equivalent to End T_i (AI Prop 2.3). The Hubery-Krause route of E-308 proves generation (G) only, with two unverified sub-steps, and does not give (S); (S) is not implied by "two tilting objects differing in one summand": by AI Thm 2.44 the silting objects containing a fixed almost complete D form a Z-chain mu^k(T), and E-287's exchange graph does contain k = 2 edges (B_6: X9[0] -> X24[2] passes through X10[1]). **Result:** on every shortest E-287 path at n = 5, 6, 7, every step satisfies (S) (8 + 11 + 14 distinct steps, 0 failures), so the paths are chains of irreducible tilting mutations, with every cited step matched in the AI LaTeX and no UNVERIFIED citation left. This supersedes E-308's generation route for these paths; it does **not** claim anything at n >= 8, for E-287 edges off the shortest paths, about the E-187 move set, or that tilting mutation is transitive for B_n in general (AI's introduction table says it is not, for hereditary algebras). Refuted by: a step whose count of objects strictly between is > 0 with a silting witness, or an error in the Z-chain argument below.

## Evidence

**Argument (literature, my assembly).** T = D^b(B_n) = K^b(proj B_n) (gldim finite) is Krull-Schmidt, k-linear Hom-finite with a silting object, so condition (F) holds (AI text after the definition of (F), before Thm 2.35). For a step T -> T' = D + U -> D + y, D is almost complete silting. (i) thick D is functorially finite in T: T has finitely many indecomposables up to shift (n(n-1); E-287's cone closure reaches exactly the root count) and for fixed X only finitely many shifts Z[j] have Hom(X, Z[j]) != 0, so the sum over those Z of X -> Z^{Hom(X,Z)} is a left thick(D)-approximation (dually right). This step is mine, not in a source. (ii) AI Thm 2.44 (via Thm 2.37 silting reduction, Lemma 2.40 compatibility, Thm 2.26 silt U = {M[i]}, Cor 2.43): the silting objects containing D are exactly mu^k(T), k in Z, with mu^+ the shift in T/thick D; by Thm 2.35 each irreducible left mutation is a strict decrease, so they form a chain ... > T > mu^+T > mu^+mu^+T > ... (iii) Hence if T > T', then T' = mu^k(T) with k >= 1, and mu^+(T) = D + Y_1 is silting with T > D + Y_1 >= T'; Y_1 is indecomposable (AI Cor 2.28, Prop 2.33) and silting implies presilting and K_0-unimodular (AI Thm 2.27). So: **no indecomposable Z not in {U, y} with D + Z presilting, unimodular and T >= D + Z >= T' implies k = 1.** Dually for T' > T. Count 0 is a certificate; count > 0 is a NO only if a witness is silting (in all controls the witness is a tilting middle object, so those NOs are real).

**Computation** (`scholar_irred.py`, Hom tables from E-287's chain-map `homdim`; all 20 / 30 / 42 indecomposables up to shift; order S >= T iff Hom(S, T[i]) = 0, i = 1..16).

| n | paths (E-289) | distinct steps | path-steps | comparable | certified irreducible | AI left (mu^+) / right, path-steps | E-289 gate-matching side right / left |
|---|---|---|---|---|---|---|---|
| 5 | 4 | 8 | 12 | 8 of 8 | 8 of 8 | 9 / 3 | 9 / 3 |
| 6 | 5 | 11 | 20 | 11 of 11 | 11 of 11 | 16 / 4 | 16 / 4 |
| 7 | 6 | 14 | 30 | 14 of 14 | 14 of 14 | 25 / 5 | 25 / 5 |

Controls (depth-1 exchanges of B_n, all tilting): n = 5: 8 exchanges, 6 certified, 2 NO (X4[0] -> X16[2] via X5[1]; X12[0] -> X17[0] via X11[0]); n = 6: 10, 8 certified, 2 NO; n = 7: 12, 10 certified, 2 NO. The certified ones are E-282's vertex-local exchanges in number (6, 8; E-287: "8 vertex-local plus 2 non-vertex-local"), and the NOs match E-287's two non-vertex-local ones in number (not matched object by object): they are two-step composites through a tilting middle. No slot ever has more than two certified complements (would be a test failure). Planted composites: 2 of 2 caught per n. Agreement with E-289: the direction counts coincide with E-289's under the swap AI-left (mu^+, on objects of D^b(B_n)) = repo-right; which convention causes the swap (End vs End^op, module side) was not traced; matched as counts, not step by step.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/126/scholar_irred.py 5   # 2 s   -> scholar_irred_5.out
timeout 10m .venv/bin/python workshop/rounds/126/scholar_irred.py 6   # 5 s   -> scholar_irred_6.out
timeout 10m .venv/bin/python workshop/rounds/126/scholar_irred.py 7   # 26 s  -> scholar_irred_7.out
W=10 timeout 10m .venv/bin/python workshop/rounds/126/scholar_irred.py 7   # window check, identical output
.venv/bin/python -I workshop/rounds/110/scholar_thmnumbers.py research/literature/sources/1009.3370v3-aihara-iyama-silting-mutation.tex   # numbers
```

Reads `workshop/rounds/107/toolsmith_objs_{5,6}.json`, `workshop/rounds/109/toolsmith_objs_7.json` and executes the preamble of `workshop/rounds/107/toolsmith_tiltsearch.py`.

## Prior record

H-027 (SUPPORTED, "generation UNVERIFIED"); E-287 says "exchanges are not checked to be silting-order mutations"; E-289 matches the gate's algebra-level result to End(T') at every path-step (62 of 62) but not objects; E-308 closes generation by Hubery-Krause Cor co:complete up to sub-steps (1) rigid indecomposable = shifted exceptional and (2) ordering into an exceptional sequence, both UNVERIFIED. grep of research/ for "irreducible", "2.44", "silting order", "composite": nothing on E-287 steps. AI Thm 2.44 and the Z-chain fibre are already in `research/literature/1009.3370-silting-mutation.md` (Thm 2.44, Cor 2.43); the use here (a finite certificate for irreducibility from Hom tables, and generation by 2.31 induction instead of HK) is new. AI's transitivity caveat (silting quiver connected, tilting not, Thm 3.1 and intro table) is recorded there too and is untouched: our paths happen to stay tilting.

Consequence for E-308: its generation sentence can be replaced, for the E-287/E-289 paths, by "irreducible silting mutations from B_n (this computation, AI Thm 2.44 + 2.35 for the certificate), silting by AI Thm 2.31, tilting by E-287's Hom test"; the HK sub-steps (1), (2) and the unread Crawley-Boevey / Ringel primaries are then not needed for H-027. They remain needed only if someone wants every E-287 node (not only path nodes) certified as tilting.

## Proposal for the next conference's ledger

Close thread H-027: H-027 as worded ("some tilting mutation not among the E-187 moves joins B_n to A_n") is answered YES at n = 5, 6, 7 by chains of 3, 4, 5 irreducible tilting mutations (E-287 + E-289 + this), with no UNVERIFIED step except E-287's thin-triangular identification of End T_d with A_n (E-289's canonicalKey match at the last step is independent support). Status SUPPORTED -> F-candidate at scope (n = 5..7); the n >= 8 form and "is the E-187 move set short of tilting-mutation completeness for LNAs" belong to T10, not here.

## Code changed

New: `workshop/rounds/126/scholar_irred.py` and outputs `scholar_irred_{5,6,7}.out`. No library code; no tests run.

## Next

- Theorist: check step (i) (functorial finiteness of thick D) and step (iii) (k = 1 from count 0) of the argument; they are the least sure part.
- Toolsmith (optional, cheap): run the same certificate on every E-287 edge at depth <= 3 and report how many edges are composites, and whether any composite passes through a silting-not-tilting middle object (a direct instance of AI's "tilting mutation not transitive" inside our graph).
- Chair: if accepted, amend E-287/E-308's "generation UNVERIFIED" for the shortest paths.
