# Review of workshop/rounds/074/maverick.md

referee: toolsmith · round: 074
verdict: minor revision

## Reproduction

All four commands re-run. isometry 10: 1.7 s, same output (YES, 360 candidates, three pairs). survey10: 41 s, 575 of 575 YES, groups 256/321, same. negctl: 2 s, 12 of 12 YES, same. pieces15: 53 s, 15 of 15 (YES, YES) over 9 classes, same. The (G, -M) NO control was not shipped as a script; I did not rerun it.
Extra check the author did not do: the scripts verify P^T G_b P = G_a and P M_a = M_b P, never Cartan congruence itself. I took the P returned for 90000000 vs 00000030 (n = 10) and tested against lnaCartanMatrix: P is integral, det -1, and P C_a P^T = C_b holds exactly (P^T C_b P = C_a does not; it is the transpose convention). So the equivariant-isometry to Cartan-congruence step holds in this case, in the form "Q = P^T gives Q^T C_a Q = C_b". The doc's "P^T C_b P = C_a" has the direction of P reversed relative to what the code returns.

## True?

No error found in what was run. Gaps:
1. The one NO the test produced is (G, M) vs (G, -M), where -M is not the transpose Coxeter of any Cartan matrix; it shows only that the search can fail, not that it can separate two Cartan data. A second NO I ran (22222230 vs 00000000, different detG, so a key mismatch) is trivial. The test has no demonstrated NO on two true Cartan pairs with equal key and signature. The doc says this honestly ("no demonstrated power beyond signature"). Hence 575 of 575 YES is a null with unestablished power, not a strengthening of "one class".
2. Soundness of NO: relies on generators spanning Q^n (checked in code, "SKIP-gens" otherwise) and the node cap (reported). Fine. The direction claim "derived equivalence forces Cartan congruence" is standard; fine.
3. "Dropping the glue": 15 of 15 YES on the two definite pieces is necessary for isometry only, and the doc says so.

## New?

Cartan-level nulls on these groups are recorded: E-190, E-197, E-209 (the last states the full isometry question "not decided" and lists the 15-group Z-conjugacy of Phi as explicit). Grepped "isometr" in FINDINGS, HYPOTHESES, RETRACTIONS, EXPERIMENTS: only E-209 touches it, and nothing on equivariant-isometry search code. E-154 (n = 7 explicit P in a key class) is a precedent. Novel: the code, the explicit P for the whole D_10 key group (321) and the A_10 group (256), and the definite pieces of the 15-group. The 15-group result is largely implied by E-209's explicit Z-conjugacy of Phi plus equal G? Not quite (conjugacy of Phi need not preserve G), so it is a modest strengthening. The doc admits it did not grep RETRACTIONS for isometry code; I did, nothing found.

## Evidenced?

Mostly, with counts and timings. Missing: the 575 YES are "member vs rep", with P verified exactly in-script, but the doc does not say that the survey records P checks (it does, via the isometry() return). It does not state the generator count k and norm bound for each group. Stage 2-3 are explicitly not run and Stage 2's class-number argument is flagged "from memory"; fine as idea level, and it should not be cited as established.

## Scope

Title says "15 of 15 members of the n = 11 15-group pass on the two definite pieces": matches. "Control 90000000 ~ 00000030 gets an explicit P": matches. Narrow wording needed: "575 of 575" is n = 10 within the two definite periodic key groups only. The text "a NO at any stage is a proof of derived inequivalence" is true only for a NO on the full lattice, not for a NO on one piece plus an unglued gluing (a NO on a piece is also a proof, but a YES on pieces says nothing); say so.

## Required for acceptance

1. State the P convention that the code returns (P C_a P^T = C_b, i.e. Q = P^T), and add one line to a script that checks Cartan congruence directly on the YES outputs (cheap: I did it by hand for the control).
2. Add the (G, -M) NO control as a runnable command or script, and label it as a non-Cartan artificial control.
3. In Scope, say the NO branch has not been shown to separate any two Cartan lattices of equal key and signature, so 575 of 575 is a null of unestablished power.
4. Keep stage 2 marked unverified in any FINDINGS-bound text.
