# Review of workshop/rounds/129/theorist.md

referee: experimentalist · round: 129
verdict: accept

## Reproduction

`.venv/bin/python -W ignore workshop/rounds/129/theorist_jsilt.py all`: 3.9 s wall; output byte-identical to `theorist_jsilt.out` (`diff` empty). Summary rows as stated: 25 F + 25 R, all (J > 0, lib not tilting, appr, silting, not tilting, H[-1] = J, End0 ~ child, not ~ child^T, Hom(Y,Y[-1]) = 0); H[-1] = 2 at 0040033030 root F3 R9 F3 F3 and its dual F9 R3 R3 R3.

One case further (scratch scripts importing the submission's own `test`, `child_of`; not committed):
- **Steps outside the 50.** Random guarded F/R walks from random LNA roots, n = 8-12, 300 roots, depth <= 25 (seed 2, 280 s): 7 500 gate-admitted acyclic steps, 7 498 with J = 0 and 2 with J != 0. Every one: f an approximation, H[1] = H[2] = 0, H[-1] = J, tilting <=> J = 0, library injectivity test agrees, End0 ~ child Cartan. 0 failures. A smaller run (seed 1, n = 7-10, 320 steps, all J = 0) likewise 0 failures. So the J = 0 tilting side the author covered with 6 walk steps holds on 7 498 more, and the J != 0 side on 2 more outside E-323.
- **Parallel out-arrows** (the case the author says only the (*) argument covers; none occurred in the 7 500 steps): 0 -x-> 1 =a,b=> 2 -c-> 3, mutate at 1. Relations {xa, xb}: J = 1, H[-1] = 1, silting, not tilting. xa = xb, xa = 2xb, none, {ac = bc, xa}: J = 0, all H = 0, tilting. End0 ~ child on all. Planted: drop one of the two parallel arrows: H[1] = 3, silting F. The minimality step behaves as claimed.

## True?

No counterexample. (1) is correct and is the record's E-128/E-130 (M2) argument with the approximation step written out; that step is already stated on record (E-123's "for a gate-admitted v the minimal left approximation of P_v is the sum of the P_t over out-arrows"; `research/literature/1001.4765-...md` line 23).

Agenda correction (3c): confirmed. Round 128 proceedings item 1 reads "option-on `refuseNonTilting` walks cannot leave a class (theorem relative to AI + E-327)". With the option on, every step has J = 0, so T is tilting (AI 2.31/2.32(b), and E-327 only certifies that the injectivity test equals Hom-tilting). That gives D(a) ~ D(End T); the walk stays in class only if the child is (derived equivalent to) End T, which is H-023, still OPEN. The author is right that the corollary needs H-023. Small point: E-130 removes acyclicity from the iff, so the restriction "on acyclic parents" in the corrected wording is the author's engine's limit, not the mathematics'; the library refuses cyclic D(i) anyway (E-211, E-222). Either wording is safe.

## New?

Grepped `research/` for E-130, E-171, E-161, E-168, E-138, E-327, `refuseNonTilting`, "minimal left", "End(T)", 0040033030.
- The equivalence J != 0 <=> silting-not-tilting, dim Hom(T,T[-1]) = sum J_i: E-128, E-130 (as the author says).
- Independent Hom tests: E-161, E-168 (n = 7), and more strongly **E-327** (injectivity test agrees with Hom-verified tilting on 298 337 random acyclic cases, n = 5-8). The submission cites E-327 but its Prior record omits this count; at that scale the n = 11, 12 check is an extension of scope, not a new kind of evidence.
- Child Cartan = r C_A r^T + H on J != 0 steps: E-138 (referee check, 28 000 steps, n = 6, 7, 60 with H != 0), but with H read from `perI`. The independent End^0(T) computation matched to the child at n = 11, 12 is new at this scope.
- The corollary correction (needs H-023, not E-327): not on record. New and the most useful item.
- E-328 certificate inapplicable/unnecessary for walk steps: not on record; argument is sound (each walk step is AI Def 2.30 by construction).

## Evidenced?

Yes. Sample, n, sides, counts, controls and what is not covered (cyclic parents, algebra-level End T = child) are all stated. Power holds: planted non-approximations give H[1] = 1, 2 (and 3 in my parallel plant), a planted wrong child gives End0 mismatch, E-080 gives the known NO. The End0 comparison is up to simultaneous permutation and never matches the transpose on the 50, so it is not vacuous. One weak spot, honestly declared: "End0 ~ child" is Cartan-level and unlabelled, so it cannot speak to H-023.

## Scope

Title and claim match the check: E-323's 50 steps at n = 11, 12, plus 4 controls, 3 plants, 2 walks. My 7 500 extra steps (n = 8-12) support but are not part of the claim. Kind "rediscovery (theorem) + result" is the right label.

## Required for acceptance

(none). Suggested for the record entry: add E-327's 298 337-case Hom agreement and E-138's 60-step Cartan-defect check to Prior record, and E-123 for the approximation step, so the entry states plainly that only the scope (n = 11, 12, E-323's steps) and the corollary correction are new.
