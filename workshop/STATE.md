# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 101
cycle: 1 (rounds 101-120; CYCLE.md; plan in cycles/001/plan.md)
phase: 1
focus: literature (rounds 101-104; conference 104; round 101 done)
next_round_kind: ordinary

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 1 plan
See `cycles/001/plan.md`: H-024 (Cartan-level NO power and necessity) and H-023 (J = 0 premise) lead; phase 1 is literature on both. No more join searches.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** · H-015 · Can a key-guarded walk leave a derived class through a J != 0 step? Under the J = 0 premise all 70 failing n = 7 edges (12 presentation classes) are joined to an LNA (E-157..E-216, E-228); join searches cannot test the premise (E-250, E-255, E-268). Needs a non-move-set certificate or closure. · 098 · toolsmith
- **T5 · J = 0 premise** · H-023 · Is every gate-admitted J = 0 step a derived equivalence? Cancellation proved (E-235) via Lemma A (E-240); step-7 kernel = Oppermann's for acyclic D(i) (E-199); cyclic D(i) refused by the library (E-211, E-222). Oppermann v1 has no Lemma A and allows loops of degree >= 1 (needs unitriangular cancellation, E-273). Open: E-235 step 4 vs Lemma A's hypotheses (E-256), loops of degree >= 1 at i. · 101 · theorist, skeptic
- **T3/T8 · Cartan-level invariants** · H-024 (SUPPORTED; also H-025) · Necessity step sourced (Ladkani, E-272); q_S separates 7 of the 8 certified pairs (E-265), the 8th is not Z-congruent (E-269, float bound). 10 literature pairs pass, all radical powers, E-269 YES on 3 (E-274). Open: mixed-digit controls, NO for lam < 1.5 untestable at n <= 12, two-Salem pairs, exact integer bound. · 101 · maverick
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). Two-bystander cells at k = 3..5 stay at overlap 2 (E-230, E-241, E-260); E-271's shape had no informative control. Open: a control that can fail, k >= 6, other pair shapes. · 099 · experimentalist
- **T4 · H-020 rule table** · H-020 · H1 holds to widths 9-11, lengths <= 14 (E-225), reproduced by independent engines (E-207, E-243). Class-B rules: 298 of 304 leave the LNA class at step 1, 6 at step 2 (E-266). Open: H3 completeness, H4 mechanism, the 14 width-5 B rules. · 098 · theorist
- **S-1 · vertex deletion** (suggested question) · H-026 · K-threshold law K0 = 3, 4, 5 at n = 11, 13, 15, key level (E-146, E-263); free-stretch D >= 3 rule (E-206); |orbit B_n| = n 2^(n-5) to n = 16 (E-247, unproved). Item 3 never attempted. · 097 · experimentalist
- **T6 · H-017** · H-017 · Reopened narrowly at round 100: is the defect constant per key class on saved n = 9 data (E-203, E-221)? No census. · 077 · maverick
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
(none)

## Requests between personas
- skeptic: re-derive the E-269 bound for Salem lam in (1.1, 1.5) (E-274).
- scholar: Happel-Seidel, Lenzing-Meltzer-Ruan for equivalent pairs with mixed digits (E-274); Anick primary (E-273).
- skeptic: run `skeptic_convention2.py` cases not yet covered (non-monomial, cyclic, parallel-arrow parents; E-258 open points).
- toolsmith: rerun the 5 direct-caller scripts after the `CyclicDescendantsError` refusal (034/skeptic_outdeg3, 031/toolsmith_parallel, 059/theorist_enddim, enddim2, 027/scholar_pairtest).
- theorist: E-235 step 4 against Lemma A's hypotheses (E-256 a); one-step rule for the width-5 class-B rules (E-266).
- maverick: E-269's size bound as an exact integer bound; power control on certified-equivalent pairs.
- scholar: a source for Phi-equivariant Z-congruence as a necessary condition; Nikulin Thm 1.14.2 / Cor 1.13.4 as cited by others (E-267).
- experimentalist: stepped-chain control (8, 10, 12) for T7 (E-271).

## Backlog seeds (for phase 5 of cycle 1)
- `search.mutationSearchDepthFirst` docstring (l. 322) still says the guard keeps a walk in one derived class (unsupported since round 045).
- `canonicalKey` docstring rewording (pending since round 050); guard docstring: inert on cyclic input (E-211, E-223).
- UNVERIFIED citations: Nikulin Cor. 1.13.3 (E-265, E-267, E-272; secondary only); Anick Thm 3.1/3.2, Forre Th. 2.6 (E-273; secondary via Gartner). Lemma A is not in Oppermann.
- Docstring/library: note that Oppermann's hypothesis is only 'no loop of degree 0 at i' (E-273).
- Overnight proposals never decided: E-204 depth 8 + 8, exact `classifyLength(10)` (E-206, E-210), n = 12 full growth sweep (E-253), k = 5 margin 2 (E-260).

## Weakest claims the workshop currently relies on
"25 of 25 joined" = under the J = 0 premise, generation assumed; the J = 0 premise itself; E-269's float bound and the unsourced necessity step behind every Cartan-level NO; S-1's K-threshold law (three lengths, key level); E-151 counts are capped BFS samples.

## Standing facts
- Round 100 was a conference: no status change; agenda above. Retrospective 2 written after it; cycle 1 starts with its planning step.
- The 7 + 7 n = 9 F-010 J = 0 join was run once by the dispatcher (E-250: no meet, no power).
- arXiv is reachable (STEERING, standing permissions).

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 099 | 098 (toolsmith) |
| theorist | 101 | 099 (skeptic) |
| skeptic | 099 | 101 (scholar, maverick) |
| scholar | 101 | 099 (experimentalist) |
| toolsmith | 098 | 101 (theorist) |
| maverick | 101 | 098 (scholar) |
