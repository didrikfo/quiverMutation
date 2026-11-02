# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 109
cycle: 1 (rounds 101-120; CYCLE.md; plan in cycles/001/plan.md)
phase: 3
focus: research (rounds 109-112; conference 112; phase summaries in cycles/001/phase-1.md, phase-2.md)
next_round_kind: ordinary

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 1 plan
See `cycles/001/plan.md`: H-024 and H-023 lead. Phase 1 (literature) done. Phase 2 done (cycles/001/phase-2.md). Phase 3 agenda (round 108, ranked; round 109 did items 1-3 in part): (1) H-024 power: a pair inequivalent by construction in E-269's domain (E-291 design); (2) T10/H-027 audit of the E-287 paths: intermediate End(T), gate/key verdict per step, then n = 7, generation chain; (3) T5 loop case: two-stage proof, Hom-computed End(T) for M1, degree-0 loop control; (4) S-1 item 3 K0 = 6 at n = 17, or H-025 converse. Blocking checks: F-candidate H-024 (planted-error NO), R-candidate 'guard keeps one class' (item 2). No more join searches.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** (agenda 2) · H-015 · Can a key-guarded walk leave a derived class through a J != 0 step? Round 109 (E-289): on all 62 steps of the n = 5..7 E-287 paths the gate admits the matching side through non-LNA End(T); the key column is implied by tilting, no key-refusal control yet; 38 of 124 cells are gate-True but not End(T') (gate is no tilting test). Needs a J != 0 key-refusal control and a closure/certificate. · 109 · toolsmith
- **T5 · J = 0 premise** (agenda 3) · H-023 · Is every gate-admitted J = 0 step a derived equivalence? Cancellation proved (E-235); loop case: M1 proved by change of generators + Kunneth, Hom control tilting (E-290); Lemma A's three hypotheses fail for loop letters; M2, M4, M5 (shared/product pivots) and degree-0 loops open; N ~ A in all degrees, Anick numbers UNVERIFIED. · 109 · theorist
- **T3/T8 · Cartan-level invariants** (agenda 1) · H-024 (SUPPORTED; also H-025) · Necessity sourced (E-272); mixed-digit positive controls pass (E-275, E-281); E-269 bound closed-form (E-276). Round 109 (E-291): YES on 4 n = 12 low-lam pairs, planted +-1 control gave no NO (157 congruent): power at lam < 1.5 not demonstrated. Open: a pair inequivalent by construction (different q_S or Z[Phi]-class), Nikulin 1.13.3 unverified. · 109 · maverick, skeptic
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). Two-bystander cells at k = 3..5 stay at overlap 2 (E-230, E-241, E-260), but no start-gap-2 pair at k <= 4 can lower, so they have no control (E-288, no power). Gate part published for gentle LNAs only (E-280). Open: a gap-2 shape that can lower, k >= 6. · 107 · experimentalist
- **T4 · H-020 rule table** · H-020 · H1 holds to widths 9-11, lengths <= 14 (E-225), independent engines (E-207, E-243). Class-B rules: 298 of 304 leave the LNA class at step 1, 6 at step 2 (E-266); first step adds one arrow at widths 5-7, but no width-independent relation-count rule (E-286). H3/H4 not closed; last item "narrowed only". · 106 · skeptic
- **H-027 · A_n vs B_n reachability** (agenda 2) · H-027 (SUPPORTED) · A_n ~ B_n (Fosse Cor. 4.6) but A_n not in orbit(B_n) (E-280). Tilting exchanges reach A_n from B_n at depth 3, 4, 5 (n = 5, 6, 7; E-287, E-289); intermediate End(T) not LNAs, each step gate-admitted; generation (Bruning) UNVERIFIED; not shown to be silting mutations. · 109 · toolsmith
- **S-1 · vertex deletion** (suggested question; agenda 4, round 108) · H-026 · K-threshold law K0 = 3, 4, 5 at n = 11, 13, 15, key level (E-146, E-263); free-stretch D >= 3 rule (E-206); |orbit B_n| = n 2^(n-5) to n = 17 is a move-orbit size, not the class: A_n ~ B_n for all n >= 5 (Fosse Cor. 4.6) yet A_n is not in the orbit (E-280); none of it is published (E-280). Item 3 never attempted. · 103 · experimentalist
- **T6 · H-017** · H-017 · Defect is NOT constant per key class at n = 9 depth 5 (E-283); only 2 key classes occur; constant only on the defect<=1 filter. Open: depth 6, n >= 10, a statistic with power. · 105 · experimentalist
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
(none)

## Requests between personas
- scholar: Anick primary (E-273); Happel-Seidel / BGS / Nikulin 1979 primaries (not on arXiv: need PDFs).
- maverick: E-269 for singular S (2 of the 6 n = 12 low-lam groups); planted controls on groups 2, 4; a pair inequivalent by construction (E-291).
- toolsmith: key-refusal control (a J != 0 step) for T10; theorist: loop models M2, M4, M5. Done in 109: intermediate End of E-287 paths are not LNAs, none gate-refused (E-289). Bruning generation and silting-order check (scholar); n = 7..9 depth to first hit (toolsmith); rebuild the E-239 depth-5 edge (E-282).
- skeptic: run `skeptic_convention2.py` cases not yet covered (non-monomial, cyclic, parallel-arrow parents; E-258 open points).
- toolsmith: (backlog) fix F-002 / `allowCyclic` docstring to the KeyError behaviour (E-277); rerun the 5 direct-caller scripts after the `CyclicDescendantsError` refusal (034/skeptic_outdeg3, 031/toolsmith_parallel, 059/theorist_enddim, enddim2, 027/scholar_pairtest).
- theorist: E-235 step 4 against Lemma A's hypotheses (E-256 a); one-step rule for the width-5 class-B rules (E-266).
- maverick: E-269's bound as an exact integer bound (closed form now known, E-276); power control on certified-equivalent pairs.
- scholar: a source for Phi-equivariant Z-congruence as a necessary condition; Nikulin Thm 1.14.2 / Cor 1.13.4 as cited by others (E-267).
- experimentalist: a gap-2 shape that can lower overlap for T7 (E-288); gate check on rows with a length-3 relation via Ladkani 1001.4765 (E-280).
- theorist/skeptic: is A_n reachable from B_n by any tilting mutation outside the E-187 moves (T10 at n = 6)? Cartan-data search on small silting-not-tilting pairs (E-278).

## Backlog seeds (for phase 5 of cycle 1)
- `search.mutationSearchDepthFirst` docstring (l. 322) still says the guard keeps a walk in one derived class (unsupported since round 045).
- `canonicalKey` docstring rewording (pending since round 050); guard docstring: inert on cyclic input (E-211, E-223).
- Cor. 4.6 of Fosse 2305.06642 misprint (E-275) to check against the PDF; `allowCyclic` docstring and F-002 text stale (E-277).
- UNVERIFIED citations: Nikulin Cor. 1.13.3 (E-265, E-267, E-272; secondary only); Anick Thm 3.1/3.2, Forre Th. 2.6 (E-273; secondary via Gartner). Lemma A is not in Oppermann.
- Docstring/library: note that Oppermann's hypothesis is only 'no loop of degree 0 at i' (E-273).
- Overnight proposals never decided: E-204 depth 8 + 8, exact `classifyLength(10)` (E-206, E-210), n = 12 full growth sweep (E-253), k = 5 margin 2 (E-260).

## Weakest claims the workshop currently relies on
"25 of 25 joined" = under the J = 0 premise, generation assumed; the J = 0 premise itself; E-269's float bound and the unsourced necessity step behind every Cartan-level NO; S-1's K-threshold law (three lengths, key level); E-151 counts are capped BFS samples.

## Standing facts
- Round 108 (conference): no F/R change; F-candidate H-024 and R-candidate 'guard keeps one class' have named blocking checks (above). Skeptic's caveat: A_n ~ B_n is Fosse Cor. 4.6, so E-287 shows a short path, not a new equivalence.
- Round 107: H-027 SUPPORTED by E-287 (n = 5, 6). Round 104 (conference): no F/R change then; E-269 on mixed pairs done (E-281: YES). Round 106: Fosse Cor. 4.6 checked against the LaTeX (E-285), R-014 written for E-218's title.
- The 7 + 7 n = 9 F-010 J = 0 join was run once by the dispatcher (E-250: no meet, no power).
- arXiv is reachable (STEERING, standing permissions).

- Round 109: E-289 (paths pass the gate; key column no evidence), E-290 (M1 loop case), E-291 (power not demonstrated at lam < 1.5). No F/R change.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 107 | 109 (theorist) |
| theorist | 109 | 107 (experimentalist) |
| skeptic | 106 | 109 (maverick) |
| scholar | 106 | 109 (toolsmith) |
| toolsmith | 109 | 103 (theorist) |
| maverick | 109 | 103 (scholar) |
