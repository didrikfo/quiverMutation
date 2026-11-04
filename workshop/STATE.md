# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 113
cycle: 1 (rounds 101-120; CYCLE.md; plan in cycles/001/plan.md)
phase: 4
focus: research (rounds 114-116; conference 116; phase summaries in cycles/001/phase-1.md .. phase-3.md)
next_round_kind: ordinary

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 1 plan
See `cycles/001/plan.md`. Phases 1-3 done (phase-1.md .. phase-3.md). Phase 4 agenda (round 112, ranked; round 113 did 1 and 2 in part): (1) T10 step 8 and a NO-capable invariant (R-candidate); (2) S-1/H-026 tail ends directly and a non-lone class (F-candidate); (3) T5: lemma for finite-dimensional non-loop parents, Lemma A at degree -1 to weight 12, M4/M5 (theorist); (4) T6 statistic with power (shuffle control) or T7 gap-2 at k = 6; (5) H-025 periodic-not-FCY LNA (scholar). Blocking checks: items 1 and 2 (remaining parts).

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** (agenda 1) · H-015 · Can a key-guarded walk leave a derived class through a J != 0 step? Hom-verified chains join 12 of E-201's 13 steps at n = 7 class 1 (E-297; meet up to coefficients, generation assumed); step 8 undecided. No NO case at n <= 7 (E-298, E-295): a NO needs an invariant with a power control (n = 9 F-010 pair hit the cap, E-212). · 113 · toolsmith
- **T5 · J = 0 premise** (agenda, round 112: 3) · H-023 · Is every gate-admitted J = 0 step a derived equivalence? Cancellation proved (E-235); loop models M1 (E-290) and M2 (E-294, formal dg algebra, g-filtration to Lemma A) proved; degree-0 loop cancels at w = 3 to weight 12 (one model). M4, M5 (shared/product pivots), Hom/End(T) for M2, Anick numbers UNVERIFIED. · 111 · theorist
- **T3/T8 · Cartan-level invariants** (H-025 is agenda 5) · H-024 (SUPPORTED; also H-025) · Necessity sourced (E-272; round 110 recheck, Phi-equivariance only inside the proof of Ladkani Lemma 3.14). E-292: NO on 19 of 19 q_S-separated pairs (lam >= 1.661), YES on 19 controls; NO is never pure lattice search. Open: lam < 1.5 NO (no LNA pair exists there, E-276), Z[Phi]-class route, Nikulin 1.13.3 unverified. · 110 · skeptic
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). Two-bystander cells at k = 3..5 stay at overlap 2 (E-230, E-241, E-260), but no start-gap-2 pair at k <= 4 can lower, so they have no control (E-288, no power). Gate part published for gentle LNAs only (E-280). Open: a gap-2 shape that can lower, k >= 6. · 107 · experimentalist
- **T4 · H-020 rule table** · H-020 · H1 holds to widths 9-11, lengths <= 14 (E-225), independent engines (E-207, E-243). Class-B rules: 298 of 304 leave the LNA class at step 1, 6 at step 2 (E-266); first step adds one arrow at widths 5-7, but no width-independent relation-count rule (E-286). H3/H4 not closed; last item "narrowed only". · 106 · skeptic
- **H-027 · A_n vs B_n reachability** · H-027 (SUPPORTED) · A_n ~ B_n (Fosse Cor. 4.6) but A_n not in orbit(B_n) (E-280). Tilting exchanges reach A_n from B_n at depth 3, 4, 5 (n = 5, 6, 7; E-287, E-289); intermediate End(T) not LNAs, each step gate-admitted; generation (Bruning) UNVERIFIED; not shown to be silting mutations. · 109 · toolsmith
- **S-1 · vertex deletion** (suggested question; agenda 2) · H-026 (SUPPORTED) · K-threshold K0 = (n-1)/2 - 2 at key level, lone-3 class, head ends, n = 11..21 (E-296, E-299: independent code, n = 21 prediction held, lower edge splits too). Not class level, not non-lone, tail only by mirror test; no derivation. F-blockers: tail directly at one n, a non-lone class. · 113 · experimentalist
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
- The 7 + 7 n = 9 F-010 J = 0 join was run once by the dispatcher (E-250: no meet, no power).
- arXiv is reachable (STEERING, standing permissions).

- Round 113: E-297 (12 of 13 chains, step 8 open), E-298 (controls), E-299 (K-law n = 21). No F/R change; R-candidate still blocked, F-candidate H-026 blocked on tail/non-lone checks.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 113 | 109 (theorist) |
| theorist | 111 | 113 (toolsmith) |
| skeptic | 113 | 113 (experimentalist) |
| scholar | 110 | 111 (toolsmith) |
| toolsmith | 113 | 110 (skeptic) |
| maverick | 111 | 113 (skeptic) |
