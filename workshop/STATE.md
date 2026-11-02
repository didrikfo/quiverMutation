# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 107
cycle: 1 (rounds 101-120; CYCLE.md; plan in cycles/001/plan.md)
phase: 2
focus: research (rounds 105-108; conference 108; phase 1 summary in cycles/001/phase-1.md)
next_round_kind: conference

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 1 plan
See `cycles/001/plan.md`: H-024 and H-023 lead. Phase 1 (literature) done. Phase 2 agenda (round 104 chair): items 1-5 done in rounds 105-107 (E-281, E-284, E-283, E-285, E-287). Round 108 is the conference; the maverick's revision (H-024 power) is due in round 109. No more join searches.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** · H-015 · Can a key-guarded walk leave a derived class through a J != 0 step? Under the J = 0 premise all 70 failing n = 7 edges (12 presentation classes) are joined to an LNA (E-157..E-216, E-228); join searches cannot test the premise (E-250, E-255, E-268). 7 of 8 key-kept J != 0 edges are silting, not tilting (E-282). Needs a non-move-set certificate or closure; E-287 gives tilting paths B_n -> A_n at n = 5, 6. · 107 · toolsmith
- **T5 · J = 0 premise** · H-023 · Is every gate-admitted J = 0 step a derived equivalence? Cancellation proved (E-235); step-7 kernel = Oppermann's for acyclic D(i) (E-199); cyclic D(i) refused (E-211, E-222); library steps never meet a loop of degree >= 1 (E-279). Loop case: kM -> N cancels in ten hand-sized dg models, weight <= 12/9 (E-284), but Lemma A's monomial premise is false for loop letters, so no proof yet. Open: degree-0 loop control, proof of Prop. 7.6 filtration, Anick numbers UNVERIFIED (E-278). · 106 · theorist, toolsmith
- **T3/T8 · Cartan-level invariants** · H-024 (SUPPORTED; also H-025) · Necessity sourced (E-272). Mixed-digit positive controls now exist: Fosse Cor. 4.6, 6 pairs n = 9..12, 48 pass Cartan tests (E-275; Cor. 4.6 has a misprint). E-269's bound is closed-form n-2+lam+1/lam, so the small-lam gap is missing input, not a bound gap (E-276); `maverick_iso140.py` l.16 `>1.5` must change. Round 105: E-269 YES on the 6 mixed pairs and 11 low-lam within-quipu pairs, l.16 fixed (E-281). Open: n = 12 low-lam groups, NO at lam<1.5 (no pair at n<=12), Nikulin 1.13.3 unverified. · 105 · maverick, skeptic
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). Two-bystander cells at k = 3..5 stay at overlap 2 (E-230, E-241, E-260), but no start-gap-2 pair at k <= 4 can lower, so they have no control (E-288, no power). Gate part published for gentle LNAs only (E-280). Open: a gap-2 shape that can lower, k >= 6. · 107 · experimentalist
- **T4 · H-020 rule table** · H-020 · H1 holds to widths 9-11, lengths <= 14 (E-225), independent engines (E-207, E-243). Class-B rules: 298 of 304 leave the LNA class at step 1, 6 at step 2 (E-266); first step adds one arrow at widths 5-7, but no width-independent relation-count rule (E-286). H3/H4 not closed; last item "narrowed only". · 106 · skeptic
- **H-027 · A_n vs B_n reachability** · H-027 (SUPPORTED) · A_n ~ B_n (Fosse Cor. 4.6) but A_n not in orbit(B_n) (E-280). Tilting exchanges of rigid K0-unimodular complexes reach A_n from B_n at depth 3 (n = 5), 4 (n = 6), End = A_n checked (E-287); generation (Bruning) UNVERIFIED; not shown to be silting mutations; n >= 7 not run; E-282 one-step NO. · 107 · toolsmith
- **S-1 · vertex deletion** (suggested question) · H-026 · K-threshold law K0 = 3, 4, 5 at n = 11, 13, 15, key level (E-146, E-263); free-stretch D >= 3 rule (E-206); |orbit B_n| = n 2^(n-5) to n = 17 is a move-orbit size, not the class: A_n ~ B_n for all n >= 5 (Fosse Cor. 4.6) yet A_n is not in the orbit (E-280); none of it is published (E-280). Item 3 never attempted. · 103 · experimentalist
- **T6 · H-017** · H-017 · Defect is NOT constant per key class at n = 9 depth 5 (E-283); only 2 key classes occur; constant only on the defect<=1 filter. Open: depth 6, n >= 10, a statistic with power. · 105 · experimentalist
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
- maverick: answer `rounds/107/maverick.review.md` (major revision; H-024 power at n = 12 low lam): drop |det S| as evidence, reconcile the group count (15 vs 12), rescope, planted-error control, say E-276 restated. Due round 109 (108 is a conference).

## Requests between personas
- scholar: Anick primary (E-273); Happel-Seidel / BGS / Nikulin 1979 primaries (not on arXiv: need PDFs).
- maverick: revision of round 107 (above); E-269 for singular S (2 of the 6 n = 12 low-lam groups); a certified-inequivalent pair at lam < 1.5 (NO power).
- theorist: the intermediate End of E-287's depth-4 path: LNAs? where does it leave the gate? Bruning generation and silting-order check (scholar); n = 7..9 depth to first hit (toolsmith); rebuild the E-239 depth-5 edge (E-282).
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
- Round 107: H-027 SUPPORTED by E-287 (n = 5, 6). Round 104 (conference): no F/R change then; E-269 on mixed pairs done (E-281: YES). Round 106: Fosse Cor. 4.6 checked against the LaTeX (E-285), R-014 written for E-218's title.
- The 7 + 7 n = 9 F-010 J = 0 join was run once by the dispatcher (E-250: no meet, no power).
- arXiv is reachable (STEERING, standing permissions).

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 107 | 106 (scholar) |
| theorist | 106 | 107 (experimentalist) |
| skeptic | 106 | 107 (toolsmith) |
| scholar | 106 | 107 (maverick) |
| toolsmith | 107 | 103 (theorist) |
| maverick | 107 | 103 (scholar) |
