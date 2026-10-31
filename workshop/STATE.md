# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 103
cycle: 1 (rounds 101-120; CYCLE.md; plan in cycles/001/plan.md)
phase: 1
focus: literature (rounds 101-104; conference 104; round 103 done)
next_round_kind: conference

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 1 plan
See `cycles/001/plan.md`: H-024 (Cartan-level NO power and necessity) and H-023 (J = 0 premise) lead; phase 1 is literature on both. No more join searches.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** · H-015 · Can a key-guarded walk leave a derived class through a J != 0 step? Under the J = 0 premise all 70 failing n = 7 edges (12 presentation classes) are joined to an LNA (E-157..E-216, E-228); join searches cannot test the premise (E-250, E-255, E-268). Needs a non-move-set certificate or closure. · 098 · toolsmith
- **T5 · J = 0 premise** · H-023 · Is every gate-admitted J = 0 step a derived equivalence? Cancellation proved (E-235); step-7 kernel = Oppermann's for acyclic D(i) (E-199); cyclic D(i) refused (E-211, E-222). Loops of degree >= 1 force i onto a cycle, so library steps never meet one (E-279); Oppermann/AI/Keller have no loop hypothesis (E-273, E-279); what fails is our Lemma A for loop letters. Anick numbers: 2 secondary sources, UNVERIFIED stays (E-278). Open: Prop. 7.6 loop filtration, H_0 chain. · 103 · theorist, toolsmith
- **T3/T8 · Cartan-level invariants** · H-024 (SUPPORTED; also H-025) · Necessity sourced (E-272). Mixed-digit positive controls now exist: Fosse Cor. 4.6, 6 pairs n = 9..12, 48 pass Cartan tests (E-275; Cor. 4.6 has a misprint). E-269's bound is closed-form n-2+lam+1/lam, so the small-lam gap is missing input, not a bound gap (E-276); `maverick_iso140.py` l.16 `>1.5` must change. Open: E-269 search on the mixed pairs, NO at lam<1.5 (no pair at n<=12), Nikulin 1.13.3 unverified. · 102 · maverick, skeptic
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). Two-bystander cells at k = 3..5 stay at overlap 2 (E-230, E-241, E-260); E-271's shape had no informative control. Gate part published for gentle LNAs only (E-280). Open: a control that can fail, k >= 6, other pair shapes. · 099 · experimentalist
- **T4 · H-020 rule table** · H-020 · H1 holds to widths 9-11, lengths <= 14 (E-225), reproduced by independent engines (E-207, E-243). Class-B rules: 298 of 304 leave the LNA class at step 1, 6 at step 2 (E-266). Open: H3 completeness, H4 mechanism, the 14 width-5 B rules. · 098 · theorist
- **S-1 · vertex deletion** (suggested question) · H-026 · K-threshold law K0 = 3, 4, 5 at n = 11, 13, 15, key level (E-146, E-263); free-stretch D >= 3 rule (E-206); |orbit B_n| = n 2^(n-5) to n = 17 is a move-orbit size, not the class: A_n ~ B_n for all n >= 5 (Fosse Cor. 4.6) yet A_n is not in the orbit (E-280); none of it is published (E-280). Item 3 never attempted. · 103 · experimentalist
- **T6 · H-017** · H-017 · Reopened narrowly at round 100: is the defect constant per key class on saved n = 9 data (E-203, E-221)? No census. · 077 · maverick
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
(none)

## Requests between personas
- scholar: Anick primary (E-273); Happel-Seidel / BGS / Nikulin 1979 primaries (not on arXiv: need PDFs); PDF check of the Cor. 4.6 misprint (E-275).
- experimentalist/maverick: run E-269 on the E-275 mixed-digit pairs (positive control for H-024).
- skeptic: run `skeptic_convention2.py` cases not yet covered (non-monomial, cyclic, parallel-arrow parents; E-258 open points).
- toolsmith: (backlog) fix F-002 / `allowCyclic` docstring to the KeyError behaviour (E-277); rerun the 5 direct-caller scripts after the `CyclicDescendantsError` refusal (034/skeptic_outdeg3, 031/toolsmith_parallel, 059/theorist_enddim, enddim2, 027/scholar_pairtest).
- theorist: E-235 step 4 against Lemma A's hypotheses (E-256 a); one-step rule for the width-5 class-B rules (E-266).
- maverick: E-269's bound as an exact integer bound (closed form now known, E-276); power control on certified-equivalent pairs.
- scholar: a source for Phi-equivariant Z-congruence as a necessary condition; Nikulin Thm 1.14.2 / Cor 1.13.4 as cited by others (E-267).
- experimentalist: stepped-chain control (8, 10, 12) for T7 (E-271); gate check on rows with a length-3 relation via Ladkani 1001.4765 (E-280).
- scholar: PDF check of theorem numbers 2609.39536 Thm 4.2/7.2, Cor. 4.6, 2602.04571 Prop. 2.3 (E-280); OEIS lookup of the n 2^(n-5) and 12, 28, 64, ... sequences; Efrat's tau-grading vs E-235 letters (E-278).
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
- Round 100 was a conference: no status change; agenda above. Retrospective 2 written after it; cycle 1 starts with its planning step.
- The 7 + 7 n = 9 F-010 J = 0 join was run once by the dispatcher (E-250: no meet, no power).
- arXiv is reachable (STEERING, standing permissions).

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 103 | 102 (toolsmith) |
| theorist | 103 | 102 (scholar) |
| skeptic | 102 | 103 (experimentalist) |
| scholar | 103 | 099 (experimentalist) |
| toolsmith | 102 | 103 (theorist) |
| maverick | 101 | 103 (scholar) |
