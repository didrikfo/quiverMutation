# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 140
cycle: 3 (rounds 141-160; CYCLE.md; plan in cycles/003/plan.md once the planning step has run). Cycle 2 closed: rounds/140/proceedings.md, cycles/002/backlog.md, cycles/002/retrospective.md
phase: none yet (planning step next)
focus: research by default; the cycle 3 plan may give a phase or personas a focus (literature, backlog, replication, problem session; CYCLE.md). Owed code from cycles/002/backlog.md ("Carried", Code) goes in the first ordinary round.
next_round_kind: planning

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151, E-332), and whether a guarded walk leaves a derived class through a J != 0 step is settled only at n = 11 (F-057); n = 12 and depth >= 8 are open. Cycle 2 retrospective: cycles/002/retrospective.md; cycle 1: cycles/001/retrospective.md; retrospective 2: retrospectives/002-rounds-045-100.md.

## Cycle 3 agenda (proposed at round 140; the planning step revises it)
1. H-023 at algebra level: End(T) vs child on E-357's 14 steps with a planted wrong-ideal control (also F-057's blocker); KY/DWZ/Oppermann hypothesis match.
2. F-055 NO side without q_S: planted same-charpoly different-class pair; certificate for the group of 140.
3. F-058: why two tau-orbit families (parity split first; N_16(8); N_15(7) control).
4. T10 at n = 12 (group of 47, depth 8), E-316 depth 5, after the two small code refusals (toolsmith).
5. H-010 proof; S-1 n = 17 h >= 4.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** · H-015 (R-017), H-032 · F-057. `refuseNonTilting` default OFF (decided 140; E-350: loses 93-180 in-class nodes per depth-4 walk at E-316, nothing measured cross-class). E-354, E-357: n = 11 rebuild and Hom check agree. Open: n = 12 group of 47 at depth 8, E-316 depth 5, 4 duals (overnight-sized, not approved). · 139 · toolsmith
- **T5 / H-023 · J = 0 premise, child = End(T)** · H-023, H-031 · E-345: child ~ End(T) at 587 of 587 decided J = 0 steps (n <= 6 acyclic); E-348: 341/341 planted wrong children NO. Open: algebra-level End(T) on the 14 E-357 steps, n >= 7, cyclic children, 5 undecided, Anick (H3)(ii) UNVERIFIED, H-031 (carried). · 135 · skeptic
- **T3/T8 · Cartan-level invariants** · H-024 (F-055), H-025 (F-058) · E-353, E-355, E-356. Open: F-055 NO side rests on q_S and the E-308/E-272 chain (Nikulin 1.13.3 UNVERIFIED); group of 140, singular S; why two nu-orbit families; m = 132 control; non-uniserial orbits. · 139 · theorist, maverick
- **T7 · H-010** · E-343: no clean interior (8,4)(10,3) cell lowers overlap; 9 capped cells, k = 5 overnight-sized; theorist proof owed. · 133 · experimentalist
- **S-1 · vertex deletion** · H-026 (SUPPORTED, odd n only) · K0 = min(h, n-4-h) at key level (E-301, E-304, E-349); even n fails (E-307). Open: tail end, non-lone law, n = 17 h >= 4 (prediction logged). · 135 · experimentalist
- **Closed:** T11 (H-030), T6 (H-017b), H-027 (F-056), H-029 (retired), T1, T2; **T9 long runs** parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
(none)

## Requests between personas
- (round 140) toolsmith: End(T) builder on the E-357 steps, nu over non-uniserial objects, exact/interval certificate (group of 140), the two small code refusals; skeptic: planted same-Cartan wrong ideal at n = 11, same-root-count control for scholar_isometry.py; theorist: Z[Phi] planted pair, (p, r) formula, H-010 proof, isometry for all n >= 5; scholar: KY/DWZ/Oppermann match, source for N_15(7), parity split; maverick: N_16(8); experimentalist: n = 12 group of 47.
- scholar: PDFs wanted (human, standing): Nikulin 1979, Anick 1982, BGS, Happel I.5, Chekanov; owed files listed in cycles/002/backlog.md. Chair housekeeping: E-308 'co:complete' -> A.7, README row for 1310.1907, F-044 N_9(5) erratum.

## Backlog
Closed at 140: see cycles/002/backlog.md (settled, retired, converted, carried).

## Weakest claims the workshop currently relies on
The J = 0 premise / child = End(T) (n <= 6 acyclic monomial; E-357 tests code); the Cartan-level NO (q_S, E-308/E-272 chain, Nikulin 1.13.3 UNVERIFIED); "two families = two orbits" (uniserials only); S-1's K-law (key level, near-definitional zero); E-151 counts are capped BFS samples.

## Standing facts
- arXiv is reachable (STEERING, standing permissions).
- Cycle 2 findings: F-055 (Cartan-level test, scoped), F-056 (B_n to A_n, n = 5-7), F-057 (guarded walks leave a class at n = 11, 12), F-058 (periodic not FCY), R-017, R-018.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 138 | 139 (scholar) |
| theorist | 139 | 139 (skeptic) |
| skeptic | 138 | 139 (theorist) |
| scholar | 139 | 138 (experimentalist) |
| toolsmith | 139 | 138 (maverick) |
| maverick | 138 | 137 (scholar) |
