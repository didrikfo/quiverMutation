# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 137
cycle: 2 (rounds 121-140; CYCLE.md; plan in cycles/002/plan.md, written 2026-11-07)
phase: 5 (round 2 of 4; rounds 138-139 backlog and replication sweep, 140 the closing conference)
focus: backlog and replication (cycles/002/plan.md, Phase focus: independent rebuild of this cycle's F-entries, backlog closures, T10 decision); 140 is the closing conference. Cycle 1 retrospective: cycles/001/retrospective.md.
next_round_kind: ordinary

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 2 plan
See cycles/002/plan.md and cycles/002/phase-1.md (what phase 1 did). Agenda for phase 5 = round-136 proceedings (rounds/136/proceedings.md); phase 4 history in cycles/002/phase-4.md (E-341 to E-349; F-058 written at conference 136), phase 3 in phase-3.md.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** (agenda 1) · H-015 (R-017), H-032 · F-057. E-350 (137): option-on loses 0 in-class nodes at E-332's n = 11 roots (10 depth-7 walks) only because no in-class refusal occurs there; at E-316's 19 n = 9 in-class J != 0 steps it loses 93-180 in-class nodes per depth-4 walk (child lost 19 of 19). Default stays off; decide formally at 140. Open: depth 8 (sharded, overnight-sized), 4 duals, depth 5+ at E-316. · 137 · toolsmith
- **T5 / H-023 · J = 0 premise, child = End(T)** (agenda, round 136: 3) · H-023, H-031 · E-345 (134): child ≅ End(T) at 587 of 587 decided J = 0 steps, n <= 6 acyclic; E-348 (135): the test says NO on 341 of 341 planted same-Cartan wrong children, ISO on 978 of 978 controls (NOs by Groebner by construction). E-329/E-340/R-018 stand. Open: 5 undecided, cyclic children, n >= 7, non-monomial cyclic Hom engine, Anick (H3)(ii). · 135 · skeptic
- **T3/T8 · Cartan-level invariants** · H-024 (F-055), H-025 (CONFIRMED -> F-058, round 136) · E-346/E-347. Open: second prime for the 4 uniserial-only groups, FCY size control at n = 14-15, the 3 one-shift groups, why two nu-orbit families (maverick: tau-orbits); Nikulin 1.13.3 UNVERIFIED. · 135 · maverick
- **T7 · H-010** (agenda, round 136: 5) · H-010 · E-343 (133): no clean interior (8,4)(10,3) cell lowers overlap (engine k = 4: 0 of 15 decided; moves at W = 13, 14 bands: 0 of 5 390 / 17 738); only a non-clean (3,3)-pair control lowers (F-022 mechanism), so E-330's null there has harness power but no clean control. Open: theorist proof that clean orbits avoid lowering rows; 9 capped cells, k = 5 overnight-sized. · 133 · experimentalist
- **S-1 · vertex deletion** (agenda, round 120: 4) · H-026 (SUPPORTED, odd n only) · K0 = min(h, n-4-h) for lone-3 classes at key level (E-301, E-304); even-n fails at the self-mirror class (E-307). E-325, E-349 (135): key changes only through K = 0 on eight classes (n = 13, 15; zero is near-definitional); boundary changes are directed K = 0 -> K = 1 `rules` rewrites. F-blockers: tail end, non-lone law, self-mirror derivation; class-level form. · 135 · experimentalist
- **T11 · gate vs Hom test** · H-030 (SUPPORTED, closed 137). E-351: `isTiltingAtVertex` = Hom criterion on acyclic quivers with one-source/one-target relations; over-admissions are E-080/E-112 shapes. Gate does not call it by default. Closed except an optional docstring/refusal on illegal relations. · 137 · theorist
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
(none)

## Requests between personas
- (round 135) experimentalist: S-1 tail end (n = 13, 15) and n = 17 h >= 4 as a test of the boundary pattern; skeptic: wrong non-monomial ideals other than KILL/SCALE for E-345; maverick: nu2 at other primes/m for the uniserial groups (GF(2^31-1) only).
- (round 134) toolsmith: faster ISO for >= 2-arrow slots (5 undecided, E-345); second nu implementation at n = 14-15 (E-346 blocker); option-on depth-9 walk and the other 5 E-332 roots (E-344). Theorist: why two nu-orbit families (E-346); cyclic-child End(T) steps. Skeptic: planted same-Cartan wrong monomial ideals / E-186-type modulus child (E-345).
- (round 133) theorist: why clean (4,3) cells never reach a lowering run-of-three containing a (3,3) pair (E-343); skeptic/toolsmith: planted in-class-loss control for the option (E-341). Next ordinary rounds: theorist child = End(T) (H-023), scholar H-025 object-level check.
- (round 131) scholar: correct literature/2112.15587 (Thm 6.3(3) is proved, not conditional); toolsmith: make `isMutable`/`isTiltingAtVertex` refuse on a cycle upstream of i; skeptic: non-monomial cyclic Hom engine; theorist/toolsmith: a non-thin J != 0 test of `refuseNonTilting`.
- (round 130) toolsmith: depth-4 refusal count on the unjoined F-041 pairs; theorist: why dim e_i A e_v <= 1 plus the gate forces J = 0 (maverick's thinness argument); experimentalist: overlap-2 control. (round 129) toolsmith: depth-7 n = 11 13-class, pure-F joins from the crossing children to the 59-class, E-149's 154 lost edges with `sameIdeal`, H-023 End(T)-presentation builder; skeptic: a dormant-thread note (T4/T6) next round; theorist: amend the corollary (rests on H-023).
- (round 127) skeptic/toolsmith: decide whether the gate should refuse when a cycle lies upstream of i (E-329 (d)); experimentalist: an interior gap-2 control (E-330); theorist: is the F step the silting mutation (maverick's reframing)?
- scholar: Happel-Seidel / BGS / Nikulin 1979 / Crawley-Boevey 1992 / Ringel 1994 primaries (not on arXiv: need PDFs; UNVERIFIED stays); owed literature file for 1310.1907.
- (round 126) toolsmith: a gate-admitted loop/cyclic parent for T5; theorist: check the Anick facts in E-326 (H3)(ii); skeptic: decide whether the gate should call `isTiltingAtVertex` by default (E-327).
- (round 125) theorist: is X_3 a function of (key, SNF(S)) on LNAs (a lemma, E-324)?; toolsmith: n = 12 other roots, pure-F route to a crossing; skeptic: decide default for refuseNonTilting with the cost of the 19 F-010 steps.
- maverick: E-269 for singular S (2 of the 6 n = 12 low-lam groups); planted controls on groups 2, 4; a pair inequivalent by construction (E-291).
- toolsmith: key-refusal control (a J != 0 step) for T10; theorist: loop models M2, M4, M5. Done in 109: intermediate End of E-287 paths are not LNAs, none gate-refused (E-289). Bruning generation and silting-order check (scholar); n = 7..9 depth to first hit (toolsmith); rebuild the E-239 depth-5 edge (E-282).
- theorist: E-235 step 4 against Lemma A's hypotheses (E-256 a); one-step rule for the width-5 class-B rules (E-266).
- maverick: E-269's bound as an exact integer bound (closed form now known, E-276); power control on certified-equivalent pairs.
- experimentalist: a gap-2 shape that can lower overlap for T7 (E-288); gate check on rows with a length-3 relation via Ladkani 1001.4765 (E-280).
- toolsmith: independent re-verification of the 78 of 78 A-side corruptions (E-306; needs a pickle rebuild, ~40 min on 4 cores).
- theorist/skeptic: is A_n reachable from B_n by any tilting mutation outside the E-187 moves (T10 at n = 6)? Cartan-data search on small silting-not-tilting pairs (E-278).

## Backlog carried to cycle 2 (list: cycles/001/backlog.md; rest settled/retired/converted at round 120)
- H-027 general n (tilting intermediates); item 8 (H-025 search); H-029 widths 8-9; E-080 shape of the 8 anomalies; T10 depth-6 ball; second-group planted power for E-269 (H-024 F-candidate).

## Old backlog seeds (done unless listed above)
- `mutationSearchDepthFirst` docstring (l. 322) 'guard keeps a walk in one derived class' is false (E-332); `canonicalKey`/guard docstrings; Cor. 4.6 misprint (E-275); UNVERIFIED citations (Nikulin 1.13.3, Anick, Forre); overnight proposals never decided (E-204, E-206/E-210, E-253, E-260).

## Weakest claims the workshop currently relies on
the pure-F crossing's certificate is SNF(S) (necessity chain E-308/E-272, sourced); "25 of 25 joined" = under the J = 0 premise, generation assumed; the J = 0 premise itself; E-269's float bound and the unsourced necessity step behind every Cartan-level NO; S-1's K-threshold law (three lengths, key level); E-151 counts are capped BFS samples.

## Standing facts
- The 7 + 7 n = 9 F-010 J = 0 join was run once by the dispatcher (E-250: no meet, no power).
- arXiv is reachable (STEERING, standing permissions).

- Round 136 conference: F-058 (H-025 CONFIRMED, scoped); T6 (H-017b key form) and H-027 closed; no R-candidate (T10 option harmlessness needs the depth-7/9 option-on walks, round 137). H-027's open piece (tilting intermediates at general n) is in the backlog. Round 132: F-057, R-018.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 135 | 137 (theorist) |
| theorist | 137 | 135 (skeptic) |
| skeptic | 135 | 137 (toolsmith) |
| scholar | 137 | 135 (maverick) |
| toolsmith | 137 | 131 (maverick) |
| maverick | 135 | 137 (scholar) |

- (round 137) toolsmith: depth-8 sharded T10 run (overnight-sized, not approved), E-316 depth 5; add an illegal-relation refusal to `isTiltingAtVertex`. scholar: owed literature files (1310.1907, 1101.3233, math/0702060, 2606.15595, 2604.25789); KLM 1203.5505 draft is in rounds/137/scholar.md. Next round: T5 (theorist/skeptic), F-058 rebuild (skeptic/experimentalist), F-057/F-056 rebuilds.
- (round 136) all: phase 5 sweep, see rounds/136/proceedings.md agenda 1-5; scholar: check the 2112.15587 summary at l. 127, 235. theorist: M3/M3b (two loops); is H-010's end-only argument enough to exclude gap-2 lowering (E-313)? skeptic: which of the 8 minimal anomalies have the E-080 shape (E-312); a second invariant on the E-318 pair (with maverick). toolsmith: J != 0 classifier and the 20-root denominator; Construction-2.2 builder for loop parents; replayable n >= 11 merges (experimentalist lists them). scholar: PDFs wanted (standing): Nikulin 1979, Happel-Seidel, BGS, Anick 1982, Chekanov 2002.
- Round 124 conf.: F-055 (Cartan-level test, scoped), R-017 (key group = one class), H-015 -> R-017, H-024 -> F-055, H-032 narrowed; docstring of `mutationSearchDepthFirst`. Phase 1 history: cycles/002/phase-1.md. Rounds 121-123: E-315 to E-322.
