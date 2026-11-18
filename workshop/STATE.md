# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 131
cycle: 2 (rounds 121-140; CYCLE.md; plan in cycles/002/plan.md, written 2026-11-07)
phase: 3 (round 4 of 4)
focus: research; phase 1 has a literature part (scholar 121-123, theorist 121-122, maverick 121); phase 5 rounds 137-139 are a backlog and replication sweep, 140 the closing conference (cycles/002/plan.md, Phase focus). Cycle 1 retrospective: cycles/001/retrospective.md.
next_round_kind: conference

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 2 plan
See cycles/002/plan.md and cycles/002/phase-1.md (what phase 1 did). Agenda for phase 3 = round-128 proceedings (rounds/128/proceedings.md); phase 2 agenda was round 124's; round 125 did T10 (E-323, E-324) and S-1 (E-325); round 126 did H-023 (E-326), H-030 (E-327), H-027 (E-328). Round 131: E-338 (thin => J = 0), E-339 (H-025 at n = 15), E-340 (gate vs Hom on cyclic parents). Round 129 (T10): E-332 pure-F crossing at n = 11, E-333 Hom check of E-323's 50 steps, E-334 Hom + invariants at the E-323 site. Phase 5 (137-139) is the backlog and replication sweep, 140 the closing conference.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** (agenda 1+2) · H-015 (REFUTED → R-017), H-032 (narrowed) · E-332: pure-F key-guarded crossing at n = 11 (SNF (1,2,4)→(1,2,2), J != 0 step); `refuseNonTilting=True` removes it. E-333/E-334: J != 0 ⇔ non-tilting silting, Hom-checked. E-337: option on loses no F-041 join at n = 10, 11 depth 3 (balls thin: no step with dim e_i A e_v >= 2); default stays off (decided 130). E-338: thin balls (dim e_i A e_v <= 1) cannot hold a J != 0 step, so E-337 is a theorem, not a test; a non-thin test is open. F-candidate E-323/E-332 (conf. 132). · 131 · theorist
- **T5 · J = 0 premise** (agenda, round 128: 4) · H-023, H-031 · E-326, E-329: no gate-admitted loop-at-i parent (vacuous); (LF), (H1w), (H3), H(I) = 0 hold on 12 136 admitted steps from cyclic monomial parents, lemma D(i)-sourced letters = D(i) model. Library gate and `isTiltingAtVertex` err on cyclic parents (E-340: 39 over-admissions in 18 145 rows, n = 5-7; no LNA walk reaches a cyclic node at n <= 8 depth <= 2-5); R-candidate for conf. 132. Open: non-monomial cyclic Hom engine, refuse-not-warn fix, Anick (H3)(ii). · 131 · maverick
- **T3/T8 · Cartan-level invariants** (H-025 periodic-not-FCY) · H-024 (SUPPORTED, F-055), H-025 (SUPPORTED) · E-339 (131): kA_15/rad^3, /rad^8 not FCY via Fuchsian R(2,3,10) (unpublished argument, referee-checked); no blind search. Necessity: Ladkani 1001.4765 + Rickard + AI Thm 2.27 (E-308); Nikulin 1.13.3 UNVERIFIED. E-292 NO on 19 of 19 pairs, E-314 power. Open: lam < 1.5 NO, Z[Phi] route, n = 10-12 groups' type. · 131 · scholar
- **T7 · H-010** (agenda, round 128: 5) · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). E-336: interior gap-2 control exists for (8,5)(10,4) (bystander (13,2) or (12,3) lowers 3 -> 2), so E-330's null has power there; none for (8,4)(10,3) (0 of 91 + 1242 screened, double mutations only); capped k = 5 cell overnight-sized. Open: overlap-2 control or argument, theorist end-only argument. · 130 · experimentalist
- **T4 · H-020 rule table** CLOSED round 130 · H-020, H-029 (SUPPORTED) · E-335: per-rule first-step delta, 106 of 106 class-B rules at widths 8-9 (lengths <= w+4), widths 5-7 <= w+3; cross-width falsifier no power; locality unproved. Close in conference-132 ledger. · 130 · skeptic
- **H-027 · A_n vs B_n reachability** · H-027 (CONFIRMED, F-056) · E-328: every step of the E-287/E-289 paths at n = 5-7 is an irreducible silting mutation (any-slot certificate), so the paths are tilting-mutation chains; one secondary input (D^b(B_n) ~ kD_n, Happel-Seidel). CONFIRMED -> F-056 at n = 5-7 (round 128). General-n generation still rests on Hubery-Krause. · 126 · scholar
- **S-1 · vertex deletion** (agenda, round 120: 4) · H-026 (SUPPORTED, odd n only) · K0 = min(h, n-4-h) for lone-3 classes at key level (ten odd-n fits, E-301, E-304); even-n fails at the self-mirror class (E-307). E-325: one orbit carries both image keys at K0 (n = 13 h = 2, 3; n = 15 h = 3); key changes only through K = 0 bridge components; K >= 1 constancy and E-307 (3) automatic. Class-level form = the bridge question. F-blockers: K = 0 transit, tail end, non-lone law, self-mirror derivation. · 125 · experimentalist
- **T6 · H-017** · H-017 · Defect not constant per key class (E-283); key shift weak, confounded with source LNA at n = 9 (E-303). Closed for cycle 1; cycle-2 form: H-017b (amendment text in H-017, E-314): per-LNA statistic, delta* fixed first, n = 10/11, overnight proposal, not run. · 119 (text) · maverick
- **T11 · presentation dependence** · H-030 (narrowed) · E-327: of the 8 over-admissions 5 have exactly E-080's shape, 2 E-080 plus an extra path, 1 is E-112's nn circuit; the injectivity test (= tiltingPlus, `isTiltingAtVertex`) refuses all 8 and agrees with Hom on 298 337 of 298 337 random cases n = 5-8. Left: decision whether the gate should use it by default (E-080 option, phase 5); docstring note. · 126 · skeptic
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
(none)

## Requests between personas
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
- item 8 (H-025 search); H-029 widths 8-9; E-080 shape of the 8 anomalies; T10 depth-6 ball; second-group planted power for E-269 (H-024 F-candidate).

## Old backlog seeds (done unless listed above)
- `mutationSearchDepthFirst` docstring (l. 322) 'guard keeps a walk in one derived class' is false (E-332); `canonicalKey`/guard docstrings; Cor. 4.6 misprint (E-275); UNVERIFIED citations (Nikulin 1.13.3, Anick, Forre); overnight proposals never decided (E-204, E-206/E-210, E-253, E-260).

## Weakest claims the workshop currently relies on
the pure-F crossing's certificate is SNF(S) (necessity chain E-308/E-272, sourced); "25 of 25 joined" = under the J = 0 premise, generation assumed; the J = 0 premise itself; E-269's float bound and the unsourced necessity step behind every Cartan-level NO; S-1's K-threshold law (three lengths, key level); E-151 counts are capped BFS samples.

## Standing facts
- The 7 + 7 n = 9 F-010 J = 0 join was run once by the dispatcher (E-250: no meet, no power).
- arXiv is reachable (STEERING, standing permissions).

- Round 128 conference: F-056 (H-027, n = 5-7); F-candidate E-323 (blocker: non-SNF/Hom check); R-candidate "gate decides tilting" (blocker: gate vs Hom on cycle-upstream parents). Four of six personas ask: is J != 0 exactly non-tilting silting?

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 130 | 131 (theorist) |
| theorist | 131 | 130 (experimentalist) |
| skeptic | 130 | 131 (scholar) |
| scholar | 131 | 130 (skeptic) |
| toolsmith | 130 | 131 (maverick) |
| maverick | 131 | 130 (toolsmith) |

- theorist: M3/M3b (two loops); is H-010's end-only argument enough to exclude gap-2 lowering (E-313)? skeptic: which of the 8 minimal anomalies have the E-080 shape (E-312); a second invariant on the E-318 pair (with maverick). toolsmith: J != 0 classifier and the 20-root denominator; Construction-2.2 builder for loop parents; replayable n >= 11 merges (experimentalist lists them). scholar: PDFs wanted (standing): Nikulin 1979, Happel-Seidel, BGS, Anick 1982, Chekanov 2002.
- Round 124 conf.: F-055 (Cartan-level test, scoped), R-017 (key group = one class), H-015 -> R-017, H-024 -> F-055, H-032 narrowed; docstring of `mutationSearchDepthFirst`. Phase 1 history: cycles/002/phase-1.md. Rounds 121-123: E-315 to E-322.
