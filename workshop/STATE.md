# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 127
cycle: 2 (rounds 121-140; CYCLE.md; plan in cycles/002/plan.md, written 2026-11-07)
phase: 2 (round 4 of 4, the conference)
focus: research; phase 1 has a literature part (scholar 121-123, theorist 121-122, maverick 121); phase 5 rounds 137-139 are a backlog and replication sweep, 140 the closing conference (cycles/002/plan.md, Phase focus). Cycle 1 retrospective: cycles/001/retrospective.md.
next_round_kind: conference

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 2 plan
See cycles/002/plan.md and cycles/002/phase-1.md (what phase 1 did). Agenda for phase 2 = round-124 proceedings (rounds/124/proceedings.md); round 125 did T10 (E-323, E-324) and S-1 (E-325); round 126 did H-023 (E-326), H-030 (E-327), H-027 (E-328). Round 127 did E-329 (T5), E-330 (T7), E-331 (T10); 128 is the phase conference. Phase 5 (137-139) is the backlog and replication sweep, 140 the closing conference.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** (agenda, round 124: 1) · H-015 (REFUTED → R-017), H-032 (narrowed) · E-318, E-323: every key-kept J != 0 step found crosses and needs an R move; E-331: n = 12 F-only depth-5 balls (20 roots) hold none. Backward claim (no LNA reaches the 13 crossing parents by pure F) not promoted: 3 of 397 J = 0 F edges have no R inverse (revision due). F-candidate at scope at conference 128; `refuseNonTilting` off. · 127 · maverick
- **T5 · J = 0 premise** (agenda, round 124: 2) · H-023, H-031 · E-326, E-329: no gate-admitted loop-at-i parent (vacuous); (LF), (H1w), (H3), H(I) = 0 hold on 12 136 admitted steps from cyclic monomial parents, lemma D(i)-sourced letters = D(i) model. Library gate and `isTiltingAtVertex` err when a cycle is upstream of i (gate not fixed; warning added). Open: Hom/End(T), non-monomial, A^op presentation-free (C), Anick facts (H3)(ii). · 127 · toolsmith
- **T3/T8 · Cartan-level invariants** (H-025 in the backlog) · H-024 (SUPPORTED, F-055 at scope; also H-025) · Necessity now has a second source: Ladkani 1001.4765 Lemma 3.5 + Rickard + Aihara-Iyama Thm 2.27, P = r^-T (E-308); Nikulin 1.13.3 still UNVERIFIED (secondary only). E-292: NO on 19 of 19 q_S-separated pairs, YES on 19 controls. Open: lam < 1.5 NO (no LNA pair exists there, E-276), Z[Phi]-class route, H-025 periodic-not-FCY search. · 117 (scholar) · skeptic E-314: planted NO power at one n = 11 group (21 of 21) and singular n = 12 group of 47 splits 14 + 33 by SNF(S) (E-291 amended).
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). E-330: E-313's ValueError was a `_copy` bug (fixed); no interior start-gap-2 pair lowers overlap in 3 non-vacuous cells (k <= 5), but no interior control exists (no power); end placement lowers 4 of 5. Open: an interior control, capped k = 5 cells, theorist argument that end-only reasoning excludes gap-2 lowering. · 127 · experimentalist
- **T4 · H-020 rule table** (converted, 118) · H-020, H-029 · H1 holds to widths 9-11, lengths <= 14 (E-225). Class-B first-step delta is a per-rule value, 82 of 82 length transfers at widths 5-6 (E-310), no formula; test is widths 8-9 (H-029). · 118 · skeptic
- **H-027 · A_n vs B_n reachability** · H-027 (SUPPORTED) · E-328: every step of the E-287/E-289 paths at n = 5-7 is an irreducible silting mutation (any-slot certificate), so the paths are tilting-mutation chains; one secondary input (D^b(B_n) ~ kD_n, Happel-Seidel). Scholar proposes closing; ledger decision at conference 128 (F-candidate at n = 5-7). General-n generation still rests on Hubery-Krause. · 126 · scholar
- **S-1 · vertex deletion** (agenda, round 120: 4) · H-026 (SUPPORTED, odd n only) · K0 = min(h, n-4-h) for lone-3 classes at key level (ten odd-n fits, E-301, E-304); even-n fails at the self-mirror class (E-307). E-325: one orbit carries both image keys at K0 (n = 13 h = 2, 3; n = 15 h = 3); key changes only through K = 0 bridge components; K >= 1 constancy and E-307 (3) automatic. Class-level form = the bridge question. F-blockers: K = 0 transit, tail end, non-lone law, self-mirror derivation. · 125 · experimentalist
- **T6 · H-017** · H-017 · Defect not constant per key class (E-283); key shift weak, confounded with source LNA at n = 9 (E-303). Closed for cycle 1; cycle-2 form: H-017b (amendment text in H-017, E-314): per-LNA statistic, delta* fixed first, n = 10/11, overnight proposal, not run. · 119 (text) · maverick
- **T11 · presentation dependence** · H-030 (narrowed) · E-327: of the 8 over-admissions 5 have exactly E-080's shape, 2 E-080 plus an extra path, 1 is E-112's nn circuit; the injectivity test (= tiltingPlus, `isTiltingAtVertex`) refuses all 8 and agrees with Hom on 298 337 of 298 337 random cases n = 5-8. Left: decision whether the gate should use it by default (E-080 option, phase 5); docstring note. · 126 · skeptic
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
- maverick: answer rounds/127/maverick.review.md (forward claim promoted as E-331; backward claim needs forward F-only balls from the parents' SNF class and iso-tests, an exact edge-inversion count, and the 3 non-invertible edges recorded).

## Requests between personas
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

- Round 120 conference: no F/R entry, no status change; F-candidate H-024 (lam >= 1.5, blocker: planted power at a second group), R-candidate "key group = one derived class" (blocker: SNF on C of the 47-group). Docstring sentence added to mutateAtVertex. Overnight proposals E-204, E-206/E-210, E-253, E-260 retired.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 127 | 122 (theorist) |
| theorist | 126 | 127 (toolsmith) |
| skeptic | 126 | 127 (experimentalist, maverick) |
| scholar | 126 | 123 (theorist) |
| toolsmith | 127 | 123 (skeptic) |
| maverick | 127 | 125 (experimentalist) |

- theorist: M3/M3b (two loops); is H-010's end-only argument enough to exclude gap-2 lowering (E-313)? skeptic: which of the 8 minimal anomalies have the E-080 shape (E-312); a second invariant on the E-318 pair (with maverick). toolsmith: J != 0 classifier and the 20-root denominator; Construction-2.2 builder for loop parents; replayable n >= 11 merges (experimentalist lists them). scholar: PDFs wanted (standing): Nikulin 1979, Happel-Seidel, BGS, Anick 1982, Chekanov 2002.
- Round 124 conference: F-055 (Cartan-level test, scoped), R-017 (key group = one class), H-015 -> R-017, H-024 -> F-055, H-032 narrowed; docstring of `mutationSearchDepthFirst`. Phase 1 history: cycles/002/phase-1.md. Rounds 121-123: E-315 to E-322.
