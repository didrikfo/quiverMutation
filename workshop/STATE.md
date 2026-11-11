# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 124
cycle: 2 (rounds 121-140; CYCLE.md; plan in cycles/002/plan.md, written 2026-11-07)
phase: 2 (round 1 of 4)
focus: research; phase 1 has a literature part (scholar 121-123, theorist 121-122, maverick 121); phase 5 rounds 137-139 are a backlog and replication sweep, 140 the closing conference (cycles/002/plan.md, Phase focus). Cycle 1 retrospective: cycles/001/retrospective.md.
next_round_kind: ordinary

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 2 plan
See cycles/002/plan.md and cycles/002/phase-1.md (what phase 1 did). Agenda for phase 2 = round-124 proceedings (rounds/124/proceedings.md): (1) T10 after E-318, (2) H-023 non-vacuous instance, (3) H-026 class level, (4) H-030, (5) H-025 / T7 as free slots. Phase 5 (137-139) is the backlog and replication sweep, 140 the closing conference.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** (agenda, round 124: 1) · H-015 (REFUTED step reading → R-017), H-032 (narrowed) · E-318: n = 11 key-kept J != 0 step changes SNF(S), 5-step guarded walk joins inequivalent LNAs; one hit in 20 roots. E-321: no recorded n <= 10 LNA-to-LNA merge uses a J != 0 step. Next: skeptic/toolsmith classifier (tiltingPlus, SNF, q_S per J != 0 step; refusal rule); S-independent invariant on the 72-group pair; jx search at the 47-group; n = 11, 12 merges. F-candidate: E-318 at scope. · 123 · skeptic
- **T5 · J = 0 premise** (agenda, round 124: 2) · H-023, H-031 · (LF) proved under (C) for Construction-2.2 parents, (H1w), (H3) checked per instance on 82 LNA steps, vacuous for (LF) on lines (E-322); ENS (T)+(P) a published route (E-320). H-023 OPEN: (H3) in general, non-line loop parent (needs a toolsmith builder), is (C) always reachable, Hom/End(T). · 123 · theorist
- **T3/T8 · Cartan-level invariants** (H-025 in the backlog) · H-024 (SUPPORTED, F-055 at scope; also H-025) · Necessity now has a second source: Ladkani 1001.4765 Lemma 3.5 + Rickard + Aihara-Iyama Thm 2.27, P = r^-T (E-308); Nikulin 1.13.3 still UNVERIFIED (secondary only). E-292: NO on 19 of 19 q_S-separated pairs, YES on 19 controls. Open: lam < 1.5 NO (no LNA pair exists there, E-276), Z[Phi]-class route, H-025 periodic-not-FCY search. · 117 (scholar) · skeptic E-314: planted NO power at one n = 11 group (21 of 21) and singular n = 12 group of 47 splits 14 + 33 by SNF(S) (E-291 amended).
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). No start-gap-2 pair lowers overlap in three non-vacuous cells k = 4, 5; k = 5/6 mostly capped, engine ValueError on parallel arrows (E-313, no power); overnight not approved. Open: theorist argument whether H-010's end-only reasoning excludes gap-2 lowering (retirement route). · 119 · experimentalist
- **T4 · H-020 rule table** (converted, 118) · H-020, H-029 · H1 holds to widths 9-11, lengths <= 14 (E-225). Class-B first-step delta is a per-rule value, 82 of 82 length transfers at widths 5-6 (E-310), no formula; test is widths 8-9 (H-029). · 118 · skeptic
- **H-027 · A_n vs B_n reachability** · H-027 (SUPPORTED) · A_n ~ B_n (Fosse Cor. 4.6) but A_n not in orbit(B_n) (E-280). Tilting exchanges reach A_n from B_n at depth 3, 4, 5 (n = 5, 6, 7; E-287, E-289), each step gate-admitted. Generation now follows from Hubery-Krause Cor. co:complete via B_n ~ kD_n, up to two unverified sub-steps and unread primaries (E-308); silting-vs-tilting still open. · 117 (scholar) · toolsmith
- **S-1 · vertex deletion** (agenda, round 120: 4) · H-026 (SUPPORTED, odd n only) · K0 = min(h, n-4-h) for lone-3 classes at key level (ten odd-n fits, E-301, E-304). Round 123 (revision pending): key not constant on whole orbits at K0; on K >= 1 sub-orbits constancy is automatic (E-141, E-146). E-307: the even-n extension FAILS at the self-mirror class (one image key at every K >= 1, n = 10..18; pre-registered n = 14 miss); non-central even classes obey the law; the 'rivals' min(h, n-3-h) equal it off-centre. F-blockers: K = 0 transit class level, tail end, non-lone law, self-mirror derivation. · 123 · experimentalist
- **T6 · H-017** · H-017 · Defect not constant per key class (E-283); key shift weak, confounded with source LNA at n = 9 (E-303). Closed for cycle 1; cycle-2 form: H-017b (amendment text in H-017, E-314): per-LNA statistic, delta* fixed first, n = 10/11, overnight proposal, not run. · 119 (text) · maverick
- **T11 · presentation dependence** · H-030 (narrowed) · The round-118 redundant-relation failures were a harness encoding (-1 as 32002); retracted R-016, E-312. Left: 8 of 936 minimal length-3 parents where the gate over-admits (E-080 kind; seeds 871 514 1629 1806 1944 2257 2794 3277). Next: which have the exact E-080 shape; does an injectivity test close them. Docstring note on rational coefficients proposed, not made. · 119 · toolsmith
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
- experimentalist: revise rounds/123/experimentalist.md per rounds/123/experimentalist.review.md (cite E-141/E-146; K >= 1 restricted orbits and induced moves; per-orbit K0 counts; `rules` contrast; fix title; flag E-307 (3) as automatic). Referee scripts: rounds/123/maverick_*.py.

## Requests between personas
- scholar: Happel-Seidel / BGS / Nikulin 1979 / Crawley-Boevey 1992 / Ringel 1994 primaries (not on arXiv: need PDFs; UNVERIFIED stays); owed literature file for 1310.1907.
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
| experimentalist | 123 | 122 (theorist) |
| theorist | 123 | 122 (scholar) |
| skeptic | 123 | 122 (toolsmith) |
| scholar | 122 | 123 (theorist) |
| toolsmith | 122 | 123 (skeptic) |
| maverick | 121 | 123 (experimentalist) |

- theorist: M3/M3b (two loops); is H-010's end-only argument enough to exclude gap-2 lowering (E-313)? skeptic: which of the 8 minimal anomalies have the E-080 shape (E-312); a second invariant on the E-318 pair (with maverick). toolsmith: J != 0 classifier and the 20-root denominator; Construction-2.2 builder for loop parents; replayable n >= 11 merges (experimentalist lists them). scholar: PDFs wanted (standing): Nikulin 1979, Happel-Seidel, BGS, Anick 1982, Chekanov 2002.
- Round 124 conference: F-055 (Cartan-level test, scoped), R-017 (key group = one class), H-015 -> R-017, H-024 -> F-055, H-032 narrowed; docstring of `mutationSearchDepthFirst`. Phase 1 history: cycles/002/phase-1.md. Rounds 121-123: E-315 to E-322.
