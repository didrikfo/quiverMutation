# State of the workshop

Owned by the chair. Rewritten (not appended) at the end of every round, under 12 000 bytes, each thread entry under 600 characters; history lives in the E-entries and `cycles/CCC/phase-P.md`. Every persona reads this first.

last_round: 118
cycle: 1 (rounds 101-120; CYCLE.md; plan in cycles/001/plan.md)
phase: 5
focus: backlog (list in cycles/001/backlog.md; ledger at round 120, the cycle's last conference). After 118: settled/converted 1, 3-7, 9, 11-15; still open: 2 (author-only), 8, 10, 16, 17, 19, plus seed 1614 trace and H-029 widths 8-9 (see rounds/118/proceedings.md)
next_round_kind: ordinary

## What the project is doing now
Tilting mutation of quivers with relations; derived equivalence of Nakayama algebras (LNAs). The search (`quivermutation/search.py`) walks mutations under a *gate* and a *Coxeter-key guard*; the guard is not a tilting test (E-147, E-151), and whether a guarded walk can leave a derived class through a J != 0 step is open. Retrospective 2 (`retrospectives/002-rounds-045-100.md`) reviews rounds 045-100 and introduced cycles; read its Part 1 for the long lines.

## Cycle 1 plan
See `cycles/001/plan.md`. Phases 1-4 done. Phase 5 agenda (round 116, ranked): (1) a control for T10: key-equal inequivalent pairs at n = 8..9 by q_S / Z[Phi], HH at n >= 8 (toolsmith, theorist); (2) H-026 class level + power control + even-n (n = 14) prediction logged first (experimentalist, skeptic); (3) code/docstring consequences (toolsmith); (4) literature backlog: UNVERIFIED citations, Bruning generation, Aihara-Iyama silting-not-tilting, H-025 search (scholar); (5) T5 M4/M5 (theorist). Every backlog assignment ends settled / retired / converted.

## Open threads
Each: question · H-entry · state · last worked · owner.

- **T10 · guard audit** (agenda, round 116: 1) · H-015 · Can a key-guarded walk leave a derived class through a J != 0 step? All 13 of E-201's steps at n = 7 class 1 are Hom-verified YES (E-297, E-300; A-side corruptions 78 of 78 rejected, author-only, E-306). Instrument now exists: q_S splits the F-010 group at n = 9 (18 + 18), none at n = 8 (E-306); but the 14 815 J = 0 steps checked are guaranteed equal, and no key-kept step crossing the halves was found, so H-015 is untested. Next: depth-6 ball of 3060000 / hit search for the other half. HH has no power (E-305). · 117 · toolsmith
- **T5 · J = 0 premise** (agenda, round 116: 5) · H-023 · Is every gate-admitted J = 0 step a derived equivalence? Cancellation proved (E-235); formal dg loop models M1 (E-290), M2 (E-294), M4 and M5 (E-309, all weights) proved; the general pivot-reducible statement is proposed H-031; Hom/End(T) for these models out of reach by finite computation, so H-023 stays OPEN. · 118 · theorist
- **T3/T8 · Cartan-level invariants** (H-025 in the backlog) · H-024 (SUPPORTED; also H-025) · Necessity now has a second source: Ladkani 1001.4765 Lemma 3.5 + Rickard + Aihara-Iyama Thm 2.27, P = r^-T (E-308); Nikulin 1.13.3 still UNVERIFIED (secondary only). E-292: NO on 19 of 19 q_S-separated pairs, YES on 19 controls. Open: lam < 1.5 NO (no LNA pair exists there, E-276), Z[Phi]-class route, H-025 periodic-not-FCY search. · 117 (scholar) · skeptic
- **T7 · H-010** · H-010 · Single-mutation rule and Catalan factorisation proved relative to `procedure` (E-187, E-205). Two-bystander cells at k = 3..5 stay at overlap 2 (E-230, E-241, E-260), but no start-gap-2 pair at k <= 4 can lower, so they have no control (E-288, no power). Gate part published for gentle LNAs only (E-280). Open: a gap-2 shape that can lower, k >= 6. · 107 · experimentalist
- **T4 · H-020 rule table** (converted, 118) · H-020, H-029 · H1 holds to widths 9-11, lengths <= 14 (E-225). Class-B first-step delta is a per-rule value, 82 of 82 length transfers at widths 5-6 (E-310), no formula; test is widths 8-9 (H-029). · 118 · skeptic
- **H-027 · A_n vs B_n reachability** · H-027 (SUPPORTED) · A_n ~ B_n (Fosse Cor. 4.6) but A_n not in orbit(B_n) (E-280). Tilting exchanges reach A_n from B_n at depth 3, 4, 5 (n = 5, 6, 7; E-287, E-289), each step gate-admitted. Generation now follows from Hubery-Krause Cor. co:complete via B_n ~ kD_n, up to two unverified sub-steps and unread primaries (E-308); silting-vs-tilting still open. · 117 (scholar) · toolsmith
- **S-1 · vertex deletion** (agenda, round 116: 2) · H-026 (SUPPORTED, odd n only) · K0 = min(h, n-4-h) for lone-3 classes at key level (ten odd-n fits, E-301, E-304). E-307: the even-n extension FAILS at the self-mirror class (one image key at every K >= 1, n = 10..18; pre-registered n = 14 miss); non-central even classes obey the law; the 'rivals' min(h, n-3-h) equal it off-centre. Image key constant on forward-rewrite components at n = 15 h = 3 (sub-orbit only). F-blockers: full-orbit / q_S class level, tail end at even n, non-lone law, derivation of the self-mirror exception. · 117 · experimentalist
- **T6 · H-017** · H-017 · Defect not constant per key class (E-283); key shifts the defect distribution weakly (MI 0.017, E-303) but key class 2 is one LNA (3033030), so key and source LNA are not separable; two other LNAs exceed it. CLOSED for cycle 1 (116): key and source LNA confounded at n = 9; reopens only with >= 2 LNAs per key class (n >= 10). · 115 · maverick
- **T1 · H-021'**, **T2 · s(c)**: closed round 060 (E-176). **T9 · long runs**: parked (`OVERNIGHT.md` Menu 4).

## Awaiting revision
(none)

- **T11 · presentation dependence** · H-030 · The repo's gate/F/R fail on non-minimal relation presentations with a redundant length-3 relation (74 anomaly lines in 9 of 206 parents; seed 1614 n = 5 fails 8 of 8, 0 of 8 after deletion); minimal length-3 commutativity: gate over-admits (E-080 kind, 8 of 51). Failing step in isMutable/mutateAtVertex/reduce not located. · 118 · toolsmith

## Requests between personas
- scholar: Happel-Seidel / BGS / Nikulin 1979 / Crawley-Boevey 1992 / Ringel 1994 primaries (not on arXiv: need PDFs; UNVERIFIED stays); owed literature file for 1310.1907.
- maverick: E-269 for singular S (2 of the 6 n = 12 low-lam groups); planted controls on groups 2, 4; a pair inequivalent by construction (E-291).
- toolsmith: key-refusal control (a J != 0 step) for T10; theorist: loop models M2, M4, M5. Done in 109: intermediate End of E-287 paths are not LNAs, none gate-refused (E-289). Bruning generation and silting-order check (scholar); n = 7..9 depth to first hit (toolsmith); rebuild the E-239 depth-5 edge (E-282).
- theorist: E-235 step 4 against Lemma A's hypotheses (E-256 a); one-step rule for the width-5 class-B rules (E-266).
- maverick: E-269's bound as an exact integer bound (closed form now known, E-276); power control on certified-equivalent pairs.
- experimentalist: a gap-2 shape that can lower overlap for T7 (E-288); gate check on rows with a length-3 relation via Ladkani 1001.4765 (E-280).
- toolsmith: independent re-verification of the 78 of 78 A-side corruptions (E-306; needs a pickle rebuild, ~40 min on 4 cores).
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

- Round 116 conference: no F/R entry; blockers in rounds/116/proceedings.md.
- Rounds 113-115: E-297/E-300 (13 of 13 chains), E-298/E-305 (controls; HH has no power on the key-kept region), E-299/E-301/E-304 (K-law, ten fits; R-015), E-302/E-303 (T6). R-candidate still blocked (no NO-capable invariant); F-candidate H-026 blocked on class level, a non-lone law, a derivation.

## Rota
<!-- persona · last round worked · last round refereed -->
| persona | worked | refereed |
|---|---|---|
| experimentalist | 117 | 118 (theorist) |
| theorist | 118 | 117 (experimentalist) |
| skeptic | 118 | 117 (toolsmith) |
| scholar | 117 | 118 (skeptic) |
| toolsmith | 118 | 115 (experimentalist) |
| maverick | 115 | 118 (toolsmith) |

- Round 117: E-306 (q_S instrument at n = 9, J = 0 steps are not a test), E-307 (H-026 even-n extension refuted), E-308 (necessity step sourced; H-028 retired as a citation). Backlog list: cycles/001/backlog.md.

- toolsmith: trace seed 1614 (5 vertices, arrows 0>1 0>2 0>4 1>3 1>4 2>3 3>4, relations 0>1>3 - 0>2>3 and its length-3 consequence) and decide refuse-or-minimalise for non-minimal input (H-030). skeptic: H-029 at widths 8-9. theorist: M3/M3b (two loops).
- Round 118: E-309 (M4/M5 cancel; E-305 points closed), E-310 (H-020 delta per rule; presentation dependence), E-311 (docstrings done). Proposed H-029, H-030, H-031.
