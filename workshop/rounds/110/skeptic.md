# E-269's Z-congruence search says NO on 19 of 19 LNA pairs made inequivalent by q_S, each NO coinciding with the q_S certificate (never a pure lattice-search NO; n = 9..12, lam 1.661..2.46); no such pair exists at lam < 1.5, so the low-lam NO side stays untested

author: skeptic · round: 110 · kind: result
thread: T3/T8 · bears on: H-024, E-269, E-265, E-276, E-281, E-291
scope: n = 9..12, slack 0.5 (all 19) and slack 5 (all 19, also run), every lnaKeys group with >= 2 quipus and S of signature (n-1,1) (1, 2, 4, 12 groups; the 13th n = 12 group is the group of 140, whose q_S agree, E-269's own pair); ONE cross-quipu pair per group; the test is `maverick_iso140.test` unchanged except `N` set per run. All lam >= 1.661. Not covered: lam < 1.5 (see Claim), the Z[Phi]-module-class route (ii), more than one pair per group, other slacks, groups without two quipus.
power: 14 of the 19 NOs have >= 2 candidate vectors at slack 0.5 (2..46), all failing the exact isometry check: there the NO is a tested NO. 5 have 0 or 1 vectors (n = 11 groups 36, 242, 282; n = 9; n = 12 group 372 first row): there the NO is weak ("nothing to test") at slack 0.5; at slack 5 they have 2, 2, 1, 5, 3 vectors, still NO. The YES side finds exact isometries on all 19 within-quipu controls (nsol 1..4 at slack 0.5, table below). No lam < 1.5 case.

## Response to referee (round 110, toolsmith: minor revision)

1. "this IS the control" replaced (power line): 14 of 19 NOs have >= 2 vectors, 5 do not. Slack-5 counts for the 5 weak NOs, run by me (`skeptic_ctrl.py`, `skeptic_ctrl{9,10,11,12}.out`), agree with the referee's n = 9 and n = 11 numbers: n = 11 groups 36 / 242 / 282: 1,1,0 -> 2,2,1 vectors; n = 9: 1 -> 5; n = 12 group 372 (lam 2.0149, |det S| 39): 0 -> 3 (ran in ~1 s after the 2 min group build, so it fit). All 0 isometries, all NO. Still thin for three n = 11 pairs (1-2 vectors); the strong evidence there is the q_S certificate, not the search.
2. Within-quipu control, nsol (exact isometries found), slack 0.5 [nvec/nsol], and at slack 5; cross-pair nvec in the last two columns:

| n | group size | ctrl 0.5 | ctrl 5 | cross vec 0.5 | cross vec 5 |
|---|---|---|---|---|---|
| 9 | 36 | 5/4 | 9/8 | 1 | 5 |
| 10 | 45 | 9/4 | 11/6 | 2 | 8 |
| 10 | 176 | 7/4 | 13/6 | 3 | 5 |
| 11 | 36 | 5/4 | 6/4 | 1 | 2 |
| 11 | 50 | 37/2 | 129/6 | 46 | 156 |
| 11 | 242 | 5/4 | 8/6 | 1 | 2 |
| 11 | 282 | 1/1 | 4/3 | 0 | 1 |
| 12 | 37 | 6/2 | 14/6 | 4 | 11 |
| 12 | 54 (2.3286) | 11/4 | 15/4 | 7 | 11 |
| 12 | 54 (2.3575) | 8/2 | 25/6 | 6 | 18 |
| 12 | 148 | 5/1 | 9/3 | 4 | 8 |
| 12 | 178 | 9/1 | 20/3 | 5 | 7 |
| 12 | 298 | 4/1 | 10/3 | 3 | 6 |
| 12 | 372 (2.0149) | 1/1 | 3/3 | 0 | 3 |
| 12 | 372 (2.0582) | 7/2 | 12/6 | 3 | 3 |
| 12 | 550 | 9/1 | 24/3 | 5 | 7 |
| 12 | 754 | 4/1 | 12/3 | 3 | 3 |
| 12 | 993 | 7/3 | 21/7 | 4 | 8 |
| 12 | 2516 (1.661) | 4/1 | 12/3 | 3 | 6 |

Every control has nsol >= 1 at both slacks (19/19 YES), every cross pair nsol = 0 at both (19/19 NO). The weak NO rows have a control with nsol >= 1 under the same bound (n = 11 group 282 and n = 12 group 372: only 1 isometry at slack 0.5, 3-4 at slack 5), so the bound does reach a congruent partner there. The group of 140 is excluded from the 19 (cross vec 12 / 21, ctrl 14/2, 31/6).
3. Claim now says: on all 19 pairs the NO coincides with the q_S certificate (hist_differs and disc_iso False on every one), so it is never a pure lattice-search NO here; the only candidate for a pure search NO is the group of 140 (E-269).
4. Prior record now cites E-265 (12 of 13 n = 12 groups with >= 2 quipus separated by q_S; 1 + 2 + 4 + 12 = 19) and E-291 (this fulfils its "Next": an inequivalent-by-construction pair with an independent certificate).
5. Item 5 (next round): not done; lam < 1.5 stays unreachable inside LNAs, nothing more claimed.
Not done: toolsmith's other gap (one cross pair per group) stays; group-level fact by E-265.

## Claim

Blocking check for F-candidate H-024, answered in two halves. (1) Where such pairs exist, the test has power, with the caveat that on these 19 pairs the NO coincides with the q_S certificate on every one, so it is never a pure lattice-search NO (the one candidate is the group of 140, E-269): for each LNA key group with two quipus (n = 9..12), a pair whose discriminant forms q_S differ (so P C_a P^T = C_b is impossible, by an argument that never uses the lattice search) is rejected by E-269's search, 19/19, and a same-quipu pair in the same group is accepted with an exact P, 19/19. A single YES on a q_S-separated pair would have been a method bug; there was none. (2) The assignment's regime lam < 1.5 cannot be reached with LNA pairs: in n = 8..12 no two-quipu group has lam < 1.661 (counts 0, 1, 2, 4, 13 groups, minima 1.7816, 1.8098, 1.8238, 1.661), which is E-276's census redone; the low-lam test therefore has no certified-NO control inside LNAs, and the lam > 1 fix to the Salem threshold is exercised for NO only at lam >= 1.661. Not claimed: that the method says NO at lam < 1.5 (the code path is the same, the bound formula n-2+lam+1/lam is the same, but not exercised); that the test is sensitive to a Z[Phi]-class difference with q_S equal (route (ii): untouched; this is the only way a NO could be a pure lattice-search NO, and the only such pair known is the group of 140, E-269). Refuted by: a q_S-separated pair at any n <= 12 with a YES, or a two-quipu LNA group with lam < 1.5.

## Evidence

Construction: groups with the same coxeter key (lnaKeys) have equal Phi charpoly, hence equal lam and |det S| (S = C(1-Phi^-1)); I take one almost-separate LNA from each of two different quipus (quipuForAlmostSeparateLNA). Filter, checked in the script: det C = 1, sig(S) = (n-1,1) for both, det S equal, charpoly(Phi) equal. Independent certificate of inequivalence (own code, not E-269's): histogram of x^T S x mod 2 over the discriminant group differs for all 19 (`hist_differs=True`), and maverick's `disc_iso` (finite-quadratic-form isometry search) also returns False for all 19. q_S is invariant under S -> P^T S P, so no unimodular P exists.

| n | groups run | lam range | NO / run | within-quipu YES |
|---|---|---|---|---|
| 9 | 1 | 1.7816 | 1/1 | 1/1 |
| 10 | 2 | 1.8098, 2.0306 | 2/2 | 2/2 |
| 11 | 4 | 1.8238..2.3193 | 4/4 | 4/4 |
| 12 | 12 | 1.661..2.4596 | 12/12 | 12/12 |

Per-group rows (lam, |det S|, rank of commutant, nvec, nsol) are in `skeptic_qsno.out` (n = 9..12 with the old 400-member cut) and `skeptic_qsno12.out` (n = 12, all members; adds the lam 1.661 group of 2516). Smallest-lam case: n = 12, group 2516, lam 1.6610, |det S| = 12, a = (2,2,3,0,2,2,2,3,0,0), b = (2,4,0,0,2,2,2,2,3,0): hist differs, 3 lattice vectors, 0 isometries, NO; within-quipu control YES. nvec = 0 or 1 rows: n = 11 groups of 36, 242, 282; n = 9; n = 12 group 372 (first row, lam 2.0149). The group-of-140 row in the .out files (hist equal, NO, 12 vectors) is E-269's pair, not a q_S control, and is not counted in the 19.

lam < 1.5: `skeptic_lamscan.py` over all members, n = 8..12: groups with >= 2 quipus 0, 1, 2, 4, 13; minimum lam 1.7816 (n = 9), 1.8098, 1.8238, 1.661. Matches E-265's counts (1, 2, 4, 13) and E-276's 1.661. Note `skeptic_qsno.py` initially cut each group at 400 members and missed the 2516 group (quipu 2 appears later); fixed, all members now used.

What a synthetic lam < 1.5 pair would need, and why I did not build it: a general integer C (det 1) in the Salem class with lam < 1.5 and a q_S different from a sibling's with the same Phi charpoly. `maverick_planted*.py` showed that +-1 entry perturbations stay in the congruence class or leave the domain; building two non-congruent unimodular C with one charpoly and different q_S means choosing a different genus of the commutant lattice by hand, not a matrix perturbation. Needs a number-field computation (Z[Phi] ideal classes; cypari now installed in the venv) -- proposal, not done.

## Reproduction

```
uv pip install --python .venv/bin/python cypari python-flint     # venv lacked them again
for n in 12 11 10 9; do timeout 10m .venv/bin/python -W ignore workshop/rounds/110/skeptic_qsno.py $n 12; done   # ~3 min total; n = 12 with 13: ~2 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/110/skeptic_lamscan.py                                      # ~2 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/110/skeptic_missing.py                                      # domain filter, which n = 12 groups run
```

## Prior record

E-265 (q_S separates 7 of 8 pairs; 12 of the 13 n = 12 groups with >= 2 quipus, which with 1 + 2 + 4 at n = 9..11 gives my 19; groups with >= 2 quipus 1, 2, 4, 13 at n = 9..12), E-269 (NO on 9 cross pairs of the 140 group, YES within), E-281 (4 cross-pair NO, 11 YES at low lam), E-276 (no certified pair at lam < 1.5; bound closed form), E-291 (its Next, "a pair inequivalent by construction with an independent certificate", is what this does; its "power not demonstrated" becomes "demonstrated at lam >= 1.661, not < 1.5"). Genuinely new: E-269's search had been run as a NO test only on the group of 140 (q_S equal there); I ran it on pairs separated by q_S, so this is the first case where NO is known true by a certificate that does not use the lattice search. Also new: the lam-ordered table showing the lowest lam where it was exercised. The "no certified-inequivalent pair at lam < 1.5" is E-276, rediscovered by lamscan.

## Code changed

None under src/ or tests/. New: workshop/rounds/110/skeptic_ctrl.py (+ skeptic_ctrl{9..12}.out), skeptic_qsno.py, skeptic_lamscan.py, skeptic_missing.py (+ .out files, < 10 KB). `iso.N` is set at runtime (the module-level N = 12 is what test() reads).

## Next

(a) experimentalist/maverick: the H-024 NO side is now powered at lam >= 1.661; mark E-291's "power not demonstrated" as "demonstrated at lam >= 1.661, not at < 1.5". (b) theorist: construct, outside LNAs, two det-1 Cartan-type forms with a Salem lam < 1.5, same charpoly, different q_S (genus choice in the commutant lattice); or compute Z[Phi] class numbers for the n = 12 lam 1.4886 charpoly and look for a Z[Phi]-class pair. (c) done in this round: slack 5 gives 1-5 vectors for the 5 weak NOs; slack 30 not run.
