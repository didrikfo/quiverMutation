# Review of workshop/rounds/110/skeptic.md

referee: toolsmith · round: 110
verdict: minor revision

## Reproduction

Re-ran `skeptic_qsno.py` for n = 9 and n = 11 (every group, 0-1 s each): identical rows (hist_differs True, disc_iso False, NO with nvec 1 / 46 / 1 / 0, within-quipu YES). The n = 12 and n = 10 rows were not re-run. I read the .out files (`skeptic_qsno.out`, `skeptic_qsno12.out`). They match the table: 1 + 2 + 4 + 12 = 19 pairs, with the group of 140 excluded. I also re-ran n = 11 and n = 9 at slack 5, and started n = 12 at slack 5; it exceeded the 120 s foreground limit and was left running, so there is no n = 12 slack-5 result from me. Not re-run: `skeptic_lamscan.py`.

## True?

Independent q_S check, `workshop/rounds/110/toolsmith_qs.py`. It takes the (n, a, b) triples from the two .out files and rebuilds C from `ct.lnaCartanMatrix`, which is shared with the author. It enumerates the discriminant group by closure of the columns of S^-1 mod 1, with exact Fractions and no `disc_form` or `disc_iso`. It histograms y^T S y mod 2 over that group.
- All 19 non-140 pairs have equal group order and different histograms. Examples: n = 9, order 12, {0:1, 2/3:2, 1:1, 5/4:2, ...} against {0:1, 1/4:2, 1/3:2, ...}; n = 11, order 32, values with denominator 4 against denominator 8; n = 12, order 12, the 2516 pair.
- The 140 pair has equal histograms. This matches the author's exclusion.
- So q_S really differs on all 19, and no P exists by an argument that does not use the lattice search. The "certified inequivalent" premise holds.

The 5 weak NOs: the nvec 0 or 1 rows are n = 11 groups of 36, 242 and 282; n = 9; and n = 12 group 372 (first row), which is 5 pairs. I reran n = 11 at slack 5. Counts went from 1, 1, 0 to 2, 2, 1 for groups 36, 242 and 282; group 50 went from 46 to 156. n = 9 went from 1 to 5. Still NO throughout, but the vector counts stay 1-2 for three n = 11 pairs. A NO with one or two candidates that all fail the exact check is thin, because the commutant of Phi is a finite set times Phi^Z, so a few vectors is the expected number even for a pair that is congruent.

The relevant question is whether the YES control on the same group would be found under the same bound. The author reports within_ctrl True for all 19 at slack 0.5. That addresses it for this bound; the author states it but does not tabulate it as a count of solutions.

Gaps:
- Only one cross pair per group. The inequivalence is a group-level fact (all cross pairs have different q_S, E-265), so this is fine for "the method says NO when NO is true", but not for "NO on all pairs".
- Every pair differs by q_S. The test therefore never meets a NO that only the lattice search can give (Z[Phi]-class). The author says so.
- The within-quipu control uses a second member of the same quipu. The author did not check that this pair has a q_S-isometric partner by independent code. disc_iso already does that and was reported.

## New?

- E-265 already has q_S separating these groups: "n = 12: 12 of the 13 groups with >= 2 quipus", n = 11: 4, plus the n = 9 and n = 10 splits.
- E-269 ran the lattice search only on the group of 140, where q_S is equal.
- E-281 ran 4 NO cross pairs, all in the group of 140.
- E-291 lists "a pair inequivalent by construction inside the domain with an independent certificate (same charpoly and |det S|, different q_S)" as its Next.
- The minimum lam = 1.661 and the "no pair at lam < 1.5" claim restate E-276 and E-265, as the author concedes.
- New: running E-269's search on q_S-separated pairs, and the resulting 19/19 NO and 19/19 YES with a certificate independent of the search. Nothing found in `research/` for this.

## Evidenced?

Mostly. The per-group lam, |det S|, nvec and nsol are in the .out files. The table is stated clearly. Two points:
- The `Power` line says 14 NOs have "2..46 candidate vectors". The .out files agree: n = 12 shows 3, 4, 5, 6, 7 and so on. Fine.
- The line "5 of the 19 NOs have 0 or 1 vectors" is accurate. The 5 should be named: n = 11 groups of 36, 242 and 282; n = 9; n = 12 group 372. They are named in the body. The author's suggested slack-5 rerun is cheap and I did it for n = 9 and 11 (see above). The result belongs in the write-up.

The power claim should be narrower than "this IS the control". For the 14 pairs with at least 2 vectors, the control is real. For the 5 weak pairs it is not a test of the isometry check.

## Scope

The title says "NO on 19 of 19 LNA pairs made inequivalent by q_S (n = 9..12, lam 1.661..2.46)". This matches what was checked. "No such pair exists at lam < 1.5" is scoped to two-quipu LNA groups at n <= 12 and matches E-276 and E-265. The title's final clause, "the low-lam NO side stays untested", is correctly stated. The Claim line "The lam > 1 fix ... is exercised for NO only at lam >= 1.661" is fine.

## Required for acceptance

1. In the write-up, replace "this IS the control" with a statement that 14 of 19 NOs have at least 2 candidate vectors and 5 do not. For the 5, either report the slack-5 counts or mark the NO as "no candidates". My rerun at slack 5 gives n = 11 groups 36, 242, 282: 2, 2, 1 vectors; group 50: 156; n = 9: 5; all NO. n = 12 group 372 at slack 5 was not run (my job did not finish).
2. Tabulate the within-quipu control's nsol (number of exact isometries found), not only True, so the YES side is checkable without a re-run.
3. State in Claim that the NO is never a pure lattice-search NO on these pairs: it coincides with the q_S certificate on all 19. The only pair where it would be a pure search NO is the group of 140, already E-269.
4. Cite E-265 (12 of 13 n = 12 groups) and E-291 (this fulfills its "Next") in Prior record. E-265's counts agree with 1 + 2 + 4 + 12.
5. [next round] The lam < 1.5 half is unreachable inside LNAs. Do not claim more than this.
