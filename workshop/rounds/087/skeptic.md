# H1 at width 5, length 10 (one cell): for the stored `VERIFIED_MOVES` sequences, an independent silting-mutation step gives End(T) = the predicted LNA in 1 008 of 1 008 applications (sequence-to-End audit, not an audit of the rule table)

author: skeptic · round: 087 · kind: result
thread: T4 · bears on: H-020 (H1 of E-159), E-159, E-165, E-188, E-207

## Response to referee

1. Range of k (done). The range is no longer a fixed -4..4. `verdict()` now takes, for each pair (T_a, T_b), the full range k in [min deg T_b - max deg T_a, max deg T_b - min deg T_a] outside which Hom_K(T_a, T_b[k]) = 0 for degree reasons, united with -4..4. The largest |k| that was ever needed is printed ('max|k| needed'): 2 at widths 4 and 5, 1 at width 3. So -4..4 was more than enough and the test is now complete by construction. All counts re-run and unchanged (below).
2. Generation (done, wording). The script tests Hom(T,T[k]) = 0 for k != 0 and End(T), not that T generates K^b(proj A). Generation is inherited from T being an iterated silting/tilting mutation of A; this rests on the Aihara-Iyama / Okuyama-Rickard citation, theorem numbers UNVERIFIED. Now said in the Scope line.
3. Saved outputs (done). New file `workshop/rounds/087/skeptic_homotopy_full.out`, untruncated: every command below (run, mirror, ctl_after, ctl_seq, ctl_cox, the three extra cells) and the calibration for widths 3 and 4, lengths 6, 7, 8, all four (OP, SIGN) conventions, each preceded by a `CMD` line. Calibration result: the two good conventions pass everything (width 3: 10/10, 28/28, 84/84; width 4: 8/8, 20/20, 56/56); the two bad ones reject every application (e.g. 20 of 20 at 7/4 with NOTLINE 10, NOTTILTING 10). ctl_seq counts re-confirmed: NOTTILTING 826, DIM 98, NOTLINE 84. The referee's unchecked split is now in the file. Note the old `.out` already held the extra-cell lines; the old `skeptic_homotopy_L10_W5.out` is left as it was.
4. Title and claim (done). Title and Scope now say that the sequences are the repo's `VERIFIED_MOVES`: this audits sequence-to-End, not the table. Prior-record wording narrowed as the referee suggested: supports H-020 at this cell and adds no statement about the rule table's completeness.
5. Output location of w5 L11, w4 L10, w4 L11 (done): `skeptic_homotopy_full.out` (re-run: 3 168 OK, 60 s; 528 OK, 7 s; 1 716 OK, 29 s) and the older `skeptic_homotopy_L10_W5.out`.

Script change: only the k-range and the printed 'max|k| needed'; no verdict changed.

scope: sequence-to-End audit only (the sequences and predictions are the repo's `VERIFIED_MOVES`; generation of K^b(proj A) by T is inherited from silting mutation and not tested, resting on the UNVERIFIED Aihara-Iyama citation). Floating `VERIFIED_MOVES` of width 5 (24 rules), length 10, every rule at every window and every LNA, no cap, 1 008 applications. Extra cells run the same way: w5 L11 (3 168), w4 L10 (528), w4 L11 (1 716), all OK. Not covered: anchored rules, widths 6..11, lengths above 11, the rule table itself (read, not re-derived), mutation of non-LNA intermediate algebras beyond what these sequences produce, the field (GF(32003) only).

## Claim

For all 1 008 applications of the 24 floating width-5 rules at length 10, performing the rule's mutation sequence as silting mutations in K^b(proj A), coded from scratch, gives a tilting complex T whose endomorphism algebra is, as a labelled-up-to-line-order algebra, exactly the LNA the rule predicts. This is stronger than E-207's Coxeter test (which only needs the polynomial to agree): Hom(T,T[k]) = 0 for every k != 0 in the degree span of the complexes (never more than |k| = 2 needed; tested at least -4..4), the Hom-dimension pattern is that of the predicted LNA, and consecutive-arrow composites vanish exactly where the predicted paths vanish, so End(T) is isomorphic to the predicted LNA. It does not claim that the rule table is complete, nor anything for other widths or lengths, nor that the repo's `mutateAtVertex` agrees step by step with mine (only the end result is compared). Refutation would be any application with a verdict other than OK.

## Evidence

Own engine (`skeptic_homotopy.py`, no import of `mutation`, `pathAlgebra`, `procedure`, `nakayama`, `lines`; it imports `lnaMoves` only for `VERIFIED_MOVES` and `anchorOf`, the object under audit). Objects are bounded complexes of projectives of a monomial LNA (Hom(P_a,P_b) at most 1-dim, composition = matrix product masked by path non-vanishing); Hom_K and Hom_K(-,-[k]) by linear algebra mod 32003; a mutation at v takes the minimal one-sided add(T/T_v)-approximation (generators of Hom modulo the radical composites through other summands), forms the cone, and removes contractible summands by Gaussian cancellation; d^2 = 0 is asserted after every step. Left or right is chosen by the sign of the vertex.

Convention calibration (fixed before the cell, on other cells): two choices of side (OP) times two sign maps. On width 3 and 4 at L = 6, 7, 8 (no overlap with the cell) the choices (OP 0, negative = right-approximation) and (OP 1, positive = left-approximation) pass everything (28/28, 20/20, 56/56 ...) and the other two fail (NOTLINE/NOTTILTING), so the data pick the sign map up to the mirror symmetry; frozen as OP 0, SIGN -1 and used unchanged for the cell. The mirror choice (OP 1, SIGN +1) also gives 1 008 of 1 008.

| mode | what it tests | L10 W5 result |
|---|---|---|
| run | rules as stored | OK 1 008 / 1 008 |
| run, mirror convention (OP 1, SIGN +1) | convention independence | OK 1 008 / 1 008 |
| ctl_after | wrong claimed result (next rule's `after`, true sequence) | 1 008 rejected (DIM 1 008) |
| ctl_seq | wrong sequence (every vertex shifted by one) | 1 008 rejected (NOTTILTING 826, DIM 98, NOTLINE 84) |
| ctl_cox | claimed result replaced by a different LNA with the SAME Coxeter polynomial (own Cartan, sympy) | 1 008 rejected (DIM 1 008) |
| wrong sign map (OP 0, SIGN +1), calibration cells | wrong convention | rejected 100% (e.g. 20 of 20 at 7/4) |

ctl_cox is the real NO control: the Coxeter test of E-207 accepts all of these by construction (equal polynomial), this check rejects all, so it has power where the Coxeter check has none. A caveat: the three controls are not drawn from independent equivalent/inequivalent pairs; the ctl_cox alternative is a fixed pseudo-random pick among key-equal LNAs, which need not be inequivalent to the predicted one by any independent theorem (most probably are, but 1 008 DIM rejections show only that the labelled End algebra differs, as it must when the LNA differs).

Extra cells (same engine, same frozen convention): w5 L11 3 168 OK; w4 L10 528 OK; w4 L11 1 716 OK. No failure, no `Illegal`, no non-local End in any run.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/087/skeptic_homotopy.py 10 5 run 0 -1      # 14 s, 1008 OK
timeout 10m .venv/bin/python workshop/rounds/087/skeptic_homotopy.py 10 5 ctl_after 0 -1 # 12 s
timeout 10m .venv/bin/python workshop/rounds/087/skeptic_homotopy.py 10 5 ctl_seq 0 -1   # 6 s
timeout 10m .venv/bin/python workshop/rounds/087/skeptic_homotopy.py 10 5 ctl_cox 0 -1   # 23 s
timeout 10m .venv/bin/python workshop/rounds/087/skeptic_homotopy.py 10 5 run 1 1        # 13 s mirror convention
timeout 10m .venv/bin/python workshop/rounds/087/skeptic_homotopy.py 7 4 calib 0 -1      # calibration; 6 7 3 / 8 4 likewise
```
Output, untruncated, with calibration: `workshop/rounds/087/skeptic_homotopy_full.out`; older excerpt `skeptic_homotopy_L10_W5.out`.

## Prior record

E-159, E-165, E-188 check H1 with `verifyMove` (shared mutation engine, Coxeter polynomial only as the derived check); E-207 re-implemented matching, application and Coxeter but kept the engine. E-161, E-184 used an independent Hom(T,T[+-1]) test on path edges, but generation of the tilting step came from the repo's mutation. Not recorded before: a mutation step and a full End(T) = LNA decision (tilting + algebra iso) with no repo mutation code, and a Coxeter-equal inequivalent control. H-020 status is unchanged; this supports H-020 at this cell and adds no statement about the rule table's completeness; it is one more independent confirmation of H1, not new mathematics. The engine, as a silting mutation in the homotopy category, equals the intended tilting mutation by Aihara-Iyama / Okuyama-Rickard (cited from memory, UNVERIFIED theorem numbers), so agreement here also supports that the repo's `mutateAtVertex` on these sequences equals the abstract mutation.

## Code changed

None in the package. New file `workshop/rounds/087/skeptic_homotopy.py`. No tests needed (no package file touched).

## Next

- skeptic or toolsmith: run the same script on width 6..8 cells (cost grows with window; w5 L11 took 54 s so widths 6..8 at L11/12 fit in a few shards) and on width 9..11 floating rules, which are the weakest H1 cells; anchored rules need windows at the ends (the script's loop `range(1, L - W + 1)` is the floating one).
- toolsmith: compare my Hom-verified end algebra step by step with `mutateAtVertex` (intermediate quivers) if the engine, not just the endpoints, is to be audited; the script exposes the intermediate complexes.
- the same engine can give the missing inequivalence certificate for T10 only if a non-tilting-search route is found; it decides tilting for a given sequence, it does not search.
