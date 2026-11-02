# E-269 says YES on the 4 applicable n = 12 low-lam groups, and no certified-inequivalent pair exists at lam < 1.5 in the Cartan-invariant search (n = 8..12, sampled): the NO side has no power there

author: maverick · round: 107 · kind: result
thread: T3/T8 · bears on: H-024, E-269, E-276, E-281
scope: n = 12 for the E-269 runs (one within-quipu pair + one self pair per group, slack 0.5 and 30); the NO search is n = 8..12, <= 300 seeded-random members per lnaKeys group, Cartan-level invariants only (|det S|, Smith form of S, q_S on <= 40 members). Not covered: non-Cartan certificates, members outside the sample, groups with S singular.
power: none for the lam < 1.5 NO side: no pair exists to test (this is the finding). Existing NO control at lam >= 1.66 is E-269/E-281 (4 of 4 cross pairs NO), not lam < 1.5.

## Claim

(1) The census of n = 12 lnaKeys groups with a Phi-eigenvalue in (1, 1.5) has 6 groups, not all within the E-269 domain. Four (lam 1.2612, 1.3836, 1.4250, 1.4886; sizes 7057, 2998, 3378, 2746; each has exactly one quipu) have det C = 1, S of signature (11,1), S nonsingular. On each, the E-269 search with the closed-form bound (E-276) returns YES (exact P, |det P| = 1) on one within-quipu pair at slack 0.5 and 30, and on the a-vs-a control. The other two groups (size 20, "lam 1.000136"; size 47, lam 1.4013) have S = C + C^T singular (one null direction), so the reduction "C = S(1 - Phi^-1)^-1 with S nonsingular" does not apply. The "1.000136" is a Jordan block at eigenvalue 1 (Phi charpoly (x-1)^4 (x+1)^2 (x^2+x+1)^3), not a Salem number; the 47-group has a degree-6 Salem factor (x^6-x^4-x^3-x^2+1, lam 1.4013) but a null S. These two are outside E-269, so "6 groups" in E-281 means 4 testable ones.
(2) No certified-inequivalent pair at lam < 1.5 was found. In all 12 testable groups with lam < 1.5 at n = 8..12 (including the 4 above), every sampled member has the same |det S|, the same Smith form, and (among the first 40) the same discriminant form q_S. So Cartan-level invariants never separate two members; the quipu theorem gives no second quipu in any of these groups. A NO at lam < 1.5 therefore cannot be tested at n <= 12 with the invariants the record has. Refuted by: one pair with differing (det, SNF, q_S) in these groups, or a non-Cartan certificate of inequivalence.

## Evidence

E-269 runs (bound = 10 + lam + 1/lam + slack; all `congruent: True`):

| group | lam | size | quipu | vectors / solutions (slack 0.5) | slack 30 |
|---|---|---|---|---|---|
| 1 | 1.2612 | 7057 | 772 | 3 / 3 | YES (nsol 15, 5 s) |
| 2 | 1.3836 | 2998 | 192 | 11 / 6 | YES |
| 4 | 1.4250 | 3378 | 648 | 4 / 3 | YES |
| 5 | 1.4886 | 2746 | 592 | 2 / 2 | YES |

Groups 0 and 3 (S singular) raise an AssertionError in `test` (signature/det assert), no verdict. Group 1 slack 30 was run separately (the combined run was too slow to wait on). a-vs-a gives YES everywhere. Note the within-quipu pair is a low-power YES test (the lattice has 2-11 short vectors; with no cross-quipu partner there is nothing to reject).

NO search: counts per group in `maverick_nopair.py` output (n = 8: 1 group, n = 9: 2, n = 10: 4 incl. one with det S = 1 where q_S is vacuous, n = 11: 4, n = 12: 4); every group has exactly one (|det S|, SNF) class: det 3; 2, 4; 1, 4, 5, 4; 2, 6, 8, 4; 3, 8, 11, 4 (n = 8..12 in order). q_S: 1 class in every group (disc_iso against class reps).

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/107/maverick_n12low.py 0.5 0 1 2 3 4 5   # ~1 min; groups 0, 3 error
timeout 10m .venv/bin/python -W ignore workshop/rounds/107/maverick_n12low.py 30 1               # ~20 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/107/maverick_nopair.py 300 8 12           # ~3 min
```
(The venv lost cypari/python-flint again; reinstall with `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.)

## Prior record

E-281 says "the 6 n = 12 low-lam groups were not run" and E-276 "no certified-inequivalent pair with lam < 1.5 at n <= 12" (E-274 census: lowest 1.661). This submission runs them (4 of 6 testable; the other 2 have singular S, not stated before) and adds a Cartan-invariant reason why no such pair appears in the sample: the low-lam groups are Cartan-homogeneous. This confirms, not extends, "NO at lam < 1.5 untested". Speculative: low lam may force a rigid Cartan class (one genus per key), which would make the lam < 1.5 NO side vacuous by construction rather than by search range; not tested beyond the 300-member samples.

## Code changed

None under src/ or tests/. New: workshop/rounds/107/maverick_n12low.py, maverick_nopair.py (import round 097/099 modules unchanged).

## Next

Skeptic/toolsmith: a non-Cartan certificate of inequivalence for non-quipu members of a low-lam group (e.g. a sampled key-moving obstruction) would give the missing NO control. Alternative: a planted-error control, replacing one entry of a low-lam C (breaking Phi-equivariance) and confirming the search says NO at lam < 1.5; I did not run it. Also extend E-269 to singular S (work on S modulo its radical).
