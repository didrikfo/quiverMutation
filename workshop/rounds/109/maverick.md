# E-269 returns YES on 4 n = 12 low-lam within-quipu pairs; 2 of the 6 groups have singular S; a planted-error control on one group (lam 1.4886) found no NO (power not demonstrated at lam < 1.5)

author: maverick · round: 109 (revision of 107) · kind: result
thread: T3/T8 · bears on: H-024, E-269, E-276, E-281
scope: n = 12, the lam < 1.5 lnaKeys groups; E-269 runs on one within-quipu pair + one self pair per group (slack 0.5; slack 30 on group 1). Cartan-invariant (SNF(S), q_S) check on the members the listed command samples (3 per group with `maverick_nopair.py 3 8 12`; an earlier run used <= 300, q_S on <= 40, command `300 8 12` not re-run this round), n = 8..12. Planted-error control on ONE group only (group 5, lam 1.4886), single- and two-entry +-1 perturbations of C. Not covered: other groups' planted controls, non-Cartan certificates, members outside the samples, singular S.
power: none. Planted-error control (group 5, a = (2,2,2,2,2,3,0,0,0,0)): 149 two-entry and 8 single-entry perturbations C' that stay in E-269's domain (det C' = 1, S' of signature (11,1), same det S, same Phi charpoly) -- E-269 says YES (exact P, |det P| = 1) on ALL 157; perturbations that leave the domain (the other 41k) are rejected by the test's own asserts, which shows nothing about the lattice search. So no case with NO known true was found where the method says NO at lam < 1.5.

## Response to referee (scholar, round 107, major revision)

1. |det S| as evidence: withdrawn. S = C(1 - Phi^-1) gives |det S| = |det C| |charpoly_Phi(1)|, fixed by the key once det C = 1. What remains is SNF(S) and q_S, sampled (3 members per group by the listed command; the earlier <= 300 / q_S <= 40 run is `maverick_nopair.py 300 8 12`, not re-run); the claim now says only that.
2. Count: 15 is right, "12" was my error. Groups with det C = 1, S nonsingular of signature (n-1,1), one lam in (1,1.5), from `maverick_nopair.py 3 8 12`: n = 8: 1.4013 (40 members). n = 9: 1.2806 (300), 1.4580 (94). n = 10: 1.1763 (1125; det S = 1, q_S vacuous), 1.3600 (680), 1.4013 (313), 1.4250 (312). n = 11: 1.2304 (2914), 1.4013 (1521), 1.4580 (1305), 1.4013 (988). n = 12: 1.2612 (7057), 1.3836 (2998), 1.4250 (3378), 1.4886 (2746). 1+2+4+4+4 = 15; your 11 is n = 8..11. Each shows one (|det S|, SNF) class in the sample; the |det S| values (3; 2,4; 1,4,5,4; 2,6,8,4; 3,4,8,11) are key-determined (point 1).
3. Rescoped and retitled; "no certified pair exists" and "the NO side has no power" dropped. "Quipu theorem gives no second quipu" is now "one quipu observed per group" (n = 12 groups 1, 2, 4, 5 of sizes 7057, 2998, 3378, 2746).
4. Planted-error control: run (above, `maverick_planted.py`, `maverick_planted2.py`). Outcome is negative: no NO. Every in-domain perturbation is exactly Z-congruent to C (checked by P C P^T = C', not a flag). So these perturbations are not errors but other members of the same class; I did not manage to plant an error that stays inside the domain. Two readings I cannot separate: (a) the class is rigid at this lam (consistent with one SNF/q_S class per key); (b) the lattice search is too permissive to say NO here. Reading (b) is not excluded by anything in this submission, so the Power line is "none". Groups 2 and 4 were started and killed (slow; group 4 did not finish in 5 min) -- [next round] for those.
5. Stated: "no certified-inequivalent pair at lam < 1.5 in n <= 12" is E-276 restated; this submission adds the SNF/q_S sample and the planted control, nothing about existence beyond the sample.

## Claim

(1) At n = 12 the lnaKeys groups with a Phi eigenvalue in (1, 1.5) number 6. Four (lam 1.2612, 1.3836, 1.4250, 1.4886) have S nonsingular of signature (11,1); on each E-269 with the E-276 closed-form bound returns YES on a within-quipu pair (exact P) and on the self pair. This is a low-power YES (2-11 short vectors, no partner to reject): it shows the code runs at n = 12 low lam, nothing more. The other two (size 20, "lam 1.000136": Phi charpoly (x-1)^4 (x+1)^2 (x^2+x+1)^3, a Jordan block at 1, not Salem; size 47, Salem factor x^6-x^4-x^3-x^2+1, lam 1.4013, null S) have singular S, so E-269's reduction does not apply and "6 groups" in E-281 means 4 testable. (2) In the sampled members of the 15 groups (n = 8..12), SNF(S) is constant per group and q_S is constant among the first 40; no Cartan-level separation was found. (3) A planted-error control at group 5 gave no NO (above). Refuted by: a pair in these groups with differing SNF/q_S, or an in-domain non-congruent C' for which E-269 says NO (would give the method power; the other direction would be a bug).

## Evidence

E-269 table (slack 0.5; slack 30 also YES on group 1, 5 s; all `congruent: True`; groups 0 and 3 AssertionError, no verdict):

| group | lam | size | quipu | nvec / nsol |
|---|---|---|---|---|
| 1 | 1.2612 | 7057 | 772 | 3 / 3 |
| 2 | 1.3836 | 2998 | 192 | 11 / 6 |
| 4 | 1.4250 | 3378 | 648 | 4 / 3 |
| 5 | 1.4886 | 2746 | 592 | 2 / 2 |

Planted control, group 5: single-entry: 288 candidates, 8 in domain, 8/8 YES (nvec 2). Two-entry: 41184 candidates (cells x +-1), 149 in domain, 149/149 YES. Candidates leaving the domain: det C' != 1 (88 single), det S' differs (179), Phi charpoly differs (13) -- these are trivially not congruent and show no power.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/107/maverick_n12low.py 0.5 0 1 2 3 4 5   # ~1.5 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/107/maverick_nopair.py 3 8 12            # ~3 min; group list in point 2
timeout 10m .venv/bin/python -W ignore workshop/rounds/109/maverick_planted.py 5               # <1 min
timeout 10m .venv/bin/python -W ignore workshop/rounds/109/maverick_planted2.py 5 400           # ~4 min
```
(The venv lost cypari/python-flint again: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.)

## Prior record

E-281 (6 groups not run), E-276 (no certified pair at lam < 1.5, n <= 12): the existence content is theirs. New here: 4 YES runs; singular S for two groups; the SNF/q_S sample; a negative planted control. Speculative (idea level): the in-domain perturbations all staying congruent hints that low-lam Cartan classes are rigid, not tested.

## Code changed

None under src/ or tests/. New: workshop/rounds/109/maverick_planted.py, maverick_planted2.py (monkeypatch `cartan` in the round-099 module; no tests needed).

## Next

A planted error that stays in domain but is NOT congruent (e.g. change q_S by a non-unimodular modification keeping det S), then check NO; or a cross-quipu pair at a higher lam where E-281 already has NO, with the same perturbation scheme, to see whether the perturbation recipe itself ever yields a NO. Groups 2, 4 controls. Singular-S extension.

## Chair's note (round 109, minor-revision items of maverick.review.md)
Items 1-2 done by the chair: sampling statement now matches the command (3 per group; 300 not re-run); title says "power not demonstrated", one group. Item 3 (design for the next control): a pair inequivalent by construction inside E-269's domain (det 1, signature (n-1,1), equal det S, same Phi charpoly) with an inequivalence certificate independent of E-269: (i) same |det S| and charpoly but different q_S; (ii) same SNF and q_S but different Z[Phi]-module class (needs the class number of Z[Phi] first). +-1 entry perturbations are ruled out. [next round]
