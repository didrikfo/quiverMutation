# Maverick's notebook

## What I now believe (after round 138)
- **T3/T8 F-058 second prime (r138): done.** nu2 on all uniserials of 040404000 (m 72), 45050400 (m 84; split into 3 object chunks to fit 10 min), N_14(4), N_14(6) (m 132) at GF(2) and GF(3): same returns, shifts and max complex sizes as GF(2^31-1). One code only at these primes (nu_bimod not re-run).
- **Two families = two periodic tau-orbits (r138, tested on small cases).** On 050404400, 040404000, 45050400, N_14(4), N_14(6), N_15(3) the returning uniserials form exactly 2 nu-orbits up to shift, one per shift; first returns tau^p X ~ X[r]: (8,2)/(6,1); (9,1)/(4,1); (7,2)/(6,1) on 45050400; (11,2)/(6,1) on both N_14(4) and N_14(6); (8,1)/(9,1) on N_15(3). All other uniserial orbits (6-23) do not return by m. Control N_11(4): 7 orbits, all (15,1). The shift at m is (m/p)(p+r).
- "Family closed under tau" is automatic (nu = tau[1]); the content is "one orbit per family".
- Module-level AR sequences are AR triangles of D^b only when pd M <= 1 and Hom(M,A) = 0 (or dual); nu2 agrees at k = 1 (11/11 yes, 0/81 otherwise). Such M: 3-4 per n = 10-11 algebra, 0 at n = 14-15, none touching a returning orbit. So that route is vacuous.
- H-025 (r135): two independent nu's agree object by object at n = 14-15. The per-step Serre check does NOT pin nu (it misses a dimension-preserving iterate error); the cross-code agreement does. A coherent wrong nu (identity) passes every assert at every prime.
- T5 (r131): the gate's cycle failure cannot be reached from LNAs but is real for random cyclic monomial algebras. T10 (r129): R inverts F only up to the ideal (compare with `sameIdeal`).

## What I tried
- r138: `138/maverick_tau.py REL m [--prime p]`: iterates nu2 on every uniserial with a stalk test at each k; prints first returns, nu-orbits among uniserials (union-find), and the AR-sequence test plus naive control. Costs: N_14(4) 132: 123 s; N_14(6) 132: 272 s; N_15(3) 144: 70 s; 040404000 72: 238 s; 45050400 84: 489 s (near the cap).
- r135: maverick_nu2.py (REL m objs [--prime p] [--serre] [--plant=...]); maverick_compare.py (vs nu_bimod). Objects can be given as a list (U4-5,...) to split a long run into chunks.
- Hom engine for cyclic MONOMIAL algebras: `skeptic_conv3.AlgG(n, arrows, [{path: 1}], False)`.

## Watch for
- numpy int64 mod 2^31-1: matrix products overflow; use split multiply (`mulmod`).
- Do not edit a script while a background loop uses it. `grep | tail` in background buffers.
- Exec-importing scripts: split on "\nif __name__ == '__main__':".
- 4 cores: four 10-min jobs in parallel are fine; time roughly 1.0-1.3x the solo time.
- A NO needs a control that says NO; a "0 found" needs the detector shown firing elsewhere.

## Next
- Why exactly two periodic tau-orbits? Get a prediction of (p, r) from a weighted-projective-line / Fuchsian model (2112.15587) and test it on N_16(8) <2,3,11] (not yet run by anyone; also run there with both codes).
- Are the non-returning orbits genuinely aperiodic (growth of total H), or periodic at larger m?
- Two-code agreement at GF(2)/GF(3) for the four uniserial groups (nu_bimod with a prime flag, if it has one).
- Idea (open): check the Serre identity against independently computed Hom(X, Y) for orbit pairs, for a nu-free certificate checker.
- Still open: gate grading idea (r131); n = 11 depth 7 13-class (r129).
