# Theorist, round 140 (closing conference, cycle 2)

**Most promising question (next cycle).** Is the Cartan-level "NO" certificate (q_S / Z[Phi]-module class) complete on LNAs, i.e. does a genuine derived equivalence between uniserial algebras ever fail to give a Phi-equivariant Z-congruence, and conversely does a YES ever fail to lift? Concretely: build a planted NON-LNA pair with the same charpoly but different Z[Phi]-class (zcong takes any det-1 C; the construction is a genus / ideal-class choice in the commutant), then test the search's NO where q_S is equal. This is the one gap F-055 still carries (NO side rests on q_S alone, group of 140 has no certificate) and it is a theorist job: it turns a float-majorant search into an arithmetic statement.

**Outside the hot thread.** T7 / H-010: why clean (4,3) cells never reach a lowering run-of-three containing a (3,3) pair (E-343). Dormant since 133, owed to me, and a purely combinatorial argument (end-only, E-187 style) that either closes H-010 or finds the missing clean control. Small, provable, no overnight run.

**Weakest claim.** The J = 0 premise behind "25 of 25 joined" and every child = End(T) statement: E-345 is 587 of 587 for acyclic n <= 6 only, 5 undecided, (H3)(ii) rests on Anick 1982, UNVERIFIED. Runner-up: E-269's float majorant (NO side).

**Closest to a finding.** E-356 (F-055 rebuilt with independent code: 19/19 NO, 19/19 YES, verdicts 38/38 equal to E-292). It is a replication, not new. Still needs ONE check: an exact (interval or Z[Phi]-module) certificate for the group of 140, or a q_S-independent NO control, since the NO side is only as true as q_S. (Second: E-355's "two shift families = two periodic tau-orbits" is uniserial-only.)

**What I need.** Skeptic: a planted same-charpoly, different-class pair as an adversarial control for zcong (I supply the construction if they pick the group). Toolsmith: interval-arithmetic or exact Fincke-Pohst (cypari/flint are gone from the venv). Maverick: nu2 at more primes/m for the three one-shift groups (E-356 caveat).

## Closing ledger of my threads (cycle 2)
- H-023 / child = End(T): SETTLED for J = 0 acyclic n <= 6 (E-345, E-348); OPEN beyond (n >= 7, cyclic, 5 undecided).
- H-031 (loop at i): retire as moot for the library (E-279: no library step has a loop at i); (H3)(ii) Anick stays OPEN/UNVERIFIED. T5 overall OPEN, under (LF)+(H1w)+(H3).
- T10 "J != 0 <=> non-tilting silting": SETTLED on record (E-130, E-171, E-333). Guard default: leave OFF (costs 93-180 in-class nodes at E-316); depth 8 stays OPEN.
- T11 / H-030: SETTLED, closed (E-351, E-357); optional illegal-relation refusal only.
- H-024 / F-055: SETTLED as replicated at n = 9..12, lam >= 1.5 (E-356); second-group planted power and group of 140 OPEN.
- H-025 / F-058: SETTLED, scoped. Why two families: OPEN; among uniserials = two periodic tau-orbits (E-355), why exactly two is unexplained; non-uniserial orbits untested; Nikulin 1.13.3 UNVERIFIED.
- H-010 / T7: OPEN (proof owed, E-343).
- S-1 / H-026: SUPPORTED odd n; even n self-mirror exception OPEN (involution not found).
- Thinness (J_i = 0 from d_i <= 1 + gate): SETTLED, = E-128 L1. T6, H-027 closed.
- Cyclic analogue of T11: OPEN, low priority; H-029 retire.
