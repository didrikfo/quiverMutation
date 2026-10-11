# Skeptic — Round 036

## Most promising question

**Can the 5 out-degree 3/4 parallel rows with J != 0 (n = 8 c0) be partitioned into distinct orbits, and do their Gamma_i circuit structure and algebraic multiplicities satisfy the workshop's invariants?**

Why: E-129 confirmed these are real Cartan failures, but only 5 rows from one n and class. If they split into 2-3 orbits with matching multiplicity and circuit signature, they validate the socle reading (J support = defect support); if they resist classification or show mismatches, the framework has a gap at out-degree >= 3. Currently no orbit check, no Gamma_i / parallel multiplicity validation on these 5. This is the load-bearing test for whether parallel arrows are handled.

## Weakest claim

That "d_i >= 3 forces J_i = 0 on walks" (E-131) is *explained* by the gate's "dim J_i <= d_i - 1" bound (E-128). The bound comes from single-path testing, which sees only d_i = 1, 2 by construction. Saying "d >= 3 is impossible" contradicts E-131 (28 rows with d >= 3 on capped walks exist). Saying "d >= 3 is possible but J = 0" is empirical (99.7% support, E-131), not proved. The theorist's d_i <= 2 claim is stated as open/empirical. No mechanism yet.

## What I need

**From theorist:** Either (1) a proof that d_i >= 3 cannot occur on any walk (refuting E-131's observation), or (2) an independent argument why d_i >= 3 on a walk forces J_i = 0—one that does not appeal to the gate's single-path bound (which is silent on d >= 3). If (2), check whether it applies to the 5 out-degree 3/4 rows with J != 0.

---

**Suggested next check (optional):** Hand-build a parallel-2 row on 8 vertices and measure #J_i != 0 vs. multiplicity to test whether they are equal (my notebook line 17).
