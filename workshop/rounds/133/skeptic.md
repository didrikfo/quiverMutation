# H-017b's per-LNA key test has power only for shifts of about 1.5 SD or more (n = 10), while the one key effect on record is 0.17 SD at the unit level (n = 9), and a YES would not single out the key: proposal to close T6 without a run

author: skeptic · round: 133 · kind: proposal
thread: T6 · bears on: H-017, H-017b (E-314 amendment), E-283, E-302, E-303
scope: power simulation of the registered statistic on the real key and LNA-feature tables at n = 10 (140 duality units) and n = 11 (1361 units), 200 replicates x 199 permutations per cell, planted effects synthetic (Gaussian per-unit noise); real defect data only at n = 9 (E-283's 645 members, 6 units) and 3 probe units at n = 10 (depth 5, 432 members). No n = 10/11 census was run. Nothing for depth >= 6.
power: planted key shift on the real n = 10 tables (23 units, 5 classes): YES 0.95 at 1.5 SD, 0.77 at 1.0 SD. Planted key-independent per-unit shift (truth NO): NO, rejection 0.01-0.06 across all four designs. But the hypothesis-level NO control fails: the source LNA's relation-length sum, which has no key term planted, reads as a "key effect" in 0.11 (n = 10, N = 15) to 0.84 (n = 11, N = 82) of samples.

## Claim

(1) **Unit.** `reachedQuipuAlgebras` walks A and A^op, so an LNA and its opposite give identical member sets. The unit is the duality orbit. At n = 9, E-283's 9 LNAs are 6 units: 3345000/5040330, 3455000/5504030 and 3505000/5050030 are pairs, and 3033030, 4444400 and 4550400 are singletons. At n = 10 the 262 LNAs make 140 units, in key classes of 1, 2, 3, 5, 5, 58 and 66 units. At n = 11 the 2647 LNAs make 1361 units in 20 classes, 17 of them with at least 3 units. Maverick's "about 15 LNAs" therefore has to be counted in units.

(2) **The motivating effect is not visible per LNA.** At n = 9 the per-unit mean defect of key class 2 (3033030, 2.711) sits 0.17 SD below the class-1 unit means (2.375-2.917, SD 0.230). E-303's shape effect does not stand out at the unit level either. Class 2 has no defect-1 members, and neither does class-1 unit 4444400. It has no defect-5 members, where every class-1 unit has exactly one. So per LNA, nothing at n = 9 separates the key from the LNA. This agrees with E-303's own conclusion and adds an effect size.

(3) **delta\*.** Take the registered statistic: the between-key/within-key sum-of-squares ratio of the residuals of the unit mean defect after OLS on [1, max C, nnz C, members]. Rank C is dropped because it is constant n for LNAs. Key labels are permuted among units. The planted shift is delta·f(key), with f standardised, in units of the within-key SD. Then delta\* (smallest grid value with power >= 0.8) is:
- 2.0 SD at n = 10 with 3 units per class (N = 15);
- 1.5 SD at n = 10 with up to 5 per class (N = 23);
- 1.0 SD at n = 11 with 3 per class (N = 51; 17 classes).

Maverick's condition for retiring the key form needs delta\* below E-303's shift, about 0.17 SD. Scaling N·delta² from the n = 11 row, that takes about 1800 units. n = 10 has 140 units and n = 11 has 1361. The condition cannot be met at n <= 11 even with every unit walked.

(4) **A YES would not be about the key.** The key is a function of the LNA. Any LNA quantity left out of the feature list but correlated with the key passes as a "key effect", and the rate grows with N. With the relation-length sum as the response and no key term: 0.11 (n = 10, N = 15), 0.27 (N = 23), 0.41 (n = 11, N = 51), 0.84 (N = 82). With the relation count: 0.10, 0.19, 0.34, 0.62. With the minimum relation length: 0.03-0.10, close to nominal. Adding those features to the regression only moves the problem to the next unlisted feature. Maverick's text already says a YES reads "key or an unmodelled LNA feature". The numbers show that at the design sizes, the "unmodelled feature" reading is the likely one.

**Proposal.** Close T6 for H-017b without a run. A NO retires nothing, by H-017b's own text, and the condition that would let it retire something cannot be reached at n <= 11 (3). A YES is confounded (4). The only effect on record is 9 to 12 times below delta\* (2). H-017 stays OPEN on its main statement (relations > cords, i.e. minimum defect >= 1). That statement has a direct falsifier and is the part worth an overnight run (Next).

Not claimed: that the key has no effect on the defect, or that a large key effect at n = 10 is excluded. The run was not made, so neither is shown. The planted noise is Gaussian, while the real per-unit spread at n = 10 is unknown (3 probe units have means 2.78, 3.30 and 3.34). Refutation of (3)/(4): a design with N <= 140 at n = 10 whose power at 0.2 SD is at least 0.8, or a feature list that keeps the relation-sum control at nominal size for N >= 23 without naming that feature.

## Evidence

Power and controls (rejection rates at alpha 0.05; 200 samples per cell, each drawn without looking; `skeptic_power10.txt`, `skeptic_power11.txt`):

| design | N, classes | d=0 | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 | LNA-only shift | maxC+noise | conf. nrel | conf. sumrel | conf. minrel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n=10, 3/class | 15, 5 | 0.02 | 0.14 | 0.42 | 0.69 | 0.82 | 0.92 | 0.01 | 0.01 | 0.10 | 0.11 | 0.03 |
| n=10, <=5/class | 23, 5 | 0.04 | 0.26 | 0.77 | 0.95 | 0.99 | 0.99 | 0.03 | 0.01 | 0.19 | 0.27 | 0.04 |
| n=11, 3/class | 51, 17 | 0.02 | 0.35 | 0.96 | 1.00 | 1.00 | 1.00 | 0.06 | 0.04 | 0.34 | 0.41 | 0.06 |
| n=11, <=5/class | 82, 17 | 0.04 | 0.66 | 1.00 | 1.00 | 1.00 | 1.00 | 0.03 | 0.04 | 0.62 | 0.84 | 0.10 |

n = 9 units (E-283 data, round-105 json): mean defect 2.711 for 3033030 (class 2, 45 members); 2.915, 2.857, 2.375, 2.917 and 2.683 for the class-1 units. The defect-1 fraction is 0 for 3033030 and 4444400, and 0.04-0.27 for the others. The defect-5 count is 0 for 3033030 and 1 for each class-1 unit.

n = 10 probe walks (E-283 walk, depth 5; these are cost sizing, not a census):

| unit (key class size) | members | defects | min | time |
|---|---|---|---|---|
| 30033030 (3) | 120 | 1:16 2:32 3:45 4:16 5:11 | 1 | 155 s |
| 00334400 (66) | 174 | 1:6 2:24 3:77 4:49 5:15 6:3 | 1 | 140 s |
| 03345000 (58) | 138 | 2:29 3:52 4:39 5:17 6:1 | 2 | 90 s |

No member has defect <= 0, so H-017's main statement holds on these three units (432 members). The cost is about 2 min per unit at n = 10, so the full census would be 140 units, about 4.7 CPU-hours. That fits sharded 10-minute commands, or one night.

One bug found and fixed in my own script: the first n = 11 run merged two key classes because `hash(-1) == hash(-2)` in Python. Keys are now indexed, not hashed. The table above is the rerun.

## Reproduction

```
.venv/bin/python workshop/rounds/133/skeptic_h017b.py units 10        # 7 s; 'units 11' about 1 min
.venv/bin/python workshop/rounds/133/skeptic_h017b.py n9              # 1 s
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_h017b.py power 10 1   # about 2.5 min
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_h017b.py power 11 1   # about 6 min
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_h017b.py collect 10 30033030   # 90-155 s per unit
```

## Prior record

- E-283: the defect is not a function of the key (n = 9).
- E-302: "key irrelevant" is withdrawn.
- E-303: the key shifts the member-level distribution (MI) but cannot be separated from the source LNA.
- H-017 amendment / E-314 item 19: maverick's H-017b text (rounds/119/maverick.md).

New here: the duality-unit count; the unit-level effect size of E-303 (0.17 SD); delta\* on the real tables; the unreachability of the retire condition; and the confound rates. No E-entry gives n = 10 depth-5 defect data. The 3 probe units are the first, and they are not a census. Nothing in RETRACTIONS bears on this.

## Code changed

None in the library. New script: `workshop/rounds/133/skeptic_h017b.py`. Outputs: `skeptic_power10.txt`, `skeptic_power11.txt`, `skeptic_defect_10_*.json` (3 files, about 35 KB in total) and `skeptic_collect_probe*.txt`. No tests were run, since no library file was touched.

## Next

- Chair / next conference ledger: close T6's H-017b form, for the reasons in (2)-(4). Record delta\* as 1.5 SD (n = 10, N = 23) and 1.0 SD (n = 11, N = 51), so that a later reopening starts from these numbers.
- If the conference wants a T6 run anyway, it should be the H-017 main-statement census, not H-017b. The OVERNIGHT proposal text, not run:
  - **Command:** at n = 10, `skeptic_h017b.py collect 10 NAME` for all 140 duality units (list from `units 10`), sharded one unit per command, about 4.7 CPU-hours.
  - **Registered:** a member with relations <= cords refutes H-017 at n = 10. "None found" reads "no member within depth 5", and no more (E-071).
  - **Positive control first:** the same walk from 3 LNAs in a quipu class must reach a defect <= 0 member. UNVERIFIED: whether `quipuRelations.certificate` admits relation-free quipus; if it does not, build the control from a quipu-with-cords algebra of defect <= 0 first.
  - **Secondary, descriptive only:** the per-class unit-mean table, with no key claim.
- Maverick: if H-017b is kept, the feature list must be fixed before the run and must include the relation-length sum. Even then (4) applies to the next unlisted feature.
