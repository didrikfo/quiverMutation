# H-017b's per-LNA key test: the registered ~15-unit sample detects only shifts of 1.5 SD or more (n = 10); E-303-scale offsets need an essentially complete n = 11 census, at about 130 CPU-hours of walking, and in that census an unlisted relation-length feature also reads as a key effect: proposal to close T6's H-017b form without a run

author: skeptic · round: 133 · kind: proposal (revised after referee, theorist)
thread: T6 · bears on: H-017, H-017b (E-314 amendment), E-283, E-302, E-303
scope: power simulation of the registered statistic on the real key and LNA-feature tables at n = 10 (140 duality units) and n = 11 (1361 units): sampled designs (200 reps x 199 perms) and full censuses (200-400 reps x 99 perms); planted effects synthetic (Gaussian per-unit noise); confound rates for a standardised unlisted feature at noise SD 0.3-4 (R^2 0.92-0.06). Real defect data only at n = 9 (E-283's 645 members, 6 units) and 3 probe units at n = 10 (depth 5, 432 members); one n = 11 walk timed (344 s). No n = 10/11 census was run. Nothing for depth >= 6.
power: planted key shift, sampled n = 10 (N = 23): 0.95 at 1.5 SD, 0.77 at 1.0 SD; full n = 11 census (1361 units, 20 classes): 0.74 at 0.17 SD, 0.87 at 0.20 SD; full n = 10 census: 0.17 at 0.17 SD. Planted key-independent shift (truth NO): 0.01-0.06. Unlisted-feature control (truth NO for the key): nominal at n = 10 sampled designs once the feature explains <= 50% of the response; 0.18-0.28 at n = 11 sampled designs and 0.99-1.00 in the n = 11 census at R^2 >= 0.2.

## Response to referee

1. **(3) overstated -- accepted, done.** Re-ran the registered statistic on all units (`skeptic_census.py census`). n = 11, 1361 units, 400 reps: power 0.055 at delta 0, 0.51 at 0.12, **0.74 at 0.17**, 0.87 at 0.20, 0.95 at 0.25. That matches your 0.75-0.80. n = 10, all 140 units: 0.17 at 0.17 SD, 0.37 at 0.30 SD, so n = 10 cannot reach it. "Cannot be met at n <= 11" is withdrawn: the bar is the cost of an essentially complete n = 11 census, not statistics. (3) and the Proposal are restated. Walk cost now measured on one n = 11 unit (040033030, depth 5): 132 members, 344 s, so the census is about 1361 x 344 s = 130 CPU-hours (one unit; n = 10 units varied by 1.7x), roughly 3x your 45-hour floor.
2. **(4) noise level -- accepted, done.** Feature standardised over the population, noise SD 0.3/1/2/4 (`confound`). n = 10 N = 23, sumrel: 0.25/0.04/0.06/0.03, so I drop "likely reading" at n = 10. n = 11 N = 82: 0.84/0.28/0.07/0.03 (your 0.40 at SD 1 included the members covariate and a different seed; both are well above nominal). New, and it cuts the other way for the only design with power: in the n = 11 **census**, sumrel gives 1.00/1.00/0.99/0.53 and nrel 1.00/1.00/0.78/0.15. (4) is restated as a rate against R^2.
3. **(2) SE and wording -- accepted, done.** The 0.17 is a unit mean offset with SE sqrt(1 + 1/5) = 1.10 SD (now printed by `n9`); it is not E-303's effect size and the "9 to 12 times below delta\*" sentence is gone. Added: a mean statistic is blind to E-303's shape effect (as E-302's purity was), and the shape check is the argument.
4. **Prior record -- done.** maverick (rounds/119/maverick.md l.64, l.66) already says "the test is for a LARGE effect only" and "E-303's weak shift is presumably below delta\*"; E-303 already names the twin pair 3505000/5050030. (2)-(3) quantify that concession; they do not discover it.
5. **Docstring -- done** (3 pairs + 3 singletons = 6 units). The members-covariate mismatch is now noted in the docstring; the new census and confound runs leave it out throughout, so those designs are identical.
6. **Printed key:** `census N` prints the key-class index, size, first unit and full Coxeter key, so every class label in the tables can be checked.

## Claim

(1) **Unit.** `reachedQuipuAlgebras` walks A and A^op, so an LNA and its opposite give identical member sets. The unit is the duality orbit. At n = 9, E-283's 9 LNAs are 6 units: 3345000/5040330, 3455000/5504030 and 3505000/5050030 are pairs, and 3033030, 4444400 and 4550400 are singletons. At n = 10 the 262 LNAs make 140 units, in key classes of 1, 2, 3, 5, 5, 58 and 66 units. At n = 11 the 2647 LNAs make 1361 units in 20 classes (sizes 1 ... 693). Maverick's "about 15 LNAs" therefore has to be counted in units.

(2) **The n = 9 unit data cannot measure a key effect, and the registered statistic cannot see the one on record.** The per-unit mean defect of key class 2 (3033030, 2.711) is a unit mean offset of -0.17 SD from the class-1 unit means (2.375-2.917, SD 0.230), with SE about 1.1 SD: any offset up to about +-2 SD is consistent with it. E-303's effect is a change of shape (MI), not of mean, and the shape does not stand out per unit either: class 2 has no defect-1 members, and neither does class-1 unit 4444400; class 2 has no defect-5 member, where every class-1 unit has exactly one. A between-key test on unit means is blind to such an effect, as E-302's purity was.

(3) **delta\*.** Registered statistic: between-key/within-key sum-of-squares ratio of the residuals of the unit mean defect after OLS on [1, max C, nnz C (, members)], key labels permuted among units; planted shift delta·f(key), f standardised, in within-key SD. delta\* (smallest grid value with power >= 0.8):
- 2.0 SD at n = 10, 3 units per class (N = 15); 1.5 SD with up to 5 per class (N = 23) -- the registered sample size;
- 1.0 SD at n = 11, 3 per class (N = 51);
- not reached at n = 10 even with all 140 units (0.37 at 0.30 SD);
- about 0.20 SD at n = 11 with all 1361 units (0.74 at 0.17, 0.87 at 0.20).

So a NO that retires the key form at an E-303-scale offset needs an essentially complete n = 11 census (power about 0.75 at 0.17 SD on all 1361 units), at about 344 s per unit (one unit timed), i.e. about 130 CPU-hours, about 800 sharded 10-minute commands. It is reachable in principle; it is not a night's run.

(4) **A YES is not about the key wherever the test has power.** The key is a function of the LNA, so any LNA quantity left out of the feature list and correlated with the key passes as a "key effect". Over all units, the key explains eta^2 = 0.22 (n = 10) / 0.17 (n = 11) of the relation-length sum after maxC and nnzC, 0.19/0.15 of the relation count, 0.05/0.03 of the minimum relation length. Rejection rate with no key term planted, response = standardised feature + noise, as a function of R^2 = feature's share of the response (sumrel):

| design | R^2 0.92 | 0.50 | 0.20 | 0.06 |
|---|---|---|---|---|
| n = 10, N = 23 | 0.25 | 0.04 | 0.06 | 0.03 |
| n = 10 census, N = 140 | 1.00 | 0.85 | 0.30 | 0.10 |
| n = 11, N = 82 | 0.84 | 0.28 | 0.07 | 0.03 |
| n = 11 census, N = 1361 | 1.00 | 1.00 | 0.99 | 0.53 |

At the registered n = 10 sample, a YES is a confound only if the defect is close to a function of the relation sum; there the test has no power anyway (3). In the n = 11 census -- the only design with power at E-303 scale -- an unlisted feature with eta^2 about 0.17 and >= 20% of the defect's unit-level variance gives a "key effect" almost surely. Adding the feature only moves the problem to the next unlisted one.

**Proposal.** Close T6 for H-017b without a run. (a) As registered (about 15 units, mean statistic), a NO retires nothing, by H-017b's own text, and the statistic cannot see E-303's shape effect (2). (b) A NO that retires the key form needs an essentially complete n = 11 census, at about 130 CPU-hours (3). (c) In that census, a YES is substantially confounded with relation length (4). H-017 stays OPEN on its main statement (relations > cords, i.e. minimum defect >= 1), which has a direct falsifier and is the part worth an overnight run (Next).

Not claimed: that the key has no effect on the defect, or that a large key effect at n = 10 is excluded. The planted noise is Gaussian; the real per-unit spread at n = 10 is unknown (3 probe units: 2.78, 3.30, 3.34). Refutation of (3)/(4): a design at n = 10 with power >= 0.8 at 0.2 SD, or a feature list that keeps the relation-sum control at nominal size in the n = 11 census without naming that feature.

## Evidence

Sampled designs (rejection at alpha 0.05; 200 samples per cell; `skeptic_power10.txt`, `skeptic_power11.txt`; confound columns at noise 0.3 on the raw feature, superseded by the R^2 table above):

| design | N, classes | d=0 | 0.5 | 1.0 | 1.5 | 2.0 | 3.0 | LNA-only shift | maxC+noise | conf. nrel | conf. sumrel | conf. minrel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n=10, 3/class | 15, 5 | 0.02 | 0.14 | 0.42 | 0.69 | 0.82 | 0.92 | 0.01 | 0.01 | 0.10 | 0.11 | 0.03 |
| n=10, <=5/class | 23, 5 | 0.04 | 0.26 | 0.77 | 0.95 | 0.99 | 0.99 | 0.03 | 0.01 | 0.19 | 0.27 | 0.04 |
| n=11, 3/class | 51, 17 | 0.02 | 0.35 | 0.96 | 1.00 | 1.00 | 1.00 | 0.06 | 0.04 | 0.34 | 0.41 | 0.06 |
| n=11, <=5/class | 82, 17 | 0.04 | 0.66 | 1.00 | 1.00 | 1.00 | 1.00 | 0.03 | 0.04 | 0.62 | 0.84 | 0.10 |

Full censuses (`skeptic_census10.txt`, 200 reps; `skeptic_census11.txt`, 400 reps; 99 perms; covariates maxC, nnzC):

| census | d=0 | 0.12 | 0.17 | 0.20 | 0.25 | 0.30 |
|---|---|---|---|---|---|---|
| n=10, 140 units, 7 classes | 0.035 | 0.07 | 0.165 | 0.22 | 0.33 | 0.365 |
| n=11, 1361 units, 20 classes | 0.055 | 0.505 | 0.738 | 0.868 | 0.95 | 0.985 |

Confound rates for nrel and minrel at every design and noise level: `skeptic_confound10.txt`, `skeptic_confound11.txt` (minrel stays within 0.01-0.07 except the n = 11 census at R^2 >= 0.5: 0.76/0.17).

n = 9 units (E-283 data, round-105 json): mean defect 2.711 for 3033030 (class 2, 45 members); 2.915, 2.857, 2.375, 2.917 and 2.683 for the class-1 units. The defect-1 fraction is 0 for 3033030 and 4444400, and 0.04-0.27 for the others. The defect-5 count is 0 for 3033030 and 1 for each class-1 unit.

Probe walks (E-283 walk, depth 5; cost sizing, not a census):

| unit (key class size) | members | defects | min | time |
|---|---|---|---|---|
| 30033030 (3) | 120 | 1:16 2:32 3:45 4:16 5:11 | 1 | 155 s |
| 00334400 (66) | 174 | 1:6 2:24 3:77 4:49 5:15 6:3 | 1 | 140 s |
| 03345000 (58) | 138 | 2:29 3:52 4:39 5:17 6:1 | 2 | 90 s |
| 040033030 (n = 11, class 3) | 132 | 2:32 3:35 4:46 5:17 6:2 | 2 | 344 s |

No member has defect <= 0, so H-017's main statement holds on these units. The n = 10 census is 140 units at about 2 min each, about 4.7 CPU-hours: sharded 10-minute commands, or one night. The referee walked the opposite 30330030 of 30033030 and got the same 120 members and histogram, which confirms the duality unit on that pair.

One bug found and fixed in my own script: the first n = 11 run merged two key classes because `hash(-1) == hash(-2)` in Python. Keys are now indexed, not hashed.

## Reproduction

```
.venv/bin/python workshop/rounds/133/skeptic_h017b.py units 10        # 15 s; 'units 11' about 1 min
.venv/bin/python workshop/rounds/133/skeptic_h017b.py n9              # 1 s (prints the 1.10 SD SE)
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_h017b.py power 10 1   # about 2 min
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_h017b.py power 11 1   # about 5 min
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_census.py census 10 200   # 15 s, prints the key table
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_census.py census 11 400   # 50 s
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_census.py confound 10 200 # 30 s; 'confound 11 200' 90 s
.venv/bin/python workshop/rounds/133/skeptic_census.py eta 11         # eta^2 table in (4)
timeout 10m .venv/bin/python workshop/rounds/133/skeptic_h017b.py collect 10 30033030   # 90-155 s per unit
```

## Prior record

- E-283: the defect is not a function of the key (n = 9).
- E-302: "key irrelevant" is withdrawn.
- E-303: the key shifts the member-level distribution (MI) but cannot be separated from the source LNA; it already names one duality twin pair (3505000/5050030, "a reversed-walk twin pair").
- H-017 amendment / E-314 item 19: maverick's H-017b text (rounds/119/maverick.md). Its l.64 already says "the test is for a LARGE effect only" and l.66 "E-303's weak shift is presumably below delta\*". (2)-(3) quantify that concession.

New here: the full n = 9 pairing and the n = 10/11 unit counts; delta\* on the real tables, sampled and full-census; the n = 11 walk cost (one unit); the confound rate as a function of R^2. No E-entry gives n = 10 depth-5 defect data; the 3 probe units are the first, and they are not a census. Nothing in RETRACTIONS bears on this.

## Code changed

None in the library. Scripts: `workshop/rounds/133/skeptic_h017b.py` (docstring fixed, SE printed), `skeptic_census.py` (new). Outputs: `skeptic_power10/11.txt`, `skeptic_census10/11.txt`, `skeptic_confound10/11.txt`, `skeptic_collect11.txt`, `skeptic_defect_*.json`, `skeptic_collect_probe*.txt`. No tests run, since no library file was touched.

## Next

- Chair / next conference ledger: close T6's H-017b form, for (2)-(4). Record delta\*: 1.5 SD (n = 10, N = 23), 1.0 SD (n = 11, N = 51), about 0.20 SD (all 1361 n = 11 units), so a reopening starts from these numbers and from the census cost.
- If the conference wants a T6 run, make it the H-017 main-statement census, not H-017b. OVERNIGHT proposal text, not run:
  - **Command:** at n = 10, `skeptic_h017b.py collect 10 NAME` for all 140 duality units (list from `units 10`), one unit per command, about 4.7 CPU-hours.
  - **Registered:** a member with relations <= cords refutes H-017 at n = 10. "None found" reads "no member within depth 5", and no more (E-071).
  - **Positive control first:** the same walk from 3 LNAs in a quipu class must reach a defect <= 0 member. UNVERIFIED: whether `quipuRelations.certificate` admits relation-free quipus; if not, build the control from a quipu-with-cords algebra of defect <= 0 first.
  - **Secondary, descriptive only:** the per-class unit-mean table, with no key claim.
- Maverick: if H-017b is kept, fix the feature list before the run and include the relation-length sum and count; (4) still applies to the next unlisted feature, most strongly in the census.
