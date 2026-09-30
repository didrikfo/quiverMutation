# Most n = 13 "reflection fits" are vacuous or chance-level; the centre formula s = first + last outside offset on the 13 interior cores survives the null (naive product 1e-7, effective about 1e-3 to 1e-4 once correlated cores are collapsed), the rest of the claimed fits mostly do not

author: skeptic · round: 007 · kind: negative (partial: the core claim survives)
thread: T2/T6 · bears on: H-021, E-060, E-061, E-062, F-053

## Claim

Null = keep each core's offsets and its orbit block sizes (from the committed census), reassign offsets to blocks uniformly at random, and run the project's own `fit()` (`d <= 6`, nontrivial pair). Under it:

1. **Existence of a fit carries little information at n = 13.** Of the 109 cores with a fit, 39 have a single orbit holding every offset: any centre closes, the fit is vacuous (`d = 0`, P(null fit) = 1.00). Of the 70 informative fits (>= 2 orbits), a shuffled partition fits with mean probability 0.73 (51 of 70 expected by chance); only 7 of 70 have P(null reaches `d <= d_obs`) < 0.05, 32 of 70 < 0.20. So at n = 13, at most about 10 percent of claimed fits individually survive at 5 percent, about 45 percent at 20 percent.
2. **The same test has power at n = 15, 16 for the 12 E-060 cores.** P(null fit) is 0.42 (n = 15) and 0.22 (n = 16) on average; 10/12 (n = 15) and 9/9 (n = 16, the 3 unmerged cores excluded by `fit`) have P(null `d <= d_obs`) < 0.05. The n-independence of `k, d` is therefore not a low-power artefact at 15 and 16. Caveat: these 12 cores were chosen after fitting at 13 (E-060), so this shows the fit is real for them, not that a random core has one.
3. **The centre formula itself survives.** For the 13 interior-block cores (E-060: 13/13), a random outside set of the same size gives `first + last = s_fitted` with probability 0.19-0.49 per core (expected 4.1 hits of 13); the naive product of the 13 P values is about 1e-7 but assumes independent cores; the 13 include one-parameter families of the same shape (`505/555/605`, `566/606/666`, `504/6004`, `45/46/56`), so the effective number of tests is about 5-6 and the joint p is about 1e-3 to 1e-4 (theorist's review). The interior class was itself selected after fitting at n = 13 (E-060); the null holds the class fixed and randomises the outside set. Redrawing the orbit partition (null A, same slide) and refitting hits `pred` with probability 0.01-0.08 (one core 0.33), expected 0.9 of 13, naive product about 1e-18 (same independence caveat). Classes where the formula fails or is vacuous:

| class (n = 13) | informative cores | formula hits | hits expected by chance (random outside set) | fails |
|---|---|---|---|---|
| int | 13 | 13 | 4.1 | none |
| endhi | 11 | 10 | 1.6 | 4506 |
| end0 | 10 | 7 | 1.0 | 4045 3556 4556 |
| allO | 10 | 7 | null B: 7.0 (prediction deterministic there); null A: about 2.3, so the 7/10 is a null-A statement | 4046 5046 5056 |
| allI | 26 | 9 | (vacuous: pred = hi) | 17 |
| allI / allO, one orbit | 39 | 39 | vacuous | - |

The E-061 count "allI 45/62, allO 10/13" includes the 36 + 3 = 39 one-orbit cores, which pass for free: the informative counts are 9/26 and 7/10. The E-061 totals "10/13 allO, 45/62 allI" should not be read as support. The end-touching classes (end0, endhi) beat chance (17/21 against 2.6 expected), but with `m` only 2-4 the per-core P is 0.2-0.5.

Does not claim: anything about `k(c)` (the formula `s = n - k`; only the slide link and existence of a fit were tested); anything for the 127 cores not at 15/16; that the null is the right one (limits below).

## Evidence

- n = 13 census: 139 cores, 109 fits (70 informative, 39 one-orbit), `skeptic_null.py 300` (shuffle trials per core; 300 for existence/`d`, per-core P resolution 0.003). Full per-core table for n = 13, 15, 16: `skeptic_null_out.txt`; `skeptic_null2_n13.txt` holds the class table and the 13 interior P values (1000 trials).
- n = 15, 16: same script on `experimentalist_census_12cores_n15/n16.jsonl`: table above (2).
- Joint P for interior cores: product of the per-core values (`int` row; 13 values listed in `skeptic_null2_n13.txt`).

## Limits (mine, stated first)

- The uniform-label null ignores that real orbits are structured (neighbouring offsets tend to share orbits), so for real cores chance closure could be higher than the null says; it is the simplest null, not the only one. A contiguous-block null was tried informally (P about 0.3-0.5 per core) but no script is committed, so that figure is not reproducible and is not relied on; the neighbour-aware null is not implemented.
- Shuffling within the given partition sizes drops the orbit sizes (2386 etc.), which are not used by `fit`.
- Null B removes the slide-orbit link but keeps `s`; null A keeps the slide and breaks the orbits. Neither tests n-independence of `k` across n at random cores (only the 12 chosen).
- n = 13 power is low because `|R| <= 4` fits of a single adjacent pair (E-061 says the same): this is now quantified, 73 percent chance fit.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/007/skeptic_null.py 300        # under 1 s (output kept in `skeptic_null_out.txt`); existence and d, n = 13, 15, 16
timeout 10m .venv/bin/python workshop/rounds/007/skeptic_null2.py 1000      # a few seconds; centre formula, n = 13
```

## Prior record

E-061 says "`|R| <= 4` fits (one adjacent pair) have almost no power" and E-061/E-060 list "a null test" as open (STATE T2). grep of `research/` for "null" finds no recorded test. Not a rediscovery. Nothing in RETRACTIONS.md concerns it.

## Code changed

None in the library. Two new scripts in `workshop/rounds/007/`. `skeptic_null.py` imports `fit` from `workshop/rounds/002/experimentalist_fit.py`. No tests touched.

## Next

- theorist: restate the E-060 result as "13/13 interior informative" and the E-061 class counts without the 39 one-orbit cores; they are tautologies.
- experimentalist: test `s = first + last outside` on a random sample of cores at 14-15 (not the pre-selected 12); these 13 were found after looking.
- skeptic (me): null for `k(33x) = 2x` (E-061) and the `k` agreement across n.
