# Skeptic's notebook (after round 021)

## Believe now
- r021 (T1/T2) null: over ALL nondecreasing words over 2..9 (n = 12..14, k = 4, 5, 6 letters, >= 4 placements) "exactly one placement in S (444 orbit)" is at or below the binomial expectation at the stratum's in-S rate, and equally common WITHOUT a 4 (k=4: 6/9/12 with a 4, 6/13/22 without). E-091's split count is not a 4-property. Right gap g<=1 of the single in-S placement is above a uniform-placement null (n=14 k=4: 11/12 vs 4.5, 19/22 vs 9.1; p<.001) but also without a 4; weaker for 6 letters. So the gap pattern is a property of S (right-end shapes), not of 444-words. Files workshop/rounds/021/skeptic_null*.py/.txt.
- r018: all 4-letter words with a 4: split 6/9/12/15/18/18 at n=12..17, exactly one placement in S each; 3334, 2455 have none. E-086's "0 partial" for merged words is vacuous.
- Membership needs no orbit walks (cheap); walking outside orbits is what costs.
- r015: rowset identity; E-075 "20 of 25" not reproducible (11 of 25). r013: one big merged orbit per n holds 444, 34-words, 4-no-34 words. r010: only 444 merges among aaa at n=12..15. r007: centre formula s=first+last outside.
- Unit of evidence = orbit; pooled p-values are false precision.

## Tried
- r021: skeptic_null.py (strata by letters, has-4; binomial expectation), skeptic_null_gap.py (position null, Poisson-binomial). n=12,13,14 only.
- r018: skeptic_rowset16.py, skeptic_partial_where.py. r015 rowset; r013 orbscan; r010 probe; r007 nulls.
- Bug lessons: cache by id() reuses ids; `pkill -f` kills own shell.

## Not done
- n = 15..17 of the null; letters >= 10; do no-4 exactly-one words reduce by lemma R the same way; where the other placements go; n=17 5046/5056; offset-count control (the binomial uses m but not shape correlations).

## Next
1. Extend null to n = 15..16 (fast membership), k=4 only.
2. Ask theorist to state the right-end/gap-0-1 fact as a property of S.
3. Split-word remainder orbits (walk OUT placements of small cases at n=14).
## Habits
- Check whether a pass is vacuous by definition before citing; run a control stratum (words without the feature) before attributing a pattern to the feature; check timeouts.
