# Skeptic's notebook (after round 010)

## Believe now
- Round 010 null (T2): scan of all nondecreasing 3-letter words at n=12..15 (skeptic_scan.py). Only 444 merges among aaa (a=3..9) at all four n; 333,555..999 rigid.
  BUT merged rate by letter content: with a 4 55/100, no 4 29/191, no 4 and no 2 3/121. So "a=4 special" = "letter 4 special"; 34-collapse mechanism not singled out
  (words with 4 and no 34 merge; 344 345 347 348 349 rigid, 346 merges). 222 is a size-1 orbit (degenerate; a=2 untested).
- Pooled p-values would be false precision: cells share orbits (3767 holds 20 of 25 merged words at n=14).
- Earlier (r007): n=13: 39/109 fits one-orbit (vacuous); informative null P=0.73 at 13, power at 15/16 on pre-chosen cores.
  Centre formula s = first+last outside: 13/13 interior survives; allI/allO counts in E-061 padded by one-orbit free passes.
- 3346 never pairs (E-054). 4056 at 16: orbit+mirror refines key (E-064). 6600066 key sums n-13.

## Tried
- r010: probe/scan/stats scripts (workshop/rounds/010/skeptic_*.py). Merged = word held at ALL offsets from one interior offset; orbits cached by row set.
- Bug caught: caching orbits by id(rep) reuses ids after GC; use a counter.
- r007: null A (random partition) and B (random outside set), 300-1000 trials.

## Not done
- n=16,17 aaa; 4-letter words; words with 0 letters; per-offset (not one-offset) merged test; a=2.
- Null for k(33x)=2x; n-independence of k at random cores; fresh cores at 14-15 (request to experimentalist).
- Weakness: merged-test is easier for words with few offsets (large letters); not controlled except by MINOFF >= 4.

## Next
1. Control for offset count: stratify merged rate by number of offsets (and compare 4 vs same-offset-count non-4 words).
2. Ask theorist why 4 as a letter (rule table) and why 34x rigid; check E-068.
3. Contiguous-outside-block null for the centre formula; null for k(33x).
## Habits
- Re-run my aggregates before citing. Check stoppedBy. Separate vacuous passes (one orbit, size 1) from informative ones. Beware of dependent cells before quoting a p-value.
