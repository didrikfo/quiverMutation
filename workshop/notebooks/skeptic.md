# Skeptic's notebook (after round 018)

## Believe now
- r018 (T1/T2): over ALL 4-letter words with a 4 (>= 4 placements) the "all in S or none" claim is false: split words 6/9/12/15/18/18 at n = 12..17 (n = 16: 18 of 120), each with exactly one placement in the 444 orbit S (last or second-to-last offset mostly; 3344 at 1). All-IN 5/10/13/16/19/19. 3334, 2455 have 0 placements in S at n = 12..17. Files workshop/rounds/018/skeptic_rowset16*.txt, skeptic_partial_where_n*.txt.
- r018 meta: E-086's "0 partial" for MERGED words is vacuous (orbits are disjoint, merged = one orbit). Only the IN count has content. My r015 line repeated it; do not cite as a test.
- Membership test needs no orbit walks (12 s at n=16); walking outside orbits is what costs (>10 min at n=16, big orbits 1e4-3e5).
- r015: rowset identity; E-075 "20 of 25" not reproducible (11 of 25).
- r013: one big merged orbit per n holds 444, 34-words, 4-no-34 words; letter 4 vs collapse-to-34 inseparable. Other merged orbits small, 2-driven, or 568/679.
- r010 null: only 444 merges among aaa at n=12..15. r007: centre formula s = first+last outside, 13/13 interior (one-orbit fits vacuous at n=13).
- Unit of evidence = orbit; pooled p-values are false precision.

## Tried
- r018: skeptic_rowset16.py (membership in S for all 4-letter words containing a 4, n=12..17, LIMIT 0 fast mode), skeptic_partial_where.py. Walked n=16 run timed out at word 3466 (tags agree with fast run).
- r015 rowset; r013 orbscan/orbstats; r010 probe/scan; r007 nulls. Bug lessons: cache by id() reuses ids; `pkill -f` pattern kills own shell.

## Not done
- Where the other placements of split words go (one second orbit or many); why exactly one in S; is all-IN = "collapse path to 34".
- n=17 5046/5056 row sets; words with zeros; a=2; offset-count control; why {2aa,4aa}; why 457 alone; letters >= 6 for 3334-like words.

## Next
1. Split-word remainder orbits (walk only the OUT placements of the 18 split words at n = 14, small orbits first).
2. Offset-count control for "merged is easier with few offsets".
3. Ask theorist: single in-S placement of rigid words via lemma R.
## Habits
- Check whether a pass is vacuous by definition before citing; re-run aggregates; check stoppedBy/timeouts; separate vacuous from informative passes.
