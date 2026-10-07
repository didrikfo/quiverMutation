# Experimentalist notebook (rewritten each round)

## What I now believe (after round 038)
- n = 8 c0 walk to 8 465 expansions (560 s cap, level 9 of 10 partial): J_i != 0 rows 243 = (d,J) (2,1) 234, (3,1) 2, (4,1) 4, (5,1) 3; dim J_i = 1 throughout (L1 holds). J_i != 0 with d >= 3: 9 rows, 5 algebras (dim A 47, 47, 64, 64, 75), all out-degree(i) = 3, all at exp >= 7798 (cap tail). (`rounds/038/experimentalist_d3table.py`, `_d3tab.txt`)
- The two (3,1) rows (exp 7798, 7810) have out(i) targets {1,3,4} distinct, out(v) = 2: no parallel arrows needed; E-131 missed them by filtering out(v) >= 3.
- out(i) = 3 is the only non-trivial separator of d = 2 from d >= 3 at J != 0; necessary in the sample, not sufficient (8 of 234 d = 2 rows have out(i) = 3). dim J, dim A, Cartan row sum, out(v) do not separate (dominated by d and depth).
- d >= 3 with J = 0: 135 rows, 59 algebras, out(i) 1..7, so J != 0 is what picks out 3.
- Earlier: a both-die square is gate-admitted J != 0 at n = 6, LNA key 0/42 (n=6), 48 (n=7), 1408 (n=8); capped n = 7 BFS found none of 44 (033). n = 8 has 11 classes, n = 9 has 19. J_i != 0 at n = 9 c0 had d = 2 only (035).

## What I tried
- 038: d3table (full rows saved, c0 only). 035: dhist. 033: bothdie, dimji. 030: BFS depth tables; 027: W walks.

## What I would do next
1. Same script on c1, c2 (`8 1`, `8 2`) and n = 9 c0, one command each; then a deeper/closed n = 8 c0 (overnight): does out(i) = 3 survive; does d >= 3 with out(i) = 2 ever show.
2. Kernel (skeptic_kernel.py) of the (3,1) rows 7798/7810; why out(i) = 3.
3. Close the n = 7 both-die reach question (reverse mutation of the 44).
4. n = 13 lone-3 orbit check (request, S-1).
- Watch: a time cap is not a verdict and counts are load-dependent; d >= 3 lives only in the last 8 % of the walk; algebra ids are keys or hashes (not isomorphism); a filter on out(v) hid rows once; pickle is 3 MB.
