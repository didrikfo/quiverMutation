# Experimentalist notebook (rewritten each round)

## What I now believe (after round 035)
- On capped walks (500 s, 4 parallel) at n = 8 c0/c1/c2 and n = 9 c0, every gate-admitted (v, i) with J_i != 0 has d_i = dim e_iAe_v = 2 and dim J_i = 1 (285 of 285). d_i = 3 (27 rows) and d_i = 4 (1 row) occur at admitted vertices but with J_i = 0. d_i >= 3 first shows at BFS level 7-8, the cap's edge, so "d_i <= 2 at J_i != 0" is supported but thinly where it matters. (`rounds/035/experimentalist_dhist.py`)
- A both-die square (p1 b1 = p2 b1, p1 b2, p2 b2 killed) is gate-admitted J != 0 at n = 6; 0 of 42 have an LNA key at n = 6, 48 at n = 7, 1 408 at n = 8; capped n = 7 BFS found none of 44 (033). Squares with an arrow side never get an LNA key (necessary test only).
- dim J_i <= d_i - 1 at admitted v (E-126 L1); dim J_i = 1 on walk rows n = 6..9.
- n = 8 c0 out-degree 3 parallel rows with J != 0 are real Cartan failures (E-127).
- Class sizes n = 8: 11 classes (2, 8, 18, 20, 26, 52, 80, 128, 128, 130, 266); n = 9: 19.

## What I tried
- 035: dhist script, 3 classes at n = 8, one at n = 9. 033: bothdie, dimji. 030: BFS depth tables of max dim; 027: W test walks.

## What I would do next
1. Look at the 28 (d >= 3, J = 0) rows: why the out-arrows jointly separate a 3-dim e_iAe_v.
2. Overnight: closed/deeper BFS histogram, n = 8 idx 0, 1 and n = 9 idx 1..3, run sequentially (parallel halves depth).
3. Close the n = 7 both-die reach question (BFS closure or reverse mutation of the 44).
4. Enumerate both-die cores more widely; positive control for W (parallel b1).
- Watch: a cap is not a verdict; an LNA key is a necessary test only; parallel runs inflate timings; seen != expanded; rows counted as (alg, v, i) vs (alg, v) differ between rounds.
