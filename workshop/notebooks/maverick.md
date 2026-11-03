# Maverick's notebook

## What I now believe (after round 111)
- S-1 / H-026 K-threshold law now holds at n = 17 key level on the WHOLE lone-3 key class (6,7)/(7,6) restricted to head ends: K=5 two image keys, K=6 two, K>=7 one (65 ends). K0 = 6 as predicted. Method: a row with a head end of K >= K1 free vertices is 0^K1 + (row of length n-K1), so enumerate allRelationLengths(n-K1) (58786 rows for K1 = 5), no 34M-row scan. Control: reproduces E-263 at n = 15 exactly (halved counts).
- Tail trailing-zero rows are NOT a shift (n=12 check: 13260 vs 1430); use head ends + class mirror closure.
- Still key level only; class/orbit level of image keys unknown. Law now 4 lengths (11,13,15,17), no derivation.
- Earlier (r105-109): E-269 within-quipu YES on 4 n=12 low-lam groups, Power none (planted perturbations all congruent); cross-quipu NO controls at lam >= 1.66 only. Convention Phi = -C^-T C, n = len + 2.

## What I tried
- r111: rounds/111/maverick_k17.py (check, run N K1). r109 planted, r107 n12low/nopair, r105 mixed/lowlam, r099 iso140, r097-089 cartanZ/genusscan/quipupairs.

## Watch for
- Venv loses cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX. Per-command 2 min tool wall: use background + until-loop for 6 min jobs.

## Next
- n = 19 K1 = 7 ran (518 s): K=7 two keys, K=8 one, so K0 = 7 as predicted. Next n = 21 needs n'=14 (too big this way).
- Even n (12,14,16): is the failing K0 floor(n/2)-2 (n=12 K0=4 per E-146)? n = 14 K1 = 5 costs n'=9, cheap.
- Planted-NO control for the key test: not found yet.
- Older: in-domain non-congruent C' for E-269; singular S (47-group).
