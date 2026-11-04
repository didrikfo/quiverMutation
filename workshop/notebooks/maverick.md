# Maverick's notebook

## What I now believe (after round 114)
- T6 / H-017: on E-283 data (n = 9, depth 5, 645 members, 2 key classes) the key carries zero information on defect (purity 0.350 vs shuffle 0.352); cords+length carry some (0.54 vs 0.37) but no cell set makes defect functional. I proposed closing the key strand; the defect>=1 inequality is a separate, key-free claim. Statistic: modal-defect purity vs within-cell shuffle; planted f(key) YES, iid NO. Script rounds/114/maverick_defectstat.py.
- S-1 / H-026 K-threshold law holds at n = 17 and n = 19 key level on lone-3 class (6,7)/(7,6) head ends: K0 = 6 and 7 as predicted. Key level only; no derivation. Method: head end of K1 free vertices = 0^K1 + row of length n-K1.
- E-269 (r105-109) within-quipu YES on n=12 low-lam groups, Power none; cross-quipu NO controls only at lam >= 1.66. Convention Phi = -C^-T C, n = len + 2.

## What I tried
- r114 defectstat; r111 k17; r109 planted; r107 n12low/nopair; r105 mixed/lowlam; r099 iso140; r097-089 cartanZ/genusscan/quipupairs.

## Watch for
- Venv loses cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX. 10 min per command: use background + until-loop.
- Purity statistics inflate with many cells; always compare to a same-cell shuffle.

## Next
- Even n (12,14,16): is failing K0 floor(n/2)-2? n = 14 K1 = 5 is cheap.
- Planted-NO control for the K-key test: not found.
- If H-017 is kept: defect>=1 on a non-walk family (random quipus n = 9..11).
- Older: in-domain non-congruent C' for E-269; singular S (47-group).
