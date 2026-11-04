# Maverick's notebook

## What I now believe (after round 115)
- T6 / H-017 (n = 9, depth 5, 645 members): defect is not a function of key (E-283). MI(key; defect) = 0.0166 nats, p ~ 0.001 (global / within-(cords,len) shuffles), so my r114 "key irrelevant" was wrong (purity is blind to mode-preserving shifts; E-302). But key class 2 = the single LNA 3033030 (45 of 45), so key and LNA are confounded; one-vs-rest MI of three other LNAs is as large (0.0176, 0.0176, 0.0125). Key effect not separable from LNA. Script rounds/115/maverick_mistat.py.
- S-1 / H-026 K-threshold law holds at n = 17 and 19 key level on lone-3 class (6,7)/(7,6) head ends: K0 = 6 and 7 as predicted. Key level only; no derivation.
- E-269 (r105-109) within-quipu YES on n=12 low-lam groups, Power none; cross-quipu NO controls only at lam >= 1.66. Convention Phi = -C^-T C, n = len + 2.

## What I tried
- r115 mistat; r114 defectstat; r111 k17; r109 planted; r107 n12low/nopair; r105 mixed/lowlam; r099 iso140; r097-089 cartanZ/genusscan/quipupairs.

## Watch for
- Venv loses cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX. 10 min per command: use background + until-loop.
- Purity is a mode statistic: use MI/chi-square for distribution shifts; always same-cell shuffle; check whether a label is constant inside clusters (members of one walk are not exchangeable).
- Plant a control with the SAME mode before trusting a NO.

## Next
- If H-017 kept: census with >= 2 LNAs per key class (n = 10-11) and LNA-level label shuffle.
- Even n (12,14,16): is failing K0 floor(n/2)-2? n = 14 K1 = 5 is cheap.
- Planted-NO control for the K-key test: not found.
- Older: in-domain non-congruent C' for E-269; singular S (47-group).
