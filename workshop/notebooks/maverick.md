# Maverick's notebook

## What I now believe (after round 097)
- NEW (r097): the discriminant quadratic form q_S of S = C + C^T (x -> x^T S x mod 2 on S^-1 Z^n / Z^n; same genus as G = C^-1 + C^-T) separates 7 of the 8 profile- and probe-equal certified pairs of E-259 (n = 11: 3/3; n = 12: 4/5). First Cartan-level NO with power. Open: n = 12 key group of 140 (reps (2,2,2,3,0,2,3,0,3,0) vs (2,2,3,0,3,0,2,2,3,0)), q_S equal, pieces ker f(Phi) equal, local module tests pass mod 2,3,5,7, global module iso / isometry undecided.
- q_S splits whole key groups: n = 9: 1/19 (the F-010 pair, calibration), n = 10: 3/40, n = 11: 4/79, n = 12: 15/176 (incl. 3 groups with unnamed members). Constant on verified-move orbits (control, small).
- Class = genus step (Nikulin 1.13.3) is from memory, UNVERIFIED; C-congruence = S-isometry + Phi-equivariance (C = -C^T Phi).
- Older: Cartan level exhausted for periodic key groups at n = 11, 12 (E-232, E-253, E-242); no periodic key group holds two quipu names (E-259). Convention: Phi = -C^-T C, n = relation string length + 2.

## What I tried
- r097: workshop/rounds/097/maverick_cartanZ.py (pairs), _genusscan.py (whole table), _control.py, _orbitcontrol.py, _pair140.py, _pieces.py.
- r094: maverick_quipupairs/quipuprofile/quipucartan. r091 groups/sweep/shift; r089 R_k; earlier rounds T4/H1.

## Watch for
- Venv lacked numpy/sympy: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint` (cypari has matkerint, bnfinit).
- Foreground >120 s goes to background; `pgrep -f <script>` inside an until-loop matches itself. Use `ps aux | grep "[s]cript"`.
- E-245/E-246 string labels: n = len + 2.

## Next
- Open pair (140 group): hermitian-lattice genus over Z[Phi] with PARI; or search the commutant lattice for det +-1 with a smarter (LLL/mod-p lifting) method.
- Skeptic/scholar check of q_S and Nikulin; test q_S on E-117/F-010 merged classes (a split would be a bug).
- Put q_S into the key; recheck which E-259 style "profile-equal" groups remain; n = 13 groups.
