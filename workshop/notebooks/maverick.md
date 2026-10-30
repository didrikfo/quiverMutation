# Maverick's notebook

## What I now believe (after round 099)
- NEW (r099): the n = 12 key group of 140 pair is NOT Z-congruent as Cartan matrices (C_a vs C_b and vs C_b^T), by an exhaustive finite search: C = S(1-Phi^-1)^-1, so P C_a P^T = C_b <=> R = P^T in the rank-14 commutant lattice {R Phi_b = Phi_a R} with R^T S_a R = S_b. S has signature (11,1), indefinite part = Salem plane, so isometries = finite x Phi^Z; majorant H = S + (Sd)(Sd)^T, bound 13.1 (checked to 42.6). 9/9 cross pairs NO, 6/6 within-quipu pairs YES, random-P control YES. So all 8 E-259 pairs are now Cartan-separated (7 by q_S, this one by Phi-equivariance). Least sure: float bound (skeptic should re-derive).
- Not claimed: S_a, S_b non-isometric (same genus; probably isometric, class number argument unverified). The obstruction is Phi-compatibility.
- r097: q_S (discriminant form of S = C + C^T) separates 7 of 8 pairs; splits key groups n = 9: 1/19, n = 10: 3/40, n = 11: 4/79, n = 12: 15/176.
- Class = genus (Nikulin 1.13.3) still UNVERIFIED (E-267). C-congruence = S-isometry + Phi-equivariance (C = -C^T Phi).
- Older: Cartan level exhausted for periodic key groups n = 11, 12 (E-232, E-253, E-242); no periodic key group holds two quipu names (E-259). Convention Phi = -C^-T C, n = relation string length + 2.

## What I tried
- r099: workshop/rounds/099/maverick_iso140.py (test(a,b), table mode), _ctl.py, _transp.py, _setup.py. First version had the commutant equation transposed (kron order) and gave false NOs on controls: ALWAYS run a within-class control with a nontrivial P before trusting a NO. qfminim: count res[0] = 2 x vectors returned; set pari.allocatemem(3e9).
- r097: maverick_cartanZ, _genusscan, _control, _orbitcontrol, _pair140, _pieces. r094 quipupairs/profile/cartan. r091, r089, earlier T4/H1.

## Watch for
- Venv lacked cypari/python-flint again: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- Foreground >120 s goes to background; `pgrep -f` in until-loop matches itself, use `ps aux | grep "[s]cript"`.
- E-245/E-246 string labels: n = len + 2.

## Next
- Run `test` over all q_S-equal pairs / key groups with signature (n-1,1) at n = 11..13 (the method needs exactly one Salem pair; for more indefinite directions the finite-times-Z^r argument needs a rank-r balance).
- Ask: does Phi-equivariance (module class of Z^12 over Z[Phi] with its S) give a cheap invariant that predicts the 9/0 split, e.g. the hermitian class on the K8-piece? Then it could join the key.
- Skeptic/scholar: verify the bound; verify Nikulin; test q_S + this on E-117/F-010 merged classes (a split would be a bug).
