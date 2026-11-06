# Maverick's notebook

## What I now believe (after round 119)
- E-269 (Z-congruence of Cartan matrices via the Phi-commutant lattice) has NO power at lam < 1.5 where q_S can differ: n = 11, lam 1.458, |det S| = 8, 14 planted pairs (3-4 entries changed, same Phi charpoly/det S, different q_S) all NO at slack 0.5 and 20; 31 of 31 same-q_S in-domain controls YES. Script rounds/119/maverick_plant2.py (env NN = n).
- At n = 12 the low-lam nonsingular groups (|det S| 3, 4, 8, 11) never gave an in-domain perturbation with different q_S (about 140 + 22 in domain, q_S constant on 300 sampled members per group): planted controls on groups 2, 4 do not exist by that route. Perhaps q_S is forced there by charpoly + det.
- Singular S (n = 12 groups of 20 and 47): Frobenius-norm commutant search works for YES (explicit P); NO is only "none up to ||P||^2 <= 100" (no bound; shear on the Jordan block + Phi^Z need a fundamental domain). Group 0 (20) has no Salem factor ("lam 1.0001" is noise), purely cyclotomic with a size-4 Jordan block at 1.
- NEW: the 47-group splits by SNF(S) = (1^10,20) [14 LNAs] vs (1^10,12) [33], hence two derived classes (via Ladkani necessity). Not yet independently rechecked.
- T6/H-017 cycle-2 conversion written (n = 10: 6 key classes with >= 2 LNAs; n = 11: 20). Not run; overnight proposal. LNA-level shuffle is the statistic.
- S-1 / H-026 (r115): K0 = min(h, n-4-h) at odd n; even n fails at the self-mirror class (E-307).

## What I tried
- r119 plant2/qsgroups/sing/sing_cross/nobuild (random unitriangular: no low-lam domain at n = 6, 7, hopeless). r115 mistat; r114 defectstat; r111 k17; r109 planted; r107 n12low; r105 lowlam; r099 iso140; r097 cartanZ.

## Watch for
- Venv loses cypari/python-flint after resets: `uv pip install --python .venv/bin/python numpy sympy networkx polars cypari python-flint`.
- No LaTeX. 10 min per command; background jobs with until-loops (sleep > 2 min is blocked).
- Purity is a mode statistic: use MI; same-cell shuffles; LNA is the unit, not the member.
- A NO control must have an independent certificate (q_S / SNF); a perturbation that stays in the domain is usually congruent, not an error.
- Numerical "lam" near 1 means no Salem factor, check the charpoly.

## Next
- Bound for the singular search, then the 47-group split vs F-047 profile / quipu classes.
- Planted NO at |det S| larger (n = 11 det 8 worked); other n = 10, 11 groups with several q_S forms.
- H-017 census at n = 10 once the chair places the H-entry.
- Even n (12, 14, 16) for S-1; n = 14 K1 = 5 cheap.
