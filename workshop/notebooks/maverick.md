# Maverick's notebook

## What I now believe (after round 135)
- **T3/T8 H-025 (r135): E-346's two-shift certificate survives an independent nu.** `135/maverick_nu2.py` works by termwise nu + cone replacement (no bimodule resolution). It agrees with `089/maverick_bimod.nu_bimod` on every dim H^t(nu^k X)_v: N_15(3), N_15(8) for k <= 144 and N_14(5) for k <= 60, on all S and P, with 0 disagreements. It also reproduces every uniserial row of E-346 (N_14(4), N_14(6) at m = 132; 050404400, 040404000, 45050400). These runs are over GF(2^31-1). GF(2) is done for N_15(3) and N_14(5) on S and P; GF(3) and GF(1000003) only on the n = 15 certificate objects. The Serre-duality identity holds at every step.
- Wrong-nu lesson: termwise object errors die on the module-map or d^2 asserts. Scalar errors are either a gauge change (rad^3: tree-like complexes) or break d^2 (rad^5). The only wrong nu that passes the asserts is a coherent autoequivalence (the identity, say). It gives a false "consistent with FCY" and only the Serre check catches it. So the Serre check is the real acceptance test.
- T5 (r131): the gate's cycle failure cannot be reached from LNAs but is real for random cyclic monomial algebras (39 over-admits / 18 145). `isTiltingAtVertex` over-admits 2 951 times. Recommended an R at library scope.
- T10, pure F (r129, E-332): the F-only key-guarded walk crosses at n = 11 depth 6; `refuseNonTilting` removes it. R inverts F only up to the ideal; compare with `sameIdeal` (129/maverick_inv.py), not with struct iso.
- E-318 pair: X_a = coker(Phi - a); K0 = min(h, n-4-h) at odd n (S-1).

## What I tried
- r135: maverick_nu2.py (REL m objs [--prime p] [--serre] [--plant=scale|kup|trunc|long|id]); maverick_compare.py (vs nu_bimod). Costs: N_15(3) 144 S,P is 43 s (nu_bimod 5 s); 45050400 84 U is 556 s (near the cap); N_14(6) 132 U is 275 s.
- r131: 131/maverick_gatehom.py. r129: maverick_inv.py, maverick_libwalk.py.
- Hom engine for cyclic MONOMIAL algebras: `skeptic_conv3.AlgG(n, arrows, [{path: 1}], False)`.

## Watch for
- **numpy int64 mod 2^31-1: matrix products overflow** (a sum of products of two residues > 2^63). It gave me a false "d^2 != 0" at k of about 30. Split-multiply (`mulmod` in maverick_nu2). Check other workshop numpy code for the same bug.
- A background job re-reads the script for each command: do not edit a script while a loop of runs uses it.
- `grep | tail` in a background loop buffers until each command ends.
- Venv may lose python-flint/cypari (`uv pip install --python .venv/bin/python cypari python-flint`); flint was absent in r135.
- Exec-importing scripts: split on "\nif __name__ == '__main__':". `/usr/bin/time` is absent.
- A NO needs a control that says NO; a "0 found" needs the detector shown firing elsewhere.

## Next
- H-025: the referee/chair may lift E-346's blocker. Optional: run uniserials at GF(2) with maverick_nu2. Also N_16(8) ⟨2,3,11] (not yet run by anyone) with both codes.
- Idea: does the Serre identity on Hom-dims, checked on all S_l and P_l at every step, pin nu up to isomorphism on these orbits? If so it is a nu-free certificate checker.
- Theorist question (r134): why two nu-orbit families with different local CY dimensions. Reframe as tau-orbits in the AR quiver of D^b: which uniserials sit in the minority family?
- Still open: gate grading idea (r131); n = 11 depth 7 13-class (r129).
