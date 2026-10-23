"""Round 082 scholar: for the n = 11 periodic-Phi indefinite LNAs (E-190: groups 535/64/24/15) factor the Coxeter
polynomial into cyclotomics (1310.1557 Thm 1, 2.10), check period = lcm{m} (not lcm{phi(m)}; the paper's misprint),
list the (a, parity b) allowed by Phi^a = (-1)^(a+b) I (S = -Phi on K0, Lemma 2.9 refined), and flag radical powers.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/082/scholar_cyclo.py [N]   (about 1 min)"""
import sys, collections, math
import numpy as np, sympy as sp
from quivermutation import coxeterTables as ct, nakayama as nk
N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
T = sp.symbols('T')
def order(Phi, cap=400):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return m
        M = M @ Phi
    return 0
groups = collections.defaultdict(list)
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float)
    Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    Phi = np.rint(-Ci.T @ C).astype(np.int64)
    o = order(Phi)
    if o == 0: continue
    chi = sp.Poly(sp.Matrix(Phi.tolist()).charpoly(T).as_expr(), T)
    groups[tuple(chi.all_coeffs())].append((rl, o, Phi))
print(f"n={N}: {len(groups)} Coxeter polynomials among periodic indefinite LNAs")
for co, mem in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    chi = sp.Poly(list(co), T)
    fac = sp.factor_list(chi.as_expr())[1]
    idx = []
    for f, e in fac:
        for m in range(1, 400):
            if sp.Poly(sp.cyclotomic_poly(m, T), T) == sp.Poly(f, T): idx.append((m, e)); break
        else: idx.append(('NOT CYCLOTOMIC', e))
    L = math.lcm(*[m for m, _ in idx if isinstance(m, int)])
    phi_l = math.lcm(*[int(sp.totient(m)) for m, _ in idx if isinstance(m, int)])
    o = mem[0][1]; I = np.eye(N, dtype=np.int64)
    M = np.eye(N, dtype=np.int64); allowed = []
    for a in range(1, o + 1):
        M = M @ mem[0][2]
        if (M == I).all(): allowed.append((a, 'a+b even'))
        elif (M == -I).all(): allowed.append((a, 'a+b odd'))
    print(f"group size {len(mem)}: chi = prod Phi_m^e {idx}; ord Phi = {o}; lcm m = {L}; lcm phi(m) = {phi_l}; "
          f"Phi^a = +-I at {allowed}")
# radical powers kA_N/rad^L (digit string str(L)*(N-L)+"0"*(L-2)) among the periodic indefinite rows
rows = {tuple(rl) if not isinstance(rl, str) else rl: co for co, mem in groups.items() for rl, _, _ in mem}
print("radical powers among periodic indefinite rows:",
      [L for L in range(3, N + 1) if any(ct.className(r) == str(L if L < 10 else '?') * (N - L) + '0' * (L - 2) for r in rows)] if N else None)
for co, mem in groups.items():
    hit = [L for L in range(3, min(N, 9) + 1) for r, _, _ in mem if ct.className(r) == str(L) * (N - L) + '0' * (L - 2)]
    if hit: print("group of", len(mem), "contains radical powers L =", hit)
