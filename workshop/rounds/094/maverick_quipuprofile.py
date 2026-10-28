"""Round 094 maverick: for each key group with >= 2 distinct quipus, F-047 profile (SNF of g(Phi) per irreducible factor g)
of up to 3 members per quipu; does the profile separate the quipus (=> Cartan-incongruent) or not (=> Cartan-congruent-so-far)?
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/094/maverick_quipuprofile.py N"""
import sys, collections
import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
N = int(sys.argv[1]); keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
for rl, k in keys.items(): groups[k].append(rl)
x = sp.symbols('x')
def profile(rl):
    C = sp.Matrix(ct.lnaCartanMatrix(N, rl).tolist() if hasattr(ct.lnaCartanMatrix(N, rl),'tolist') else ct.lnaCartanMatrix(N, rl))
    Phi = -(C.T).inv() * C   # Phi = -C^-T C ; Z-conj class same as the repo convention up to transpose/inverse
    cp = sp.Poly(Phi.charpoly(x).as_expr(), x)
    out = []
    for g, e in sp.factor_list(cp.as_expr())[1]:
        gm = sp.zeros(N, N); P = sp.eye(N)
        for c in reversed(sp.Poly(g, x).all_coeffs()): gm += c * P; P = P * Phi
        s = smith_normal_form(gm, domain=sp.ZZ)
        out.append((str(g), tuple(sorted(abs(int(s[i, i])) for i in range(N)))))
    return tuple(sorted(out))
for k, v in groups.items():
    q = collections.defaultdict(list)
    for rl in v:
        r = quipuForAlmostSeparateLNA(N, list(rl))
        if r is not None: q[r].append(rl)
    if len(q) < 2: continue
    prof = {r: {profile(rl) for rl in m[:3]} for r, m in q.items()}
    allp = list(prof.values())
    split = len({frozenset(p) for p in allp}) > 1 or any(len(p) > 1 for p in allp)
    print("group size", len(v), "quipus", list(q), "profiles/quipu", [len(p) for p in allp],
          "distinct profile sets", len({frozenset(p) for p in allp}), "-> SEPARATED" if len({frozenset(p) for p in allp}) == len(allp) and not any(len(p)>1 for p in allp) else "-> NOT fully separated")
