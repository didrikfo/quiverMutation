"""H-017 Euler signature power check (round 050).
Signature (pos,neg,zero) of C+C^T. Since C^-1+C^-T = C^-1 (C+C^T) C^-T (Sylvester), the signature of the Euler
form symmetrisation equals that of C+C^T: a function of the Cartan matrix, and a derived invariant.
Tables, per length n: LNA status (QUIPU / UNPLACED / NOT_QUIPU) x [pos <= n-2]; Coxeter-polynomial groups
split by signature; quipus (hereditary, known class) as controls.
usage: python workshop/rounds/050/maverick_sigpower.py N [N ...]"""
import sys, collections
import numpy as np
from quivermutation import coxeterTables as ct, nakayama as nk, invariants as inv, quipuForms as qf

def sig(C):
    C = np.array(C, dtype=float); ev = np.linalg.eigvalsh(C + C.T)
    return (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()))

for n in map(int, sys.argv[1:]):
    st = ct.lnaStatus(n); keys = ct.lnaKeys(n)
    tab = collections.Counter(); sigs = {}
    for rl in st:
        s = sig(ct.lnaCartanMatrix(n, rl)); sigs[rl] = s
        tab[(st[rl], s[0] <= n - 2)] += 1
    print(f"n={n}: LNAs {len(st)}")
    for k in sorted(tab): print("  status", k[0], "pos<=n-2:", k[1], tab[k])
    # poly groups split by signature
    idx = ct.lnaKeyIndex(n); split = 0; mixed_status = 0; groups = 0
    ex = None
    for key, mem in idx.items():
        if len(mem) < 2: continue
        groups += 1
        ss = {sigs[m] for m in mem}
        if len(ss) > 1:
            split += 1
            if ex is None: ex = [(ct.className(m), sigs[m], st[m]) for m in mem][:4]
    print(f"  poly groups with >=2 LNAs: {groups}; split by signature: {split}; example {ex}")
    # hereditary quipus
    qs = collections.Counter(); qsig = {}
    for p in qf.allQuipusOfOrder(n):
        s = sig(inv.cartanMatrix(nk.QuipuAlgebra(*p)))
        qs[(s[0] >= n - 1)] += 1; qsig[inv.coxeterKey(nk.QuipuAlgebra(*p))] = s
    print("  quipus pos>=n-1:", dict(qs))
    # cross: LNA sharing Coxeter poly with a quipu, signature different from that quipu's?
    diff = same = 0; exs = []
    for rl, key in keys.items():
        if key in qsig:
            if sigs[rl] == qsig[key]: same += 1
            else:
                diff += 1
                if len(exs) < 3: exs.append((ct.className(rl), sigs[rl], qsig[key], st[rl]))
    print(f"  LNA with a quipu's poly: signature equal {same}, different {diff}; e.g. {exs}")
    # unplaced split
    un = [rl for rl in st if st[rl] == ct.UNPLACED]
    print("  UNPLACED:", len(un), "with pos<=n-2:", sum(sigs[r][0] <= n-2 for r in un))

# Failure-to-separate control: quipus sharing a Coxeter polynomial (F-010 cospectral pairs), hereditary of different
# trees; report signature and Smith form of C+C^T
if len(sys.argv) > 1 and sys.argv[-1] == "9":
    import sympy
    from sympy.matrices.normalforms import smith_normal_form
    idx = ct.quipuKeyIndex(9)
    for key, names in idx.items():
        if len(names) > 1:
            out = []
            for p in qf.allQuipusOfOrder(9):
                if qf.formatQuipu(p) in names:
                    C = np.array(inv.cartanMatrix(nk.QuipuAlgebra(*p)), dtype=int); G = C + C.T
                    S = smith_normal_form(sympy.Matrix(G.tolist()), domain=sympy.ZZ)
                    out.append((qf.formatQuipu(p), sig(C), tuple(sorted(abs(int(S[i, i])) for i in range(9)))))
            print("cospectral quipus", out)
