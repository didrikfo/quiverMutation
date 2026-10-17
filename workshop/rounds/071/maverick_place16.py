"""Round 071: place the 16 unplaced n = 11 LNAs (E-162) against the four periodic-indefinite groups (E-190): Coxeter key
equal? signature equal? usage (repo root): timeout 10m .venv/bin/python workshop/rounds/071/maverick_place16.py"""
import numpy as np, collections
from quivermutation import coxeterTables as ct, nakayama as nk
N = 11
names = """334500030 345500030 350500030 303345000 303455000 303505000 304444400 304550400 305504030 305040330 305050030 444440030 455040030 550403030 504033030 505003030""".split()
keys = ct.lnaKeys(N)
def sig(C): ev = np.linalg.eigvalsh(C + C.T); return (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()))
def periodic(Phi, cap=200):
    I = np.eye(len(Phi), dtype=np.int64); M = Phi.copy()
    for m in range(1, cap + 1):
        if (M == I).all(): return m
        M = M @ Phi
    return 0
gk = {}
for rl in nk.allRelationLengths(N):
    C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
    if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
    if periodic(np.rint(-Ci.T @ C).astype(np.int64)): gk.setdefault(keys[rl], [rl, 0, sig(C)])[1] += 1
print("group keys (size, signature (pos,neg,zero) of C+C^T):", sorted((v[1], v[2]) for v in gk.values()))
sixteen = [tuple(int(c) for c in s) for s in names]
ks = {keys[r] for r in sixteen}
print("16 unplaced: distinct Coxeter keys", len(ks), "; any equal to a group key:", any(k in gk for k in ks))
C = np.array(ct.lnaCartanMatrix(N, sixteen[0]), float); Ci = np.linalg.inv(C)
print("signatures of the 16:", collections.Counter(sig(np.array(ct.lnaCartanMatrix(N, r), float)) for r in sixteen))
print("periodic Phi among the 16:", [periodic(np.rint(-np.linalg.inv(np.array(ct.lnaCartanMatrix(N, r), float)).T @ np.array(ct.lnaCartanMatrix(N, r), float)).astype(np.int64)) for r in sixteen][:4], "...")
