"""Round 050 response: list the LNAs whose Coxeter key equals a quipu's but whose signature differs from EVERY
quipu sharing that key; check uniqueness of quipu per key / signature per key; test F-048 (periodic Phi + indefinite
Euler form) on them and count how many F-048 rows exist. usage: python maverick_sixteen.py N [N ...]"""
import sys, collections
import numpy as np
from quivermutation import coxeterTables as ct, nakayama as nk, invariants as inv, quipuForms as qf

def sig(C):
    ev = np.linalg.eigvalsh(np.array(C, float) + np.array(C, float).T)
    return (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()), int((abs(ev) <= 1e-9).sum()))

def periodic(C, cap=600):
    C = np.array(C, float); Phi = np.rint(-C.T @ np.linalg.inv(C)).astype(object)
    I = np.eye(len(C), dtype=int).astype(object); M = Phi
    for m in range(1, cap + 1):
        if (M == I).all(): return m
        M = M.dot(Phi)
        if max(abs(x) for x in M.flat) > 10**6: return 0
    return 0

def indefinite(C):
    C = np.array(C, float); X = np.linalg.inv(C); return np.linalg.eigvalsh(X + X.T).min() < -1e-9

for n in map(int, sys.argv[1:]):
    st = ct.lnaStatus(n); keys = ct.lnaKeys(n)
    qk = collections.defaultdict(list)
    for p in qf.allQuipusOfOrder(n):
        qk[inv.coxeterKey(nk.QuipuAlgebra(*p))].append((qf.formatQuipu(p), sig(inv.cartanMatrix(nk.QuipuAlgebra(*p)))))
    multi = {k: v for k, v in qk.items() if len({s for _, s in v}) > 1}
    print(f"n={n}: quipu keys {len(qk)}; keys with >1 quipu {sum(len(v)>1 for v in qk.values())}; keys whose quipus differ in signature {len(multi)}")
    f48 = set(); sigs = {}
    for rl in st:
        C = ct.lnaCartanMatrix(n, rl); sigs[rl] = sig(C)
        if n >= 10 and sigs[rl][1] > 0 and periodic(C) and indefinite(C): f48.add(rl)
    diff_all = []; diff_any = 0; shared = 0
    for rl, key in keys.items():
        if key in qk:
            ss = {s for _, s in qk[key]}
            if len(qk[key]) > 1: shared += 1
            if sigs[rl] not in ss: diff_all.append(rl)
            if sigs[rl] != qk[key][0][1] or len(ss) > 1: diff_any += 1
    print(f"  LNAs with a quipu key {sum(k in qk for k in keys.values())} (sharing key with >1 quipu: {shared}); signature differs from ALL sharers: {len(diff_all)}; differs from some sharer: {diff_any}")
    print(f"  F-048 rows (periodic, indefinite): {len(f48)}; of the {len(diff_all)} listed, in F-048: {sum(r in f48 for r in diff_all)}")
    for rl in diff_all:
        print("   ", ct.className(rl), st[rl], sigs[rl], [(a, s) for a, s in qk[keys[rl]]], "F048" if rl in f48 else "not-F048")
