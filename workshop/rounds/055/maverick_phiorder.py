"""Exact order of Phi (Phi^k = I) for each cyclotomic key group with >= 2 classes.  usage: python workshop/rounds/055/maverick_phiorder.py N"""
import sys, collections
sys.path.insert(0, "workshop/rounds/054"); sys.path.insert(0, "workshop/rounds/029")
import sympy, maverick_pq as mp
from quivermutation import coxeterTables as ct
n = int(sys.argv[1]); keys, comps, *_ = mp.classesOf(n)
by = collections.defaultdict(list)
for c, mem in comps.items(): by[keys[mem[0]]].append(mem)
for k, cl in by.items():
    if len(cl) < 2: continue
    C = sympy.Matrix(ct.lnaCartanMatrix(n, cl[0][0])); Phi = -(C.inv().T) * C
    P = sympy.eye(n)
    for e in range(1, 61):
        P = P * Phi
        if P == sympy.eye(n): print("classes", len(cl), "Phi^%d = I" % e); break
    else: print("classes", len(cl), "no Phi^k = I, k <= 60")
