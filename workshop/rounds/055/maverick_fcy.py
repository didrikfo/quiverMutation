"""Powering pre-check for the candidate invariant 'Serre functor periodicity on objects' (round 055).  The candidate is non-vacuous only
if the Coxeter polynomial is a product of cyclotomic polynomials (Phi has finite order, K_0-level fractional CY).  For each key group
(= Coxeter polynomial) holding >= 2 orbit+mirror classes, report: cyclotomic?, Phi order, F-047 profile separates the classes?
usage: python workshop/rounds/055/maverick_fcy.py N [N ...]"""
import sys, collections, math
sys.path.insert(0, "workshop/rounds/054"); sys.path.insert(0, "workshop/rounds/029")
import sympy
import maverick_pq as mp
import toolsmith_snfresolve as sr
from quivermutation import coxeterTables as ct
T = sympy.symbols('T')
def cyclo(poly):
    """order of Phi if poly is a product of cyclotomics (sympy factor list), else None"""
    L = 1
    for g, e in sympy.factor_list(poly)[1]:
        for k in range(1, 400):
            if sympy.Poly(sympy.cyclotomic_poly(k, T), T) == sympy.Poly(g, T) or sympy.Poly(sympy.cyclotomic_poly(k, T), T) == -sympy.Poly(g, T):
                L = L * k // math.gcd(L, k); break
        else: return None
    return L
for n in [int(a) for a in sys.argv[1:]]:
    keys, comps, oid, norb, _ = mp.classesOf(n)
    bykey = collections.defaultdict(list)
    for c, mem in comps.items(): bykey[keys[mem[0]]].append(mem)
    rows = []
    for k, cl in bykey.items():
        if len(cl) < 2: continue
        sets = [{sr.profile(n, m) for m in mem} for mem in cl]
        sep = any(not (s & t) for s in sets for t in sets if s is not t)
        m0 = cl[0][0]
        C = sympy.Matrix(ct.lnaCartanMatrix(n, m0)); Phi = -(C.inv().T) * C
        cp = Phi.charpoly(T).as_expr()
        rows.append((len(cl), sep, cyclo(cp), str(sympy.factor(cp))[:60]))
    print(f"n={n}: key groups with >=2 classes {len(rows)}; profile-separated {sum(r[1] for r in rows)}; cyclotomic (Phi finite order) {sum(r[2] is not None for r in rows)}")
    for r in rows: print("   classes", r[0], "profile-separates", r[1], "Phi order", r[2], r[3])
