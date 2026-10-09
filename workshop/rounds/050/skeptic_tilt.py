"""Round 050 (skeptic), T10 (i): INDEPENDENT tilting test for one vertex step, not using tiltingPlus / perI / the gate.
For an acyclic bound quiver algebra A (repo PathAlgebra) and a vertex k, build the complex T = (+_{i!=k} P_i) + T_k,
T_k = ( P_k --(right mult by b)--> +_{b: k->h} P_h ), placed in degrees (S, S+1) (S = 0 or -1: both tried), P_x = A e_x (left modules),
Hom(P_x,P_y) = paths x~>y mod I; composition g o f = path f then path g.  (The orientation 'P_h -> P_k, Hom = paths y~>x' is NOT tilting for the
repo's F step already on A3, see selftest; this left-module convention is the one tiltingPlus describes: p -> (p beta) injective.)
Computes, by linear algebra over Q (rank mod the prime 2^31-1), dim Hom_K(T_i, T_j[m]) for m = -1, 0, 1 (chain maps mod homotopy).
T is a tilting complex iff Hom(T,T[m]) = 0 for m != 0 (generation of K^b(proj A) is Okuyama-Rickard, NOT tested here).
Then End(T) has Cartan matrix H[i][j] = dim Hom(T_i,T_j); compare with the Cartan matrix of the claimed child (labelled).
Only dimension of Hom spaces is compared with the child: End(T) ~ child as algebras is NOT tested."""
import sys
sys.path.insert(0, '.')
from fractions import Fraction
from quivermutation import arrowPaths as ap, procedure

P = 2147483647
class Alg:
    def __init__(self, alg):
        self.q = alg.quiver; self.rels = procedure.relationsFrom(alg); self._piv = {}; self._bas = {}
    def piv(self, i, j):
        if (i, j) not in self._piv: self._piv[(i, j)] = ap.idealBasis(self.q, self.rels, i, j)
        return self._piv[(i, j)]
    def basis(self, i, j):          # basis of paths i ~> j modulo I (non-pivot paths)
        if (i, j) not in self._bas:
            pv = self.piv(i, j) if i != j else {}
            self._bas[(i, j)] = [p for p in ap.allPathsBetween(self.q, i, j) if p not in pv]
        return self._bas[(i, j)]
    def nf(self, i, j, comb):
        if not comb: return {}
        if i == j: return {k: v for k, v in comb.items() if v != 0}
        return ap.reduceAgainstPivots(comb, self.piv(i, j))

def rank(rows):
    piv = {}; rk = 0
    for r in rows:
        r = {k: (v.numerator * pow(v.denominator, -1, P)) % P for k, v in r.items()}
        r = {k: v for k, v in r.items() if v}
        while r:
            h = min(r)
            if h in piv:
                f = r[h]; pr = piv[h]
                for k, v in pr.items(): r[k] = (r.get(k, 0) - f * v) % P
                r = {k: v for k, v in r.items() if v}
            else:
                inv = pow(r[h], -1, P); piv[h] = {k: v * inv % P for k, v in r.items()}; rk += 1; break
    return rk

def complexT(A, k, V, S=0):
    """list over summands i in V of (degs, d): degs {deg: [vertices]}, d {deg: {(a,b): comb}} (a index in deg, b in deg+1)."""
    out = {}
    for i in V:
        if i != k: out[i] = ({0: [i]}, {})
        else:
            outs = ap.arrowsOutOf(A.q, k)
            out[i] = ({S: [k], S + 1: [b[1] for b in outs]}, {S: {(0, a): {(outs[a],): Fraction(1)} for a in range(len(outs))}})
    return out

def compose_mat(A, g, f, Xv, Yv, Zv):
    """g o f; f: X->Y entries f[(a,b)] comb of paths Yv[b]~>Xv[a]; g[(b,c)] paths Zv[c]~>Yv[b]; result paths Zv[c]~>Xv[a]."""
    res = {}
    for (a, b), fc in f.items():
        for (b2, c), gc in g.items():
            if b2 != b: continue
            acc = res.setdefault((a, c), {})
            for pg, cg in gc.items():
                for pf, cf in fc.items():
                    p = pg + pf; acc[p] = acc.get(p, 0) + cg * cf
    return res

def homdim(A, X, Y, m):
    """dim Hom_K(X, Y[m]) for complexes X=(degs,d), Y=(degs,d)."""
    Xd, Xdi = X; Yd, Ydi = Y
    def hom_basis(dx, dy):   # maps X^dx -> Y^dy: list of (a,b,path)
        out = []
        for a, va in enumerate(Xd.get(dx, [])):
            for b, vb in enumerate(Yd.get(dy, [])):
                for p in A.basis(va, vb): out.append((a, b, p))
        return out
    degs = sorted(set(Xd) | {d - m for d in Yd} | {d for d in Xd})
    C = {d: hom_basis(d, d + m) for d in degs}
    dimC = sum(len(v) for v in C.values())
    def apply_chain(d, a, b, p):   # basis elt f in Hom(X^d,Y^{d+m}); returns coords of (d_Y f - f d_X) in Hom(X^d, Y^{d+m+1}) and Hom(X^{d-1}, Y^{d+m})
        out = {}
        # d_Y f : X^d -> Y^{d+m+1}
        vb = Yd[d + m][b]
        for (b1, c), cm in Ydi.get(d + m, {}).items():
            if b1 != b: continue
            vc = Yd[d + m + 1][c]
            comb = {}
            for pc, cc in cm.items(): comb[p + pc] = comb.get(p + pc, 0) + cc
            va = Xd[d][a]
            for pp, cf in A.nf(va, vc, comb).items(): out[('A', d, a, c, pp)] = cf
        # f d_X : X^{d-1} -> Y^{d+m}  (f in degree d precomposed with d_X^{d-1}: X^{d-1}->X^d)
        for (a0, a1), cm in Xdi.get(d - 1, {}).items():
            if a1 != a: continue
            v0 = Xd[d - 1][a0]; vb_ = Yd[d + m][b]
            comb = {}
            for pc, cc in cm.items(): comb[pc + p] = comb.get(pc + p, 0) + cc
            for pp, cf in A.nf(v0, vb_, comb).items(): out[('B', d - 1, a0, b, pp)] = out.get(('B', d - 1, a0, b, pp), 0) - cf
        return out
    # Phi: chain map defect; combine keys
    rowsPhi = []
    for d in degs:
        for (a, b, p) in C[d]:
            rowsPhi.append({(k): Fraction(v) for k, v in apply_chain(d, a, b, p).items()})
    # rows are images of basis vectors; rank of the image = rank of these rows. The 'A'-part of degree d and 'B'-part of degree d+1 land in the same space: rename
    def norm(r):
        o = {}
        for k, v in r.items():
            if k[0] == 'A': kk = (k[1], k[2], k[3], k[4])     # map X^d -> Y^{d+m+1}, component (a,c,path)
            else: kk = (k[1] + 1, k[2], k[3], k[4])           # f d_X : X^{d-1} -> Y^{d+m}  is a map X^{d-1}->Y^{(d-1)+m+1}; deg index d-1
            o[kk] = o.get(kk, 0) + v
        return {k: v for k, v in o.items() if v}
    # careful: B-part is a map from X^{d-1}, so its source degree index is d-1 (k[1] already d-1). A-part source degree is d.
    def norm(r):
        o = {}
        for k, v in r.items():
            kk = (k[1], k[2], k[3], k[4]); o[kk] = o.get(kk, 0) + v
        return {k: v for k, v in o.items() if v}
    rkPhi = rank([norm(r) for r in rowsPhi])
    # Psi: homotopies s^d : X^d -> Y^{d+m-1}; image  d_Y s + s d_X  in Hom^m
    rowsPsi = []
    for d in degs:
        for (a, b, p) in hom_basis(d, d + m - 1):
            out = {}
            vb = Yd[d + m - 1][b]; va = Xd[d][a]
            # d_Y s : X^d -> Y^{d+m}
            for (b1, c), cm in Ydi.get(d + m - 1, {}).items():
                if b1 != b: continue
                vc = Yd[d + m][c]; comb = {}
                for pc, cc in cm.items(): comb[p + pc] = comb.get(p + pc, 0) + cc
                for pp, cf in A.nf(va, vc, comb).items(): out[(d, a, c, pp)] = out.get((d, a, c, pp), 0) + cf
            # s d_X : X^{d-1} -> Y^{d+m-1}
            for (a0, a1), cm in Xdi.get(d - 1, {}).items():
                if a1 != a: continue
                v0 = Xd[d - 1][a0]; comb = {}
                for pc, cc in cm.items(): comb[pc + p] = comb.get(pc + p, 0) + cc
                for pp, cf in A.nf(v0, vb, comb).items(): out[(d - 1, a0, b, pp)] = out.get((d - 1, a0, b, pp), 0) + cf
            rowsPsi.append({k: Fraction(v) for k, v in out.items() if v})
    rkPsi = rank(rowsPsi)
    return dimC - rkPhi - rkPsi

def stepTest(alg, k, S=0):
    """alg: PathAlgebra (acyclic).  Returns dict(H=matrix dim Hom(T_i,T_j), Hm1, Hp1 = max dims at m=-1,+1, A=Hom(A,A) check)."""
    A = Alg(alg); V = sorted(alg.quiver.nodes); T = complexT(A, k, V, S)
    res = {}
    for m in (-1, 0, 1):
        res[m] = [[homdim(A, T[i], T[j], m) for j in V] for i in V]
    return V, res
