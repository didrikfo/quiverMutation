"""Round 053 (toolsmith), T10 agenda item 1: End(T) of the tilting complex as a QUIVER WITH RELATIONS, compared with the mutated algebra.

Library part (this file, no search): given a PathAlgebra `a` (acyclic) and a vertex k, build T = (+_{i != k} P_i) + T_k,
T_k = (P_k -> +_{k->h} P_h) in degrees (S, S+1), S = -1 (the convention of skeptic_replay13 / E-159), P_x with Hom(P_x,P_y) = paths x ~> y mod I,
composition `f then g` = path concatenation f + g (rounds/050/skeptic_tilt.py).  Then, over Q (exact Fractions):
  * Hom_K(T_i, T_j) = chain maps / homotopies, with an explicit basis and coordinates;  composition of basis elements;
  * radical filtration rad^m(i,j) (i != j; End(T_i) = K is checked), arrows of End(T) = rad/rad^2, lifts chosen by a fixed rule;
  * the algebra map  K[Q_End] -> End(T)  on all paths, hence the image dimensions (= dim Hom) and the relation ideal of End(T).
Comparison with the claimed child c (PathAlgebra, c = reducePathAlgebra(quiverMutationAtVertex(a, v))):
  (1) arrows: number of arrows i -> j of End(T) equals number of arrows of c between the same labels, in one of the two directions;
  (2) dims: dim Hom(T_i,T_j) = dim e_.. c e_.. (paths mod I);
  (3) relations: every relation of c, with its arrows sent to lifted morphisms (times scalars when needed), is zero in End(T).
      (1)+(2)+(3) => K Q_c / I_c -> End(T) is a surjective algebra map between equal finite dimensions, hence an isomorphism (decisive);
      if (3) cannot be settled the verdict is only 'invariants agree'.
Usage: imported by toolsmith_endt_run.py."""
import sys
sys.path.insert(0, '.')
from fractions import Fraction
from itertools import product
from quivermutation import arrowPaths as ap, procedure
_src = open('workshop/rounds/050/skeptic_tilt.py').read()
exec(_src[_src.index('class Alg'):_src.index('def rank')])
exec(_src[_src.index('def complexT'):_src.index('def compose_mat')])


class Echelon:
    """Exact echelon form over Q of sparse vectors (dicts), with optional tags (coefficient dicts) carried along."""
    def __init__(self):
        self.rows = {}           # pivot key -> (vec, tag)
    def reduce(self, v, tag=None):
        v = {k: Fraction(x) for k, x in v.items() if x != 0}; tag = dict(tag or {})
        while v:
            h = min(v, key=repr)
            if h not in self.rows: return v, tag, h
            pv, pt = self.rows[h]; f = v[h]
            for k, x in pv.items():
                nv = v.get(k, 0) - f * x
                if nv == 0: v.pop(k, None)
                else: v[k] = nv
            for k, x in pt.items():
                nt = tag.get(k, 0) - f * x
                if nt == 0: tag.pop(k, None)
                else: tag[k] = nt
        return v, tag, None
    def add(self, v, tag=None):
        v, tag, h = self.reduce(v, tag)
        if h is None: return False
        inv = 1 / v[h]
        self.rows[h] = ({k: x * inv for k, x in v.items()}, {k: x * inv for k, x in tag.items()}); return True


def nullspace(cols):
    """cols: list of sparse image vectors of the basis e_0..e_{m-1}; returns basis of {c : sum c_i cols_i = 0} as dicts index->coef."""
    E = Echelon(); out = []
    for i, c in enumerate(cols):
        v, tag, h = E.reduce(c, {i: Fraction(1)})
        if h is None: out.append(tag)
        else:
            inv = 1 / v[h]; E.rows[h] = ({k: x * inv for k, x in v.items()}, {k: x * inv for k, x in tag.items()})
    return out


class TiltEnd:
    def __init__(self, alg, k, S=-1):
        self.A = Alg(alg); self.k = k; self.V = sorted(alg.quiver.nodes); self.T = complexT(self.A, k, self.V, S)
        self._hom = {}; self._comp = {}
        for i in self.V:
            for j in self.V: self.hom(i, j)
    # --- chain maps from X to Y (degree 0), keys (d, a, b, path)
    def _chain_space(self, X, Y):
        Xd, Xdi = X; Yd, Ydi = Y; A = self.A
        unk = []
        for d in sorted(set(Xd) & set(Yd)):
            for a, va in enumerate(Xd[d]):
                for b, vb in enumerate(Yd[d]):
                    for p in A.basis(va, vb): unk.append((d, a, b, p))
        # condition: d_Y f^d = f^{d+1} d_X, components X^d[a] -> Y^{d+1}[c]
        def image(u):
            d, a, b, p = u; out = {}
            for (b1, c), cm in Ydi.get(d, {}).items():           # d_Y after f
                if b1 != b: continue
                comb = {p + pc: cc for pc, cc in cm.items()}
                for pp, cf in A.nf(Xd[d][a], Yd[d + 1][c], comb).items(): out[('L', d, a, c, pp)] = out.get(('L', d, a, c, pp), 0) + cf
            for (a0, a1), cm in Xdi.get(d - 1, {}).items():      # f^d after d_X^{d-1}: X^{d-1}[a0] -> X^d[a1=a] -> Y^d[b]
                if a1 != a: continue
                comb = {pc + p: cc for pc, cc in cm.items()}
                for pp, cf in A.nf(Xd[d - 1][a0], Yd[d][b], comb).items(): out[('R', d - 1, a0, b, pp)] = out.get(('R', d - 1, a0, b, pp), 0) - cf
            # rename so that L and R parts of the same target component coincide: target = (source degree, source idx, target idx, path)
            o = {}
            for key, cf in out.items(): o[key[1:]] = o.get(key[1:], 0) + cf
            return {kk: v for kk, v in o.items() if v != 0}
        ker = nullspace([image(u) for u in unk])
        # homotopies s^d : X^d -> Y^{d-1};  h = d_Y s + s d_X  : X^d -> Y^d
        hom = []
        for d in sorted(set(Xd)):
            if d - 1 not in Yd: continue
            for a, va in enumerate(Xd[d]):
                for b, vb in enumerate(Yd[d - 1]):
                    for p in A.basis(va, vb):
                        out = {}
                        for (b1, c), cm in Ydi.get(d - 1, {}).items():
                            if b1 != b: continue
                            comb = {p + pc: cc for pc, cc in cm.items()}
                            for pp, cf in A.nf(va, Yd[d][c], comb).items(): out[(d, a, c, pp)] = out.get((d, a, c, pp), 0) + cf
                        for (a0, a1), cm in Xdi.get(d - 1, {}).items():
                            if a1 != a: continue
                            comb = {pc + p: cc for pc, cc in cm.items()}
                            for pp, cf in A.nf(Xd[d - 1][a0], vb, comb).items(): out[(d - 1, a0, b, pp)] = out.get((d - 1, a0, b, pp), 0) + cf
                        out = {kk: v for kk, v in out.items() if v != 0}
                        if out: hom.append(out)
        # chain-map vectors as dicts over unknown keys
        kervecs = [{unk[i]: c for i, c in kv.items()} for kv in ker]
        E = Echelon()
        for h in hom: E.add(h)
        picks = []; Q = Echelon()
        for h in hom: Q.add(h, {})
        for kv in kervecs:
            v, tag, h = Q.reduce(kv, {len(picks): Fraction(1)})
            if h is None: continue
            # kv independent of im + earlier picks: it becomes pick number len(picks)
            inv = 1 / v[h]; Q.rows[h] = ({kk: x * inv for kk, x in v.items()}, {kk: x * inv for kk, x in tag.items()}); picks.append(kv)
        return picks, Q
    def hom(self, i, j):
        if (i, j) not in self._hom: self._hom[(i, j)] = self._chain_space(self.T[i], self.T[j])
        return self._hom[(i, j)]
    def dim(self, i, j): return len(self.hom(i, j)[0])
    def coords(self, i, j, f):
        """coordinates of chain map f (dict key->coef, keys already normal-form paths) in the basis picks of Hom(i,j)."""
        picks, Q = self.hom(i, j); v, tag, h = Q.reduce(f, {})
        # Q.reduce tracks tags only as accumulated 'tag' of rows; rows of the im-part carry empty tags, picks carry unit tags
        assert h is None, 'f is not a combination of homotopies and picks: not a chain map'
        # tag accumulates  -sum(f_h * tag_h);  coefficient of the pick = minus that
        return {kk: -c for kk, c in tag.items() if c != 0}
    def compose(self, i, j, l, f, g):
        """f in Hom(T_i, T_j) (chain map dict), g in Hom(T_j, T_l): g after f, chain map i -> l."""
        A = self.A; Xd = self.T[i][0]; Yd = self.T[j][0]; Zd = self.T[l][0]; out = {}
        gi = {}
        for (d, b, c, q), cg in g.items(): gi.setdefault((d, b), []).append((c, q, cg))
        for (d, a, b, p), cf in f.items():
            for c, q, cg in gi.get((d, b), []):
                comb = {p + q: cf * cg}
                for pp, cc in A.nf(Xd[d][a], Zd[d][c], comb).items(): out[(d, a, c, pp)] = out.get((d, a, c, pp), 0) + cc
        return {kk: v for kk, v in out.items() if v != 0}
    def basis_maps(self, i, j): return self.hom(i, j)[0]
    def ident(self, i):
        Xd = self.T[i][0]; return {(d, a, a, ()): Fraction(1) for d in Xd for a in range(len(Xd[d]))}
