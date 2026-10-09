"""Minimal non-monomial test of CHZ Cor 3.6 (path-wise) against Prop 3.5 (socle form).

Q: a:1->2, b:2->4, c:1->3, d:3->4, e:4->5;  I = <(ab+cd)e>  (admissible, non-monomial).
Paths compose left to right (p q = p then q), right modules P_i = e_i A.
For S = {v}: Prop 3.5(2) says (filt S_v, S_v^perp) is a derived equivalence iff
v not in supp soc P_i for all i != v.  Cor 3.6(2) path-wise uses tail-maximal paths.
Also J_i = ker(x -> (x b)_b on e_i A e_v) of Aihara-Iyama 2.32(b) (= Hom(S_v, e_iA)).
Run: .venv/bin/python workshop/rounds/055/scholar_chz_nonmonomial.py
"""
import itertools
import numpy as np

arrows = {'a': (1, 2), 'b': (2, 4), 'c': (1, 3), 'd': (3, 4), 'e': (4, 5)}
V = [1, 2, 3, 4, 5]
# all paths as tuples of arrow names; trivial path = ('e',i)
paths = [('#', i) for i in V]  # trivial
def ext(p):
    out = []
    end = p[-1][1] if p[0] != '#' else None
    return out
allp = []
def endv(p):
    return p[1] if p[0] == '#' else arrows[p[-1]][1]
def startv(p):
    return p[1] if p[0] == '#' else arrows[p[0]][0]
frontier = [('#', i) for i in V]
allp += frontier
while frontier:
    new = []
    for p in frontier:
        for n, (s, t) in arrows.items():
            if s == endv(p):
                q = (n,) if p[0] == '#' else p + (n,)
                new.append(q)
    allp += new
    frontier = new
idx = {p: k for k, p in enumerate(allp)}
def mul(p, q):
    if endv(p) != startv(q):
        return None
    if p[0] == '#':
        return q
    if q[0] == '#':
        return p
    return p + q
# generator z = (ab+cd)e as dict path->coef
z = {('a', 'b', 'e'): 1, ('c', 'd', 'e'): 1}
def vec(d):
    v = np.zeros(len(allp))
    for p, c in d.items():
        v[idx[p]] += c
    return v
gens = []
for u in allp:
    for w in allp:
        d = {}
        for p, c in z.items():
            up = mul(u, p)
            if up is None: continue
            upw = mul(up, w)
            if upw is None: continue
            d[upw] = d.get(upw, 0) + c
        if d: gens.append(vec(d))
I = np.array(gens)
def rank(M):
    return 0 if len(M) == 0 else np.linalg.matrix_rank(np.array(M))
rI = rank(I)
def in_I_rank(vs):  # rank of I + span(vs) - rank I  (dimension of vs mod I)
    return rank(np.vstack([I] + list(vs))) - rI
def e_i_A_e_j(i, j):
    return [p for p in allp if startv(p) == i and endv(p) == j]
def soc_dim(i, j):
    """dim of {x in e_iAe_j : x*a in I for all arrows a out of j}, modulo I."""
    B = e_i_A_e_j(i, j)
    outs = [n for n, (s, t) in arrows.items() if s == j]
    # kernel of x -> (x a mod I)_a on span(B) / (I cap span B)
    if not B: return 0
    # work in quotient: dim(span B + I)/I minus dim of image
    dimB = in_I_rank([vec({p: 1}) for p in B])
    if not outs: return dimB
    # build map from span(B) to (A/I)^outs: compute dim of kernel via ranks
    # kernel of composite span(B) -> prod (A/I): K = {x: x a in I all a}
    # dim(K + I-part) computed by brute: nullspace of [x a rows | I-basis columns]
    # Use: dim image = rank of map into quotient; dim(B mod I) - dim image = dim soc part
    # image dimension = rank of stacked [ (x a)_a ; I⊗outs ] - rank(I⊗outs) over concatenated spaces
    n = len(allp)
    Ibig = np.zeros((len(I) * len(outs), n * len(outs)))
    for k, a in enumerate(outs):
        Ibig[k * len(I):(k + 1) * len(I), k * n:(k + 1) * n] = I
    rows = []
    for p in B:
        r = np.zeros(n * len(outs))
        for k, a in enumerate(outs):
            r[k * n + idx[mul(p, (a,))]] = 1
        rows.append(r)
    # contributions of I∩span(B) map to I-part automatically (ideal), so image rank:
    img = rank(np.vstack([Ibig] + rows)) - rank(Ibig)
    return dimB - img
def tailmax(i, j):
    out = []
    for p in e_i_A_e_j(i, j):
        if in_I_rank([vec({p: 1})]) == 0: continue
        outs = [n for n, (s, t) in arrows.items() if s == j]
        if all(in_I_rank([vec({mul(p, (a,)): 1})]) == 0 for a in outs):
            out.append(p)
    return out
print("dim I =", rI)
print("i | supp soc P_i (dims)      | ends of tail-maximal paths")
for i in V:
    ss = {j: soc_dim(i, j) for j in V if soc_dim(i, j) > 0}
    tm = sorted({endv(p) for j in V for p in tailmax(i, j)})
    print(i, "|", ss, "|", tm)
v = 4
print("J_i = Hom(S_4, e_iA) dims (i != 4):", {i: soc_dim(i, v) for i in V if i != v})
Sc = [i for i in V if i != v]
print("Prop 3.5(2) closure Phi+(S^c) subset S^c, S={4}: ",
      all(soc_dim(i, v) == 0 for i in Sc), "(socle form)")
print("Cor 3.6(2) path-wise closure, S={4}:", all(v not in {endv(p) for j in V for p in tailmax(i, j)} for i in Sc))
