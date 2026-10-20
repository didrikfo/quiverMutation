"""Round 077 (scholar), T5: J_j = Hom(S_v, e_j A) in the repo's convention (workshop/rounds/033/scholar_socle.py:
J_j = {x in e_j A e_v : x a = 0 in A for every arrow a out of v}), for MONOMIAL parents with true paths.
By AI Thm 2.32(b) / 2.50(c) (arXiv:1009.3370v3), T (right-mutation complex at loopless v) is tilting iff J_j = 0 for all j != v.
Cases: E-211 parent at v=1 (cyclic D(1)); A_3/rad^2 at v=1 (Oppermann Ex., tilting) and v=2 (silting not tilting, E-217 side note)
as controls that the count can say both yes and no.  Run: .venv/bin/python workshop/rounds/077/scholar_jcheck.py   (<1 s)"""
def zero(p, rels): return any(any(p[k:k+len(r)] == r for k in range(len(p)-len(r)+1)) for r in rels)
def J(arrows, rels, v, maxlen=12):
    out = {}
    verts = sorted({x for t in arrows.values() for x in t})
    for j in verts:
        if j == v: continue
        # nonzero paths j ~> v (length >= 1), then test x*a = 0 for every arrow a with tail v
        paths, fr = [], [((a,), arrows[a][1]) for a in arrows if arrows[a][0] == j]
        for _ in range(maxlen):
            fr = [(p, h) for p, h in fr if not zero(p, rels)]
            paths += [p for p, h in fr if h == v]
            fr = [(p+(a,), arrows[a][1]) for p, h in fr for a in arrows if arrows[a][0] == h]
        assert not fr or all(zero(p, rels) for p, h in fr), "not finite-dimensional?"
        outs = [a for a in arrows if arrows[a][0] == v]
        out[j] = sum(1 for p in paths if all(zero(p+(a,), rels) for a in outs))  # monomial: x a in I iff the path is zero
    return out
cases = {
 "E-211 parent, v=1": (dict(a=(1,2), s=(2,3), b=(3,4), c=(4,3)), [("a","s"),("s","b","c"),("c","b")], 1),
 "A3/rad^2, v=1":     (dict(f=(1,2), g=(2,3)), [("f","g")], 1),
 "A3/rad^2, v=2":     (dict(f=(1,2), g=(2,3)), [("f","g")], 2),
}
for k,(A,R,v) in cases.items(): print(k, "J_j =", J(A,R,v))

# R077 response to referee (item 3): DUAL count J^dual_j = #{nonzero paths v ~> j (length >= 1) with a x = 0 for every arrow a INTO v}
# (= Hom(S_v, e_v A-side of the opposite algebra); the tilting criterion for the other-side mutation).  Same cases.
def Jdual(arrows, rels, v, maxlen=12):
    out, ins = {}, [a for a in arrows if arrows[a][1] == v]
    verts = sorted({x for t in arrows.values() for x in t})
    paths, fr = [], [((a,), arrows[a][1]) for a in arrows if arrows[a][0] == v]
    for _ in range(maxlen):
        fr = [(p, h) for p, h in fr if not zero(p, rels)]
        paths += [(p, h) for p, h in fr]
        fr = [(p+(a,), arrows[a][1]) for p, h in fr for a in arrows if arrows[a][0] == h]
    for j in verts:
        if j == v: continue
        out[j] = sum(1 for p, h in paths if h == j and all(zero((a,)+p, rels) for a in ins))
    return out
print("--- dual count")
for k,(A,R,v) in cases.items(): print(k, "J^dual_j =", Jdual(A,R,v))
