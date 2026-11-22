"""Round 135 maverick: object-by-object comparison of the two nu implementations.
nu_A = rounds/089 maverick_bimod.nu_bimod (bimodule resolution of D(A), used by E-346);
nu_B = rounds/135 maverick_nu2.nu (termwise nu + cone replacement, independent code).
For each object X (S_l, P_l, or uniserials) compares dim H^t(nu^k X)_v for every k <= m, v, t.
Usage (repo root): timeout 10m .venv/bin/python workshop/rounds/135/maverick_compare.py REL m [objs] [--plant=scale|kup|trunc|long|id]
"""
import sys, time
R = "/home/user/quiverMutation/workshop/rounds"
sys.path.insert(0, R + "/135"); sys.path.insert(0, R + "/134"); sys.path.insert(0, R + "/089"); sys.path.insert(0, R + "/087")
import maverick_nu2 as B
from maverick_bimod import resolve, nu_bimod
from scholar_rankgrowth import Alg as AlgA, Cx as CxA
from scholar_simple_nu import cohomology

def resA(A, s, e):          # minimal projective resolution of [s, e] (as in rounds/134 scholar_uniserial.res)
    ivs = [(s, e)]
    while True:
        a, b = ivs[-1]; nb = a + A.p[a] - 1
        if b + 1 > nb: break
        ivs.append((b + 1, nb))
    X = CxA()
    for k, (a, _) in enumerate(ivs): X.add(k, -k, a)
    for k in range(1, len(ivs)): X.setd(k, k - 1, 1)
    return X

if __name__ == "__main__":
    a = sys.argv[1:]
    for md in ("scale", "kup", "id", "trunc", "long"):
        if "--plant=" + md in a: B.PLANT = md; a.remove("--plant=" + md)
    rel, m = a[0], int(a[1]); spec = a[2] if len(a) > 2 else "S,P"
    A = AlgA([int(c) for c in rel]); L = resolve(A)
    B.ALG = B.Alg(rel); B.NV = B.ALG.n
    assert [A.p[i] for i in range(1, A.n + 1)] == [B.ALG.U[i] - i + 1 for i in range(1, A.n + 1)], "Kupisch series differ"
    t0 = time.time(); agree = dis = 0; first = []
    for name, iv in B.objects(B.ALG, spec):
        XA = resA(A, *iv); XB = B.minimise(B.replace(B.stalk(iv), B.ALG))
        for k in range(1, m + 1):
            XA = nu_bimod(A, L, XA)[0]; XB = B.nu(XB)
            HA = cohomology(A, XA); HB = B.cohom(XB)
            if HA == HB: agree += 1
            else:
                dis += 1
                if len(first) < 5: first.append((name, k, sum(HA.values()), sum(HB.values())))
    print(rel, "m", m, "objs", spec, "plant", B.PLANT, "(object, k) pairs agreeing", agree, "disagreeing", dis,
          "first disagreements (obj, k, totalH_A, totalH_B)", first, f"{time.time()-t0:.0f}s", flush=True)
