"""Round 119 maverick: E-269's Z-congruence search extended to SINGULAR S (n = 12 low-lam groups 0 and 3 of rounds/107/maverick_n12low.py: det S = 0, ker S rank 1,
Phi has a Jordan block at eigenvalue 1 so the isometry group is infinite beyond Phi^Z and E-276's closed-form bound does not apply).
Method: R = P^T ranges over the commutant lattice {R Pb = Pa R} (pari matkerint, exact integers); enumerate ALL its vectors with Frobenius norm tr(R^T R) <= B (qfminim, B escalated
over BS), test R^T Ca R = Cb and |det R| = 1 exactly.  YES is certified by the explicit P; NO means only "no P with ||P||_F^2 <= B" -- NOT a proof (no bound argument).
Invariant for NO controls: q of the nondegenerate quotient lattice Z^n / ker S (unimodular U with last column = saturated kernel; S' = U^T S U; disc form of its 11x11 block).
usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/119/maverick_sing.py group(0|3) [B list comma] [maxpairs]"""
import sys, collections, time, itertools; sys.path.insert(0, "workshop/rounds/097")
import numpy as np, sympy as sp
from fractions import Fraction
from maverick_cartanZ import cartan, disc_form, phi, pari
from quivermutation import coxeterTables as ct
from quivermutation.quipuForms import quipuForAlmostSeparateLNA
pari.allocatemem(2 * 10 ** 9)
N = 12
def sing_groups():
    keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
    for rl, kk in keys.items(): groups[kk].append(tuple(rl))
    out = []
    for m in groups.values():
        C = np.array(ct.lnaCartanMatrix(N, list(m[0]))).astype(float)
        if round(np.linalg.det(C)) != 1 or abs(np.linalg.det(C + C.T)) > 1e-6: continue
        ev = np.linalg.eigvalsh(C + C.T)
        if (ev < -1e-9).sum() != 1: continue
        w = np.linalg.eigvals(-np.linalg.inv(C.T) @ C); big = [z.real for z in w if abs(z.imag) < 1e-9 and z.real > 1 + 1e-3]
        lam = big[0] if big else 1.0
        if lam < 1.5: out.append((lam, m))
    return sorted(out, key=lambda t: t[0])
def qbar_fp(Cm):
    """fingerprint of the discriminant form of S on Z^n/ker S."""
    S = Cm + Cm.T; ker = S.nullspace()[0]; ker = ker * sp.ilcm(*[x.q for x in ker]); g = sp.igcd(*[int(x) for x in ker]); ker = (ker / g)
    # unimodular U with last column ker: complete via pari matsnf-free approach: Hermite on [ker] -> use pari mathnf/ matkerint trick
    k = [int(x) for x in ker]; n = N
    # find unimodular completion: extended gcd row operations
    U = sp.eye(n); v = sp.Matrix(k)
    # reduce v to e_1 * 1 by unimodular row ops E (E v = e_n) then U = E^-1
    E = sp.eye(n); w = list(k)
    for i in range(n - 1):
        a, b = w[i], w[i + 1]
        if a == 0: continue
        if b == 0: Op = sp.eye(n); Op[i, i], Op[i, i + 1], Op[i + 1, i], Op[i + 1, i + 1] = 0, 1, 1, 0
        else:
            g_ = int(sp.gcd(a, b)); s, t, _ = sp.gcdex(a, b); s = int(s); t = int(t)
            Op = sp.eye(n); Op[i, i], Op[i, i + 1], Op[i + 1, i], Op[i + 1, i + 1] = -b // g_, a // g_, s, t
        E = Op * E; w = list(Op * sp.Matrix(w)); 
    assert abs(int(w[-1])) == 1 and all(x == 0 for x in w[:-1]), w
    Ui = E.inv(); Sp = Ui.T * S * Ui
    assert all(Sp[i, n - 1] == 0 for i in range(n)), "kernel completion"
    Sb = Sp[:n - 1, :n - 1]; d, El, G, bil, q, add = disc_form(Sb)
    Sf = [[Fraction(int(Sb[i, j])) for j in range(n - 1)] for i in range(n - 1)]
    return abs(d), tuple(sorted(collections.Counter(sum(x[i] * Sf[i][j] * x[j] for i in range(n - 1) for j in range(n - 1)) % 2 for x in El).items()))
def search(Ca, Cb, Bs_, cap=400000):
    n = N; Pa, Pb = phi(Ca), phi(Cb)
    A = sp.kronecker_product(sp.eye(n), Pb.T) - sp.kronecker_product(Pa, sp.eye(n))
    K = pari.matkerint(pari.matrix(n * n, n * n, [int(v) for v in A])); k = int(pari.matsize(K)[1])
    Bm = [np.array([int(K[i, j]) for i in range(n * n)], dtype=object).reshape(n, n) for j in range(k)]
    G = pari.matrix(k, k, [int(np.trace(Bm[i].T.dot(Bm[j]))) for i in range(k) for j in range(k)])
    U = pari.qflllgram(G); Gr = U.mattranspose() * G * U
    Ui = np.array([[int(U[i, j]) for j in range(k)] for i in range(k)], dtype=object)
    Cai = np.array(Ca.tolist(), dtype=object); Cbi = np.array(Cb.tolist(), dtype=object); res_all = []
    for B in Bs_:
        res = pari.qfminim(Gr, B, cap, 0); nv = int(pari.matsize(res[2])[1]); over = 2 * nv != int(res[0])
        sols = 0; first = None
        if not over:
            for c in range(nv):
                x = Ui.dot(np.array([int(res[2][i, c]) for i in range(k)], dtype=object))
                R = sum((int(x[j]) * Bm[j] for j in range(k)), np.zeros((n, n), dtype=object)); P = R.T
                if (P.dot(Cai).dot(P.T) == Cbi).all() and abs(int(pari.matdet(pari.matrix(n, n, [int(v) for v in P.flatten()])))) == 1:
                    sols += 1; first = first or P.tolist()
        res_all.append(dict(B=B, rank=k, nvec=nv, overflow=over, nsol=sols, P=first))
        if sols or over: break
    return res_all
if __name__ == "__main__":
    gi = int(sys.argv[1]); Bs_ = [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "40,100,250").split(",")]; maxp = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    keys = ct.lnaKeys(N); tgt = (2, 4, 5, 0, 0, 0, 4, 0, 0, 0) if gi == 0 else (2, 2, 4, 0, 0, 3, 3, 0, 3, 0); kk = keys[tgt]; mem = [tuple(r) for r, k in keys.items() if k == kk]; lam = 0.0; print("group", gi, "lam %.4f" % lam, "size", len(mem), flush=True)
    fp = collections.defaultdict(list); q = collections.defaultdict(list)
    for rl in mem:
        fp[qbar_fp(cartan(N, list(rl)))].append(rl)
        try: r = quipuForAlmostSeparateLNA(N, list(rl))
        except Exception: r = None
        if r is not None: q[r].append(rl)
    print("q-bar fingerprints:", {str(k): len(v) for k, v in fp.items()}, "| quipus:", {str(k): len(v) for k, v in q.items()}, flush=True)
    pairs = [("self", mem[0], mem[0])]
    for v in q.values():
        if len(v) >= 2: pairs.append(("within-quipu", v[0], v[1])); pairs.append(("within-quipu", v[0], v[-1]))
    if len(q) >= 2:
        vs = list(q.values()); pairs.append(("cross-quipu", vs[0][0], vs[1][0]))
    keysf = list(fp)
    for kf in keysf:
        if len(fp[kf]) >= 2: pairs.append(("same-qbar", fp[kf][0], fp[kf][1])); pairs.append(("same-qbar", fp[kf][0], fp[kf][-1]))
    if len(keysf) >= 2: pairs.append(("qbar-DIFFERENT", fp[keysf[0]][0], fp[keysf[1]][0]))
    for tag, a, b in pairs[:maxp]:
        t = time.time(); r = search(cartan(N, list(a)), cartan(N, list(b)), Bs_)
        print(tag, a, b, [{k: v for k, v in x.items() if k != "P"} for x in r], "%.0fs" % (time.time() - t), flush=True)
