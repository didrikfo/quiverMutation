"""Round 138 maverick: tau-orbit / AR-triangle test of E-346's two nu-orbit families (uses rounds/135 maverick_nu2).

For a linear Nakayama algebra A (REL digits as in the record) and every uniserial X = [a, e] (top a, socle e):
  1. iterate nu up to m (maverick_nu2) and record every k <= m with nu^k X ~ X[s] (stalk return), so the first
     return (p, q) and the return at m;
  2. AR test.  tau_D = nu[-1] on D^b(A).  If M = [a, e] is non-projective, pd M <= 1 and Hom(M, A) = 0, then
     tau_D M = tau_A M = [a+1, e+1] (a module), and the AR sequence 0 -> [a+1,e+1] -> [a,e+1] + [a+1,e] -> M -> 0
     is the AR triangle of D^b.  F = nu^m[-s] is an autoequivalence; if F M ~ M then F tau M ~ tau M and F maps
     the AR triangle to itself, so F permutes the summands of the middle term E.  Prediction: when M returns at m
     with shift s, every summand of E returns with the same s (or two summands are swapped by F).  Conversely a
     returning summand forces nothing; we report both directions.
  Control (NO where NO is true): the same test applied to the *naive* mod-A AR sequences of all non-projective M,
  ignoring the condition (there tau_D M != tau_A M, the sequence is not an AR triangle and nothing is predicted).
  Check of the condition itself: for every M satisfying it, H^*(nu M) must be tau_A M in degree -1.
Usage (repo root): timeout 10m .venv/bin/python workshop/rounds/138/maverick_tau.py REL m [--prime p]
"""
import sys, time

src = open("workshop/rounds/135/maverick_nu2.py").read().split("\nif __name__ == '__main__':")[0]
NU = {}
exec(compile(src, "maverick_nu2", "exec"), NU)


def main():
    a = sys.argv[1:]
    if "--prime" in a:
        i = a.index("--prime"); NU["P"] = int(a[i + 1]); del a[i:i + 2]
    rel, m = a[0], int(a[1])
    alg = NU["Alg"](rel); NU["ALG"] = alg; NU["NV"] = alg.n
    n, U = alg.n, alg.U
    unis = [(x, e) for x in range(1, n + 1) for e in range(x, U[x] + 1)]
    print("REL", rel, "n", n, "P", NU["P"], "m", m, "U", U[1:], flush=True)
    t0 = time.time()
    ret = {}      # iv -> list of (k, s)
    taucheck = {}
    links = {}    # iv -> [(k, Y, s)] with H(nu^k iv) = H(Y) in degree -s, Y another uniserial (same nu-orbit up to shift,
                  # as far as cohomology tells: a complex with H a uniserial stalk IS that stalk)
    for iv in unis:
        X = NU["minimise"](NU["replace"](NU["stalk"](iv), alg)); NU["check"](X)
        rs = []
        for k in range(1, m + 1):
            X = NU["nu"](X)
            H = NU["cohom"](X)
            if k == 1: taucheck[iv] = H
            ts = {t for (_, t) in H}
            if len(ts) == 1 and H == {(v, list(ts)[0]): 1 for v in range(iv[0], iv[1] + 1)}:
                rs.append((k, -list(ts)[0]))
            if len(ts) == 1 and all(c == 1 for c in H.values()):
                vs = sorted(v for (v, _) in H)
                if vs == list(range(vs[0], vs[-1] + 1)) and (vs[0], vs[-1]) != iv and (vs[0], vs[-1]) in set(unis):
                    links.setdefault(iv, []).append((k, (vs[0], vs[-1]), -list(ts)[0]))
        ret[iv] = rs
    print("iterated", len(unis), "uniserials", f"{time.time()-t0:.0f}s", flush=True)

    def at_m(iv):
        d = dict(ret.get(iv, []))
        return d.get(m)

    def first(iv):
        return ret[iv][0] if ret.get(iv) else None

    for iv in unis:
        if at_m(iv) is not None:
            p, q = first(iv)
            ln = links.get(iv, [])
            print("  returns", iv, "first (p,q) =", (p, q), "tau^p ~ [q-p] with q-p =", q - p, "at m:", at_m(iv))
    # nu-orbits (up to shift) among uniserials: union-find over the links
    par = {iv: iv for iv in unis}
    def fd(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for iv, ln in links.items():
        for _, Y, _ in ln: par[fd(iv)] = fd(Y)
    orb = {}
    for iv in unis: orb.setdefault(fd(iv), []).append(iv)
    print("  nu-orbits up to shift among uniserials (k <= m):", len(orb))
    for o in sorted(orb.values(), key=len, reverse=True):
        fr = sorted({first(iv) for iv in o if at_m(iv) is not None})
        print("    size", len(o), "returning at m", sum(1 for iv in o if at_m(iv) is not None), "first returns", fr,
              "members", " ".join("%d-%d" % iv for iv in o))

    def hom_to_A(M):
        return any(NU["std_ok"](M, alg.proj(i)) for i in range(1, n + 1))

    def cond(M):
        x, e = M
        if e >= U[x]: return False                       # projective
        pd1 = (U[e + 1] == U[x])                          # kernel [e+1, U_x] = P_{e+1}
        return pd1 and not hom_to_A(M)

    def cond_dual(N):
        x, e = N
        if x <= alg.L[e]: return False                     # injective
        id1 = (alg.L[x - 1] == alg.L[e])                   # cokernel [L_e, x-1] = I_{x-1}
        return id1 and not any(NU["std_ok"](alg.inj(i), N) for i in range(1, n + 1))

    def test(Ms, label):
        """Ms: list of right ends M (AR sequence tau M -> E -> M); all members of the sequence are checked."""
        agree = dis = swap = 0; bad = []
        for M in Ms:
            x, e = M
            tM = (x + 1, e + 1)
            E = [(x, e + 1)] + ([(x + 1, e)] if x < e else [])
            members = [M, tM] + E
            if all(at_m(Y) is None for Y in members): continue
            sM = at_m(M) if at_m(M) is not None else at_m(tM)
            if sM is None: sM = [at_m(Y) for Y in E if at_m(Y) is not None][0]
            for Y in members:
                if at_m(Y) == sM: agree += 1
                else: dis += 1; bad.append((M, Y, sM, at_m(Y)))
        print(f"  {label}: {len(Ms)} sequences, members agree {agree}, disagree {dis} (sequences with no returning member skipped)")
        for b in bad[:12]: print("     disagree: M", b[0], "neighbour", b[1], "shift(M)", b[2], "shift(nb)", b[3])
        return dis

    nonproj = [M for M in unis if M[1] < U[M[0]]]
    good = [M for M in nonproj if cond(M)]
    # condition check: H(nu M) = tau_A M in degree -1
    okc = sum(1 for M in good if taucheck[M] == {(v, -1): 1 for v in range(M[0] + 1, M[1] + 2)})
    print(f"  condition (pd<=1, Hom(M,A)=0): {len(good)} of {len(nonproj)} non-projective; H(nu M) = tau_A M[1] on {okc} of {len(good)}")
    naive_ok = sum(1 for M in nonproj if M not in good and taucheck[M] == {(v, -1): 1 for v in range(M[0] + 1, M[1] + 2)})
    print(f"  (outside the condition H(nu M) = tau_A M[1] anyway on {naive_ok} of {len(nonproj) - len(good)})")
    gd = [(N[0] - 1, N[1] - 1) for N in unis if cond_dual(N)]     # right ends tau^-1 N of sequences starting at N
    okd = sum(1 for M in gd if taucheck[M] == {(v, -1): 1 for v in range(M[0] + 1, M[1] + 2)})
    print(f"  dual condition (id<=1, Hom(DA,N)=0): {len(gd)}; H(nu tau_A^-1 N) = N[1] on {okd} of {len(gd)}")
    good = sorted(set(good) | set(gd))
    test(good, "AR-triangle test (condition or dual condition)")
    test(nonproj, "naive control (all non-projective)")
    print("done", f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
