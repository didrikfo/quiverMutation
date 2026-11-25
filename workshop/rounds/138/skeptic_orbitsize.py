"""Round 138 skeptic: largest complex met along the nu-orbit (k = 0..m) of chosen objects, using
rounds/135/maverick_nu2.py UNCHANGED (imported, globals set here). Used to compare the complex sizes the FCY size
controls reach with those the F-058 certificate objects reach.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/138/skeptic_orbitsize.py REL m OBJS
  e.g. 3333333333330 144 S5,S11,P2,P14 ; OBJS as in maverick_nu2 (S, P, U, S5, P2, U4-5)."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "135"))
import maverick_nu2 as N

MYARGV = sys.argv[1:]
rel, m, spec = MYARGV[0], int(MYARGV[1]), MYARGV[2]
N.ALG = N.Alg(rel); N.NV = N.ALG.n
t0 = time.time()
for name, iv in N.objects(N.ALG, spec):
    X = N.minimise(N.replace(N.stalk(iv), N.ALG)); N.check(X)
    sizes = [sum(len(v) for v in X.T.values())]
    for k in range(m):
        X = N.nu(X); sizes.append(sum(len(v) for v in X.T.values()))
    H = N.cohom(X); ts = {t for (_, t) in H}
    ret = len(ts) == 1 and H == {(v, list(ts)[0]): 1 for v in range(iv[0], iv[1] + 1)}
    print(f"{rel} {name} m={m}: max summands {max(sizes)} at k={sizes.index(max(sizes))}, end {sizes[-1]}, "
          f"returns {'[%d]' % -list(ts)[0] if ret else 'no'}  ({time.time()-t0:.0f}s)", flush=True)
