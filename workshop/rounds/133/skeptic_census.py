"""Round 133 referee response (skeptic).  Full-census power of the registered H-017b statistic, and the confound
rate as a function of how much of the response the unmodelled feature explains.

usage (from the repository root):
  .venv/bin/python workshop/rounds/133/skeptic_census.py census N [REPS] [SEED]   -> power on ALL units at n = N
  .venv/bin/python workshop/rounds/133/skeptic_census.py confound N [REPS] [SEED] -> standardised-feature confound rates
  .venv/bin/python workshop/rounds/133/skeptic_census.py eta N                    -> eta^2 of key on each feature (after maxC, nnzC)

Statistic and permutation exactly as skeptic_h017b.ftest (between/within SS ratio of OLS residuals, key labels
permuted among units, class sizes kept), vectorised.  Covariates [1, maxC, nnzC] (the members covariate is unknown
before a walk and is left out, as in the null and confound runs of skeptic_h017b)."""
import sys, os, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
MYARGV = list(sys.argv)
from skeptic_h017b import units, features, design, name


def prep(n):
    us = units(n)
    kid = {k: i for i, k in enumerate(sorted({u["key"] for u in us}))}   # indexed, never hashed
    return us, kid


def resid(y, X):
    A = np.column_stack([np.ones(len(y)), X])
    return y - A @ np.linalg.lstsq(A, y, rcond = None)[0]


def ptest(res, lab, K, perms, rng):
    """p-value of T = between/within SS of res over labels lab (0..K-1), by permutation of lab."""
    cnt = np.bincount(lab, minlength = K).astype(float); ok = cnt > 0
    tot = ((res - res.mean()) ** 2).sum()
    def T(l):
        s = np.bincount(l, weights = res, minlength = K)
        b = (s[ok] ** 2 / cnt[ok]).sum() - res.sum() ** 2 / len(res)
        return b / max(1e-12, tot - b)
    t0 = T(lab); c = sum(T(rng.permutation(lab)) >= t0 - 1e-12 for _ in range(perms))
    return (c + 1) / (perms + 1)


def census(n, reps, seed):
    us, kid = prep(n); rng = np.random.default_rng(seed)
    lab0 = np.array([kid[u["key"]] for u in us]); K = len(kid)
    print("n=%d full census: %d units, %d key classes, sizes %s" % (n, len(us), K, sorted(collections.Counter(lab0).values())))
    print("key index -> (class size, first unit, key):")
    for k, i in sorted(kid.items(), key = lambda kv: kv[1]):
        m = [u for u in us if u["key"] == k]
        print("  k%02d  %4d  %s  %s" % (i, len(m), name(m[0]["rep"]), k))
    X = np.array([[features(u)[c] for c in ("maxC", "nnzC")] for u in us], float)
    for delta in (0, 0.12, 0.17, 0.2, 0.25, 0.3):
        rej = 0
        for _ in range(reps):
            f = rng.standard_normal(K); f = (f - f.mean()) / f.std()
            y = rng.standard_normal(len(us)) + delta * f[lab0]
            rej += ptest(resid(y, X), lab0, K, 99, rng) < 0.05
        p = rej / reps
        print("  delta %.2f SD: power %.3f (+- %.3f, %d reps, 99 perms)" % (delta, p, (p * (1 - p) / reps) ** .5, reps), flush = True)


def confound(n, reps, seed):
    us, kid = prep(n); rng = np.random.default_rng(seed)
    mu = {c: np.mean([features(u)[c] for u in us]) for c in ("nrel", "sumrel", "minrel")}
    sd = {c: np.std([features(u)[c] for u in us]) for c in ("nrel", "sumrel", "minrel")}
    print("n=%d population SD over units: %s" % (n, {c: round(float(v), 3) for c, v in sd.items()}))
    designs = [("per 3", lambda: design(us, rng, 3)), ("per 5", lambda: design(us, rng, 5)), ("census", lambda: us)]
    for dl, dg in designs:
        for c in ("sumrel", "nrel", "minrel"):
            row = []
            for s_noise in (0.3, 1, 2, 4):
                rej = 0
                for _ in range(reps):
                    s = dg(); lab = np.array([kid[u["key"]] for u in s])
                    X = np.array([[features(u)[d] for d in ("maxC", "nnzC")] for u in s], float)
                    z = (np.array([features(u)[c] for u in s], float) - mu[c]) / sd[c]
                    y = z + s_noise * rng.standard_normal(len(s))
                    rej += ptest(resid(y, X), lab, len(kid), 99, rng) < 0.05
                row.append("noise %.1f (R2 %.2f): %.2f" % (s_noise, 1 / (1 + s_noise ** 2), rej / reps))
            print("  %s N=%d %-6s | %s" % (dl, len(s), c, " | ".join(row)), flush = True)


def eta(n):
    us, kid = prep(n)
    lab = np.array([kid[u["key"]] for u in us]); K = len(kid)
    X = np.array([[features(u)[d] for d in ("maxC", "nnzC")] for u in us], float)
    for c in ("sumrel", "nrel", "minrel"):
        r = resid(np.array([features(u)[c] for u in us], float), X)
        cnt = np.bincount(lab, minlength = K); ok = cnt > 0; s = np.bincount(lab, weights = r, minlength = K)
        b = (s[ok] ** 2 / cnt[ok]).sum() - r.sum() ** 2 / len(r)
        print("n=%d all %d units: eta^2(key | maxC, nnzC) of %s = %.3f" % (n, len(us), c, b / ((r - r.mean()) ** 2).sum()))


if __name__ == "__main__":
    a = MYARGV[1:]
    reps = int(a[2]) if len(a) > 2 else 200; seed = int(a[3]) if len(a) > 3 else 1
    if a[0] == "census": census(int(a[1]), reps, seed)
    elif a[0] == "confound": confound(int(a[1]), reps, seed)
    elif a[0] == "eta": eta(int(a[1]))
