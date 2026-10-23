"""Round 083 theorist response. Part C: omega-monotonicity of Def. 3.2 (merges experimentalist_omega.py, larger).
Part D: Lemma A toy with TWO phi and TWO alpha (4 thetas, 2 star letters), exact ranks mod p, weight <= WMAX.
Usage: .venv/bin/python workshop/rounds/083/theorist_omega.py NTRIALS WMAX"""
import sys, random, importlib.util, os
HERE = os.path.dirname(os.path.abspath(__file__))
s = importlib.util.spec_from_file_location('sf', os.path.join(HERE, 'theorist_starfilt.py'))
sf = importlib.util.module_from_spec(s); s.loader.exec_module(sf); opp = sf.opp

def partC(ntr, seed0):
    tot = dict(parents=0, letters=0, neg_norm=0, up=0, eq_nonstar=0, star_terms=0, multiA=0, theta=0, W=0)
    for t in range(ntr):
        rng = random.Random(seed0 + t)
        q = sf.random_parent(rng, False); tot['parents'] += 1
        if len(q.A) > 1: tot['multiA'] += 1
        m = opp.M(q, {})
        for L in opp.all_letters(q):
            tot['letters'] += 1
            w0 = (L,)
            nrm = lambda w: sum(m.deg(x) for x in w)
            om = lambda w: nrm(w) + sf.nstar(w)
            if nrm(w0) < 0: tot['neg_norm'] += 1
            if L[0] == 'm' and L[1] is None and L[3] is not None: tot['theta'] += 1
            for w, c in m.partial({w0: 1}).items():
                if nrm(w) < 0: tot['neg_norm'] += 1
                if om(w) > om(w0): tot['up'] += 1
                if sf.nstar(w) == 0 and om(w) == om(w0): tot['eq_nonstar'] += 1
                if sf.nstar(w) == 1: tot['star_terms'] += 1
    return tot

# ---------- Part D ----------
LET = {  # name: (src, tgt, deg, weight); i=0, j1=1, j2=2, k=3, m1=4, m2=5
 'as1': (1, 0, 0, 1), 'as2': (2, 0, 0, 1),
 'psi': (3, 0, 1, 2), 'c1': (3, 1, 0, 1), 'c2': (3, 2, 0, 1),
 'ph1': (0, 4, 0, 1), 'ph2': (0, 5, 0, 1),
 't11': (1, 4, 1, 2), 't12': (2, 4, 1, 2), 't21': (1, 5, 1, 2), 't22': (2, 5, 1, 2),
 'z1': (4, 3, 0, 1), 'z2': (5, 3, 0, 1)}   # cycles i -> m_p -> k -> i
D1 = {'psi': [(('c1', 'as1'), 1), (('c2', 'as2'), 1)]}
for p in (1, 2):
    for q in (1, 2): D1['t%d%d' % (p, q)] = [(('as%d' % q, 'ph%d' % p), 1)]
TH = [k for k in LET if k.startswith('t')]
def words(w, s, t):
    out = []
    def rec(cur, v, wt):
        if wt == w:
            if v == t: out.append(tuple(cur))
            return
        for n, (sv, tv, dg, wl) in LET.items():
            if sv == v and wt + wl <= w:
                cur.append(n); rec(cur, tv, wt + wl); cur.pop()
    rec([], s, 0); return out
def deg(w): return sum(LET[x][2] for x in w)
def bd(w, ctrl):
    out = {}
    for k, x in enumerate(w):
        if x in D1 and not (ctrl and x in TH):
            sg = (-1) ** deg(w[:k])
            for rep, c in D1[x]:
                nw = w[:k] + rep + w[k+1:]; out[nw] = out.get(nw, 0) + sg * c
    return {a: b for a, b in out.items() if b}
def inK(w): return any(x in TH for x in w) or any(w[k].startswith('as') and w[k+1].startswith('ph') for k in range(len(w)-1))
def partD(WMAX, ctrl=False):
    bad = cells = nw = 0
    for wt in range(1, WMAX + 1):
        for s_ in range(6):
            for t_ in range(6):
                K = {}
                for w in words(wt, s_, t_):
                    if inK(w): K.setdefault(deg(w), []).append(w)
                idx = {d: {w: n for n, w in enumerate(L)} for d, L in K.items()}
                rk = {}
                for d, L in K.items():
                    rows = []
                    for w in L:
                        b = bd(w, ctrl)
                        assert all(inK(u) for u in b)
                        tg = idx.setdefault(d - 1, {})
                        rows.append({tg.setdefault(u, len(tg)): c for u, c in b.items()})
                    rk[d] = sf.rank_mod(rows, 0)
                for d, L in K.items():
                    h = len(L) - rk.get(d, 0) - rk.get(d + 1, 0); nw += len(L); cells += 1
                    if h: bad += 1
    return bad, cells, nw
if __name__ == '__main__':
    n = int(sys.argv[1]); W = int(sys.argv[2])
    print('Part C seeds 5000..:', partC(n, 5000))
    print('Part D two phi two alpha (bad, cells, words):', partD(W))
    print('Part D control theta not killed (bad cells):', partD(W, True)[0])
