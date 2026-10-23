"""Round 083 theorist, thread T5: the cancellation step of Prop. 7.6 of arXiv:1504.02617v1.
Part A  star-count structure of the printed differential (Def. 3.2): for every M-letter L of a random
        parent, partial(L) = d0(L) + d1(L), d0 star-free, d1 has exactly one star; d1 is
        sum_a a*.(a psi) for psi: j->i, +-phi a* for phi: i->m (no loop at i), 0 otherwise.
        Parents: loops and cycles allowed anywhere except LOOPS AT i (cycles through i allowed).
        Control: loops at i allowed -> the shape of d1 on a phi b^-1 changes.
Part B  the associated graded complex E = (T(U,W,c,Y,Theta), d1) on a small quiver (cycle through i):
        kernel K of E -> T(U,W,c,Y)/(phi a*) is acyclic in every weight <= WMAX (exact rank computation mod p).
Usage: .venv/bin/python workshop/rounds/083/theorist_starfilt.py [ntrials] [WMAX]
"""
import sys, random, itertools, importlib.util, os
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('opp', os.path.join(HERE, '..', '069', 'scholar_oppermann_sign.py'))
opp = importlib.util.module_from_spec(spec); spec.loader.exec_module(opp)

def nstar(w): return sum(1 for L in w if L[0] == 's')

def random_parent(rng, loops_at_i):
    q = opp.random_instance(rng)
    if not loops_at_i:
        # drop every loop at i (vertex 0) together with its differential; remove from d-values
        gone = [a for a, (s, t, d) in q.arr.items() if s == 0 and t == 0 and a not in q.A]
        for a in gone:
            del q.arr[a]; q.d.pop(a, None)
        for a in list(q.d):
            q.d[a] = {w: c for w, c in q.d[a].items() if not any(x in gone for x in w)}
    return q

def partA(ntr, loops_at_i, seed):
    rng = random.Random(seed)
    stats = dict(letters=0, bad_s=0, bad_shape=0, bad_sq=0, theta=0, deg2=0)
    for _ in range(ntr):
        q = random_parent(rng, loops_at_i)
        m = opp.M(q, {})
        for L in opp.all_letters(q):
            stats['letters'] += 1
            dl = m.partial({(L,): 1})
            d1 = {w: c for w, c in dl.items() if nstar(w) == 1}
            if any(nstar(w) > 1 for w in dl): stats['bad_s'] += 1
            _, a, p, b = L; sp, tp, dp = q.arr[p]
            exp = {}
            isloop = (sp == q.i and tp == q.i)
            if a is None and b is None and tp == q.i:
                for al in q.A: exp[(('m', al, p, None), ('s', al))] = 1
            elif a is None and b is not None and not isloop:
                stats['theta'] += 1
                exp[(('s', b), ('m', None, p, None))] = (-1) ** m.deg(('m', None, p, None))
            if isloop:
                stats['loopletters'] = stats.get('loopletters', 0) + 1
                if a is not None or b is not None: stats['loop_theta'] = stats.get('loop_theta', 0) + 1
                if d1 != exp: stats['loop_d1_differs_from_noloop_shape'] = stats.get('loop_d1_differs_from_noloop_shape', 0) + 1
            elif d1 != exp: stats['bad_shape'] += 1
            # d1^2 = 0 : the star-raising part of partial applied twice, star count 2 part
            dd = m.partial(dl)
            if any(nstar(w) == 2 for w in dd): stats['bad_sq'] += 1
            # (partial^2 = 0 overall is E-202)
    return stats

# ---------------- Part B ----------------
P = 1000003
def rank_mod(rows, ncols):
    # rows: list of dict col->val ; gaussian elimination mod P (sparse)
    piv = {}
    r = 0
    for row in rows:
        row = {c: v % P for c, v in row.items() if v % P}
        while row:
            c = min(row)
            if c in piv:
                f = row[c]; pr = piv[c]
                for cc, vv in pr.items():
                    nv = (row.get(cc, 0) - f * vv) % P
                    if nv: row[cc] = nv
                    else: row.pop(cc, None)
            else:
                inv = pow(row[c], P - 2, P)
                row = {cc: vv * inv % P for cc, vv in row.items()}
                piv[c] = row; r += 1; break
    return r

# letters: name -> (src, tgt, deg, weight); vertices i=0, j=1 (target of A), k=2, m=3
LET = {
 'as': (1, 0, 0, 1),   # a* : j -> i   (Y)
 'psi': (2, 0, 1, 2),  # psi: k -> i   (W, degree 1)
 'c':  (2, 1, 0, 1),   # c = a psi : k -> j
 'phi': (0, 3, 0, 1),  # phi: i -> m   (U_i)
 'th': (1, 3, 1, 2),   # theta = phi a^-1 : j -> m
 'z':  (3, 3, 0, 1),   # loop elsewhere
 'z2': (3, 2, 0, 1),   # m -> k : closes the cycle i -> m -> k -> i  (through i)
 'u':  (3, 1, 1, 1),   # extra U arrow m -> j
}
# differential d1 (applied-first tuple order): psi -> (c, as) [c applied first then a*], th -> sign*(as, phi)
D1 = {'psi': [(('c', 'as'), 1)], 'th': [(('as', 'phi'), 1)]}

def words(w, s=None, t=None):
    """all composable words (applied-first) of total weight w from s to t; returns list of tuples"""
    out = []
    def rec(cur, v, wt):
        if wt == w:
            if t is None or v == t: out.append(tuple(cur))
            return
        for n, (sv, tv, dg, wl) in LET.items():
            if sv == v and wt + wl <= w:
                cur.append(n); rec(cur, tv, wt + wl); cur.pop()
    for v0 in ([s] if s is not None else range(4)): rec([], v0, 0)
    return out

def deg(w): return sum(LET[x][2] for x in w)
def boundary(w, ctrl=False):
    out = {}
    for k, x in enumerate(w):
        if x in D1 and not (ctrl and x == 'th'):
            sg = (-1) ** deg(w[:k])
            for rep, c in D1[x]:
                nw = w[:k] + rep + w[k+1:]
                out[nw] = out.get(nw, 0) + sg * c
    return {a: b for a, b in out.items() if b}

def forbidden(w):   # a* applied first then phi
    return any(w[k] == 'as' and w[k+1] == 'phi' for k in range(len(w) - 1))

def partB(WMAX, ctrl=False):
    bad = 0; checked = 0; tot_kernel = 0
    for wt in range(1, WMAX + 1):
        for s in range(4):
            for t in range(4):
                ws = words(wt, s, t)
                bydeg = {}
                for w in ws: bydeg.setdefault(deg(w), []).append(w)
                # kernel K = words containing th or a forbidden pair
                Kw = {d: [w for w in L if ('th' in w) or forbidden(w)] for d, L in bydeg.items()}
                idx = {d: {w: n for n, w in enumerate(L)} for d, L in Kw.items()}
                ranks = {}
                for d, L in Kw.items():
                    if d - 1 in idx or True:
                        rows = []
                        for w in L:
                            b = boundary(w, ctrl)
                            # boundary of a kernel word lies in kernel: check
                            for u in b: assert ('th' in u) or forbidden(u), (w, u)
                            rows.append({idx.setdefault(d - 1, {}).setdefault(u, len(idx[d - 1])): c for u, c in b.items()})
                        ranks[d] = rank_mod(rows, 0)
                for d, L in Kw.items():
                    n = len(L)
                    h = n - ranks.get(d, 0) - ranks.get(d + 1, 0)
                    tot_kernel += n; checked += 1
                    if h != 0: bad += 1; print('NONZERO H(K)', wt, s, t, d, h)
    return bad, checked, tot_kernel

if __name__ == '__main__':
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    WMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    print('Part A, no loop at i (loops/cycles elsewhere, cycles through i allowed):', partA(ntr, False, 1))
    print('Part A, control with loops at i allowed:', partA(ntr, True, 2))
    print('Part B (bad, cells checked, kernel words):', partB(WMAX))
    print('Part B control, theta not killed (expect nonzero H(K)):', partB(WMAX, True)[0] > 0)
