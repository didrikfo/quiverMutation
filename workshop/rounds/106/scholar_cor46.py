"""Independent check of Fosse 2305.06642 Cor. 4.6 (LaTeX read, v1) against Thm 4.2.
Each of the 8 algebras of Cor. 4.6 is sent to its quipu (Thm 4.2: k_0=n_0,
k_i=n_i+1-n_{i-1}-l_{i-1}, k_{r+1}=n+1-n_r-l_r, m_i=l_i-2) and the quipu TREES are
compared up to graph isomorphism. Two variants of the printed formulas are run:
  printed : right column rows 3,4 last length n_0+l_0  (as printed)
  fixed   : last length n_0+2 (E-275's correction)
Controls: the 6 forms that are not suspect must agree (YES); `printed` must fail
admissibility (n_r+l_r<=n) or disagree on some inputs (it can say NO).
Also A_n = rad^(n-1) vs B_n (r=0, n_0=n-3, l=3): same quipu for n>=5.
usage (repo root): python -I workshop/rounds/106/scholar_cor46.py"""
import itertools, random, sys
import networkx as nx

def quipu_tree(n, st, ln):
    r = len(st) - 1
    k = [st[0]] + [st[i] + 1 - st[i-1] - ln[i-1] for i in range(1, r+1)] + [n + 1 - st[r] - ln[r]]
    m = [l - 2 for l in ln]
    if min(k) < 0 or min(m) < 0: return None
    G = nx.Graph(); c = 0; prev = None
    def new():
        nonlocal c; c += 1; G.add_node(c); return c
    def chain(v, cnt):
        for _ in range(cnt):
            w = new()
            if v is not None: G.add_edge(v, w)
            v = w
        return v
    v = chain(None, k[0]) if k[0] else None
    for i in range(r+1):
        b = new()
        if v is not None: G.add_edge(v, b)
        chain(b, m[i])
        v = chain(b, k[i+1])
    assert G.number_of_nodes() == n, (G.number_of_nodes(), n)
    return G

def admissible(n, st, ln):
    r = len(st) - 1
    if min(st) < 1 or any(l < 2 for l in ln): return False
    if st[-1] + ln[-1] > n: return False
    for i in range(r):
        if not (st[i] < st[i+1] and st[i]+ln[i] < st[i+1]+ln[i+1]): return False
        if st[i+1] < st[i] + ln[i] - 1: return False
    return True

def forms(n, st, ln, fixed):
    r = len(st) - 1; n0, l0, nr, lr = st[0], ln[0], st[-1], ln[-1]
    mir_st = [n - st[i] - ln[i] + 1 for i in range(r, -1, -1)]   # n-n_i-l_i+1, reversed
    rev_ln = ln[::-1]
    lastE = n0 + 2 if fixed else n0 + l0                          # the suspect entry
    tail = n - nr - lr + 3
    return {
     'L1': (st, ln),
     'R1': (mir_st, rev_ln),
     'L2': (st, ln[:-1] + [tail]),
     'R2': ([lr - 2] + mir_st[1:], [tail] + rev_ln[1:]),
     'L3': ([l0 - 2] + st[1:], [n0 + 2] + ln[1:]),
     'R3': (mir_st, rev_ln[:-1] + [lastE]),
     'L4': ([l0 - 2] + st[1:], [n0 + 2] + ln[1:-1] + [tail]),
     'R4': ([lr - 2] + mir_st[1:], [tail] + rev_ln[1:-1] + [lastE]),
    }

def run(fixed, trials=4000, seed=1):
    rnd = random.Random(seed); tot = bad_adm = bad_iso = 0; per = {}
    for _ in range(trials):
        r = rnd.randint(0, 3); n = rnd.randint(8, 22)
        st, ln = [], []; pos = rnd.randint(1, 4)
        ok = True
        for i in range(r + 1):
            l = rnd.randint(3, 6)
            st.append(pos); ln.append(l); pos = pos + l - 1 + rnd.randint(0, 3)
        if not admissible(n, st, ln): continue
        tot += 1
        base = quipu_tree(n, st, ln)
        for name, (s2, l2) in forms(n, st, ln, fixed).items():
            if not admissible(n, s2, l2):
                bad_adm += 1; per[name] = per.get(name, 0) + 1; continue
            T = quipu_tree(n, s2, l2)
            if T is None or not nx.is_isomorphic(base, T):
                bad_iso += 1; per[name+'-iso'] = per.get(name+'-iso', 0) + 1
    return tot, bad_adm, bad_iso, per

for fixed in (False, True):
    print('fixed' if fixed else 'printed', run(fixed))
# A_n vs B_n
for n in range(4, 13):
    A = quipu_tree(n, [1], [n-1]); B = quipu_tree(n, [n-3], [3])
    print('n', n, 'quipu A_n ~ B_n:', nx.is_isomorphic(A, B))
