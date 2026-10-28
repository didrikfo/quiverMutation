"""Round 093 theorist, T5: does gr_G of E-235 step 4 land in the hypotheses of Lemma A (E-240)?
Hypotheses checked on the printed d_1 (Def. 3.2 via workshop/rounds/069 executable) of random parents with
no loop at i (generator of round 083):
 H1 theta letters = exactly U_i x A (every composable pair (phi, a*) of X1 x Y has a theta, no extras)
 H2 d_1(theta_{phi,a}) = +-phi.a*  (a single monomial in X1.Y), with phi in U_i letters, a* in Y
 H3 X1 = {U_i letters (m,None,p,None), p: i->v, v!=i} and Y = {a*} are disjoint letter sets
 H4 W-count G: every non-theta letter has gr_G D = 0: d_1(L)=0, or L = psi (W) and every term of d_1(L)
    has W-count one lower (a* . c with c = a psi counting 0), theta keeps its W-count (0 -> 0)
Usage: .venv/bin/python workshop/rounds/093/theorist_step4.py NTRIALS"""
import sys, random, importlib.util, os
HERE = os.path.dirname(os.path.abspath(__file__))
s = importlib.util.spec_from_file_location('sf', os.path.join(HERE, '..', '083', 'theorist_starfilt.py'))
sf = importlib.util.module_from_spec(s); s.loader.exec_module(sf); opp = sf.opp

def isW(L, q):   # psi: not-A arrow ending at i, no star/inverse decoration (loops excluded by generator)
    _, a, p, b = L
    return a is None and b is None and q.arr[p][1] == q.i
def wc(w, q): return sum(1 for x in w if x[0] == 'm' and isW(x, q))
def nstar(w): return sum(1 for x in w if x[0] == 's')

def run(ntr, seed0=0):
    st = dict(parents=0, letters=0, theta=0, H1_bad=0, H2_bad=0, H3_bad=0, H4_bad=0,
              psi=0, psi_terms=0, spectators_d1_nonzero_nonpsi=0, multiA=0, noncomposable_pairs=0)
    for t in range(ntr):
        rng = random.Random(seed0 + t)
        q = sf.random_parent(rng, False); st['parents'] += 1
        if len(q.A) > 1: st['multiA'] += 1
        m = opp.M(q, {})
        Ui = [p for p, (sv, tv, d) in q.arr.items() if p not in q.A and sv == q.i and tv != q.i]
        letters = opp.all_letters(q)
        thetas = {(L[2], L[3]) for L in letters if L[1] is None and L[3] is not None}
        if thetas != {(p, a) for p in Ui for a in q.A}: st['H1_bad'] += 1
        X1 = {('m', None, p, None) for p in Ui}; Y = {('s', a) for a in q.A}
        if X1 & Y: st['H3_bad'] += 1
        for L in letters:
            st['letters'] += 1
            dl = m.partial({(L,): 1})
            d1 = {w: c for w, c in dl.items() if nstar(w) == 1}
            _, a, p, b = L
            if a is None and b is not None:
                st['theta'] += 1
                (w, c), = d1.items() if len(d1) == 1 else [((None,), 0)]
                ok = (len(w) == 2 and w[0] in Y and w[1] in X1 and w[0][1] == b and w[1][2] == p and abs(c) == 1)
                if not ok: st['H2_bad'] += 1
                if any(wc(w, q) != wc((L,), q) for w in d1): st['H4_bad'] += 1
            elif isW(L, q):
                st['psi'] += 1
                for w in d1:
                    st['psi_terms'] += 1
                    if wc(w, q) != wc((L,), q) - 1: st['H4_bad'] += 1
            elif d1:
                st['spectators_d1_nonzero_nonpsi'] += 1; st['H4_bad'] += 1
    return st

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    print(run(n))
