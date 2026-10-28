"""Round 093 theorist, response to referee: (1) degree range of sf.random_parent(rng, False); (2) sign of d_1 theta vs
(-1)^deg(phi); (3) the FULL D (all star counts), not only the star==1 part: classify every term by
(star count, W-count change). gr_G (W-count filtration) keeps delta = 0 terms; expected: only theta -> +-phi a*.
Usage: .venv/bin/python workshop/rounds/093/theorist_graded.py NTRIALS"""
import sys, random, importlib.util, os, collections
HERE = os.path.dirname(os.path.abspath(__file__))
def load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
s4 = load('s4', os.path.join(HERE, 'theorist_step4.py')); sf = s4.sf; opp = sf.opp
def run(ntr):
    arrdeg = collections.Counter(); letdeg = collections.Counter(); parent_maxdeg = collections.Counter()
    thetadeg = collections.Counter(); sign_bad = 0; nth = 0
    terms = collections.Counter()       # (kind, nstar, delta) -> count over terms of the FULL D
    grbad = 0; negdelta_nonpsi = 0
    for t in range(ntr):
        q = sf.random_parent(random.Random(t), False); m = opp.M(q, {})
        for a, (sv, tv, d) in q.arr.items(): arrdeg[d] += 1
        parent_maxdeg[max(d for _, _, d in q.arr.values())] += 1
        Ui = {p for p, (sv, tv, d) in q.arr.items() if p not in q.A and sv == q.i and tv != q.i}
        for L in opp.all_letters(q):
            dg = m.deg(L); letdeg[dg] += 1
            dl = m.partial({(L,): 1})
            _, a, p, b = L
            istheta = a is None and b is not None
            kind = 'theta' if istheta else ('psi' if s4.isW(L, q) else 'other')
            gr = {}
            for w, c in dl.items():
                dlt = s4.wc(w, q) - s4.wc((L,), q)
                terms[(kind, s4.nstar(w), dlt)] += 1
                if dlt == 0: gr[w] = c
            if istheta:
                nth += 1; thetadeg[m.deg(('m', None, p, None))] += 1
                exp = {(('s', b), ('m', None, p, None)): (-1) ** m.deg(('m', None, p, None))}
                if {w: c for w, c in dl.items() if s4.nstar(w) == 1} != exp: sign_bad += 1
                if gr != exp: grbad += 1
            elif gr: grbad += 1
    print('arrow degrees (all arrows of the 600 parents):', dict(sorted(arrdeg.items())))
    print('letter degrees:', dict(sorted(letdeg.items())))
    print('max arrow degree per parent:', dict(sorted(parent_maxdeg.items())))
    print('theta letters', nth, ' deg(phi) over thetas:', dict(sorted(thetadeg.items())))
    print('d_1 theta != (-1)^deg(phi) phi a* exactly (sign included):', sign_bad)
    print('delta-0 (gr_G) part of the FULL D != expected (theta: that monomial; others: 0):', grbad)
    for k in sorted(terms): print('terms (kind, nstar, Wcount change) =', k, terms[k])
if __name__ == '__main__': run(int(sys.argv[1]) if len(sys.argv) > 1 else 600)
