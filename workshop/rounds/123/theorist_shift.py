"""Round 123 response (theorist): the referee's vertex-shift form of the (LF) Proposition, checked on letters.
Under (C) (max w(A) < min w(p), p a non-A parent letter from i) put phi(i) = c strictly between, phi = 0 elsewhere; every kM
letter L: s -> t gets weight w(L) + phi(t) - phi(s).  Claim: all > 0.  Run on the 82 LNA instances (path-length weights),
the two-A-arrow toy kQ/(a1 c) (w(a2) = 5: (C) fails, no shift; w(a2) = 1: shift), and the referee's toy
kQ/(a1 c d - a2 d, a1 e) with w(e) = 1 ((C) fails, cycles still positive) and w(e) = 2 ((C) holds).
Control: c set to max w(A) (outside the open interval) must give a letter of weight <= 0 whenever a star exists.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/123/theorist_shift.py"""
import sys, importlib.util
from fractions import Fraction as Fr
sys.path.insert(0, '.')
s_ = importlib.util.spec_from_file_location('lnap', 'workshop/rounds/123/theorist_lnaparent.py')
P = importlib.util.module_from_spec(s_); s_.loader.exec_module(P)
lf, opp, h = P.lf, P.opp, P.h

def shift_check(m, letters, wt, arr, w, i, c=None):
    A = [L for L, (s, t, g) in arr.items() if s == i and g == 0]
    nonA = [L for L, (s, t, g) in arr.items() if s == i and L not in A]
    mA, mN = max(w[a] for a in A), min((w[p] for p in nonA), default=float('inf'))
    C = mA < mN
    if c is None:
        if not C: return C, None, None
        c = Fr(mA + mN, 2) if mN != float('inf') else mA + 1
    phi = lambda v: c if v == i else 0
    sw = {L: wt(L) + phi(m.tgt(L)) - phi(m.src(L)) for L in letters}
    return C, min(sw.values()), [h.nm(L) for L in letters if sw[L] <= 0]

def toy(arr, A, d, w, i=0):
    q = lf.mk(arr, A, d, w); m = opp.M(q, {'last': -1})
    letters = opp.all_letters(q) + [('s', al) for al in sorted(q.A)]
    wt = lambda L: -q.w[L[1]] if L[0] == 's' else q.w[L[2]] + (q.w[L[1]] if L[1] else 0) - (q.w[L[3]] if L[3] else 0)
    ok, mins = P.lf_line(m, letters, wt)
    return m, letters, wt, ok

if __name__ == '__main__':
    tot = good = ctrl = 0
    for n, row in P.ROWS:
        for op in (False, True):
            nn, rels, _ = P.lna(n, row, op)
            arr, d = P.minimal_model(n, rels, op)
            w = {L: abs(t - s) for L, (s, t, g) in arr.items()}
            for i in range(1, n + 1):
                if not any(s == i and g == 0 for s, t, g in arr.values()): continue
                q, m, letters, wt, dval, theta = P.kM(arr, d, i)
                C, mn, bad = shift_check(m, letters, wt, arr, w, i)
                tot += 1; good += bool(C and mn > 0)
                mA = max(w[L] for L, (s, t, g) in arr.items() if s == i and g == 0)
                _, mn0, bad0 = shift_check(m, letters, wt, arr, w, i, c=mA)
                ctrl += bool(bad0)
                print(f'A{n} {"".join(map(str, row))} {"op" if op else "  "} i={i}: (C) {C}, min shifted letter weight {mn}, '
                      f'nonpositive {bad}; control c = max w(A): min {mn0}, nonpositive {len(bad0)}')
    print(f'[SHIFT] LNA instances {tot}: (C) and all shifted letters > 0 on {good}; control caught on {ctrl}')
    for w2 in (5, 1):
        arr = {'a1': (0, 1, 0), 'c': (1, 2, 0), 'a2': (0, 2, 0), 'f': (0, 2, 1)}
        w = {'a1': 1, 'c': 1, 'a2': w2, 'f': 2}
        m, letters, wt, ok = toy(arr, ['a1', 'a2'], {'f': {('a1', 'c'): 1}}, w)
        print(f'[TOY1] kQ/(a1 c), w(a2) = {w2}: (LF) {ok}; shift {shift_check(m, letters, wt, arr, w, 0)}')
    for we in (1, 2):
        arr = {'a1': (0, 1, 0), 'c': (1, 2, 0), 'a2': (0, 2, 0), 'd': (2, 3, 0), 'e': (1, 4, 0), 'f': (0, 3, 1), 'g': (0, 4, 1)}
        w = {'a1': 1, 'c': 1, 'a2': 2, 'd': 1, 'e': we, 'f': 3, 'g': 1 + we}
        dd = {'f': {('a1', 'c', 'd'): 1, ('a2', 'd'): -1}, 'g': {('a1', 'e'): 1}}
        m, letters, wt, ok = toy(arr, ['a1', 'a2'], dd, w)
        print(f'[TOY2] kQ/(a1 c d - a2 d, a1 e), w(e) = {we}: (LF) {ok}; letter weights {sorted((h.nm(L), wt(L)) for L in letters)}; '
              f'shift ((C), min, nonpositive) {shift_check(m, letters, wt, arr, w, 0)}')
