"""Round 075, E-210 request: are the n = 9 images of the E-190 pair, kA_9/rad^8 and the length-3 relation on 6..9,
joined by orbit + mirror moves?  Row convention: a row of an n-vertex LNA has n - 2 entries, so the n = 9 rows are
8000000 and 0000030 (E-210 wrote 80000000 / 00000300, one trailing zero too many: as 8-entry rows those are n = 10).
Both readings are run.  Tools: freeMoves.orbitReport (forward orbit incl. edge/double moves, which change n),
freeMoves.movesJoin (bidirectional), freeMoves.derivedOrbits (union-find over all rows of one n, i.e. also reverse
rewrites) + mirrorRow.   usage: timeout 10m .venv/bin/python workshop/rounds/075/experimentalist_join9.py"""
import time
from quivermutation import freeMoves as fm
def row(s): return tuple(int(c) for c in s)
def show(r): return ''.join(map(str, r))
def n_of(r): return len(r) + 2
def mir(r): return fm.mirrorRow(n_of(r), r)
def study(sa, sb):
    a, b = row(sa), row(sb); N = n_of(a); assert n_of(b) == N
    print('\n== n = %d: %s vs %s   mirror: %s %s' % (N, sa, sb, show(mir(a)), show(mir(b))))
    for free in (False, True):
        for name, x in (('a', a), ('b', b)):
            t = time.time(); w = fm.orbitReport(N, x, free=free, limit=500000)
            ns = sorted({n_of(r) for r in w.rows})
            print(' free=%d orbit(%s): %d rows, %s, n in %s, mirror-closed %s, has a/mirror a/b/mirror b: %s %s %s %s (%.1fs)' % (
                free, name, len(w.rows), w.stoppedBy, ns, {mir(r) for r in w.rows} == w.rows,
                a in w.rows, mir(a) in w.rows, b in w.rows, mir(b) in w.rows, time.time() - t))
            if len(w.rows) <= 6: print('    rows', sorted(show(r) for r in w.rows))
        j1 = fm.movesJoin(N, a, b, free=bool(free)); j2 = fm.movesJoin(N, a, mir(b), free=bool(free))
        print(' free=%d movesJoin(a,b)=%s movesJoin(a,mirror b)=%s' % (free, j1 and show(j1), j2 and show(j2)))
    for free, edges, doubles in ((False, False, False), (False, True, True), (True, True, True)):
        lnas, orbs = fm.derivedOrbits(N, free=free, edges=edges, doubles=doubles)
        comp = {r: k for k, rows in orbs.items() for r in rows}
        par = {k: k for k in orbs}
        def f(k):
            while par[k] != k: k = par[k]
            return k
        for r in lnas:
            x, y = f(comp[r]), f(comp[mir(r)])
            if x != y: par[x] = y
        sz = {}
        for r in lnas: sz[f(comp[r])] = sz.get(f(comp[r]), 0) + 1
        print(' partition free=%d edges=%d doubles=%d: %d orbits, %d with mirror; same orbit %s, orbit+mirror %s; sizes %d, %d; largest %d' % (
            free, edges, doubles, len(orbs), len(sz), comp[a] == comp[b], f(comp[a]) == f(comp[b]), sz[f(comp[a])], sz[f(comp[b])], max(sz.values())))
study('8000000', '0000030')     # n = 9
study('80000000', '00000300')   # the literal E-210 strings, n = 10
study('90000000', '00000030')   # E-190 pair, n = 10 (control: known NOT joined, E-190)
