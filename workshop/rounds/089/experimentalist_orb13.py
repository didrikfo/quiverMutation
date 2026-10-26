"""Round 089, S-1: forward orbit size of B_n (row 0..0 3 0) under orbit+free+mirror moves, n = 6..N; test n 2^(n-5).
Also: A_n orbit size and whether movesJoin(A_n, B_n or its mirror) finds a path.
usage: timeout 10m .venv/bin/python workshop/rounds/089/experimentalist_orb13.py [Nmax=13]"""
import sys, time
from quivermutation import freeMoves as fm
N = int(sys.argv[1]) if len(sys.argv) > 1 else 13
A = lambda n: tuple([n - 1] + [0] * (n - 3))
B = lambda n: tuple([0] * (n - 4) + [3, 0])
for n in range(6, N + 1):
    t = time.time()
    wa = fm.orbitReport(n, A(n), free=True, limit=200000)
    wb = fm.orbitReport(n, B(n), free=True, limit=200000)
    m = fm.mirrorRow(n, B(n))
    j = A(n) in wb.rows or A(n) in set(fm.mirrorRow(n, r) for r in wb.rows)
    print(n, 'orbitA', len(wa.rows), wa.stoppedBy, 'orbitB', len(wb.rows), wb.stoppedBy, 'pred', n * 2 ** (n - 5), 'A in orbit(B)+mirror', j, 'mirrorclosed', m in wb.rows, '%.1fs' % (time.time() - t), flush=True)
