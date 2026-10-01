"""Shortest labelled move sequences in the reduced (derived) walk.
  .venv/bin/python workshop/rounds/013/theorist_path.py N WORD1 OFF1 WORD2 OFF2
Prints each step: row, kind of move (rule description / edge / double sequence), whether a 2-arrow spectator was added first.
"""
import sys
sys.path.insert(0, '.')
import batch
from collections import deque
from quivermutation import freeMoves as fm, lnaMoves, edgeMoves, doubleMutation

def labelled(n, cur):
    cur = fm.stripLengthTwo(cur)
    starts = [(cur, '')]
    for v in fm.addableLengthTwo(n, cur):
        d = list(cur); d[v - 1] = 2
        starts.append((tuple(d), '+2@%d ' % v))
    out = {}
    for st, pre in starts:
        for desc, w in lnaMoves.matchingWindows(n, list(st), lnaMoves.ALL_MOVES):
            ap = lnaMoves.applyAt(n, list(st), desc, w)
            if ap is None: continue
            r = fm.stripLengthTwo(tuple(ap[0]))
            if r != cur and r not in out:
                out[r] = pre + 'rule w=%d before=%s after=%s seq=%s' % (w, desc[1], desc[2], ap[1])
        for r0, seq in edgeMoves.rewritesOf(n, st):
            r = fm.stripLengthTwo(tuple(r0))
            if r != cur and r not in out: out[r] = pre + 'edge seq=%s' % (seq,)
        for r0, seq in doubleMutation.rewritesOf(n, st):
            r = fm.stripLengthTwo(tuple(r0))
            if r != cur and r not in out: out[r] = pre + 'double seq=%s' % (seq,)
    return out

def bfs(n, a, b, limit=200000):
    a = fm.stripLengthTwo(a); b = fm.stripLengthTwo(b)
    par = {a: None}; q = deque([a])
    while q:
        c = q.popleft()
        if c == b: break
        if len(par) > limit: return None
        for r, lab in labelled(n, c).items():
            if r not in par:
                par[r] = (c, lab); q.append(r)
    if b not in par:
        return 'closed orbit of %d rows, target absent' % len(par) if not q else None
    path = []; c = b
    while par[c]: path.append((c, par[c][1])); c = par[c][0]
    return a, path[::-1]

if __name__ == '__main__':
    n = int(sys.argv[1]); w1, o1, w2, o2 = sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
    a = tuple(batch._rowFor(n, w1, o1)); b = tuple(batch._rowFor(n, w2, o2))
    res = bfs(n, a, b)
    if isinstance(res, str) or res is None: print(n, w1, o1, '->', w2, o2, ':', res); sys.exit()
    a, path = res
    print('n=%d %s@%d -> %s@%d : %d steps' % (n, w1, o1, w2, o2, len(path)))
    print(' ', ''.join(map(str, a)))
    for r, lab in path: print(' ', ''.join(map(str, r)), ' <-', lab)
