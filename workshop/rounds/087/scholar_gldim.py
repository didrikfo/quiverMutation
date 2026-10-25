"""Round 087 scholar: global dimension of the linear Nakayama algebra with relation string RELSTRING (digits as in scholar_rankgrowth).
Module [x,y] (top x, socle y) has projective cover P_x = [x, x+p_x-1]; syzygy [y+1, x+p_x-1] (0 if empty).
usage: .venv/bin/python workshop/rounds/087/scholar_gldim.py RELSTRING ..."""
import sys
sys.path.insert(0, "workshop/rounds/087")
import scholar_rankgrowth as rg
for s in sys.argv[1:]:
    A = rg.Alg([int(c) for c in s]); n = A.n; best = 0
    for i in range(1, n + 1):
        x, y, d = i, i, 0
        while True:
            top = x + A.p[x] - 1
            if y == top: break
            x, y, d = y + 1, top, d + 1
        best = max(best, d)
    print(s, "n", n, "gl.dim", best)
