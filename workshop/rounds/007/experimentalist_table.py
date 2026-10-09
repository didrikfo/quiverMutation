"""Tabulate experimentalist_kd_*.txt: per word and n: s, k = n - s, raw d, d after pairing equal-size singletons o, s-o (mirror pairs, E-066) and dropping the centre.
  .venv/bin/python workshop/rounds/007/experimentalist_table.py"""
import glob, re, ast
rows = {}
for f in sorted(glob.glob("workshop/rounds/007/experimentalist_kd_*.txt")):
    for line in open(f):
        m = re.match(r"(\d+) (\d+) hi (\d+) classes (\[\[.*?\]\]) sizes (\[.*?\]) \| sums (\[.*?\]) \| k (\S+) \| d (\d+) \| big-classes (\d+) \| (\w+)", line)
        if not m: continue
        n, w, hi = int(m[1]), m[2], int(m[3])
        cl, sz, sums = ast.literal_eval(m[4]), ast.literal_eval(m[5]), ast.literal_eval(m[6])
        size = {c[0]: z for c, z in zip(cl, sz)}
        s = sums[0] if len(sums) == 1 else None
        deff = None
        if s is not None:
            sing = [c[0] for c in cl if len(c) == 1]
            left = [o for o in sing if 2 * o != s and not (s - o in sing and size[s - o] == size[o])]
            deff = len(left)
        rows[(w, n)] = (s, n - s if s is not None else None, int(m[8]), deff, m[10], len(cl), len(sums))
for w in sorted({k[0] for k in rows}):
    print(w, " ".join("n%d:%s" % (n, ("s=%s k=%s d=%s d_eff=%s %s" % rows[(w, n)][:5]) if (w, n) in rows else "-") for n in range(14, 18)))
