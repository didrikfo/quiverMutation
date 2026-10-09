"""T6 (round 023, maverick): 'mirror chain' peeling formula, tested over ALL LNAs n = 8, 9, 10 (depth data as in maverick_variants.py).
Relation = (start s, end e = s + len). For a big relation (s, e): the left chain is the longest sequence of relations R1, R2, ... with e(R1) = s + 1,
e(R_{k+1}) = s(R_k) + 1; the right chain: s(R1) = e - 1, s(R_{k+1}) = e(R_k) - 1. depth = 1 + min over big relations of min(left, right).
Variants: 'chain' (all lengths, ambiguous choice resolved by shortest-chain-max: take the longest), 'chain2' (only length-2 links, as E-101 / theorist_peel.py).
usage: maverick_chain.py   (repository root)"""
from quivermutation import coxeterTables as ct
data = {}
for line in open('workshop/rounds/022/theorist_blocked_depths.txt'):
    n, d, _, dep = line.split(); data[(int(n), d)] = int(dep)
def formula(d, links):
    R = [(i + 1, i + 1 + int(c)) for i, c in enumerate(d) if int(c)]
    best = None
    def left(s, memo={}):  # longest chain of relations ending at s+1, then at start+1 ...
        out = 0
        for (s2, e2) in R:
            if e2 == s + 1 and links(e2 - s2): out = max(out, 1 + left(s2))
        return out
    def right(e):
        out = 0
        for (s2, e2) in R:
            if s2 == e - 1 and links(e2 - s2): out = max(out, 1 + right(e2))
        return out
    for (s, e) in R:
        if e - s < 3: continue
        v = 1 + min(left(s), right(e)); best = v if best is None else min(best, v)
    return best
if __name__ == '__main__':
    for n in (8, 9, 10):
        lst = ["".join(map(str, r)) for r in sorted(ct.lnaStatus(n))]
        for name, links in (('chain2', lambda l: l == 2), ('chain', lambda l: True)):
            for cls in ('one', 'two+'):
                tot = bad = 0; ex = []
                for d in lst:
                    nb = sum(1 for c in d if int(c) >= 3)
                    if nb == 0 or (nb == 1) != (cls == 'one'): continue
                    tot += 1; f = formula(d, links)
                    if f != data.get((n, d), 1): bad += 1; ex.append((d, data.get((n, d), 1), f))
                print(n, name, cls, 'LNAs', tot, 'mismatch', bad, ex[:8])
