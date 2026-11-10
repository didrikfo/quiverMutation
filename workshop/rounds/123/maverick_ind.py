import sys, collections
sys.path.insert(0, 'workshop/rounds/123'); sys.path.insert(0, 'workshop/rounds/117')
import experimentalist_classlevel as E
n, h, moves = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
cls, _ = E.keyclass(n, h); S = {m for m in cls if E.Kof(m, n) >= 1}
c = collections.Counter()
for m in S:
    for r in E.neighbours(n, m, moves):
        if r not in S: continue
        a, b = m[1:], r[1:]
        if a == b: c['same image'] += 1
        elif b in E.neighbours(n-1, a, moves) or a in E.neighbours(n-1, b, moves): c['image move'] += 1
        else: c['neither'] += 1
print(n, h, moves, dict(c))
