"""Per k=4 cell: do bystanders plus pair (8:3)(9:3) jointly contain a run of three arrows?
Bystander (t,m) covers arrows t..t+m-1; pair covers 8..11."""
PAIR = set(range(8, 12))
def runs(s, n=3):
    s = sorted(s); out = []; cur = []
    for a in s:
        if cur and a == cur[-1] + 1: cur.append(a)
        else:
            if cur: out.append(cur)
            cur = [a]
    if cur: out.append(cur)
    return out
cells = [((11,2),(12,2)), ((6,2),(7,2)), ((7,2),(11,2)), ((6,2),(11,2)), ((7,2),(12,2)),
         ((7,3),(12,2)), ((10,3),(5,2)), ((6,4),(12,2))]
print("cell | bystander-only arrows, longest run | bystanders+pair longest run | bystander-only run>=3 | any 3 consecutive covered arrows containing a non-pair arrow (joint)")
for c in cells:
    B = set(a for t, m in c for a in range(t, t + m))
    Bo = B - PAIR
    rb = max(len(r) for r in runs(B)); rj = max(len(r) for r in runs(B | PAIR))
    joint3 = any(all(x in B | PAIR for x in (s, s+1, s+2)) and any(x not in PAIR for x in (s, s+1, s+2)) for s in range(0, 20))
    only3 = any(all(x in B for x in (s, s+1, s+2)) for s in range(0, 20))
    print(c, "| B", sorted(B), "longest", rb, "| B+pair longest", rj, "| bystanders-alone run3 (union of their arrows):", only3, "| joint run3 using a bystander arrow:", joint3)
