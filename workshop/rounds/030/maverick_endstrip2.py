"""S-1 rerun of E-112's free-end deletion with E-115-corrected class labels (unresolved cospectral LNAs placed by the F-047 profile
where it gives a single class). K0 = min free-run length at the end that is deleted from. Reports, per K0, source classes whose
image classes (over all LNAs, both ends) are one class. n=11: sources/images unresolved after the profile are dropped (counted).
usage: python workshop/rounds/030/maverick_endstrip2.py N"""
import sys, os, collections
sys.path.insert(0, "workshop/rounds/027"); sys.path.insert(0, "workshop/rounds/029")
import maverick_classes as mc, toolsmith_snfresolve as sr
from quivermutation import piecewiseHereditary as pwh
def corrected(n):
    lab, bad = mc.classes(n); lab = dict(lab)
    keyset = {lab[m][0] for m in bad}; p2c = collections.defaultdict(set)
    for m, l in lab.items():
        if l[0] in keyset and l[1] != '?' and isinstance(l, tuple) and isinstance(l[1], str): p2c[(l[0], sr.profile(n, m))].add(l)
    left = 0
    for m in bad:
        cs = p2c.get((lab[m][0], sr.profile(n, m)))
        if cs and len(cs) == 1: lab[m] = next(iter(cs))
        else: left += 1
    return lab, left
if __name__ == '__main__':
    n = int(sys.argv[1]); lab, l1 = corrected(n); lab2, l2 = corrected(n - 1)
    print('n', n, 'unresolved after profile', l1, '; n-1', l2)
    unres = lambda l: isinstance(l, tuple) and l[1] == '?'
    def sp(m): return [(i + 1, i + 1 + a) for i, a in enumerate(m) if a]
    members = collections.defaultdict(list)
    for m, l in lab.items():
        if not unres(l): members[l].append(m)
    for K0 in (1, 2, 3, 4, 5):
        ok = bad = used = dropped = 0; detail = []
        for c, mem in sorted(members.items(), key=lambda t: -len(t[1])):
            S = set(); k = 0
            for m in mem:
                R = sp(m)
                if not R: continue
                first = min(s for s, e in R); last = max(e for s, e in R)
                for cond, v in ((first - 1 >= K0, 1), (n - last >= K0, n)):
                    if cond:
                        t = lab2[tuple(pwh.removeVertex(n, list(m), v)[1])]
                        if unres(t): dropped += 1
                        else: S.add(t); k += 1
            if k == 0: continue
            used += k
            if len(S) == 1: ok += 1
            else: bad += 1; detail.append((str(c)[-24:], len(mem), k, len(S)))
        print('K0 =', K0, 'classes', ok + bad, 'transported', ok, 'not', bad, 'ends used', used, 'dropped(unresolved image)', dropped, detail)
