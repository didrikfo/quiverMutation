"""Sharp test of H-017: do below-diagonal (relations < cords) candidates with the class's Coxeter polynomial
AND Smith profile actually reach an LNA of that class by mutation search?  Whole cells (rels-cords = -2,-1),
up to K candidates per (poly, cords, rels) cell, search depth D.
Round 010 (toolsmith) extension of rounds/004/maverick_verify.py; round 013 adds node counts (visits and distinct
canonical states) to every "reached" line; the reached set is computed exactly as families.verify does.
usage: toolsmith_verify.py N DEPTH K [maxdiag] [--cand I[,J..]] [--list] [--budget-hours H]
  --list           print the numbered candidates (after the Smith filter and the K cap) and stop, no search
  --cand I,J       search only candidates with these 0-based indices (numbering as printed by --list)
  --budget-hours H exit 2 (before starting the next candidate) once H hours have elapsed; the exit message
                   names the candidates not yet searched. A candidate already running is not interrupted."""
import sys, os, collections, time
sys.path.insert(0, os.getcwd())   # families.py lives in the repository root
flags = {}
_pos = []
_it = iter(sys.argv[1:])
for a in _it:
    if a in ("--cand", "--budget-hours"): flags[a] = next(_it)
    elif a == "--list": flags[a] = True
    else: _pos.append(a)
sys.argv = [sys.argv[0]] + _pos
t0 = time.time()
import families as fm
import numpy as np, sympy
import copy
from quivermutation import search as se, pathAlgebra, lines, fingerprint
def verifyCounted(match, order, depth):
    """families.verify with a visitor: returns (reached names, visits, distinct canonical states).
    Same walk as se.linesReachedFrom(algebra, depth) (both directions, Coxeter guard on, no dedup);
    visits = nodes the undeduped walk shows the visitor (start points included), distinct = fingerprint.canonicalKey classes."""
    wanted = {name for name, _status in match['lnas']}
    algebra = qr.algebraFromMatch(match)
    visits = [0]; keys = set()
    def visit(alg, path):
        visits[0] += 1; keys.add(fingerprint.canonicalKey(alg))
    reached = {}
    for index, start in enumerate([algebra, pathAlgebra.dualPathAlgebra(algebra)]):
        collected = []
        se.mutationSearchDepthFirst(copy.deepcopy(start), depth, [], 'lines', printOutput=False, collected=collected, visitor=visit)
        for found in lines.mutationListLineCleanup(collected, printOutput=False):
            a = found[0] if index == 0 else pathAlgebra.dualPathAlgebra(found[0])
            reached[lines.relSetToString(sorted(se._renumberedLine(a).rels))] = found[1]
    names = {"".join(str(x) for x in lines.relationStringToLineRelLengths(order, r)) for r in reached}
    return names & wanted, visits[0], len(keys)
from sympy.matrices.normalforms import smith_normal_form
from quivermutation import nakayama as nk, invariants as inv
def cart(m):
    C = np.eye(n, dtype=int); succ = collections.defaultdict(list)
    for (a, b), bit in zip(m['edges'], m['orientation']):
        (t, h) = (a, b) if bit else (b, a); succ[t].append(h)
    rels = [tuple(p) for p in m['relations']]
    def bad(p): return any(p[i:i+len(r)] == r for r in rels for i in range(len(p)-len(r)+1))
    def walk(p):
        for h in succ[p[-1]]:
            q = p + (h,)
            if not bad(q): C[p[0], h] += 1; walk(q)
    for v in range(n): walk((v,))
    return C
def snf(C):
    S = smith_normal_form(sympy.Matrix((C + C.T).tolist()), domain=sympy.ZZ)
    return tuple(sorted(abs(int(S[i, i])) for i in range(n)))
from quivermutation import quipuRelations as qr, coxeterTables as ct
n, depth, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
maxdiag = int(sys.argv[4]) if len(sys.argv) > 4 else -1
res = qr.search(n, minArrows=2, statuses=(ct.NOT_QUIPU, ct.UNPLACED), keepPerKey=None)
seen = collections.Counter()
lnasnf = {}
want = {int(x) for x in flags["--cand"].split(",")} if "--cand" in flags else None
budget = float(flags["--budget-hours"]) * 3600 if "--budget-hours" in flags else None
idx = -1; skipped = []

for key, group in res['examples'].items():
    for m in group:
        cords = sum(1 for x in m['parameters'][1] if x > 0); r = len(m['relations'])
        if r - cords > maxdiag: continue
        for nm, _st in m['lnas']:
            if nm not in lnasnf: lnasnf[nm] = snf(np.array(inv.cartanMatrix(nk.LinearNakayamaAlgebra(n, [int(c) for c in nm])), dtype=int))
        if snf(cart(m)) not in {lnasnf[nm] for nm, _ in m['lnas']}: continue   # Smith-form filter (same survivors as the full F-047 profile)
        cell = (key[:3], cords, r)
        if seen[cell] >= K: continue
        seen[cell] += 1
        idx += 1
        if "--list" in flags:
            print("cand", idx, "cords", cords, "rels", r, m['quipu'], "arrows", fm.arrowsOf(m), "rels", fm.relationsAsVertices(m)); continue
        if want is not None and idx not in want: continue
        if budget is not None and time.time() - t0 > budget:
            skipped.append(idx); continue
        print("cand", idx, end=" ")
        t = time.time(); reached, nvis, ndist = verifyCounted(m, n, depth)
        print("cords", cords, "rels", r, m['quipu'], "arrows", fm.arrowsOf(m), "rels", fm.relationsAsVertices(m),
              "-> reached", sorted(reached), "nodes", nvis, "distinct", ndist, "%.0fs" % (time.time() - t), flush=True)
if skipped:
    print("BUDGET SPENT after %.2f h: candidates not searched: %s" % ((time.time() - t0) / 3600, skipped), flush=True)
    sys.exit(2)
