"""Round 071 T4: independent re-check of floating width-w VERIFIED_MOVES at one length.
usage: skeptic_verify.py LENGTH WIDTH
Own code: LNA enumeration, window matching, rewrite application, Cartan/Coxeter polynomial.
Package code used only for the rule table and for the mutation engine (sequence legality + result)."""
import sys, itertools
from functools import lru_cache
import sympy as sp
from quivermutation import lnaMoves as lm, nakayama, mutation, pathAlgebra, lines
L, W = int(sys.argv[1]), int(sys.argv[2])

def admissible(n, r):  # r[i] = arrows of relation starting at vertex i+1 (0 = none), len n-2 .. we use len n-1 tolerant
    rel = [(i+1, a) for i, a in enumerate(r) if a]
    if any(a < 2 or s + a > n for s, a in rel): return False
    return all(s1 + a1 < s2 + a2 for (s1, a1), (s2, a2) in zip(rel, rel[1:]))

def enum(n):  # all admissible relation-length rows, own recursion
    out = []
    def rec(i, prev_end, row):
        if i == n - 2:
            out.append(tuple(row)); return
        rec(i + 1, prev_end, row + [0])
        for a in range(2, n - i):
            if i + 1 + a > n: break
            if i + 1 + a > prev_end: rec(i + 1, i + 1 + a, row + [a])
    rec(0, 0, [])
    return out

@lru_cache(None)
def cox(n, r):
    rel = [(i+1, a) for i, a in enumerate(r) if a]
    C = sp.zeros(n, n)
    for i in range(1, n+1):
        for j in range(i, n+1):
            if not any(i <= s and s + a <= j for s, a in rel): C[i-1, j-1] = 1
    Phi = -C.T * C.inv()
    return sp.Poly(Phi.charpoly(sp.Symbol('x')).as_expr(), sp.Symbol('x')).all_coeffs().__str__()

def match(r, rule, ws):
    w, before = rule[0], rule[1]
    inside = []
    for i, a in enumerate(r):
        if not a: continue
        s = i + 1; lo, hi = s, s + a - 1
        if hi < ws or lo > ws + w - 1: continue
        if lo < ws or hi > ws + w - 1: return False
        inside.append((s - ws, a))
    return sorted(inside) == sorted(before)

def apply(n, r, rule, ws):
    w, before, after, seq = rule[:4]
    new = list(r)
    for s, a in before: new[ws + s - 1] = 0
    for s, a in after:
        if new[ws + s - 1]: return None
        new[ws + s - 1] = a
    if not admissible(n, new): return None
    return tuple(new), [v + ws - 1 if v > 0 else v - ws + 1 for v in seq]

rules = [x for x in lm.VERIFIED_MOVES if lm.anchorOf(x) is None and x[0] == W]
import os
if os.environ.get('CONTROL') == 'seq':   # sabotage: wrong sequence (first vertex shifted by one inside window)
    rules = [(x[0], x[1], x[2], tuple([(v + 1 if v > 0 else v - 1) if abs(v) < W else v for v in x[3]])) + tuple(x[4:]) for x in rules]
if os.environ.get('CONTROL') == 'after':  # sabotage: wrong claimed result, true sequence (swap after of rule i with rule i+1 of same before-size)
    rules = [(x[0], x[1], rules[(i+1) % len(rules)][2], x[3]) for i, x in enumerate(rules)]
rows = enum(L)
ref = [tuple(x) for x in nakayama.allRelationLengths(L)]
print("rules", len(rules), "LNAs", len(rows), "package LNAs", len(ref), "same set", set(rows) == set(ref), flush=True)
tot = bad = skipped = 0; steps = {}
nocheck = os.environ.get('NOWRONG') == '1'   # control only: disable the WRONG check so COX is reached
cnt = {'ILLEGAL': 0, 'WRONG': 0, 'COX': 0}
for rule in rules:
    for r in rows:
        for ws in range(1, L - W + 1):
            if not match(r, rule, ws): continue
            ap = apply(L, r, rule, ws)
            if ap is None: skipped += 1; continue
            pred, seq = ap
            alg = nakayama.LinearNakayamaAlgebra(L, list(r))
            cur = lm._copy(alg); ok = True
            if os.environ.get('CONTROL') == 'trunc': seq = seq[:-1]   # sabotage: legal but incomplete sequence
            for v in seq:  # step-by-step: admissible?, then mutate
                t = cur if v > 0 else lm._quiet(pathAlgebra.dualPathAlgebra, cur)
                if not lm._quiet(mutation.mutationIsPossibleAtVertex, t, abs(v)): ok = False; break
                cur = lm._quiet(mutation.quiverMutationAtVertices, cur, [v])
            tot += 1
            if not ok: bad += 1; cnt['ILLEGAL'] += 1; print("ILLEGAL", rule, r, ws); continue
            curL = lm.asRelLengths(cur, L)   # shared helper; None if not an LNA
            if os.environ.get('CONTROL') == 'samecount':  # Coxeter-only sabotage: the 'mutated' LNA is replaced by a different LNA with the same relation count
                k = sum(1 for a in pred if a)
                curL = next(q for q in rows[rows.index(tuple(pred)) + 1:] + rows if q != tuple(pred) and sum(1 for a in q if a) == k)
            f = []
            if curL is None: f.append('WRONG'); curL = tuple(pred)
            curL = tuple(curL)
            if curL != tuple(pred) and not nocheck: f.append('WRONG')
            if cox(L, r) != cox(L, curL): f.append('COX')   # tests the MUTATED algebra cur, not the table's pred
            for k in f: cnt[k] += 1; print(k, rule, r, ws)
            if f: bad += 1
    print(len(steps), tot, bad, flush=True); steps[len(steps)] = 1
print("DONE length", L, "width", W, "rules", len(rules), "applications", tot, "skipped_none", skipped, "failures", bad, cnt)
