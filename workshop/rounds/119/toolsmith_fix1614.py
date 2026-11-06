"""Round 119 (toolsmith), H-030: seed 1614 with the repo given coefficients -1 instead of 32002 (= P-1).
Control: with the skeptic's 32002 the repo fails (reproduces E-310); with -1 it should agree with tilting.
usage: toolsmith_fix1614.py [KIND N SEED0 COUNT]   from repo root"""
import sys, random
sys.path.insert(0, '.'); sys.path.insert(0, 'workshop/rounds/118'); sys.path.insert(0, 'workshop/rounds/091')
import skeptic_conv3 as c3
import skeptic_convention as sc
P = sc.P

def signed(c): return c - P if c > P // 2 else c
_orig = c3.repo_input
def repo_input_signed(A):
    q, rels = _orig(A)
    return q, [{k: signed(v) for k, v in r.items()} for r in rels]

def run(A, i):
    out = []
    invR = sc.invariants(A, i, 'R'); invL = sc.invariants(A, i, 'L')
    kids, gF, gR = c3.repo_children(A, i)
    hasout = any(a == i for a, _, _ in A.arrows); hasin = any(b == i for _, b, _ in A.arrows)
    for name, inv, gate, has in (('F', invR, gF, hasout), ('R', invL, gR, hasin)):
        if not has: continue
        til = inv[0]; bad = []
        if gate != til: bad.append('gate!=tilt')
        kid = kids[name]
        if kid[0] != 'err' and til and not c3.verdict(kid, inv)[0]: bad.append('nomatch')
        out.append((name, i, bad))
    return out

def sweep(kind, n, s0, cnt):
    res = {'signed': 0, 'unsigned': 0}; parents = 0
    for seed in range(s0, s0 + cnt):
        A = c3.gen(kind, n, random.Random(seed))
        if A is None: continue
        parents += 1
        for label, fn in (('unsigned', _orig), ('signed', repo_input_signed)):
            c3.repo_input = fn
            for i in range(n):
                if any(a == b == i for a, b, _ in A.arrows): continue
                for name, v, bad in run(A, i):
                    if bad:
                        res[label] += 1
                        if seed == 1614 or label == 'signed': print(label, 'seed', seed, name, v, bad)
    c3.repo_input = _orig
    print('parents', parents, 'failing vertex/side cases', res)

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: a = ['N', '5', '1614', '1']
    sweep(a[0], int(a[1]), int(a[2]), int(a[3]))
