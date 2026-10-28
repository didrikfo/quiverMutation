"""Round 094 response: per-step counters. usage: skeptic_conv2_local.py N SEEDS
For each (parent, vertex i, side): is i ON a zero relation of arrow length>=3 (vertex-length>=4), ADJACENT to one (arrow-neighbour of a vertex on it, not on it), or FAR; gate vs tilting split by that class and by relation-length."""
import sys, os, random
sys.path.insert(0, '.'); sys.path.insert(0, '/home/user/quiverMutation/workshop/rounds/091')
import skeptic_convention as sc
from skeptic_convention2 import random_alg2, gateR
from collections import Counter
def sweep(n, seeds):
    cnt = Counter(); bad = []
    for seed in range(seeds):
        A = random_alg2(n, random.Random(seed))
        if A is None: continue
        longs = [z for z in A.zero if len(z) >= 4]
        on = set(v for z in longs for v in z)
        adj = set(b for a, b, _ in A.arrows if a in on) | set(a for a, b, _ in A.arrows if b in on)
        ML = max([len(z) for z in A.zero] + [0])
        for i in range(n):
            hasout = any(a == i for a, _, _ in A.arrows); hasin = any(b == i for _, b, _ in A.arrows)
            if not (hasout or hasin): continue
            cl = 'ON' if i in on else ('ADJ' if i in adj else ('far(long parent)' if longs else 'no long rel'))
            try: out, kids, invR, invL = sc.run_case(A, i)
            except Exception as e: cnt['crash'] += 1; continue
            for side, has, gate, tilt, key in (('F', hasout, out['gate'], out['tiltR'], 'vs T_R'), ('R', hasin, None, out['tiltL'], 'vs T_L')):
                if not has: continue
                g = gate if side == 'F' else gateR(A, i)
                cnt['%s %-16s gate=%-5s tilting=%-5s' % (side, cl, g, tilt)] += 1
                if tilt and out[side] != 'err':
                    v = out[side][key]; ok = bool(v['G=A^T'] and v['N=M^T'] is not False)
                    cnt['%s %-16s MATCH %s' % (side, cl, ok)] += 1
                    if not ok: bad.append((side, seed, i))
                if cl == 'ON': cnt['%s ON by parent max vertex-length %d: gate=%s tilting=%s' % (side, ML, g, tilt)] += 1
    for k in sorted(cnt): print(cnt[k], k)
    print('mismatches:', bad[:10])
if __name__ == '__main__': sweep(int(sys.argv[1]), int(sys.argv[2]))
