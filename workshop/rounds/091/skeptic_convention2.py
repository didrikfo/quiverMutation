"""Round 091 response: (a) zero relations of length 2, 3 and 4 (minimal), (b) gate for R = isMutable on the opposite algebra vs T_L tilting.
usage: skeptic_convention2.py sweep N SEEDS"""
import sys, os, random
sys.path.insert(0, '.'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import skeptic_convention as sc
from collections import Counter

def paths_vs(arrows, n):
    out = []
    def rec(p):
        if len(p) >= 3: out.append(tuple(p))
        for (a, b) in arrows:
            if a == p[-1]: rec(p + [b])
    for v in range(n): rec([v])
    return out

def random_alg2(n, rng):
    arrows = [(a, b) for a in range(n) for b in range(a + 1, n) if rng.random() < 0.38]
    if len(arrows) < n - 1 or set(v for a in arrows for v in a) != set(range(n)): return None
    cand = paths_vs(arrows, n)   # all vertex sequences with >= 2 arrows
    zero = [p for p in cand if rng.random() < (0.35 if len(p) == 3 else 0.2)]
    # minimality: drop any z that contains a shorter z
    zero = [z for z in zero if not any(len(w) < len(z) and any(z[s:s + len(w)] == w for s in range(len(z) - len(w) + 1)) for w in zero)]
    return sc.Alg(n, arrows, zero)

def gateR(A, i):
    from quivermutation import procedure
    from tests.helpers import path_algebra
    pa = path_algebra(arrows=[(a + 1, b + 1) for a, b, _ in A.arrows], rels=[[tuple(v + 1 for v in z)] for z in A.zero])
    return procedure.isMutable(*procedure._opposite(pa.quiver, procedure.relationsFrom(pa)), i + 1)

def sweep(n, seeds):
    cnt = Counter(); exs = []
    for seed in range(seeds):
        A = random_alg2(n, random.Random(seed))
        if A is None: continue
        L = max([len(z) for z in A.zero] + [0]); cnt['parents with max zero-rel vertex-length %d' % L] += 1
        for i in range(n):
            hasout = any(a == i for a, _, _ in A.arrows); hasin = any(b == i for _, b, _ in A.arrows)
            if not (hasout or hasin): continue
            try: out, kids, invR, invL = sc.run_case(A, i)
            except Exception as e:
                cnt['crash ' + repr(e)[:60]] += 1; exs.append(('crash', seed, i)); continue
            longrel = any(len(z) >= 4 for z in A.zero)
            if hasout:
                cnt['F: gate %s, T_R tilting %s' % (out['gate'], out['tiltR'])] += 1
                if out['tiltR'] and out['F'] != 'err':
                    v = out['F']['vs T_R']
                    cnt['F tilting: cyclic=%s G=A^T %s N=M^T %s ; long-rel parent %s' % (v['cyc'], v['G=A^T'], v['N=M^T'], longrel)] += 1
                    if out['R_vs_L_differ']: cnt['F power (T_R,T_L differ): other-side match %s' % out['F']['vs T_L']['N=M^T,G=A^T']] += 1
                    cnt['F opp-literal (his left): End %s' % out['opp']['F'][1]] += 1
                    if not (v['G=A^T'] and v['N=M^T'] is not False): exs.append(('F', seed, i))
            if hasin:
                g = gateR(A, i)
                cnt['R: gate %s, T_L tilting %s' % (g, out['tiltL'])] += 1
                if out['tiltL'] and out['R'] != 'err':
                    v = out['R']['vs T_L']
                    cnt['R tilting: cyclic=%s G=A^T %s N=M^T %s ; long-rel parent %s' % (v['cyc'], v['G=A^T'], v['N=M^T'], longrel)] += 1
                    if out['R_vs_L_differ']: cnt['R power: other-side match %s' % out['R']['vs T_R']['N=M^T,G=A^T']] += 1
                    cnt['R opp-literal (his right): End %s' % out['opp']['R'][1]] += 1
                    if not (v['G=A^T'] and v['N=M^T'] is not False): exs.append(('R', seed, i))
    for k in sorted(cnt): print(cnt[k], k)
    print('failures to match own:', exs[:10])

if __name__ == '__main__': sweep(int(sys.argv[2]), int(sys.argv[3]))
