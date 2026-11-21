"""Round 134 response (theorist): replay the `rand` samples of theorist_endt.py step for step (same rng, same gate) and print
  - per step: J_lib, Hom(T,T[-1]), Hom(T,T[+1])  -> counts of J == Hm1 and of Hm1 = Hp1 = 0 at J = 0   (referee item 3)
  - every planted KILL child whose verdict is ISO (par and gen), with parent, child relations and KILL relations (item 4)
  - every par-tagged KILL child whose verdict is NO (the referee's seed-977 counterexample type)       (item 1)
The 'dir' / 'opp' iso tests are NOT re-run here (they are in theorist_endt_rand_*.out); only KILL is.
usage (repo root):  .venv/bin/python -W ignore workshop/rounds/134/theorist_audit.py N SEED NMAX
"""
import sys, random, collections, time
sys.path.insert(0, 'workshop/rounds/134')
import theorist_endt as T
import networkx as nx
from quivermutation import mutation

def main(N, seed, nmax):
    rng = random.Random(seed); rows = 0; t0 = time.time(); cnt = collections.Counter()
    while rows < N:
        al = T.TH.randomAlg(rng)
        if al.n > nmax: continue
        pa = T.TH.toLib(al)
        for k in sorted(pa.quiver.nodes):
            if not mutation.mutationIsPossibleAtVertex(pa, k): continue
            rows += 1
            ch = T.J.child_of(pa, k)
            if ch is None: cnt['no child'] += 1; continue
            if not nx.is_directed_acyclic_graph(ch.quiver): cnt['cyclic child'] += 1; continue
            nodes, arrows, rels = T.child_data(ch)
            Jv = T.J.J_lib(pa, k)
            E = T.EndT(pa, k).table().rad(); C = E.C; TT = E.T
            Hm1 = sum(T.J.homdim(C, X[1], Y[1], -1) for X in TT for Y in TT)
            Hp1 = sum(T.J.homdim(C, X[1], Y[1], 1) for X in TT for Y in TT)
            cnt['steps with child'] += 1
            cnt['J == Hm1'] += (Jv == Hm1)
            if Jv == 0: cnt['J = 0'] += 1; cnt['J = 0 and Hm1 = Hp1 = 0'] += (Hm1 == 0 and Hp1 == 0)
            else: print('J != 0 step: n=%d k=%d J=%d Hm1=%d Hp1=%d' % (al.n, k, Jv, Hm1, Hp1), flush=True)
            R2 = T.planted(nodes, arrows, rels, 'KILL')
            if R2 is None: continue
            v = T.iso(E, nodes, arrows, R2[0]); tag = R2[1]
            cnt['KILL %s %s' % (tag, v.split(':')[0])] += 1
            if v == 'ISO' or (tag == 'par' and v.startswith('NO')):
                print('KILL %s %s (J=%d): n=%d k=%d arrows=%s rels=%s' % (tag, v, Jv, al.n, k, al.arrows, al.rels))
                print('    child arrows %s\n    child rels %s\n    KILL  rels %s' % (arrows, rels, R2[0]), flush=True)
    print('seed %d, %d gate-admitted steps, %.0f s' % (seed, rows, time.time() - t0))
    for kk, vv in sorted(cnt.items()): print('  ', vv, kk)

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
