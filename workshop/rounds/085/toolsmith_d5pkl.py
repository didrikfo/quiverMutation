"""Round 083 (toolsmith): write one skeptic_kept.py-format pickle per key-kept depth-5 parent: toolsmith_d5pkl.py PAR5.pkl OUTDIR i1 i2 ...  (needs toolsmith_d5scan.py build output)"""
import sys, pickle
A_ = sys.argv; par = pickle.load(open(A_[1], 'rb')); OUT = A_[2]; idx = [int(a) for a in A_[3:]]
sys.argv = ['x', 'ctl', 'x']
src = open('workshop/rounds/085/toolsmith_d5scan.py').read().split("if MODE == 'build'")[0]
sys.path.insert(0, '.'); exec(compile(src, 's', 'exec'))
for i in idx:
    x, p = par[i]; ks = [(k, v) for k, v, st in jcands(x) if st == 'kept']; assert len(ks) == 1
    pickle.dump([dict(parent=x, kind=ks[0][0], v=ks[0][1], depth=5, path=p, rootalg=None)], open('%s/k%d.pkl' % (OUT, i), 'wb')); print(i, ks)
