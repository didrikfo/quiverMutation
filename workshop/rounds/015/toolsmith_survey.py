"""Survey of reachedQuipuAlgebras members at n=8 by (path length, arrows, relations); pickles members with arrows >= n.
usage: toolsmith_survey.py N WALKDEPTH FIRST LAST   (repository root)"""
import sys, os, time, collections, pickle
sys.path.insert(0, os.getcwd())
from quivermutation import quipuRelations as qr, coxeterTables as ct, nakayama as nk
n, wd, lo, hi = [int(x) for x in sys.argv[1:5]]
status = ct.lnaStatus(n); out = {}
for k, rl in enumerate(sorted(status)):
    if k < lo or k > hi: continue
    t = time.time(); r = qr.reachedQuipuAlgebras(nk.LinearNakayamaAlgebra(n, list(rl)), wd)
    h = collections.Counter((len(p), len(c[0]), len(c[1])) for c, p in r.items())
    print("lna", k, "".join(map(str, rl)), "members", len(r), "(%.0fs)" % (time.time() - t))
    for key in sorted(h): print("   L %d arrows %d rels %d : %d" % (*key, h[key]))
    out[k] = [(len(p), c) for c, p in r.items() if len(c[0]) >= n]
    print("  arrows>=n members:", len(out[k]), flush=True)
pickle.dump(out, open("workshop/rounds/015/toolsmith_survey_n%d_%d_%d.pkl" % (n, lo, hi), "wb"))
