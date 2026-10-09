"""Round 051 (toolsmith): positive-control input for horizon 13.  A random forward J = 0 walk (F moves only, as rounds/049 controls) of length L from a class-1 LNA/dual,
written as a one-record pickle in the format of the collector ('fail' record with childObj), to be searched with toolsmith_depth13.py front/shard (child depth 7, target depth 6):
a path of length L = 13 exists by construction, so the search must report a hit of total <= 13.
Usage (repo root): toolsmith_ctrl13.py out.pkl L seed"""
import sys, pickle, random, copy
ARGV = sys.argv; OUT, L, SEED = ARGV[1], int(ARGV[2]), int(ARGV[3])
sys.argv = ['x', '/tmp/tsm/none.pkl', '1', '5', '6', '400000', 'none']
src = open('workshop/rounds/049/toolsmith_tiltpath.py').read().split("def ball(")[0]
exec(compile(src, 'tp', 'exec'))
rnd = random.Random(SEED); alg = copy.deepcopy(rnd.choice(classes[base])); path = []
for step in range(L):
    ms = [m for m in moves(alg) if m[0] == 'F']
    if not ms: break
    kind, v, alg = rnd.choice(ms); path.append(v)
print('walk length', len(path), 'path', path, 'relabelingCost', fingerprint.relabelingCost(alg.quiver), 'key', fingerprint.canonicalKey(alg) is not None)
pickle.dump([dict(kind='fail', childObj=alg, path=path)], open(OUT, 'wb'))
