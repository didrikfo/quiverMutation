"""Round 137 (toolsmith), review item 2: classify the in-class nodes the option loses on the E-316 walks (`walk e316-I step D off|on`,
then this).  Each lost node's recorded off-path is replayed from the start with `childOf`; we count its non-tilting steps
(`isTiltingAtVertex` False).  A node whose only non-tilting step is step 1 (the E-316 step, whose child E-316 places in the parent's
J = 0 component) is a tilting descendant of a derived-equivalent child: a certified loss.  Also: is the step-1 child itself lost?
usage: T10OUT=... .venv/bin/python -W ignore workshop/rounds/137/toolsmith_e316lost.py D"""
import sys, os, json, pickle, collections
sys.path.insert(0, 'workshop/rounds/137'); D = int(sys.argv[1]); sys.argv = ['x']
import toolsmith_t10walk as T
from quivermutation import fingerprint
tot = collections.Counter()
for i, (rec, root, path, step) in enumerate(T.E316):
    off = pickle.load(open(T.fpath('e316-%d' % i, 'step', D, 'off'), 'rb')); on = pickle.load(open(T.fpath('e316-%d' % i, 'step', D, 'on'), 'rb'))
    S = T.e316Start(i); v = int(step[1:]); _, c1 = T.childOf(S, v); k1 = fingerprint.canonicalKey(c1)
    lost = [k for k in off if k not in on and off[k][0]]; shape = collections.Counter()
    for k in lost:
        a, nt = S, []
        for j, u in enumerate(off[k][1]):
            if not T._tilt(a, u): nt.append(j + 1)
            a = T.childOf(a, u)[1]
        assert fingerprint.canonicalKey(a) == k
        shape['only step 1' if nt == [1] else 'none' if not nt else 'other'] += 1
    row = dict(i=i, rec=rec, off=len(off), on=len(on), L_in=len(lost), child1Lost=(k1 in off and k1 not in on), **shape)
    tot.update({x: row.get(x, 0) for x in ('L_in', 'only step 1', 'none', 'other')}); tot['child1Lost'] += row['child1Lost']
    print(json.dumps(row), flush=True)
print('TOTAL', dict(tot))
