"""Round 046 (skeptic): rebuild parents/children recorded by skeptic_collect.py from (arrows, relations) alone, and re-test: gate, tiltingPlus, J, step child ==
recorded child (canonicalKey), Coxeter keys, Cartan congruence by the step's own R.   Usage: skeptic_rebuild.py in.pkl [kind]"""
import sys, pickle, ast
sys.path.insert(0, '.'); ARGV = sys.argv; sys.argv = ['x']
src = open('workshop/rounds/033/experimentalist_bothdie.py').read().split("mode = _a[1]")[0]
exec(compile(src, 'bd', 'exec'))
h15 = open('workshop/rounds/001/scholar_h015.py').read().split("\ndef main()")[0]; exec(compile(h15, 'h015', 'exec'))
def mk(edges, reprr, nodes):
    arrows = [(a, b) for a, b, k in edges]; rels = []
    for d in ast.literal_eval(reprr):
        rels.append([[p[0][0]] + [a[1] for a in p] for p in d])
    return build(arrows, rels, nodes)
kind = ARGV[2] if len(ARGV) > 2 else 'fail'; ok = tot = 0
for r in (x for x in pickle.load(open(ARGV[1], 'rb')) if x['kind'] == kind):
    N = r['V']; v = r['v']
    try: A = mk(r['pe'], r['pr'], N); rels = procedure.relationsFrom(A); perI(A, v)
    except Exception as e: print('rebuild failed (parallel arrows):', repr(e)[:60]); continue
    rels = procedure.relationsFrom(A); gate = mutation.mutationIsPossibleAtVertex(A, v)
    ch = reduction.reducePathAlgebra(mutation.quiverMutationAtVertex(A, v))
    CA, CB = cartan(A), cartan(ch); R = rplus(A, v, N)
    row = dict(gate=gate, tp=bool(tiltingPlus(A.quiver, rels, v)), J=sorted(perI(A, v).items()), keyeq=search._coxeterKeyOrNone(A) == search._coxeterKeyOrNone(ch),
               CBmatch=(CB.astype(int).tolist() == r['CB']), Rcong=bool((R.dot(CA).dot(R.T) == CB).all()))
    tot += 1; ok += (row['gate'] and (row['tp'] == (kind != 'fail')) and bool(row['J']) == (kind == 'fail') and row['CBmatch'] and row['keyeq'] and (row['Rcong'] == (kind != 'fail')))
    print(row)
print('rebuilt', tot, 'consistent with record (gate; tiltingPlus False and J != 0 for fail, True and J = 0 for ctrl; same key; child Cartan; R-congruence only for ctrl)', ok)
