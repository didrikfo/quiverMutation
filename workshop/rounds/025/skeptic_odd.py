"""Skeptic: rerun scholar_longsquare walk, print parent arrows/rels for J!=0 rejects with outdeg 1 and no long square.
Usage: skeptic_odd.py n cls budget"""
import sys
n, cls, bud = sys.argv[1:4]
src = open('workshop/rounds/023/scholar_longsquare.py').read()
hook = "                raw = mutation.quiverMutationAtVertex(alg, v)\n"
assert hook in src
add = ("                if r['kerdim'] and r['outdeg'] == 1 and not r['longsq']:\n"
       "                    print('ODD v=', v, 'arrows=', sorted(alg.quiver.edges()), 'rels=', alg.rels, flush=True)\n")
src = src.replace(hook, add + hook)
sys.argv = ['x', n, '--class', cls, '--budget-sec', bud]
exec(compile(src, 'sl', 'exec'))
