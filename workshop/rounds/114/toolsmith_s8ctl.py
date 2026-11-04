"""Round 114 (toolsmith): controls for the step-8 chain. Single-edge corruptions and end swaps, split by failing test.
  PART=0|1|swap .venv/bin/python workshop/rounds/114/toolsmith_s8ctl.py   (about 5-8 min)"""
import sys, os
PART = os.environ.get('PART', 'all')   # 0, 1 (corrupt side A / B) or swap
sys.argv = ['x', 'list']
src = open('workshop/rounds/113/toolsmith_c1_chain.py').read().split("t0 = time.time()\nnexp, Z")[0]
exec(compile(src, 'c1', 'exec'))
nexp, Z = load_steps(); A, B = Z[8]['alg'], Z[8]['ch']
mk = lambda t: [(w[0], int(w[1:])) for w in t.split()]
mvA, mvB = mk("R4 F1 R7 F1 F4 F3"), mk("F1 F3 R7 R4 R3 F5")
ok = verify(A, B, mvA, mvB); print('chain', ok[0], ok[1], ok[2], flush=True)
tot = acc = 0; anyf = {}; solo = {}; accl = []
for side, mv in ((0, mvA), (1, mvB)):
  if PART in ('all', str(side)):
    for e in range(len(mv)):
        for kind in 'FR':
            for v in sorted(A.quiver.nodes):
                if (kind, v) == mv[e]: continue
                tam = (e, kind, v)
                o = verify(A, B, mvA, mvB, tam if side == 0 else None, tam if side == 1 else None)[0]
                tot += 1; acc += o
                if o: accl.append((side, e, kind, v))
                else:
                    for t in VERIFY_FAILS: anyf[t] = anyf.get(t, 0) + 1
                    if len(VERIFY_FAILS) == 1: t, = VERIFY_FAILS; solo[t] = solo.get(t, 0) + 1
print('corruptions', tot, 'accepted', acc, accl, flush=True)
print('rejected', tot - acc, 'failing tests', anyf, 'ONE test only', solo, flush=True)
sw = []
for j, rj in enumerate(Z):
    if j == 8 or PART not in ('all', 'swap'): continue
    o = verify(A, B, mvA, mvB, swapB=rj['ch']); sw.append((j, bool(o[0]), tuple(sorted(VERIFY_FAILS))))
print('B swapped for another child: YES for', [j for j, o, _ in sw if o], 'of', len(sw), sw, flush=True)
