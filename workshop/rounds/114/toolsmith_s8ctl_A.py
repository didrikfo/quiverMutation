"""Round 114 (toolsmith): A-side random sample of single-edge corruptions; for the step-8 chain. Single-edge corruptions and end swaps, split by failing test.
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
import random, time
random.seed(114); N = int(os.environ.get('N', 10)); T0 = time.time()
cands = [(e, k, v) for e in range(len(mvA)) for k in 'FR' for v in sorted(A.quiver.nodes) if (k, v) != mvA[e]]
random.shuffle(cands); tot = acc = 0; res = []
for tam in cands[:N]:
    if time.time() - T0 > 480: break
    o = verify(A, B, mvA, mvB, tam, None)[0]; tot += 1; acc += o
    print(tam, 'ACCEPT' if o else 'reject', sorted(VERIFY_FAILS), round(time.time() - T0), flush=True)
print('A-side sample', tot, 'accepted', acc, flush=True)
