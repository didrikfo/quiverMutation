"""Round 117 (toolsmith): backlog 2. The rest (cands[34:84]) of the round-114 A-side single-edge corruption control for the step-8 chain (same seed 114 and shuffle as workshop/rounds/114/toolsmith_s8ctl_A.py).
  SLICE=34:47 .venv/bin/python workshop/rounds/117/toolsmith_s8ctl_Arest.py   (4 shards: 34:47 47:60 60:73 73:84; about 9 min each)"""
import sys, os
PART = "all"
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
LO, HI = (int(x) for x in os.environ["SLICE"].split(":"))
print("ncands", len(cands), flush=True)
for tam in cands[LO:HI]:
    if time.time() - T0 > 560: break
    o = verify(A, B, mvA, mvB, tam, None)[0]; tot += 1; acc += o
    print(tam, 'ACCEPT' if o else 'reject', sorted(VERIFY_FAILS), round(time.time() - T0), flush=True)
print('A-side sample', tot, 'accepted', acc, flush=True)
