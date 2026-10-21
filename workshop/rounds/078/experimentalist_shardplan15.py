"""S-1 / round 078: shard plan for the n = 15 scan (PROPOSAL, nothing is scanned).  Shard k of K takes rows with enumeration index
i in [k*ceil(R/K), (k+1)*ceil(R/K)) of nakayama.allRelationLengths(15); R = Catalan(14) = 2 674 440.
  experimentalist_shardplan15.py plan [SAMPLE]    times enumeration, key per row, key of all 13 vertex deletions per row (samples spread over the whole index range)
  experimentalist_shardplan15.py shard k K        prints index range and first row of shard k (the scan body is the toolsmith_s1witness.py logic, not run here)
Run from repo root."""
import sys, time, itertools, math
from quivermutation import coxeterTables as ct, nakayama, piecewiseHereditary as ph
n, R = 15, math.comb(28, 14) // 15
def rng(k, K): s = -(-R // K); return k * s, min(R, (k + 1) * s)
if sys.argv[1] == 'shard':
    k, K = int(sys.argv[2]), int(sys.argv[3]); lo, hi = rng(k, K)
    r = next(itertools.islice(nakayama.allRelationLengths(n), lo, None))
    print('shard', k, 'of', K, 'rows', lo, hi, hi - lo, 'first', ''.join(map(str, r)))
else:
    S = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    t = time.time(); c = sum(1 for _ in itertools.islice(nakayama.allRelationLengths(n), 200000)); te = time.time() - t
    print('enum 200000 rows %.2fs; skip to a shard start costs <= %.1fs' % (te, te / 200000 * R))
    # sample S rows evenly spaced by index in a single pass of the first 1/2 and last: use stride over a full pass is 2.7s*13 -> do stride sample
    t = time.time(); rows = []
    for i, r in enumerate(nakayama.allRelationLengths(n)):
        if i % (R // S) == 0: rows.append(tuple(r))
    print('full enumeration pass %.1fs, %d sampled rows' % (time.time() - t, len(rows)))
    t = time.time()
    for r in rows: ct.lnaCoxeterKey(n, r)
    tk = (time.time() - t) / len(rows)
    print('key %.2f ms/row' % (tk * 1e3), flush=True)
    t = time.time(); nd = 0
    for r in rows:
        for i in range(1, n + 1):
            x = ph.removeVertex(n, list(r), i)
            if x is not None: ct.lnaCoxeterKey(x[0], tuple(x[1])); nd += 1
    td = (time.time() - t) / len(rows)
    print('all deletions: %.1f images/row, %.2f ms/row (key of every image at n-1)' % (nd / len(rows), td * 1e3))
    tot = (tk + td) * R
    print('full scan key + all deletion keys: %.0f CPU-s; 24 shards: %.1f min each; %d rows/shard' % (tot, tot / 24 / 60, -(-R // 24)))
