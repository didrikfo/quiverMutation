"""T1 census: the reduced-walk orbit of every offset of every single-cluster core.

Rebuilt from the recipe in `experimentalist.md` (the round's own script was in a
scratchpad and was never committed). For each core of `--max-word` with a
placement at length `n`, every offset `o` with a row `batch._rowFor(n, c, o)`:
walk the orbit under the reduced walk to closure or `--limit`, and record

  held    -- offsets whose (reduced) row lies in the orbit;
  mirrors -- offsets `p` whose mirror row `freeMoves.mirrorRow(n, c@p)` does.

That keeps both readings of H-021's "mirror": the loose one (the orbit holds
the mirror of *some* placement) and the strict one (the mirror of a placement
at a *different* offset from those the orbit holds). One JSON line per core, so
the output is resumable: a core already in `--out` is skipped.

    .venv/bin/python workshop/rounds/001/experimentalist_census.py 14 --shard 0/4
"""

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))))

import batch
from quivermutation import freeMoves


def census(n, word, limit):
    rows = {}
    for offset in range(0, n):
        row = batch._rowFor(n, word, offset)
        if row is not None:
            rows[offset] = tuple(row)
    if not rows:
        return None
    reduced = {o: freeMoves._startOf(r, freeMoves.REDUCED) for o, r in rows.items()}
    mirrored = {o: freeMoves._startOf(freeMoves.mirrorRow(n, r), freeMoves.REDUCED)
                for o, r in rows.items()}
    orbits, done = [], set()
    for offset in sorted(rows):
        if offset in done:
            continue
        walk = freeMoves.orbitReport(n, rows[offset], free = freeMoves.REDUCED,
                                     limit = limit)
        held = sorted(o for o in rows if reduced[o] in walk.rows)
        mirrors = sorted(o for o in rows if mirrored[o] in walk.rows)
        done.update(held)
        orbits.append(dict(held = held, mirrors = mirrors, size = len(walk.rows),
                           closed = walk.closed))
    return dict(n = n, core = word, offsets = sorted(rows), orbits = orbits)


def main():
    parser = argparse.ArgumentParser(description = __doc__.split("\n")[0])
    parser.add_argument("n", type = int)
    parser.add_argument("--max-word", type = int, default = 4)
    parser.add_argument("--max-arrows", type = int, default = 6)
    parser.add_argument("--limit", type = int, default = 1500000)
    parser.add_argument("--shard", default = "0/1", help = "k/K: every Kth core from the kth")
    parser.add_argument("--cores", default = "", help = "comma-separated words; default all")
    parser.add_argument("--out", default = None)
    # overnight.py passes this to every job: stop between cores once it is
    # spent and exit 2 ("out of budget"); rerunning resumes from --out.
    parser.add_argument("--budget-hours", type = float, default = 0, dest = "budgetHours")
    args = parser.parse_args()

    k, total = (int(part) for part in args.shard.split("/"))
    words = (args.cores.split(",") if args.cores
             else batch._singleCores(args.max_word, args.max_arrows, False))
    words = words[k::total]
    out = args.out or "logs/t1-census-n%d-w%d-shard%dof%d.jsonl" % (
        args.n, args.max_word, k, total)
    os.makedirs(os.path.dirname(out) or ".", exist_ok = True)
    finished = set()
    if os.path.exists(out):
        with open(out) as handle:
            finished = {json.loads(line)["core"] for line in handle if line.strip()}
    print("n = %d, %d cores in this shard, %d already done -> %s"
          % (args.n, len(words), len(finished & set(words)), out), flush = True)
    started = time.time()
    with open(out, "a") as handle:
        for word in words:
            if word in finished:
                continue
            if args.budgetHours and time.time() - started > args.budgetHours * 3600:
                print("stopped on the budget; rerun to resume", flush = True)
                sys.exit(2)
            began = time.time()
            record = census(args.n, word, args.limit)
            if record is None:
                continue
            record["seconds"] = round(time.time() - began, 1)
            handle.write(json.dumps(record) + "\n")
            handle.flush()
            print("%s  %s  %.1fs" % (word, " ".join(
                "{%s}%d%s" % (",".join(map(str, o["held"])), o["size"],
                              "" if o["closed"] else " CAP") for o in record["orbits"]),
                record["seconds"]), flush = True)


if __name__ == "__main__":
    main()
