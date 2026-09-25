#!/usr/bin/env python
"""Read a shape-atlas ledger: which quiver shapes the walks pass through.

    python batch.py atlas 9 --depth 4 --jobs 7          walk, and write the ledger
    python atlas.py 9 --depth 4                         read it
    python atlas.py 9 --depth 4 --validate --page logs/atlas-n9.html

Resolves every quiver in the ledger into exact keys at four levels up to
relabelling (`quivermutation/shapeKeys.py`), writes the tables as parquet beside
the ledger, and prints the hubs, the bridges, the hubs among the leftovers, the
commonest cycles through a line, and every candidate merge -- a shape two
classes share -- after replaying it.  `--validate` runs H-022's four checks.
Spec: docs/superpowers/specs/2026-09-24-shape-atlas-design.md.
"""

import argparse
import os
import sys

from quivermutation import atlasPage
from quivermutation import shapeAtlas


def main(argv = None):
    parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("length", type = int)
    parser.add_argument("--depth", type = int, default = 4)
    parser.add_argument("--sample", type = int, default = 0)
    parser.add_argument("--seed", type = int, default = 0)
    parser.add_argument("--level", type = int, default = 2, choices = (0, 1, 2, 3),
                        help = "the level hubs and cycles are counted at (default 2)")
    parser.add_argument("--top", type = int, default = 30)
    parser.add_argument("--validate", action = "store_true",
                        help = "run H-022's four checks")
    parser.add_argument("--page", default = None, help = "write the drawings here")
    args = parser.parse_args(argv)

    path = shapeAtlas.ledgerPath(args.length, args.depth, args.sample, args.seed)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print("no ledger at {0}; run `python batch.py atlas {1} --depth {2}` first".format(
            path, args.length, args.depth))
        return 1
    # Streamed, and handed to `validate` as a path: the n = 9 ledger parsed whole
    # does not fit beside the analysis (task 8a).
    tables = shapeAtlas.resolve(shapeAtlas.iterLedger(path), args.length)
    shapeAtlas.writeTables(tables, path[:-len('.jsonl')])
    shapeAtlas.report(tables, args.level, args.top, sys.stdout)

    candidates = shapeAtlas.candidateMerges(tables)
    print("\ncandidate merges: {0}".format(len(candidates)))
    for candidate in candidates:
        verdict = shapeAtlas.replay(args.length, candidate)
        print("  {0} {1} ~ {2}  via {3} / {4}: {5}".format(
            'MERGE ' if verdict['ok'] else 'FAILED', candidate['first']['cls'],
            candidate['second']['cls'], candidate['first']['path'],
            candidate['second']['path'], verdict['reason']))

    if args.validate:
        print("\nvalidation (H-022)")
        for name, outcome in shapeAtlas.validate(tables, args.length, path).items():
            print("  {0}: {1}".format(name, outcome))
    if args.page:
        title = "Shape atlas n = {0}, depth {1}".format(args.length, args.depth)
        with open(args.page, "w", encoding = "utf-8") as handle:
            handle.write(atlasPage.render(title, atlasPage.sectionsFrom(
                tables, args.level, args.top)))
        print("\nwrote {0}".format(args.page))
    return 0


if __name__ == "__main__":
    sys.exit(main())
