#!/usr/bin/env python
"""Classify the linear Nakayama algebras of a given length up to derived equivalence.

    python classify.py 8
    python classify.py 8 --depth 6 --out A_8.csv

Writes a CSV (and a parquet alongside it) with one row per LNA, giving the class
it belongs to, the mutation path from the class representative, its Coxeter
polynomial, the vertex numbering, and the quipu the class corresponds to.

The class names are the quipus themselves, in the P^(m_0,...,m_r)_(k_0,...,k_r+1)
notation of arXiv:2305.06642, canonicalised so that one quipu has one name.
"""

import argparse
import collections
import sys

import quiverMutation as qm


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("length", type=int, help="number of vertices in the line quiver")
    parser.add_argument("--depth", type=int, default=6,
                        help="starting depth of the mutation search (default 6). "
                             "It decays for later rows to keep the run affordable.")
    parser.add_argument("--resolve-depth", type=int, default=6,
                        help="depth of the extra search used to settle classes that "
                             "share a Coxeter polynomial (default 6)")
    parser.add_argument("--out", default=None, help="output CSV path")
    parser.add_argument("--quiet", action="store_true", help="only print the summary")
    args = parser.parse_args(argv)

    if args.length < 2:
        parser.error("a line quiver needs at least 2 vertices")

    table, report = qm.classifyLength(
        args.length, args.depth, args.resolve_depth, args.out, printOutput=not args.quiet)

    sizes = collections.Counter(row[1] for row in table.rows())
    print()
    print("{0} LNAs of length {1} in {2} classes".format(
        len(table), args.length, len(table.classNames())))
    for className, size in sizes.most_common():
        print("  {0:<28} {1:>6}".format(className, size))
    if report["separated"]:
        print()
        print("Classes proved distinct despite sharing a Coxeter polynomial:")
        for polynomial, classNames in report["separated"].items():
            print("  {0}: {1}".format(polynomial, sorted(classNames)))
    if report["candidate"]:
        print()
        print("NOT SETTLED -- these share a Coxeter polynomial and no hereditary")
        print("form was found for all of them. Try a larger --resolve-depth.")
        for polynomial, classNames in report["candidate"].items():
            print("  {0}: {1}".format(polynomial, sorted(classNames)))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
