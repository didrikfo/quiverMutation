#!/usr/bin/env python
"""Read a classification `classify.py` has written, without opening a spreadsheet.

    python classes.py 9                       the classes, largest first
    python classes.py 9 --collisions          where the Coxeter polynomial stops separating
    python classes.py 9 --kind "not piecewise hereditary"
    python classes.py 9 --min-size 100        the classes that hold most of the LNAs
    python classes.py 9 --quipu "P^(1,4)_(1,0,1)"
    python classes.py 9 --members "P^(1,4)_(1,0,1)"      the LNAs in one class
    python classes.py 9 --page A_9.html       the same thing as a page to browse

The table has one row per LNA and there are Catalan(n-1) of them -- 1430 at
n = 9, 16796 at n = 11, 58786 at n = 12 -- while the number of classes stays
small.  Everything here but `--members` reads the six columns it needs and
groups; nothing loads the per-LNA rows for the whole table.
"""

import argparse
import sys

from quivermutation import classview as cv


def printClasses(classification, found, verbose = False):
    """One line per class: name, size, kind, and what names it."""
    if not found:
        print("no classes match")
        return
    width = max(len(c.name) for c in found)
    for mutationClass in found:
        line = "{0:<{1}}  {2:>6}".format(mutationClass.name, width, mutationClass.size)
        if mutationClass.kind != cv.QUIPU or mutationClass.name != mutationClass.hereditaryForm:
            line += "  {0}".format(mutationClass.hereditaryForm or mutationClass.kind)
        print(line)
        if verbose:
            print("{0}  representative: {1}".format(
                " " * width, mutationClass.representative or "(no relations)"))
            print("{0}  Coxeter: {1}".format(" " * width, mutationClass.coxeterPolynomial))


def printSummary(classification, found):
    total = sum(c.size for c in found)
    print()
    print("{0} LNAs of length {1} in {2} classes".format(
        total, classification.length, len(found)))
    for kind, (classes, members) in sorted(classification.kindCounts().items()):
        print("  {0:<26} {1:>4} classes, {2:>6} LNAs".format(kind, classes, members))
    unsettled = classification.unsettled()
    if unsettled:
        print("  {0} LNAs are still unplaced -- the classification is unfinished".format(
            len(unsettled)))


def printCollisions(classification):
    """The classes no Coxeter polynomial can tell apart."""
    collisions = classification.coxeterCollisions()
    if not collisions:
        print("no two classes of length {0} share a Coxeter polynomial".format(
            classification.length))
        return
    print("{0} Coxeter polynomial(s) carried by more than one class:".format(len(collisions)))
    for polynomial, found in sorted(collisions.items()):
        print()
        print("  {0}".format(polynomial))
        for mutationClass in found:
            print("    {0}  ({1} LNAs, {2})".format(
                mutationClass.name, mutationClass.size, mutationClass.kind))


def printMembers(classification, className):
    mutationClass = classification.classNamed(className)
    if mutationClass is None:
        print("no class named {0!r}".format(className))
        return 1
    print("{0}: {1} LNAs, {2}".format(
        mutationClass.name, mutationClass.size, mutationClass.kind))
    print("Coxeter polynomial: {0}".format(mutationClass.coxeterPolynomial))
    if mutationClass.hereditaryForm:
        print("hereditary form:    {0}".format(mutationClass.hereditaryForm))
    print()
    print("{0:<40}  {1}".format("Relations", "Mutation path from the representative"))
    for relations, path, _numbering in classification.members(className):
        print("{0:<40}  {1}".format(relations or "(none)", path or "-"))
    return 0


def main(argv = None):
    parser = argparse.ArgumentParser(description = __doc__.splitlines()[0])
    parser.add_argument("length", type = int, help = "the length classify.py was run for")
    parser.add_argument("--directory", default = ".",
                        help = "where the table was written (default: here)")
    parser.add_argument("--collisions", action = "store_true",
                        help = "show only where the Coxeter polynomial fails to separate")
    parser.add_argument("--members", metavar = "CLASS",
                        help = "list the LNAs of one class, with the path to each")
    parser.add_argument("--kind", choices = cv.KINDS, help = "keep only classes of this kind")
    parser.add_argument("--quipu", metavar = "NAME",
                        help = "keep only the class with this hereditary form")
    parser.add_argument("--coxeter", metavar = "POLY",
                        help = "keep only classes with this Coxeter polynomial")
    parser.add_argument("--min-size", type = int, dest = "minSize")
    parser.add_argument("--max-size", type = int, dest = "maxSize")
    parser.add_argument("--verbose", "-v", action = "store_true",
                        help = "also print each class' representative and polynomial")
    parser.add_argument("--page", metavar = "FILE",
                        help = "write the whole classification as a page to browse")
    args = parser.parse_args(argv)

    try:
        classification = cv.Classification.forLength(args.length, args.directory)
    except FileNotFoundError as missing:
        print(missing, file = sys.stderr)
        return 2

    if args.page:
        from quivermutation import classpage
        with open(args.page, "w") as f:
            f.write(classpage.renderStandalone(classification))
        print("wrote {0}".format(args.page))
        return 0

    if args.members:
        return printMembers(classification, args.members)

    if args.collisions:
        printCollisions(classification)
        return 0

    found = classification.classes(
        kind = args.kind, quipu = args.quipu, coxeterPolynomial = args.coxeter,
        minSize = args.minSize, maxSize = args.maxSize)
    printClasses(classification, found, verbose = args.verbose)
    if not any((args.kind, args.quipu, args.coxeter, args.minSize, args.maxSize)):
        printSummary(classification, found)
    return 0


if __name__ == "__main__":
    sys.exit(main())
