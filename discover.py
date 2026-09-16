#!/usr/bin/env python
"""Search for mutation rules that break the overlap the quipu theorem cannot.

    python discover.py 9 --steps 4 --jobs 8
    python discover.py 9 --steps 5 --jobs 8 --resume
    python discover.py --verify-only

Every LNA the quipu theorem misses contains an *overlapping run*: a group of
relations chained by overlaps of two or more arrows (`unplaced.py` measures
which runs actually occur, and how many rows each one blocks).  The verified
moves of `lnaMoves` do not help with those, because they slide such a run along
the quiver without ever reducing its overlap -- which is why seeding plus move
orbits stalls at 57% of the table at n = 8 and 45% at n = 9 (H-003).

So this does not look for rewrites in general.  It plants a run that is known to
block real rows in the interior of a long quiver, walks the mutation sequences
near it, and keeps only what lands on an LNA whose **worst overlap is strictly
smaller** -- a move in the direction that matters.  Everything kept is then
verified the hard way.

Verification is the whole job, and checking that a rewrite produces the
predicted LNA is not enough: a sequence containing an inadmissible step still
returns a quiver, just not a derived equivalent one.  `lnaMoves.verifyMove`
checks the result, the admissibility of every step, and that the Coxeter
polynomial does not move, at every window position of every LNA over a range of
lengths (see R-005, which is what happens without it).

Built to be left running.  Work is one (run, embedding) unit at a time, spread
over processes, and every finished unit is appended to a JSONL file, so a run
that is killed loses at most the units in flight and `--resume` skips the rest.
"""

import argparse
import collections
import json
import os
import sys
import time
from concurrent import futures

import lnaMoves
import unplaced
import quiverMutation as qm


def maximumOverlap(relLengths):
    """The largest overlap, in arrows, between consecutive relations.

    Zero or one means almost separate, which is what the quipu theorem needs.
    Reducing this is the whole point of the search.
    """
    relations = [(start, arrows)
                 for start, arrows in enumerate(relLengths, start = 1) if arrows]
    return max((earlier[0] + earlier[1] - later[0]
                for earlier, later in zip(relations, relations[1:])),
               default = 0)


def targetRuns(length, howMany):
    """The overlapping runs blocking the most unplaced rows at this length."""
    _, rows = unplaced.unplacedRows(length)
    counts = collections.Counter()
    for relationString in rows:
        relLengths = qm.relationStringToLineRelLengths(length, relationString)
        for run in set(unplaced.overlappingRuns(relLengths)):
            counts[run] += 1
    return [(run, blocked) for run, blocked in counts.most_common(howMany)]


def unitKey(run, quiverLength, offset, steps):
    return "{0}|{1}|{2}|{3}".format(
        ",".join("{0}:{1}".format(s, a) for s, a in run), quiverLength, offset, steps)


def searchOneUnit(task):
    """Walk the mutations near one planted run and keep the links that help.

    Returns a JSON-safe record.  Run in a worker process, so it takes and returns
    plain data and imports nothing the parent has not already imported.
    """
    run, quiverLength, offset, steps, margin = task
    started = time.monotonic()
    relLengths = lnaMoves.embedPattern(quiverLength, run, offset)
    record = {
        "key": unitKey(run, quiverLength, offset, steps),
        "run": [list(r) for r in run],
        "quiverLength": quiverLength,
        "offset": offset,
        "steps": steps,
        "links": [],
        "seconds": 0.0,
    }
    if relLengths is None:
        record["seconds"] = time.monotonic() - started
        record["skipped"] = "the run does not fit at this offset"
        return record

    startOverlap = maximumOverlap(relLengths)
    centreLo = offset + min(start for start, _ in run)
    centreHi = offset + max(start + arrows - 1 for start, arrows in run) - 1
    reached = lnaMoves.localMutationSequences(
        quiverLength, relLengths, centreLo, centreHi, steps, margin)

    for name, sequence in reached.items():
        after = [int(c) for c in name]
        overlap = maximumOverlap(after)
        if overlap >= startOverlap:
            continue                      # no progress: still just as tangled
        description = lnaMoves.describeLink(quiverLength, relLengths, after, sequence)
        record["links"].append({
            "before": lnaMoves.className(relLengths),
            "after": name,
            "sequence": list(sequence),
            "overlapBefore": startOverlap,
            "overlapAfter": overlap,
            # A link that is not a local rewrite cannot become a table rule, but
            # it is still a mutation path from a blocked LNA to a less tangled
            # one, which is worth keeping.
            "description": None if description is None else [
                description[0],
                [list(r) for r in description[1]],
                [list(r) for r in description[2]],
                list(description[3]),
            ],
        })
    record["seconds"] = time.monotonic() - started
    return record


def paddedRule(rule, left, right):
    """The same rewrite with a wider window, pinning the arrows beside it empty.

    `describeLink` returns the *smallest* window containing everything a link
    touches, and that is often too small: `matchesAt` requires only that no
    relation outside reaches *into* the window, so a rewrite whose result depends
    on the arrow just beyond the edge being free will match where it should not,
    and fail with "wrong result".  Padding is the fix, and it is a strictly
    stronger precondition -- the wider rule fires in a subset of the places the
    narrow one did -- so a rule that verifies padded is sound, just less general.

    A relation at relative start `s` sits at absolute arrow `windowStart + s - 1`,
    and a signed offset `v` names vertex `windowStart + |v| - 1`, so moving the
    window `left` arrows to the left adds `left` to every start and to every
    offset's magnitude.  Padding on the right only widens the window.
    """
    width, before, after, offsets = rule
    # Negative padding is how `tightFormOf` undoes it, so nothing here may assume
    # the shift is positive.
    shift = lambda relations: tuple(
        (start + left, arrows) for start, arrows in relations)
    return (width + left + right,
            tuple(sorted(shift(before))),
            tuple(sorted(shift(after))),
            tuple(v + left if v > 0 else v - left for v in offsets))


def verifyOneCandidate(task):
    """Verify one candidate rewrite, widening its window if the tight one fails.

    Returns the first form that verifies clean, or the tight form's failure if
    none does.

    "Clean" means no failures *and* at least `minConfirmations` confirmations,
    and the second half is not a formality.  Widening a window narrows where the
    rule fires, so padding can reach a form that fires once or twice in the
    lengths being checked and passes vacuously.  That happened on the first run
    of this script: a four-mutation rewrite of the maximally overlapping pair
    passed at lengths 8 and 9 on one confirmation each, and failed at 10, 11 and
    12 as soon as there was room for it to fire properly (R-008).  A rule is only
    as good as the number of places it has been checked in.
    """
    description, lengths, maxPadding, minConfirmations = task
    width, before, after, offsets = description
    tight = (width,
             tuple(tuple(r) for r in before),
             tuple(tuple(r) for r in after),
             tuple(offsets))
    started = time.monotonic()
    attempts = [(0, 0, tight)]
    for pad in range(1, maxPadding + 1):
        for left, right in ((pad, pad), (pad, 0), (0, pad)):
            attempts.append((left, right, paddedRule(tight, left, right)))

    first = None
    for left, right, rule in attempts:
        confirmed, failures = lnaMoves.verifyMove(rule, lengths, checkCoxeter = True)
        result = {
            "rule": [rule[0], [list(r) for r in rule[1]], [list(r) for r in rule[2]],
                     list(rule[3])],
            "padding": [left, right],
            "confirmed": confirmed,
            "failures": len(failures),
            "firstFailures": [list(map(str, f)) for f in failures[:3]],
            "lengths": list(lengths),
            "seconds": time.monotonic() - started,
        }
        result["enoughConfirmations"] = confirmed >= minConfirmations
        if first is None:
            first = result
        if not failures and confirmed >= minConfirmations:
            return result
    return first


def readDone(path):
    """The unit keys already recorded, and the records themselves."""
    records = []
    if not os.path.exists(path):
        return set(), records
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except ValueError:
                continue              # a half-written last line, from a kill
    return {r["key"] for r in records if "key" in r}, records


def append(path, record):
    """Append one record and flush, so a kill cannot lose a finished unit."""
    with open(path, "a") as f:
        f.write(json.dumps(record) + "\n")
        f.flush()
        os.fsync(f.fileno())


def candidatesFrom(records):
    """The distinct local rewrites found, with where each was seen."""
    seen = {}
    for record in records:
        for link in record.get("links", []):
            if link["description"] is None:
                continue
            width, before, after, offsets = link["description"]
            key = (width,
                   tuple(tuple(r) for r in before),
                   tuple(tuple(r) for r in after),
                   tuple(offsets))
            seen.setdefault(key, []).append(record["quiverLength"])
    return seen


def runSearch(args, out):
    targets = targetRuns(args.length, args.targets)
    print("Targets: the {0} overlapping runs blocking the most unplaced rows "
          "of A_{1}".format(len(targets), args.length))
    for run, blocked in targets:
        print("  {0:<28} blocks {1} rows".format(
            " ".join("({0}:{1})".format(s, a) for s, a in run), blocked))

    embeddings = [(quiverLength, offset)
                  for quiverLength in args.quiver_lengths
                  for offset in (args.margin, args.margin + 1)]
    tasks = []
    doneKeys, _ = readDone(out)
    for run, _ in targets:
        for quiverLength, offset in embeddings:
            key = unitKey(run, quiverLength, offset, args.steps)
            if key in doneKeys:
                continue
            tasks.append((run, quiverLength, offset, args.steps, args.margin))
    print()
    print("{0} work units ({1} already done), {2} at a time, {3} mutation steps each".format(
        len(tasks), len(doneKeys), args.jobs, args.steps))
    if not tasks:
        return

    deadline = None if args.budget_hours is None else time.monotonic() + args.budget_hours * 3600
    completed = 0
    with futures.ProcessPoolExecutor(max_workers = args.jobs) as pool:
        pending = {pool.submit(searchOneUnit, task): task for task in tasks}
        try:
            for future in futures.as_completed(pending):
                record = future.result()
                append(out, record)
                completed += 1
                helpful = len(record.get("links", []))
                print("  [{0}/{1}] {2}  {3} link(s) reducing overlap  {4:.0f}s".format(
                    completed, len(tasks), record["key"], helpful, record["seconds"]),
                    flush = True)
                if deadline is not None and time.monotonic() >= deadline:
                    print("\nbudget reached; stopping. Everything finished is in {0}; "
                          "continue with --resume.".format(out))
                    for f in pending:
                        f.cancel()
                    break
        except KeyboardInterrupt:
            print("\ninterrupted; {0} finished units are in {1}".format(completed, out))
            for f in pending:
                f.cancel()
            raise


def tightFormOf(result):
    """The unpadded rule a verification record came from, as a key.

    A record stores whichever form was tried last, padded or not, so the padding
    has to be undone to recognise the candidate again on a resume.
    """
    width, before, after, offsets = result["rule"]
    left, right = result.get("padding", [0, 0])
    rule = (width,
            tuple(tuple(r) for r in before),
            tuple(tuple(r) for r in after),
            tuple(offsets))
    return paddedRule(rule, -left, -right)


def runVerification(args, out):
    _, records = readDone(out)
    candidates = candidatesFrom(records)
    known = set(lnaMoves.VERIFIED_MOVES)
    fresh = {rule: places for rule, places in candidates.items() if rule not in known}

    # Verifying is itself hours of work at these lengths, so it resumes too.
    verifiedFile = out.replace(".jsonl", ".verified.jsonl")
    _, done = readDone(verifiedFile)
    alreadyVerified = {tightFormOf(result) for result in done if "rule" in result}
    outstanding = sorted(rule for rule in fresh if rule not in alreadyVerified)

    print()
    print("{0} distinct local rewrites found, {1} not already in the move table, "
          "{2} not yet verified".format(len(candidates), len(fresh), len(outstanding)))
    if not outstanding:
        return

    lengths = list(args.verify_lengths)
    print("Verifying each over lengths {0} -- result, admissibility of every "
          "step, and the Coxeter polynomial.".format(lengths))
    deadline = (None if args.budget_hours is None
                else time.monotonic() + args.budget_hours * 3600)
    verified, rejected, tooThin = [], [], []
    with futures.ProcessPoolExecutor(max_workers = args.jobs) as pool:
        tasks = [(rule, lengths, args.max_padding, args.min_confirmations)
                 for rule in outstanding]
        pending = {pool.submit(verifyOneCandidate, task): task for task in tasks}
        for future in futures.as_completed(pending):
            result = future.result()
            append(verifiedFile, result)
            if result["failures"] == 0 and result["confirmed"] >= args.min_confirmations:
                verified.append(result)
            elif not result["failures"]:
                tooThin.append(result)
            else:
                rejected.append(result)
            if deadline is not None and time.monotonic() >= deadline:
                print("budget reached during verification; {0} candidates still "
                      "unchecked. Continue with --verify-only --resume.".format(
                          len(outstanding) - len(verified) - len(rejected) - len(tooThin)))
                for f in pending:
                    f.cancel()
                break

    print()
    print("{0} verified, {1} rejected, {2} unproven (no failures but fewer than "
          "{3} confirmations).".format(
              len(verified), len(rejected), len(tooThin), args.min_confirmations))
    if tooThin:
        print()
        print("Unproven -- re-verify these over more lengths before trusting them. A "
              "rewrite that fires once or twice can pass by accident:")
        for result in tooThin[:10]:
            print("  {0}  only {1} confirmation(s)".format(
                lnaMoves.formatMove((result["rule"][0],
                                     tuple(tuple(r) for r in result["rule"][1]),
                                     tuple(tuple(r) for r in result["rule"][2]),
                                     tuple(result["rule"][3]))),
                result["confirmed"]))
    if verified:
        print()
        print("Paste into lnaMoves.VERIFIED_MOVES:")
        for result in sorted(verified, key = lambda r: -r["confirmed"]):
            width, before, after, offsets = result["rule"]
            print("    ({0}, {1}, {2}, {3}),   # {4} confirmations, no failures{5}".format(
                width,
                tuple(tuple(r) for r in before),
                tuple(tuple(r) for r in after),
                tuple(offsets),
                result["confirmed"],
                "" if result["padding"] == [0, 0]
                else ", window widened by {0}".format(result["padding"])))
    if rejected:
        print()
        print("Rejected (this is the normal outcome and the reason for the check):")
        for result in rejected[:10]:
            print("  {0}  {1} confirmed, {2} failures  {3}".format(
                lnaMoves.formatMove((result["rule"][0],
                                     tuple(tuple(r) for r in result["rule"][1]),
                                     tuple(tuple(r) for r in result["rule"][2]),
                                     tuple(result["rule"][3]))),
                result["confirmed"], result["failures"],
                result["firstFailures"][0][-1] if result["firstFailures"] else ""))


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("length", type = int, nargs = "?", default = 9,
                        help = "the length whose unplaced rows choose the targets (default 9)")
    parser.add_argument("--targets", type = int, default = 6,
                        help = "how many of the commonest blocking runs to aim at (default 6)")
    parser.add_argument("--steps", type = int, default = 4,
                        help = "longest mutation sequence to walk (default 4). H-008 is the "
                               "suspicion that the rules that matter need four or five, which "
                               "is why the three-step searches found little.")
    parser.add_argument("--margin", type = int, default = 3,
                        help = "how far from the planted run mutations may happen (default 3)")
    parser.add_argument("--quiver-lengths", type = int, nargs = "+", default = [13, 14],
                        help = "lengths to plant the run in (default 13 14). Long enough that "
                               "the run sits clear of both ends, which is H-007.")
    parser.add_argument("--verify-lengths", type = int, nargs = "+", default = [7, 8, 9, 10],
                        help = "lengths to verify a candidate over (default 7 8 9 10)")
    parser.add_argument("--jobs", type = int, default = max(1, (os.cpu_count() or 2) - 1),
                        help = "worker processes (default: one fewer than the cores)")
    parser.add_argument("--out", default = None, help = "JSONL file of finished work units")
    parser.add_argument("--resume", action = "store_true",
                        help = "skip work units already in the JSONL file")
    parser.add_argument("--budget-hours", type = float, default = None,
                        help = "stop cleanly after this many hours")
    parser.add_argument("--min-confirmations", type = int, default = 8,
                        help = "how many places a rewrite must be confirmed in before it "
                               "counts as verified (default 8). Fewer is not evidence: "
                               "see R-008.")
    parser.add_argument("--max-padding", type = int, default = 2,
                        help = "how many arrows a candidate's window may be widened by "
                               "when the tight one fails to verify (default 2). The tight "
                               "window describeLink returns is often too small to state "
                               "the rule's real precondition.")
    parser.add_argument("--verify-only", action = "store_true",
                        help = "do not search; just verify the candidates already found")
    args = parser.parse_args(argv)

    out = args.out or "discovery_A{0}_steps{1}.jsonl".format(args.length, args.steps)
    if not args.resume and not args.verify_only and os.path.exists(out):
        parser.error("{0} exists; pass --resume to continue it, or --out for a new "
                     "file".format(out))

    if not args.verify_only:
        runSearch(args, out)
    runVerification(args, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
