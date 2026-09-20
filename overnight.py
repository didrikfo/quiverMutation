#!/usr/bin/env python
"""Run the long jobs overnight, and survive the night -- on any platform.

    python overnight.py                 the default jobs, 9 hours
    python overnight.py --hours 7
    python overnight.py --only merges10 classify10
    python overnight.py --run "batch.py cores 16 --max-word 4 --jobs 7"
    python overnight.py --dry-run       print what would run

`OVERNIGHT.md` is the menu: exact commands, what each parameter changes, and
what a run of each size costs.  `--run` is why that file can hold a month of
nights without this one being edited -- it takes any command as a job, gives it
the night's budget, and logs and restarts it like the named ones.

`overnight.sh` does this with caffeinate and bash job control, neither of which
exists on Windows. This does the same three things from Python:

* keeps the machine awake for the whole run (`SetThreadExecutionState` on
  Windows, `caffeinate` on macOS if present) -- closing a laptop's lid can still
  suspend it, so leave it open;
* restarts a job that dies, with the job's own resume, until it finishes or the
  budget is spent -- exit 0 and 1 are finished answers, 2 is "out of budget";
* logs every job unbuffered to `logs/<job>-<stamp>.log`.

The default jobs, aimed by research F-032 and H-013 (2026-09-17), and H-019
(2026-09-19):

    merges10    python merges.py 10 --depths 5 6 7 8     7 processes
    merges11    python merges.py 11 --depths 4 5 6       7 processes
    classify10  python classify.py 10 --resume           1 process, an independent
                                                          route to the n = 10 count
    sample14    python batch.py sample 14 --count 4000   2 processes, the leftover
    sample16    python batch.py sample 16 --count 4000     rate at lengths that
                                                           cannot be enumerated

The two sampling jobs are the ones to widen when there is machine time to spare:
raise `--count` and rerun, and the ledger makes the extra draws the only work
done. `python batch.py sample 14 --summary` reads the answer back.

The length-15 night, aimed at H-018 and H-019 (2026-09-19), which is three jobs
and is **not** in the default set -- ask for it by name, because it wants the
whole machine:

    python overnight.py --hours 9 --only sample15 cores15 cores13

    sample15    python batch.py sample 15 --count 6000    7 processes, the
                    --depth 4                               leftover rate at a
                                                            length that cannot
                                                            be enumerated, and a
                                                            mutation search out
                                                            of every leftover
    cores15     python batch.py cores 15 --max-word 4     6 processes, every
                                                            overlapping core at
                                                            every offset of a
                                                            line with room
    cores13     python batch.py cores 13 --max-word 4     2 processes, the same
                                                            census at a length
                                                            F-042 already
                                                            covered, as the
                                                            control

`cores13` is what makes `cores15` readable. A census at one length is a list; a
law is what changes between two, which is how F-042 was found and why it is
invisible at n = 9.

Stop it with Ctrl-C, which stops the children. Every job resumes from its own
checkpoint when rerun, so stopping loses only the work in flight.
"""

import argparse
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def jobs(hours):
    budget = ["--budget-hours", str(hours)]
    return {
        'merges10': ["merges.py", "10", "--depths", "5", "6", "7", "8", "--jobs", "7"] + budget,
        'merges11': ["merges.py", "11", "--depths", "4", "5", "6", "--jobs", "7"] + budget,
        'classify10': ["classify.py", "10", "--resume"] + budget,
        'sample14': ["batch.py", "sample", "14", "--count", "4000", "--jobs", "2"] + budget,
        'sample16': ["batch.py", "sample", "16", "--count", "4000", "--jobs", "2"] + budget,
        'sample15': ["batch.py", "sample", "15", "--count", "6000", "--depth", "4",
                     "--jobs", "7"] + budget,
        'cores15': ["batch.py", "cores", "15", "--max-word", "4",
                    "--jobs", "6"] + budget,
        'cores13': ["batch.py", "cores", "13", "--max-word", "4",
                    "--jobs", "2"] + budget,
    }


#: What to run in the morning to read each job's answer back, and what the
#: answer is aimed at.  Keyed by job so that `--only` prints only the relevant
#: ones: a note telling a reader to summarise a job that did not run is worse
#: than no note, because it is the sort of thing that gets pasted anyway.
MORNING = {
    'merges10': ("python merges.py 10 --summary", "H-013"),
    'merges11': ("python merges.py 11 --summary", "H-013"),
    'classify10': ("python classes.py 10", "the n = 10 class count"),
    'sample14': ("python batch.py sample 14 --summary", "H-019"),
    'sample16': ("python batch.py sample 16 --summary", "H-019"),
    'sample15': ("python batch.py sample 15 --count 6000 --depth 4 --summary",
                 "H-019, and H-017 from what the searches reached"),
    'cores15': ("python batch.py cores 15 --max-word 4 --summary",
                "H-018, and F-040's count of separated clusters"),
    'cores13': ("python batch.py cores 13 --max-word 4 --summary",
                "the control length for cores15 -- the pattern is the "
                "difference between the two"),
}


#: Set `ES_CONTINUOUS | ES_SYSTEM_REQUIRED` and then hold, because the flag is
#: per *thread*: it lasts exactly as long as the process that set it, which is
#: what makes killing that process the way to release it again.
#:
#: The flags are written as the single decimal 2147483649 and not as
#: `0x80000000 -bor 0x00000001`, which is the spelling the Windows branch below
#: uses and which **does not work here**: PowerShell reads `0x80000000` as a
#: signed Int32, so the `-bor` comes out as -2147483647 and the call fails to
#: convert its argument. Above Int32's range PowerShell reads a decimal literal
#: as Int64, and the cast to UInt32 is then exact. The snippet exits 1 when the
#: call returns 0 -- its failure return -- rather than holding a process open
#: that is not keeping anything awake, which is what the first version did.
_WINDOWS_STAY_AWAKE = (
    "Add-Type -Name Power -Namespace Win32 -MemberDefinition '"
    "[DllImport(\"kernel32.dll\")] public static extern uint "
    "SetThreadExecutionState(uint esFlags);';"
    "if ([Win32.Power]::SetThreadExecutionState([uint32]2147483649) -eq 0) "
    "{ exit 1 };"
    "while ($true) { Start-Sleep -Seconds 60 }"
)


def isWSL():
    """Whether this Linux is WSL, where the host that sleeps is the Windows one.

    Worth its own check rather than folding into the platform test, because the
    failure it prevents is silent and costs a whole night: under WSL
    `sys.platform` is `linux`, so the Windows branch below never runs, and
    nothing in Linux's own power management has any say over whether the machine
    the VM sits on goes to sleep at midnight.
    """
    if not sys.platform.startswith('linux'):
        return False
    try:
        with open("/proc/version") as handle:
            return "microsoft" in handle.read().lower()
    except OSError:
        return False


def keepAwake():
    if sys.platform == 'win32':
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)
        return None
    if sys.platform == 'darwin':
        try:
            return subprocess.Popen(["caffeinate", "-ims", "-w", str(os.getpid())])
        except OSError:
            pass
    elif isWSL():
        try:
            handle = subprocess.Popen(
                ["powershell.exe", "-NoProfile", "-NonInteractive",
                 "-Command", _WINDOWS_STAY_AWAKE],
                stdin = subprocess.DEVNULL, stdout = subprocess.DEVNULL,
                stderr = subprocess.DEVNULL)
        except OSError:
            print("warning: running under WSL and powershell.exe could not be "
                  "started, so nothing is stopping the Windows host sleeping")
            return None
        # The holder exits 1 the moment the call fails, so a short wait is
        # enough to tell a request that took hold from one that did not.  Worth
        # the two seconds: the failure mode this replaces was a process sitting
        # there all night having kept nothing awake.
        try:
            handle.wait(timeout = 2.0)
        except subprocess.TimeoutExpired:
            print("keeping the Windows host awake through WSL interop "
                  "(PID {0}); the lid can still suspend it, so leave it "
                  "open".format(handle.pid))
            return handle
        print("warning: the Windows keep-awake call failed (exit {0}); nothing "
              "is stopping the host sleeping, so check its power settings "
              "before leaving".format(handle.returncode))
        return None
    print("warning: nothing is stopping this machine sleeping; check its power settings")
    return None


def _adhocName(command, taken):
    """A log-file name for a `--run` job: short, readable and unique.

    Built from the words before the first flag, with the `.py` dropped, so
    `batch.py cores 16 --jobs 8` becomes `batch-cores-16` -- that is what a
    person scanning `logs/` is looking for, and keeping every bare word means
    no command has to be special-cased to keep its name legible.  A collision
    gets a number rather than overwriting, since two ad-hoc jobs in one night
    that differ only in their flags are exactly the case this is for.
    """
    words = []
    for word in shlex.split(command):
        if word.startswith("-"):
            break
        words.append(os.path.splitext(os.path.basename(word))[0])
    stem = "-".join(words) or "run"
    name = stem
    index = 2
    while name in taken:
        name = "{0}-{1}".format(stem, index)
        index += 1
    return name


def main(argv = None):
    parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hours", type = float, default = 9)
    parser.add_argument("--only", nargs = "+", default = None)
    parser.add_argument("--dry-run", action = "store_true", dest = "dryRun")
    parser.add_argument("--run", action = "append", default = [], dest = "adhoc",
                        metavar = "COMMAND",
                        help = "a command to run as a job tonight, as one "
                               "quoted string, e.g. --run \"batch.py cores "
                               "16 --max-word 3 --jobs 8\"; repeatable, and "
                               "--budget-hours is added for you")
    args = parser.parse_args(argv)

    table = jobs(args.hours)
    chosen = list(args.only or ([] if args.adhoc else list(table)))
    for command in args.adhoc:
        name = _adhocName(command, table)
        table[name] = shlex.split(command) + ["--budget-hours", str(args.hours)]
        chosen.append(name)
    if not chosen:
        print("nothing to run", file = sys.stderr)
        return 1
    python = sys.executable
    stamp = time.strftime("%Y%m%d-%H%M")
    logs = os.path.join(HERE, "logs")
    os.makedirs(logs, exist_ok = True)
    deadline = time.time() + args.hours * 3600 + 20 * 60   # grace for in-flight work

    commands = {name: [python, "-u"] + table[name] for name in chosen}
    for name, command in commands.items():
        print("[{0}] {1}".format(name, " ".join(command)))
    if args.dryRun:
        return 0

    awake = keepAwake()
    env = dict(os.environ, PYTHONUNBUFFERED = "1")
    running, attempts, logPaths = {}, {}, {}

    def start(name):
        attempts[name] = attempts.get(name, 0) + 1
        logPaths[name] = os.path.join(logs, "{0}-{1}.log".format(name, stamp))
        handle = open(logPaths[name], "a")
        handle.write("=== {0} attempt {1} started {2} ===\n".format(
            name, attempts[name], time.strftime("%Y-%m-%d %H:%M:%S")))
        handle.flush()
        running[name] = (subprocess.Popen(commands[name], cwd = HERE, env = env,
                                          stdout = handle, stderr = subprocess.STDOUT), handle)
        print("[{0}] started as PID {1}, log {2}".format(name, running[name][0].pid, logPaths[name]),
              flush = True)

    print("overnight run of {0} h started {1}; Ctrl-C to stop".format(
        args.hours, time.strftime("%Y-%m-%d %H:%M")), flush = True)
    for name in commands:
        start(name)
    try:
        while running:
            time.sleep(15)
            for name in list(running):
                process, handle = running[name]
                code = process.poll()
                if code is None:
                    if time.time() > deadline:
                        process.terminate()
                        handle.write("=== terminated at the deadline ===\n")
                    continue
                handle.write("=== {0} exited {1} at {2} ===\n".format(
                    name, code, time.strftime("%Y-%m-%d %H:%M:%S")))
                handle.close()
                del running[name]
                if code in (0, 1, 2):
                    print("[{0}] {1} (exit {2})".format(
                        name, "stopped on its budget" if code == 2 else "finished", code), flush = True)
                elif attempts[name] < 20 and time.time() < deadline - 20 * 60:
                    print("[{0}] died with {1}; resuming".format(name, code), flush = True)
                    start(name)
                else:
                    print("[{0}] died with {1}; giving up".format(name, code), flush = True)
    except KeyboardInterrupt:
        for process, handle in running.values():
            process.terminate()
        print("stopped; rerun to resume")
    finally:
        if awake is not None:
            awake.terminate()
    print("done at {0}. Logs:".format(time.strftime("%Y-%m-%d %H:%M")))
    for name, path in logPaths.items():
        print("  ", path)
    print("\nIn the morning, to read back what each job established:")
    for name in commands:
        summary, aimedAt = MORNING.get(name, (None, None))
        if summary:
            print("   {0:<12} {1}".format(name, summary))
            print("   {0:<12}   against {1}".format("", aimedAt))
        else:
            # An ad-hoc `--run` job has no entry, and the honest thing to
            # print is its own command with `--summary` on the end, which
            # is what every task here answers to.
            words = [word for word in table[name]
                     if word not in ("--budget-hours", str(args.hours))]
            print("   {0:<12} python {1} --summary".format(
                name, " ".join(words)))
    print("\nRecord the outcomes in research/EXPERIMENTS.md whichever way they "
          "went -- a run that\nfound nothing is worth recording precisely so it "
          "is not repeated.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
