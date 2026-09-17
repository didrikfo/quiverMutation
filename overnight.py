#!/usr/bin/env python
"""Run the long jobs overnight, and survive the night -- on any platform.

    python overnight.py                 the default jobs, 9 hours
    python overnight.py --hours 7
    python overnight.py --only merges10 classify10
    python overnight.py --dry-run       print what would run

`overnight.sh` does this with caffeinate and bash job control, neither of which
exists on Windows. This does the same three things from Python:

* keeps the machine awake for the whole run (`SetThreadExecutionState` on
  Windows, `caffeinate` on macOS if present) -- closing a laptop's lid can still
  suspend it, so leave it open;
* restarts a job that dies, with the job's own resume, until it finishes or the
  budget is spent -- exit 0 and 1 are finished answers, 2 is "out of budget";
* logs every job unbuffered to `logs/<job>-<stamp>.log`.

The default jobs, aimed by research F-032 and H-013 (2026-09-17):

    merges10    python merges.py 10 --depths 5 6 7 8     7 processes
    merges11    python merges.py 11 --depths 4 5 6       7 processes
    classify10  python classify.py 10 --resume           1 process, an independent
                                                          route to the n = 10 count

Stop it with Ctrl-C, which stops the children. Every job resumes from its own
checkpoint when rerun, so stopping loses only the work in flight.
"""

import argparse
import os
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
    }


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
    print("warning: nothing is stopping this machine sleeping; check its power settings")
    return None


def main(argv = None):
    parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hours", type = float, default = 9)
    parser.add_argument("--only", nargs = "+", default = None)
    parser.add_argument("--dry-run", action = "store_true", dest = "dryRun")
    args = parser.parse_args(argv)

    table = jobs(args.hours)
    chosen = args.only or list(table)
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
    print("In the morning: `python merges.py 10 --summary` and `python merges.py 11 --summary`, "
          "and record the outcome against H-013 in research/EXPERIMENTS.md.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
