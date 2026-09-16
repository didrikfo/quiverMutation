#!/usr/bin/env bash
#
# Run the two long jobs overnight, side by side, and survive the night.
#
#     ./overnight.sh                  # both jobs, 9 hours, cores split automatically
#     ./overnight.sh 7                # both jobs, 7 hours
#     ./overnight.sh 9 classify       # just the classification
#     ./overnight.sh 9 discover       # just the rule search
#
# Job A, on one core: finish the n = 10 classification.  Every step checkpoints
# now, so this is resumable at class granularity and is restarted automatically
# if it dies.
#
# Job B, on the rest: search for a mutation rule that breaks the maximally
# overlapping pair, which blocks a fifth of the rows job A has to search for
# (research H-009).  Four mutations first, then five.  Also checkpointed, also
# restarted.
#
# Both write their progress to disk continuously, so whatever the night gets
# through is kept.  Read the logs in the morning; nothing here needs watching.
#
# Two things that have killed runs before (research E-008), both handled here:
#
#   * the machine sleeping.  This script wraps itself in caffeinate on macOS and
#     systemd-inhibit on Linux, so the laptop stays awake with the lid open.
#     Closing the lid still suspends most laptops -- leave it open.
#   * an over-broad `pkill -f` typed to clean up, which matches the shell that
#     issued it.  To stop these jobs, kill the PIDs printed below, or Ctrl-C
#     this script, which stops its children with it.

set -u

HOURS="${1:-9}"
WHICH="${2:-both}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-$HERE/.venv/bin/python}"
LENGTH="${LENGTH:-10}"
TARGET_LENGTH="${TARGET_LENGTH:-9}"
LOGS="$HERE/logs"

if [ ! -x "$PYTHON" ]; then
    PYTHON="$(command -v python3 || command -v python)"
    echo "note: no .venv found, using $PYTHON"
fi

# Keep the machine awake.  Re-exec under whichever inhibitor exists, once.
if [ -z "${OVERNIGHT_AWAKE:-}" ]; then
    export OVERNIGHT_AWAKE=1
    # Try each inhibitor on a trivial command first.  systemd-inhibit is present
    # but unusable in a container with no session bus, and exec'ing into a
    # failing inhibitor would end the run before it started.
    if command -v caffeinate >/dev/null 2>&1 && caffeinate -i true >/dev/null 2>&1; then
        exec caffeinate -ims "$0" "$@"
    elif command -v systemd-inhibit >/dev/null 2>&1 \
         && systemd-inhibit --what=idle --why=probe true >/dev/null 2>&1; then
        exec systemd-inhibit --what=idle:sleep --why="quiver mutation overnight run" \
             "$0" "$@"
    else
        echo "warning: could not stop this machine sleeping -- no working caffeinate"
        echo "         or systemd-inhibit. If it sleeps, the jobs stop. Check its"
        echo "         power settings, and leave the lid open."
    fi
fi

mkdir -p "$LOGS"
STAMP="$(date +%Y%m%d-%H%M)"
CORES="$( (command -v nproc >/dev/null 2>&1 && nproc) \
        || sysctl -n hw.ncpu 2>/dev/null || echo 4 )"

# Job A gets one core; job B gets the rest, less one left for the machine.
if [ "$WHICH" = "discover" ]; then
    DISCOVER_JOBS=$(( CORES > 2 ? CORES - 1 : 1 ))
else
    DISCOVER_JOBS=$(( CORES > 3 ? CORES - 2 : 1 ))
fi

# Restart a job if it dies, until it succeeds or the budget is spent.  Exit 2 is
# classify.py saying "out of time, resume me"; 0 and 1 are both finished answers.
persist() {
    local name="$1" log="$2"; shift 2
    local attempt=1
    while :; do
        echo "=== $name attempt $attempt started $(date) ===" >> "$log"
        "$@" >> "$log" 2>&1
        local code=$?
        echo "=== $name exited $code at $(date) ===" >> "$log"
        case $code in
            0|1) echo "[$name] finished (exit $code); see $log"; return $code ;;
            2)   echo "[$name] stopped on its budget; see $log";  return 0 ;;
        esac
        if [ "$attempt" -ge 20 ]; then
            echo "[$name] died 20 times; giving up. See $log"; return $code
        fi
        echo "[$name] died with $code, resuming in 10s (attempt $attempt)"
        attempt=$(( attempt + 1 ))
        sleep 10
    done
}

echo "Quiver mutation, overnight run of $HOURS hour(s), started $(date)"
echo "  $CORES cores; logs in $LOGS"
echo

PIDS=()

if [ "$WHICH" = "both" ] || [ "$WHICH" = "classify" ]; then
    LOG_A="$LOGS/classify-$LENGTH-$STAMP.log"
    echo "[A] classifying n = $LENGTH  ->  $LOG_A"
    persist "classify n=$LENGTH" "$LOG_A" \
        "$PYTHON" "$HERE/classify.py" "$LENGTH" --resume --budget-hours "$HOURS" &
    PIDS+=($!)
fi

if [ "$WHICH" = "both" ] || [ "$WHICH" = "discover" ]; then
    LOG_B="$LOGS/discover-$STAMP.log"
    echo "[B] searching for a rule on $DISCOVER_JOBS core(s)  ->  $LOG_B"
    (
        # Four mutations first: cheaper, and H-009 may already be settled there.
        # Then five, which no search has ever reached.  The two passes share one
        # budget -- each is given what is left of it, not a fresh copy, or the
        # night would run to twice the hours asked for.
        ENDS_AT=$(( $(date +%s) + $(awk -v h="$HOURS" 'BEGIN{printf "%d", h*3600}') ))
        for steps in 4 5; do
            LEFT=$(( ENDS_AT - $(date +%s) ))
            if [ "$LEFT" -le 120 ]; then
                echo "[discover] out of budget before steps=$steps; resume with:"
                echo "  $PYTHON $HERE/discover.py $TARGET_LENGTH --steps $steps --resume"
                break
            fi
            persist "discover steps=$steps" "$LOG_B" \
                "$PYTHON" "$HERE/discover.py" "$TARGET_LENGTH" \
                --steps "$steps" --jobs "$DISCOVER_JOBS" --resume \
                --budget-hours "$(awk -v s="$LEFT" 'BEGIN{printf "%.4f", s/3600}')"
        done
    ) &
    PIDS+=($!)
fi

echo
echo "Running as PIDs: ${PIDS[*]}"
echo "Stop them with:  kill ${PIDS[*]}    (not pkill -f -- see the note above)"
echo "Watch with:      tail -f $LOGS/*-$STAMP.log"
echo

wait
echo
echo "All jobs done at $(date). In the morning:"
[ -n "${LOG_A:-}" ] && echo "  * $LOG_A"
[ -n "${LOG_A:-}" ] && echo "      tail it; rerun the same command with --resume if unfinished"
[ -n "${LOG_B:-}" ] && echo "  * $LOG_B"
[ -n "${LOG_B:-}" ] && echo "      the verified rules, if any, are at the end of each pass"
