#!/usr/bin/env bash
#
# Run the long jobs overnight, and survive the night.
#
#     ./overnight.sh                  # both jobs, 9 hours, cores split automatically
#     ./overnight.sh 7                # both jobs, 7 hours
#     ./overnight.sh 9 classify       # just the classification
#     ./overnight.sh 9 probe          # just the deep probe
#
# Job A: finish the n = 10 classification. Every step of it checkpoints, so this
# is resumable at class granularity and is restarted automatically if it dies.
#
# Job B: the deep interior probe research H-010 asks for -- whether any
# sequence, at any depth, lowers the overlap of an isolated pair. Six mutations
# took 43 minutes; seven is the next unknown. **This one is not resumable**: it
# holds its search in memory and prints at the end, so a kill loses it. It is
# run under a hard timeout for that reason, and if the night is not long enough
# the answer is to raise the hours, not to expect a partial result.
#
# Both jobs print unbuffered, so the logs are live. Whatever job A gets through
# is kept; job B either answers or does not.
#
# Two things that have killed runs before (research E-008), both handled here:
#
#   * the machine sleeping. This wraps itself in caffeinate on macOS and
#     systemd-inhibit on Linux. Closing the lid still suspends most laptops --
#     leave it open.
#   * an over-broad `pkill -f` typed to clean up, which matches the shell that
#     issued it. To stop these jobs, kill the PIDs printed below, or Ctrl-C this
#     script, which stops its children with it.

set -u

HOURS="${1:-9}"
WHICH="${2:-both}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-$HERE/.venv/bin/python}"
LENGTH="${LENGTH:-10}"
PATTERN="${PATTERN:-1:3,2:3}"
PROBE_STEPS="${PROBE_STEPS:-7}"
PROBE_CLEARANCE="${PROBE_CLEARANCE:-9}"
LOGS="$HERE/logs"

if [ ! -x "$PYTHON" ]; then
    PYTHON="$(command -v python3 || command -v python)"
    echo "note: no .venv found, using $PYTHON"
fi

# Keep the machine awake. Try each inhibitor on a trivial command first:
# systemd-inhibit is present but unusable without a session bus, and exec'ing
# into a failing inhibitor would end the run before it started.
if [ -z "${OVERNIGHT_AWAKE:-}" ]; then
    export OVERNIGHT_AWAKE=1
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
SECONDS_TOTAL="$(awk -v h="$HOURS" 'BEGIN{printf "%d", h*3600}')"

# Restart a job if it dies, until it succeeds or the budget is spent. Exit 2 is
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
echo "  logs in $LOGS (python runs unbuffered, so tail -f them)"
echo

PIDS=()

if [ "$WHICH" = "both" ] || [ "$WHICH" = "classify" ]; then
    LOG_A="$LOGS/classify-$LENGTH-$STAMP.log"
    echo "[A] classifying n = $LENGTH  ->  $LOG_A"
    persist "classify n=$LENGTH" "$LOG_A" \
        "$PYTHON" -u "$HERE/classify.py" "$LENGTH" --resume --budget-hours "$HOURS" &
    PIDS+=($!)
fi

if [ "$WHICH" = "both" ] || [ "$WHICH" = "probe" ]; then
    LOG_B="$LOGS/probe-$STAMP.log"
    echo "[B] probing $PATTERN to $PROBE_STEPS mutations  ->  $LOG_B"
    echo "    (not resumable -- under a hard $HOURS h timeout)"
    (
        echo "=== probe $PATTERN --steps $PROBE_STEPS started $(date) ===" >> "$LOG_B"
        timeout "$SECONDS_TOTAL" \
            "$PYTHON" -u "$HERE/probe.py" "$PATTERN" \
            --steps "$PROBE_STEPS" --clearance "$PROBE_CLEARANCE" >> "$LOG_B" 2>&1
        code=$?
        if [ "$code" -eq 124 ]; then
            echo "=== probe ran out of the night at $(date); no partial result ===" >> "$LOG_B"
            echo "[probe] out of time -- rerun with more hours, or fewer --steps"
        else
            echo "=== probe exited $code at $(date) ===" >> "$LOG_B"
            echo "[probe] finished (exit $code); see $LOG_B"
        fi
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
[ -n "${LOG_B:-}" ] && echo "      'nothing lower' is a real result -- record it against H-010"
