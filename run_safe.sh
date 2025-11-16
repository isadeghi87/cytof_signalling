#!/usr/bin/env bash
set -euo pipefail

# Safe pipeline runner: redirect heavy output to logs to keep terminal responsive
LOGDIR="logs"
mkdir -p "$LOGDIR"

echo "Starting safe pipeline run: $(date)" > "$LOGDIR/run_safe.log"

if [ ! -x ./run_all.sh ]; then
  echo "Warning: ./run_all.sh not executable, attempting to execute with bash" | tee -a "$LOGDIR/run_safe.log"
  bash ./run_all.sh > "$LOGDIR/run_all.log" 2>&1 || { echo "run_all.sh failed; see $LOGDIR/run_all.log" | tee -a "$LOGDIR/run_safe.log"; exit 1; }
else
  ./run_all.sh > "$LOGDIR/run_all.log" 2>&1 || { echo "run_all.sh failed; see $LOGDIR/run_all.log" | tee -a "$LOGDIR/run_safe.log"; exit 1; }
fi

echo "run_all.sh completed: $(date)" >> "$LOGDIR/run_safe.log"
echo "Logs: $LOGDIR/run_all.log (details), $LOGDIR/run_safe.log (summary)"

# show last 20 lines of log so user sees progress without flooding terminal
tail -n 20 "$LOGDIR/run_all.log" || true

exit 0
