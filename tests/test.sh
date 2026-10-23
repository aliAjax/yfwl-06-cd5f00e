#!/bin/bash
set -euo pipefail

LOG_DIR=/logs/verifier
REWARD_FILE="$LOG_DIR/reward.txt"
STATUS_FILE="$LOG_DIR/verification-status.txt"
TEST_LOG="$LOG_DIR/test-stdout.txt"
REPORT="$LOG_DIR/pytest-results.xml"
completed=0

mkdir -p "$LOG_DIR"

write_reward() {
  local value="$1"
  chmod 0644 "$REWARD_FILE" 2>/dev/null || true
  printf '%s\n' "$value" > "$REWARD_FILE"
  chmod 0444 "$REWARD_FILE"
}

write_status() {
  chmod 0644 "$STATUS_FILE" 2>/dev/null || true
  printf '%s\n' "$1" > "$STATUS_FILE"
  chmod 0444 "$STATUS_FILE"
}

on_exit() {
  if [ "$completed" -ne 1 ]; then
    write_status "INCOMPLETE: verifier exited before pytest completion"
    write_reward 0.000
  fi
}

on_signal() {
  write_status "INCOMPLETE: verifier interrupted before pytest completion"
  write_reward 0.000
  completed=1
  exit 0
}

trap on_exit EXIT
trap on_signal INT TERM HUP

write_reward 0.000
write_status "INCOMPLETE: pytest has not completed"
rm -f "$REPORT" "$TEST_LOG" "$LOG_DIR/scoring-report.json"

python3 -m pip install --no-cache-dir --disable-pip-version-check pytest==8.4.1 >/dev/null 2>&1

set +e
timeout --signal=TERM --kill-after=10s 145s \
  python3 -u -m pytest -p no:cacheprovider -rA \
  --junitxml="$REPORT" /tests/test_outputs.py > "$TEST_LOG" 2>&1
RC=$?
set -e
cat "$TEST_LOG" 2>/dev/null || true

if [ "$RC" -ne 0 ] && [ "$RC" -ne 1 ]; then
  write_status "INCOMPLETE: pytest did not complete; runner exited with code $RC"
  write_reward 0.000
elif [ ! -s "$REPORT" ]; then
  write_status "INCOMPLETE: pytest completion report missing"
  write_reward 0.000
else
  if python3 /tests/score.py \
    --config /tests/scoring.json \
    --junit "$REPORT" \
    --report-out "$LOG_DIR/scoring-report.json" \
    --reward-out "$REWARD_FILE"; then
    REWARD=$(cat "$REWARD_FILE")
    case "$REWARD" in
      0.[0-9][0-9][0-9]|1.000) ;;
      *)
        write_status "INCOMPLETE: weighted scorer produced invalid reward"
        write_reward 0.000
        completed=1
        exit 0
        ;;
    esac
    chmod 0444 "$REWARD_FILE"
    write_status "COMPLETED: weighted scoring finished with reward=$REWARD"
  else
    write_status "INCOMPLETE: weighted scorer rejected the report"
    write_reward 0.000
  fi
fi

completed=1
exit 0

