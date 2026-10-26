#!/bin/bash
# Study 2 GEARS queue (NOT started automatically). Usage: scripts/study2_run_queue.sh "seed1 seed2 seed3" "1 2 3"   [from the repo root]
# For every (split, training-seed) pair: skip if finished (predictions_test.npz exists), otherwise run/resume the checkpointed job. Runs under `caffeinate -is` and refuses to start on battery.
set -u
SPLITS=${1:-"seed1 seed2 seed3"}; SEEDS=${2:-"1 2 3"}
PY=/Users/zeynep/miniforge3/envs/hic-perturb-gears/bin/python
if ! pmset -g batt | grep -q "AC Power"; then echo "On battery: connect AC power first (caffeinate -s cannot prevent sleep on battery)."; exit 1; fi
for sp in $SPLITS; do for sd in $SEEDS; do
  out="results/study2/gears_k562_gears_${sp}_s${sd}"
  if [ -f "$out/predictions_test.npz" ]; then echo "skip $sp s$sd (done)"; continue; fi
  echo "[$(date)] start/resume $sp s$sd"
  caffeinate -is $PY -u scripts/study2_gears_run.py --label "gears_${sp}" --seed "$sd" --device mps > "results/study2_stdout_${sp}_s${sd}.log" 2>&1 \
    || { echo "FAILED $sp s$sd (checkpoint kept; rerun this script to resume)"; exit 2; }
done; done
echo "[$(date)] queue finished"
