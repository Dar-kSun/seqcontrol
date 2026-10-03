#!/usr/bin/env bash
# Unattended GPU queue: robustness re-runs and the wide-window flank sweep.
# Each step logs to logs/, is retried once on failure, and the queue moves on
# either way. Progress: logs/overnight_status.txt. Takes several hours.
#
# Usage, inside the WSL evo2 environment (see docs/model-choice.md):
#     bash scripts/run_overnight.sh
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source "${CONDA_ROOT:-/opt/miniforge}/bin/activate" evo2
mkdir -p logs
STATUS=logs/overnight_status.txt
echo "queue started $(date -Is) on $(hostname), commit $(git rev-parse --short HEAD 2>/dev/null)" > "$STATUS"

step() {  # step <name> <command...>
    local name=$1; shift
    for attempt in 1 2; do
        echo "$(date -Is) START $name (attempt $attempt): $*" >> "$STATUS"
        python "$@" > "logs/overnight_${name}.log" 2>&1
        local code=$?
        echo "$(date -Is) END   $name exit=$code" >> "$STATUS"
        [ $code -eq 0 ] && return 0
        sleep 60  # let the GPU driver release memory before retrying
    done
    return 1
}

# 0. Wait for a flank sweep that is already running, then make sure it finished.
while pgrep -f "python scripts/03_flank_sweep.py" > /dev/null; do sleep 60; done
[ -f results/flank_sweep.json ] || step flank_sweep scripts/03_flank_sweep.py

# 1. Short robustness re-runs under the other FP8 recipe (~50 min).
step swap_fp8current     scripts/02_permutation.py --control trna-swap --precision fp8-current
step rotation_fp8current scripts/02_permutation.py --control window-rotation --precision fp8-current
step sweep_fp8current    scripts/03_flank_sweep.py --precision fp8-current

# 2. Wide windows (4,097 bp) so the sweep reaches radii up to 1,900 bp (hours).
step baseline_w4097 scripts/01_baseline.py --window 4097
step sweep_w4097    scripts/03_flank_sweep.py --window 4097 --radii 0,50,100,250,500,1000,1500,1900

echo "$(date -Is) QUEUE DONE" >> "$STATUS"
