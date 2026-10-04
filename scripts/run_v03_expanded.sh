#!/usr/bin/env bash
# GPU queue for plan v0.3, Arm 2: widened tRNA labels (ClinVar benign 1+ star).
# Same pattern as run_overnight.sh. Progress: logs/v03_status.txt. ~2 hours.
set -o pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
source "${CONDA_ROOT:-/opt/miniforge}/bin/activate" evo2
set -u  # after activation: conda's activate script uses unset variables
mkdir -p logs
STATUS=logs/v03_status.txt
echo "queue started $(date -Is), commit $(git rev-parse --short HEAD 2>/dev/null)" > "$STATUS"

step() {  # step <name> <command...>
    local name=$1; shift
    for attempt in 1 2; do
        echo "$(date -Is) START $name (attempt $attempt): $*" >> "$STATUS"
        python "$@" > "logs/v03_${name}.log" 2>&1
        local code=$?
        echo "$(date -Is) END   $name exit=$code" >> "$STATUS"
        [ $code -eq 0 ] && return 0
        sleep 60
    done
    return 1
}

step baseline_expanded scripts/01_baseline.py --labels expanded
step swap_expanded     scripts/02_permutation.py --control trna-swap --labels expanded
step sweep_expanded    scripts/03_flank_sweep.py --labels expanded --radii 0,100,400

echo "$(date -Is) QUEUE DONE" >> "$STATUS"
