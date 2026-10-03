#!/usr/bin/env bash
# Build the Linux environment that runs Evo 2 (see docs/model-choice.md).
# Tested on Ubuntu 24.04 under WSL2, RTX 4060 Laptop GPU (sm_89), driver 595.97.
# Usage, from inside WSL:  bash scripts/setup_evo2_wsl.sh
set -euo pipefail

PREFIX=${CONDA_PREFIX_DIR:-$HOME/miniforge}
ENV=evo2

if [ ! -x "$PREFIX/bin/conda" ]; then
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$PREFIX"
fi
source "$PREFIX/bin/activate"

conda create -y -n "$ENV" python=3.12
conda activate "$ENV"

# Transformer Engine (FP8) is required by evo2_1b_base; it brings PyTorch 2.6 (CUDA 12.6).
conda install -y -c nvidia -c conda-forge cuda-nvcc cuda-cudart-dev transformer-engine-torch=2.3.0
pip install ninja packaging wheel
MAX_JOBS=4 pip install flash-attn==2.8.0.post2 --no-build-isolation
pip install evo2==0.6.0
pip install -e "$(dirname "$0")/.."  # this repo, so scripts can import seqcontrol

python -c "import torch, transformer_engine as te, flash_attn, evo2; \
print('torch', torch.__version__, 'te', te.__version__, 'flash_attn', flash_attn.__version__, \
'cuda_ok', torch.cuda.is_available())"
