"""Does batch size change Evo 2 scores? Compares batch sizes 1, 8 and 16.

Usage (inside the WSL evo2 environment):
    python scripts/check_batch_size.py

Scores the ref and alt 501 bp windows of 8 single-base variants at random positions
in tests/fixtures/chrM_1_1000.fa, and prints the largest per-sequence score change
between batch sizes next to the typical variant effect. Calls Evo 2 directly, not
through Evo2Adapter, because the adapter fixes batch size at 1.
Results are recorded in docs/model-choice.md.
"""

from pathlib import Path

import numpy as np
from evo2 import Evo2

from seqcontrol.variants import Variant, variant_windows

FIXTURE = Path(__file__).resolve().parent.parent / "tests" / "fixtures" / "chrM_1_1000.fa"
N_VARIANTS = 8


def main() -> None:
    lines = FIXTURE.read_text().splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    rng = np.random.default_rng(0)
    refs, alts = [], []
    for pos in rng.choice(np.arange(260, 740), size=N_VARIANTS, replace=False):
        ref_base = seq[pos - 1]
        alt_base = next(b for b in "ACGT" if b != ref_base)
        variant = Variant("MT", int(pos), ref_base, alt_base)
        ref, alt = variant_windows(seq, variant, 501, circular=True)
        refs.append(ref)
        alts.append(alt)
    seqs = refs + alts

    model = Evo2("evo2_1b_base")
    s1 = np.array(model.score_sequences(seqs, batch_size=1))
    s1_again = np.array(model.score_sequences(seqs, batch_size=1))
    s8 = np.array(model.score_sequences(seqs, batch_size=8))
    s16 = np.array(model.score_sequences(seqs, batch_size=16))

    def effects(s):
        return s[N_VARIANTS:] - s[:N_VARIANTS]

    print(f"batch 1 vs batch 1 again  max|diff| {np.abs(s1 - s1_again).max():.6f}")
    print(f"batch 1 vs batch 8        max|diff| {np.abs(s1 - s8).max():.6f}")
    print(f"batch 1 vs batch 16       max|diff| {np.abs(s1 - s16).max():.6f}")
    print(f"median |variant effect| at batch 1  {np.median(np.abs(effects(s1))):.6f}")
    flips = int(np.sum(np.sign(effects(s1)) != np.sign(effects(s16))))
    print(f"variant effects that change sign, batch 1 vs 16: {flips} of {N_VARIANTS}")


if __name__ == "__main__":
    main()
