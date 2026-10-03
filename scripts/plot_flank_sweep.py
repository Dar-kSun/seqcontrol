"""Draw the flank-sweep figure from results/flank_sweep.json (no GPU).

Usage:
    pip install -e ".[plot]"
    python scripts/plot_flank_sweep.py   # writes results/flank_sweep.png
"""

from seqcontrol import config
from seqcontrol.plots import flank_sweep_figure


def main() -> None:
    out = config.ROOT / "results" / "flank_sweep.png"
    flank_sweep_figure(config.ROOT / "results", out)
    print(f"wrote {out.relative_to(config.ROOT)}")


if __name__ == "__main__":
    main()
