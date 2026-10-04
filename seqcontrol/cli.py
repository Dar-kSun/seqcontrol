"""Command-line entry point: seqcontrol card / run / smoke."""

import subprocess
import sys

import typer

from seqcontrol import __version__, config

app = typer.Typer(help="Control experiments for genomic language models.", no_args_is_help=True)

# Each run step is a script in scripts/; GPU steps need the Evo 2 environment.
STEPS = {
    "data": ["scripts/fetch_data.py"],
    "baseline": ["scripts/01_baseline.py"],
    "swap": ["scripts/02_permutation.py", "--control", "trna-swap"],
    "rotation": ["scripts/02_permutation.py", "--control", "window-rotation"],
    "sweep": ["scripts/03_flank_sweep.py"],
}
GPU_STEPS = {"baseline", "swap", "rotation", "sweep"}
CSV_STEPS = {"swap", "rotation", "sweep"}  # can recompute summaries without a GPU


@app.callback()
def main() -> None:
    """Control experiments for genomic language models."""


@app.command()
def version() -> None:
    """Print the installed seqcontrol version."""
    typer.echo(__version__)


@app.command()
def card(no_figure: bool = typer.Option(False, help="Skip the matplotlib figure.")) -> None:
    """Write the tRNA trust card to results/ from the saved results."""
    from seqcontrol.trustcard import write_card

    out = write_card(config.ROOT / "results", config.ROOT / "results", figure=not no_figure)
    typer.echo(f"wrote {out.relative_to(config.ROOT)}")


@app.command()
def run(
    step: str = typer.Argument(..., help=f"one of: {', '.join(STEPS)}, all"),
    from_csv: bool = typer.Option(False, help="Recompute summaries from saved CSVs (no GPU)."),
) -> None:
    """Run one pipeline step (or all of them) via its script in scripts/."""
    names = list(STEPS) if step == "all" else [step]
    for name in names:
        if name not in STEPS:
            raise typer.BadParameter(f"unknown step {name!r}; choose from {', '.join(STEPS)}")
        if from_csv and name not in CSV_STEPS:
            typer.echo(f"skip {name}: no --from-csv mode")
            continue
        cmd = [sys.executable, *STEPS[name], *(["--from-csv"] if from_csv else [])]
        if name in GPU_STEPS and not from_csv:
            typer.echo(f"{name}: needs the Evo 2 GPU environment (docs/model-choice.md)")
        typer.echo("$ " + " ".join(cmd[1:]))
        code = subprocess.call(cmd, cwd=config.ROOT)
        if code != 0:
            raise typer.Exit(code)


@app.command()
def smoke(model: str = "evo2_1b_base") -> None:
    """Load the model, score 1 kb of real mtDNA, report time and memory (needs a GPU)."""
    raise typer.Exit(
        subprocess.call(
            [sys.executable, "scripts/00_smoke_test.py", "--model", model], cwd=config.ROOT
        )  # fmt: skip
    )


if __name__ == "__main__":
    app()
