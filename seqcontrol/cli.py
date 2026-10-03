"""Command-line entry point. Subcommands (run, card, smoke) arrive in later milestones."""

import typer

from seqcontrol import __version__

app = typer.Typer(help="Control experiments for genomic language models.", no_args_is_help=True)


@app.callback()
def main() -> None:
    """Control experiments for genomic language models."""


@app.command()
def version() -> None:
    """Print the installed seqcontrol version."""
    typer.echo(__version__)


if __name__ == "__main__":
    app()
