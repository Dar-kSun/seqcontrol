"""The trust card must report exactly what the committed result files say."""

import json
import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

from seqcontrol.cli import app
from seqcontrol.trustcard import render_markdown, write_card

RESULTS = Path(__file__).resolve().parent.parent / "results"
needs_results = pytest.mark.skipif(
    not (RESULTS / "permutation_trna-swap.json").exists(), reason="results not committed yet"
)


@needs_results
def test_card_numbers_match_result_files():
    card = render_markdown(RESULTS)
    base = json.loads((RESULTS / "baseline_fp8-delayed.json").read_text())
    swap = json.loads((RESULTS / "permutation_trna-swap.json").read_text())
    native = base["regions"]["Mt_tRNA"]["auroc"]
    assert f"{native['value']:.3f} [{native['ci95'][0]:.3f}, {native['ci95'][1]:.3f}]" in card
    cdi = swap["cdi"]
    assert f"{cdi['value']:.2f} [{cdi['ci95'][0]:.2f}, {cdi['ci95'][1]:.2f}]" in card
    assert base["model"] in card
    assert "clinical" in card  # the no-clinical-claims caveat is always present


@needs_results
def test_card_without_figure(tmp_path):
    md = write_card(RESULTS, tmp_path, figure=False)
    text = md.read_text(encoding="utf-8")
    assert "![" not in text
    assert "## Read with care" in text


def test_card_refuses_to_invent_missing_results(tmp_path):
    with pytest.raises(FileNotFoundError):
        render_markdown(tmp_path)


@needs_results
def test_card_skips_controls_that_were_not_run(tmp_path):
    for name in ["baseline_fp8-delayed.json", "baseline_fp8-delayed.csv",
                 "baseline_fp8-current.csv", "permutation_trna-swap.json"]:  # fmt: skip
        shutil.copy(RESULTS / name, tmp_path / name)
    card = render_markdown(tmp_path)
    assert "tRNA swap" in card
    assert "Window rotation" not in card and "Flank shuffle" not in card


def test_cli_rejects_unknown_step():
    result = CliRunner().invoke(app, ["run", "nonsense"])
    assert result.exit_code != 0
