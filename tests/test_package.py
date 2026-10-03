from typer.testing import CliRunner

import seqcontrol
from seqcontrol.cli import app


def test_version_is_set():
    assert seqcontrol.__version__


def test_cli_version():
    result = CliRunner().invoke(app, ["version"])
    assert result.exit_code == 0
    assert seqcontrol.__version__ in result.stdout
