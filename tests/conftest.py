from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def chrm_1kb() -> str:
    lines = (FIXTURES / "chrM_1_1000.fa").read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))
