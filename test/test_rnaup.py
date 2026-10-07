import shutil
from pathlib import Path

import pytest

from model.score_rnaup import anchor_interaction, parse_rnaup_output, run_rnaup


OUTPUT = "((((((&))))))  121,160 : 1,40  (-12.30 = -20.00 + 5.00 + 2.70)"


def test_parser_and_anchor_mapping() -> None:
    parsed = parse_rnaup_output(OUTPUT)
    assert parsed
    assert parsed["rnaup_dG_total"] == -12.3
    assert parsed["interaction_target_start"] == 121
    ok, overlap, start, end = anchor_interaction(parsed, 0, 120, 160)
    assert ok and overlap == 1.0 and (start, end) == (120, 160)
    off, _, _, _ = anchor_interaction(parsed, 1000, 120, 160)
    assert not off


def test_mock_rnaup_command_and_result(tmp_path: Path) -> None:
    executable = tmp_path / "RNAup"
    executable.write_text(
        "#!/bin/sh\n"
        "printf '>target\\n>query\\n((((((&))))))  121,160 : 1,40  "
        "(-12.30 = -20.00 + 5.00 + 2.70)\\n'\n"
    )
    executable.chmod(0o755)
    utr = "ACGT" * 100
    result = run_rnaup(
        utr, "ACGT" * 10, 120, 160, context=120, rnaup_exe=str(executable)
    )
    assert result["status"] == "eligible"
    assert result["anchor_overlap"] == 1.0


@pytest.mark.integration
def test_real_rnaup_if_installed() -> None:
    if not shutil.which("RNAup"):
        pytest.skip("RNAup is not installed")
    utr = "ACGT" * 80
    result = run_rnaup(utr, "ACGT" * 10, 120, 160)
    assert result["status"] in {"eligible", "failed"}
    if result["status"] == "eligible":
        assert isinstance(result["rnaup_dG_total"], float)
