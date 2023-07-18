import filecmp
import importlib.util
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).absolute().parents[1] / "examples"


@pytest.mark.parametrize(
    "test_file,expected_output",
    (
        ("reaxpro_user_case.py", "reaxpro_user_case"),
        ("reaxpro_user_case_full_semantic.py", "reaxpro_user_case"),
    ),
)
def test_reaxpro_user_case_base(test_file, expected_output):
    test_file = EXAMPLES / test_file

    spec = importlib.util.spec_from_file_location("module.name", test_file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.myJob.check()  # True if normal termination

    expected = Path(__file__).parent / test_file.stem
    assert expected.exists()

    result = Path(module.myJob.path)
    assert result.exists()

    comparison = filecmp.dircmp(result, expected)
    comparison.report_full_closure()

    input_files = {
        "energetics_input.dat",
        "lattice_input.dat",
        "mechanism_input.dat",
        "simulation_input.dat",
        "slurm.run",
    }

    # Seems like zacros is not entirely deterministic despite fixed seed
    # At least if the input files are in the common files, that's good enough
    assert input_files.issubset(comparison.common_files)
