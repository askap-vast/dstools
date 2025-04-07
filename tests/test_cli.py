import os
import warnings

import pytest
from click.testing import CliRunner

from dstools.cli import subtract_model

warnings.filterwarnings("ignore", category=DeprecationWarning, append=True)

scripts = [
    "askap-preprocess",
    "create-model",
    "insert-model",
    "selfcal",
    "subtract-model",
    "extract-ds",
    "plot-ds",
]


@pytest.mark.parametrize("script", scripts)
def test_script_help(script):
    exit_code = os.system(f"dstools-{script} --help")

    assert exit_code == 0
