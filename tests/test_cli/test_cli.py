import warnings

import click
import pytest
from click.testing import CliRunner

from dstools.cli import (
    askap_preprocess,
    create_model,
    derotate_feeds,
    extract_ds,
    insert_model,
    plot_ds,
    selfcal,
    subtract_model,
)

warnings.filterwarnings("ignore", category=DeprecationWarning, append=True)

commands = [
    ("askap-preprocess", askap_preprocess.main),
    ("derotate-feeds", derotate_feeds.main),
    ("create-model", create_model.main),
    ("insert-model", insert_model.main),
    ("selfcal", selfcal.main),
    ("subtract-model", subtract_model.main),
    ("extract-ds", extract_ds.main),
    ("plot-ds", plot_ds.main),
]


def assert_logged_error_result(result):
    output = click.unstyle(result.output)

    assert result.exit_code == 1
    assert "ERROR" in output


@pytest.mark.parametrize("script, command", commands)
def test_script_help(script, command):
    runner = CliRunner()
    result = runner.invoke(command, ["--help"])

    assert result.exit_code == 0
    assert "Usage:" in result.output
    assert "--help" in result.output


@pytest.mark.parametrize(
    "module, args",
    [
        (askap_preprocess, ["dummy.ms"]),
        (derotate_feeds, ["dummy.ms"]),
        (create_model, ["dummy.ms"]),
        (insert_model, ["missing-model", "dummy.ms"]),
        (selfcal, ["dummy.ms"]),
        (subtract_model, ["dummy.ms"]),
        (extract_ds, ["dummy.ms", "out.ds"]),
    ],
)
def test_casa_cli_unsupported_dependency_failure(monkeypatch, module, args):
    runner = CliRunner()
    monkeypatch.setattr(module, "HAS_CASA_SUPPORT", False)

    result = runner.invoke(module.main, args)

    assert_logged_error_result(result)
