import io
import logging
import os
import subprocess
import sys
import tempfile
from contextlib import contextmanager

import pytest

from dstools.logger import (
    filter_pipe_output,
    filter_stdout,
    parse_stdout_stderr,
    redirect_c_output,
    setupLogger,
)


def test_setupLogger_verbose_true(caplog):
    setupLogger(verbose=True)

    logger = logging.getLogger(__name__)
    logger.addHandler(caplog.handler)

    logger.debug("debug message")
    assert "debug message" in caplog.text


def test_setupLogger_verbose_false(caplog):
    setupLogger(verbose=False)

    logger = logging.getLogger(__name__)
    logger.addHandler(caplog.handler)

    logger.debug("debug message")
    assert "debug message" not in caplog.text

    logger.info("info message")
    assert "info message" in caplog.text


def test_setupLogger_file_capture(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("temp")
    setupLogger(verbose=False, filename=tmp_path / "test.log")

    logger = logging.getLogger(__name__)

    logger.info("test message")

    assert (tmp_path / "test.log").exists()


debug_lines = [
    "Warning:  Using post-Aug94 ATCA flux scale for 1934-638",
    "Warning:  Correlations flagged or edge-rejected:",
    "PGPLOT /png: writing new file as",
]


@pytest.mark.parametrize("msg", debug_lines)
def test_parse_stdout_stderr_debug_lines_remove_warning(msg, caplog):
    logger = logging.getLogger(__name__)
    setupLogger(verbose=True)

    logger.addHandler(caplog.handler)

    p = subprocess.Popen(
        f"echo {msg}",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        executable="/bin/bash",
    )

    parse_stdout_stderr(p, logger, print_stdout=False)

    msg = msg.replace("Warning:  ", "")

    assert msg in caplog.text


def test_parse_stdout_stderr_warning_log_strips_warning_text(caplog):
    logger = logging.getLogger(__name__)
    setupLogger(verbose=True)

    logger.addHandler(caplog.handler)

    with caplog.at_level(logging.WARNING, logger=logger.name):
        p = subprocess.Popen(
            "echo This is a legitimate warning! 1>&2",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            executable="/bin/bash",
        )

        parse_stdout_stderr(p, logger, print_stdout=False)

    assert "This is a legitimate warning!" in caplog.text


def test_parse_stdout_stderr_to_stdout(caplog, capfd):
    logger = logging.getLogger(__name__)
    setupLogger(verbose=True)

    logger.addHandler(caplog.handler)

    p = subprocess.Popen(
        "echo This is a message intended for stdout and not logs!",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        executable="/bin/bash",
    )

    parse_stdout_stderr(p, logger, print_stdout=True)

    out, _ = capfd.readouterr()

    assert "This is a message intended for stdout and not logs!" in out
    assert "This is a message intended for stdout and not logs!" not in caplog.text


def run_filter_test(input_lines, substrings, expected_output):
    """Test logic for simulating filter_pipe_output interaction with a pipe / file descriptor."""

    pipe_r, pipe_w = os.pipe()

    # Write test input to the write-end of the pipe
    with os.fdopen(pipe_w, "w") as w:
        for line in input_lines:
            w.write(line + "\n")

    # Capture original stream output via a temp file
    with tempfile.TemporaryFile(mode="w+b") as temp_output:
        orig_fd = temp_output.fileno()

        # Save a duplicate of stdout (so we can restore it)
        saved_stdout = os.dup(1)
        os.dup2(orig_fd, 1)

        try:
            filter_pipe_output(pipe_r, substrings, 1, "stdout")

            temp_output.seek(0)
            output = temp_output.read().decode().splitlines()
        finally:
            os.dup2(saved_stdout, 1)
            os.close(saved_stdout)

    return output


def test_filter_pipe_output_filters_correctly():
    input_lines = [
        "this should be shown",
        "debug: this should be filtered",
        "warning: this too",
        "normal line",
    ]
    substrings = ["debug", "warning"]

    expected_output = ["this should be shown", "normal line"]

    output = run_filter_test(input_lines, substrings, expected_output)
    assert output == expected_output


def test_filter_pipe_output_no_filtering():
    input_lines = ["debug msg", "info msg"]
    substrings = ["warning"]
    expected_output = ["debug msg", "info msg"]

    output = run_filter_test(input_lines, substrings, expected_output)
    assert output == expected_output


def test_filter_pipe_output_all_filtered():
    input_lines = ["debug msg", "warning msg"]
    substrings = ["debug", "warning"]
    expected_output = []

    output = run_filter_test(input_lines, substrings, expected_output)
    assert output == expected_output


def simulate_c_output(stdout_lines, stderr_lines):
    """Helper to simulate C-level output (stdout and stderr)."""

    for line in stdout_lines:
        os.write(1, f"{line}\n".encode())
    for line in stderr_lines:
        os.write(2, f"{line}\n".encode())

    return


def test_redirect_c_output_filtering_stdout_only(capfd):
    stdout_lines = ["debug msg", "info msg"]
    stderr_lines = ["warning msg", "critical msg"]

    with redirect_c_output(
        substrings=["debug", "warning"],
        filter_stdout=True,
        filter_stderr=False,
    ):
        simulate_c_output(stdout_lines, stderr_lines)

        out, err = capfd.readouterr()

        assert "debug msg" not in out
        assert "info msg" in out
        assert "warning msg" in err
        assert "critical" in err


def test_redirect_c_output_filtering_stderr_only(capfd):
    stdout_lines = ["debug msg", "info msg"]
    stderr_lines = ["warning msg", "critical msg"]

    with redirect_c_output(
        substrings=["debug", "warning"],
        filter_stdout=False,
        filter_stderr=True,
    ):
        simulate_c_output(stdout_lines, stderr_lines)

        out, err = capfd.readouterr()

        assert "debug" in out
        assert "info" in out
        assert "warning" not in err
        assert "critical" in err


def test_redirect_c_output_filtering_both_streams(capfd):
    stdout_lines = ["debug msg", "info msg"]
    stderr_lines = ["warning msg", "critical msg"]

    with redirect_c_output(
        substrings=["debug", "warning"],
        filter_stdout=True,
        filter_stderr=True,
    ):
        simulate_c_output(stdout_lines, stderr_lines)

        out, err = capfd.readouterr()

        assert "debug" not in out
        assert "info" in out
        assert "warning" not in err
        assert "critical" in err


def test_filter_stdout_decorator_behavior(capfd):
    @filter_stdout("nuisance", "annoying")
    def decorated_func():
        simulate_c_output(
            stdout_lines=["useful msg", "nuisance msg"],
            stderr_lines=["important error", "annoying warning"],
        )

    decorated_func()

    out, err = capfd.readouterr()

    assert "useful msg" in out
    assert "nuisance msg" not in out
    assert "important error" in err
    assert "annoying warning" not in err
