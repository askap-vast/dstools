import logging
import os
import selectors
import threading
from contextlib import contextmanager
from functools import wraps
from typing import Callable, Iterable, Optional

import colorlog


def setupLogger(verbose: bool, filename: Optional[str] = None) -> None:
    level = logging.DEBUG if verbose else logging.INFO

    # Get root logger disable any existing handlers, and set level
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.handlers = []

    # Turn off some bothersome verbose logging modules
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
    logging.getLogger("urllib").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("h5py").setLevel(logging.INFO)

    if filename:
        formatter = logging.Formatter(
            "%(levelname)-8s %(asctime)s - %(name)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler = logging.FileHandler(filename)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)

        root_logger.addHandler(file_handler)

    colorformatter = colorlog.ColoredFormatter(
        "%(log_color)s%(levelname)-8s%(reset)s %(asctime)s - %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        reset=True,
        log_colors={
            "DEBUG": "cyan",
            "INFO": "green",
            "WARNING": "yellow",
            "ERROR": "red",
            "CRITICAL": "red,bg_white",
        },
    )

    stream_handler = colorlog.StreamHandler()
    stream_handler.setFormatter(colorformatter)
    stream_handler.setLevel(level)

    root_logger.addHandler(stream_handler)

    return


def parse_stdout_stderr(process, logger, print_stdout: bool = False):
    """Parse STDOUT and STDERR from a subprocess and redirect to logger."""

    sel = selectors.DefaultSelector()
    sel.register(process.stdout, selectors.EVENT_READ)
    sel.register(process.stderr, selectors.EVENT_READ)

    # Filter uninteresting warnings and set to DEBUG level
    debug_lines = [
        "### Warning:  Using post-Aug94 ATCA flux scale for 1934-638",
        "### Warning:  Correlations flagged or edge-rejected:",
        "PGPLOT /png: writing new file as",
    ]

    lines_to_parse = True
    while lines_to_parse:
        for key, val in sel.select():
            line = key.fileobj.readline()
            if not line:
                lines_to_parse = False
                break

            line = line.decode().rstrip()
            debug_line = any(debug_str in line for debug_str in debug_lines)

            if print_stdout:
                print(line)
            elif debug_line or key.fileobj is process.stdout:
                logger.debug(line)
            else:
                logger.warning(line.replace("### Warning:  ", ""))

    return


logger = logging.getLogger(__name__)


def filter_pipe_output(
    pipe_r: int,
    substrings: Iterable[str],
    original_stream_fd: int,
    stream: str,
) -> None:
    """Read lines from a pipe file descriptor and filter to debug/warning."""

    with os.fdopen(pipe_r) as pipe:
        for line in pipe:
            if any(substr in line for substr in substrings):
                continue
            else:
                os.write(original_stream_fd, line.encode())

    return


@contextmanager
def redirect_c_output(
    substrings: Iterable[str],
    filter_stdout: bool = True,
    filter_stderr: bool = True,
):
    """A context manager to redirect STDOUT / STDERR for both python and C level output."""

    # Save original file descriptors
    original_stdout_fd = os.dup(1) if filter_stdout else None
    original_stderr_fd = os.dup(2) if filter_stderr else None

    # Create pipes for capturing output
    stdout_pipe_r, stdout_pipe_w = os.pipe() if filter_stdout else (None, None)
    stderr_pipe_r, stderr_pipe_w = os.pipe() if filter_stderr else (None, None)

    # Start threads to read and filter the output
    stdout_thread = None
    if filter_stdout:
        stdout_thread = threading.Thread(
            target=filter_pipe_output,
            args=(stdout_pipe_r, substrings, original_stdout_fd, "STDOUT"),
            daemon=True,
        )
        stdout_thread.start()

    stderr_thread = None
    if filter_stderr:
        stderr_thread = threading.Thread(
            target=filter_pipe_output,
            args=(stderr_pipe_r, substrings, original_stderr_fd, "STDERR"),
            daemon=True,
        )
        stderr_thread.start()

    try:
        # Redirect STDOUT / STDERR to pipe for filtering
        if filter_stdout:
            os.dup2(stdout_pipe_w, 1)
        if filter_stderr:
            os.dup2(stderr_pipe_w, 2)

        yield

    finally:
        # Restore original file descriptors
        if filter_stdout:
            os.dup2(original_stdout_fd, 1)
        if filter_stderr:
            os.dup2(original_stderr_fd, 2)

        # Close the write ends of the pipes to signal EOF to the threads
        if filter_stdout:
            os.close(stdout_pipe_w)
        if filter_stderr:
            os.close(stderr_pipe_w)

        # Wait for threads to finish
        if stdout_thread:
            stdout_thread.join()
        if stderr_thread:
            stderr_thread.join()

        # Close the original file descriptors and pipes
        if filter_stdout:
            os.close(original_stdout_fd)
        if filter_stderr:
            os.close(original_stderr_fd)

    return


def filter_stdout(
    *substrings: str,
    filter_stdout: bool = True,
    filter_stderr: bool = True,
):
    """A decorator to filter C-level STDOUT/STDERR from within CASA function calls."""

    def decorator(func: Callable):
        @wraps(func)
        def wrapper(*args, **kwargs):
            with redirect_c_output(
                substrings,
                filter_stdout=filter_stdout,
                filter_stderr=filter_stderr,
            ):
                return func(*args, **kwargs)

        return wrapper

    return decorator
