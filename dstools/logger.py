import logging
import selectors
from typing import Optional

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

    sel = selectors.DefaultSelector()
    sel.register(process.stdout, selectors.EVENT_READ)
    sel.register(process.stderr, selectors.EVENT_READ)

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
            debug_line = any(l in line for l in debug_lines)

            if print_stdout:
                print(line)
            elif debug_line or key.fileobj is process.stdout:
                logger.debug(line)
            else:
                logger.warning(line.replace("### Warning:  ", ""))

    return
