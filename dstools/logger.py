import logging
from typing import Optional

import colorlog


def setupLogger(verbose: bool, filename: Optional[str] = None) -> None:
    level = logging.DEBUG if verbose else logging.INFO

    # Get root logger disable any existing handlers, and set level
    root_logger = logging.getLogger()
    root_logger.handlers = []
    root_logger.setLevel(level)

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

    root_logger.addHandler(stream_handler)

    return
