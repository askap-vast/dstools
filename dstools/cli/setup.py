import logging
import os
import shutil
import subprocess
from pathlib import Path

import dstools
from dstools.logger import setupLogger

logger = logging.getLogger(__name__)
setupLogger(verbose=False)


def main():
    casa_bin = shutil.which("casa")
    casa_dir = Path(casa_bin).parent.parent

    # Attempt import of DStools and dependencies to test installation
    cmd = [
        casa_bin,
        "--nologger",
        "--nologfile",
        "-c",
        "import dstools",
        "import astropy",
        "import click",
    ]
    out = subprocess.run(cmd, capture_output=True)
    dstools_installed = "ModuleNotFoundError" not in str(out.stderr)

    if dstools_installed:
        logger.info(f"DStools already installed into CASA environment at {casa_dir}")
        exit()

    # Get installation paths of DStools and CASA python distribution
    cmd = [
        "casa",
        "--nologger",
        "--nologfile",
        "-c",
        "import numpy",
        "print(numpy.__path__[0])",
    ]
    out = str(subprocess.run(cmd, capture_output=True).stdout).split(r"\n")

    dstools_python_path = Path(dstools.__path__[0]).parent
    dstools_casa_path = Path([line for line in out if "numpy" in line][0]).parent

    # Copy DStools and dependencies into CASA python environment
    packages = [
        "dstools",
        "astropy",
        "click",
    ]
    for package in packages:
        os.system(f"cp -r {dstools_python_path}/{package}* {dstools_casa_path}")

    logger.info(f"DStools installed into CASA environment at {casa_dir}")


if __name__ == "__main__":
    main()
