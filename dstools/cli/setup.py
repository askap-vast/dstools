import glob
import logging
import shutil
import subprocess
from pathlib import Path

from dstools.logger import setupLogger

logger = logging.getLogger(__name__)
setupLogger(verbose=False)


def main():
    casa_bin = shutil.which("casa")
    casa_dir = Path(casa_bin).parent.parent
    dstools_pkgs = glob.glob(f"{casa_dir}/lib/py/lib/python*/site-packages/*dstools*")
    dstools_installed = len(dstools_pkgs) > 0
    if len(dstools_pkgs) > 0:
        logger.error(f"DStools already installed into CASA environment at {casa_dir}")
        exit(1)

    cmd = "import sys\nimport subprocess\nsubprocess.check_call([sys.executable, '-m', 'pip', 'install', 'py-dstools'])"
    call = [
        casa_bin,
        "--nologger",
        "--nologfile",
        "-c",
        cmd,
    ]
    subprocess.run(
        call,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    logger.info(f"DStools installed into CASA environment at {casa_dir}")


if __name__ == "__main__":
    main()
