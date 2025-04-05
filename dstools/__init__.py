from importlib.metadata import version

from casaconfig import config

__version__ = version("radio-dstools")

config.logfile = "/dev/null"
