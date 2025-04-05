from importlib.metadata import version

from casaconfig import config

__version__ = version("radio-dstools")

config.logfile = "/dev/null"
config.measures_auto_update = False
config.data_auto_update = False
