from importlib.metadata import version

__version__ = version("radio-dstools")

try:
    from casaconfig import config
except ImportError:
    config = None
else:
    config.logfile = "/dev/null"
