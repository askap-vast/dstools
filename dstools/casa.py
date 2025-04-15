import casatasks

from dstools.logger import filter_stdout


def applycal(*args, **kwargs):
    return casatasks.applycal(*args, **kwargs)


@filter_stdout("Forcing use of OLD VisibilityIterator.")
def clearcal(*args, **kwargs):
    return casatasks.clearcal(*args, **kwargs)


@filter_stdout("combineSpws progress")
def cvel(*args, **kwargs):
    return casatasks.cvel(*args, **kwargs)


def exportfits(*args, **kwargs):
    return casatasks.exportfits(*args, **kwargs)


def flagdata(*args, **kwargs):
    return casatasks.flagdata(*args, **kwargs)


def gaincal(*args, **kwargs):
    return casatasks.gaincal(*args, **kwargs)


@filter_stdout(
    "XYZHAND keyword not found in AN table.",
    "No systemic velocity",
    "No rest frequency",
)
def importuvfits(*args, **kwargs):
    return casatasks.importuvfits(*args, **kwargs)


def imsubimage(*args, **kwargs):
    return casatasks.imsubimage(*args, **kwargs)


def listobs(*args, **kwargs):
    return casatasks.listobs(*args, **kwargs)


@filter_stdout("There is only one selected SPW, no need to combine")
def mstransform(*args, **kwargs):
    return casatasks.mstransform(*args, **kwargs)


def phaseshift(*args, **kwargs):
    return casatasks.phaseshift(*args, **kwargs)


def split(*args, **kwargs):
    return casatasks.split(*args, **kwargs)


@filter_stdout("Restoring with an empty model image")
def tclean(*args, **kwargs):
    return casatasks.tclean(*args, **kwargs)


def uvsub(*args, **kwargs):
    return casatasks.uvsub(*args, **kwargs)
