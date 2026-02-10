import pickle
import subprocess
import sys
import tempfile
from functools import wraps

from dstools.logger import filter_stdout


def run_in_subprocess(taskname):
    """Run CASA task in subprocess to avoid binding clash with python-casacore in main process."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Load a temporary file to dump return contents
            with tempfile.NamedTemporaryFile(suffix=".pkl") as tf:
                # Get byte representation of task arguments
                args = pickle.dumps((args, kwargs)).hex()

                # Execute task via dstools.casa.main() in a subprocess
                casa_cmd = [
                    sys.executable,
                    "-m",
                    "dstools.casa",
                    taskname,
                    tf.name,
                    args,
                ]
                subprocess.run(
                    casa_cmd,
                    stdout=None,
                    stderr=None,
                )

                # Load return values from temp file
                with open(tf.name, "rb") as f:
                    result = pickle.load(f)

            return result

        return wrapper

    return decorator


def casatask_subprocess():
    """Subprocess called CLI function to isolate loading of CASA bindings to casacore."""

    # Load casatasks imports within subprocess
    from casatasks import __dict__ as import_dict

    taskname = sys.argv[1]
    return_data_path = sys.argv[2]

    # Decode task arguments
    arg_bytes = sys.argv[3]
    args, kwargs = pickle.loads(bytes.fromhex(arg_bytes))

    # Import and run the task inside subprocess
    casatask = import_dict.get(taskname)
    if casatask is None:
        raise NotImplementedError(f"Task {taskname} not importable from casatasks.")

    result = casatask(*args, **kwargs)

    # Write return value to temporary file
    with open(return_data_path, "wb") as f:
        pickle.dump(result, f)

    return


@run_in_subprocess("applycal")
def applycal(*args, **kwargs):
    pass


@run_in_subprocess("clearcal")
@filter_stdout("Forcing use of OLD VisibilityIterator.")
def clearcal(*args, **kwargs):
    pass


@run_in_subprocess("cvel")
@filter_stdout("combineSpws progress")
def cvel(*args, **kwargs):
    pass


@run_in_subprocess("exportuvfits")
def exportfits(*args, **kwargs):
    pass


@run_in_subprocess("flagdata")
def flagdata(*args, **kwargs):
    pass


@run_in_subprocess("gaincal")
def gaincal(*args, **kwargs):
    pass


@run_in_subprocess("importuvfits")
@filter_stdout(
    "XYZHAND keyword not found in AN table.",
    "No systemic velocity",
    "No rest frequency",
)
def importuvfits(*args, **kwargs):
    pass


@run_in_subprocess("imsubimage")
def imsubimage(*args, **kwargs):
    pass


@run_in_subprocess("listobs")
def listobs(*args, **kwargs):
    pass


@run_in_subprocess("mstransform")
@filter_stdout("There is only one selected SPW, no need to combine")
def mstransform(*args, **kwargs):
    pass


@run_in_subprocess("phaseshift")
def phaseshift(*args, **kwargs):
    pass


@run_in_subprocess("split")
def split(*args, **kwargs):
    pass


@run_in_subprocess("exportuvfits")
def exportuvfits(*args, **kwargs):
    pass


@run_in_subprocess("concat")
def concat(*args, **kwargs):
    pass


@run_in_subprocess("tclean")
@filter_stdout("Restoring with an empty model image")
def tclean(*args, **kwargs):
    pass


@run_in_subprocess("uvsub")
def uvsub(*args, **kwargs):
    pass


if __name__ == "__main__":
    casatask_subprocess()
