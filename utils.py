import os

def colored(msg):

    return "\033[91m{}\033[0m".format(msg)


def prompt(msg):

    msg = msg[:-1] + " (y/n)?\n"

    resp = input(colored(msg))
    if resp not in ["y", "n"]:
        resp = input(colored(msg))

    return True if resp == "y" else False


def nearest_power(number):

    # Ensure number is even
    number = number // 2 * 2

    # Get nearest multiple of 5
    number = number // 5 * 5

    if number == 4090:
        number = 4096

    return int(number)


def update_param(name, val, dtype):
    while True:
        newval = input("Enter new {} (currently {}): ".format(name, val))

        # Accept old value if nothing entered
        if not newval:
            break

        # Test for correct input type
        try:
            val = dtype(newval)
            break
        except ValueError:
            print("{} must be of type {}".format(name, type(val)))

    return val


def resolve_array_config(band, config):
    """Determine reffreq, primary beam and cell sizes from array parameters."""

    wavelengths = {
        "low": 0.4026,
        "mid": 0.2450,
        "L": 0.0967,
        "C": 0.0461,
        "X": 0.0200,
    }
    frequencies = {
        "low": "888.49",
        "mid": "1367.49",
        "L": "2100",
        "C": "5500",
        "X": "9000",
    }
    primary_beams = {
        "low": 1,
        "mid": 0.8,
        "L": 0.75,
        "C": 0.25,
        "X": 0.25,
    }
    max_baselines = {
        "6km": 6000,
        "750_no6": 750,
        "750_6": 5020,
        "H168": 185,
    }

    freq = frequencies[band]

    wavelength = wavelengths[band]
    baseline = max_baselines[config]

    # Convert to resolution in arcsec
    resolution = wavelength / baseline * 206264.8

    imradius = primary_beams[band]
    cell = round(resolution / 5, 2)

    return freq, imradius, cell


def import_data(input_file, proj_dir, msname, reimport=False):

    if reimport:
        os.system(f"mv {proj_dir}/{msname} {proj_dir}/{msname}.bak")
        os.system(f"rm -r {proj_dir}/{msname}.flagversions")


    if input_file.endswith('.calfits'):
        importuvfits(
            fitsfile=input_file,
            vis=f"{proj_dir}/{msname}",
        )
    elif input_file.endswith('.cal'):
        importmiriad(
            mirfile=input_file,
            vis=f"{proj_dir}/{msname}",
        )
    elif input_file.endswith('.ms'):
        os.system(f"cp -r {input_file} {proj_dir}/{msname}")

    return
