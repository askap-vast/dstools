import click
import os

import dstools


@click.command()
@click.argument("ms")
def main(ms):
    root = dstools.__path__[0]

    # Fix beam pointing
    os.system(f"python {root}/fix_dir.py {ms}")

    # Convert instrumental pol visibilities from average to total flux
    os.system(f"python {root}/rescale_askapsoft.py {ms}")

    return


if __name__ == "__main__":
    main()
