import logging
from pathlib import Path

import click

from dstools.imaging import WSCleanModel
from dstools.logger import setupLogger
from dstools.ms import MeasurementSet
from dstools.utils import parse_coordinates

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.option(
    "-I",
    "--interactive",
    is_flag=True,
    default=False,
    help="Disable automatic target mask and draw interactively.",
)
@click.option(
    "-p",
    "--mask-pos",
    type=str,
    nargs=2,
    default=None,
    help=(
        "Coordinates around which to apply automatic target mask."
        "(provide as separate values, e.g. -p <RA> <DEC>). Default to phasecentre."
    ),
)
@click.option(
    "-r",
    "--mask-radius",
    default=None,
    type=float,
    help="Radius of automatic target mask in arcseconds.",
)
@click.argument("model_dir", type=Path)
@click.argument("ms", type=MeasurementSet)
def main(mask_pos, mask_radius, interactive, model_dir, ms):
    setupLogger(verbose=False)

    if not model_dir.exists():
        logger.error(f"Path {model_dir} does not exist.")
        exit(1)

    # Read model images in
    model = WSCleanModel(model_dir)

    # Generate a mask array
    if interactive:
        mask = model.get_auto_mask()
    else:
        mask_pos = ms.phasecentre if mask_pos is None else parse_coordinates(mask_pos)
        mask = model.get_circular_mask(mask_pos, mask_radius)

    # Apply mask edits to output model images
    logger.info(f"Masking all model images in {model_dir}.")
    model.apply_mask(mask)

    # Insert masked visibilities into MODEL_DATA column
    logger.info(f"Inserting masked model into {ms}.")
    model.insert_into(ms)

    return


if __name__ == "__main__":
    main()
