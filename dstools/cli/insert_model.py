import logging
import os
from pathlib import Path

import click
from dstools.imaging import WSCleanModel
from dstools.logger import setupLogger
from dstools.viewer import Image, Viewer

logger = logging.getLogger(__name__)


@click.command(context_settings={"show_default": True})
@click.argument("model_dir", type=Path)
@click.argument("ms", type=Path)
def main(model_dir, model_format, ms):

    setupLogger(verbose=False)

    if not model_dir.exists():
        logger.error(f"Path {model_dir} does not exist.")
        exit(1)

    # Read model images in
    model = WSCleanModel(model_dir)

    # Set up interactive masking session
    logger.info("Launching interactive mask viewer...")
    print("  click to draw mask polygon vertices")
    print("  press 'x' to mask pixels within polygon")
    print("  press 'c' to unmask pixels within polygon")
    print("  close viewer when happy with masking.")

    images = [
        Image(name="image", path=model.image),
        Image(name="residual", path=model.residual),
        Image(name="model", path=model.model),
    ]
    viewer = Viewer(images=images)

    # Apply mask edits to output model images
    logger.info(f"Masking all model images in {model_dir}.")
    model.apply_mask(viewer.mask)

    # Insert masked visibilities into MODEL_DATA column
    logger.info(f"Inserting masked model into {ms}.")
    model.insert_into(ms)

    return


if __name__ == "__main__":
    main()
