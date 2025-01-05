import logging
import os
from dataclasses import dataclass
from pathlib import Path

import click
from dstools.logger import setupLogger
from dstools.miriad import CABBContinuumPipeline, MiriadWrapper

logger = logging.getLogger(__name__)


@dataclass
class Band:
    freq: str
    spec: str
    IF: str


BANDS = {
    "L": Band(freq="2100", spec="2.1", IF="1"),
    "C": Band(freq="5500", spec="5.5", IF="1"),
    "X": Band(freq="9000", spec="9.0", IF="2"),
    "K": Band(freq="12000", spec="12.0", IF="1"),
}


@click.command(context_settings={"show_default": True})
@click.option("-B", "--band", type=str, default="L")
@click.option("-p", "--primary-cal", type=str, default="1934-638")
@click.option("-s", "--gain-cal", type=str, default=None)
@click.option("-t", "--target", type=str, default=None)
@click.option("-l", "--leakage-cal", type=str, default=None)
@click.option(
    "-m",
    "--mfinterval",
    default="1.0",
    type=str,
    help="Time interval to solve for antenna gains in bandpass calibration.",
)
@click.option(
    "-b",
    "--bpinterval",
    default="1.0",
    type=str,
    help="Time interval to solve for bandpass in bandpass calibration.",
)
@click.option(
    "-g",
    "--gpinterval",
    default="0.1",
    type=str,
    help="Time interval to solve for antenna gains in gain calibration.",
)
@click.option(
    "-f",
    "--nfbin",
    default="4",
    type=str,
    help="Number of frequency subbands in which to solve for gains/leakage.",
)
@click.option(
    "-n",
    "--num-flag-rounds",
    default=1,
    type=int,
    help="Number of rounds in each autoflagging / calibration loop.",
)
@click.option(
    "-r",
    "--refant",
    type=click.Choice(["1", "2", "3", "4", "5", "6"]),
    default="3",
    help="Reference antenna.",
)
@click.option(
    "--shiftra",
    type=str,
    default="0",
    help="Offset between pointing and phase centre right ascension in arcsec.",
)
@click.option(
    "--shiftdec",
    type=str,
    default="0",
    help="Offset between pointing and phase centre declination in arcsec.",
)
@click.option(
    "-P",
    "--strong-pol",
    is_flag=True,
    default=False,
    help="Solve for absolute XY phase and leakage (requires good leakage calibrator)",
)
@click.option(
    "-F",
    "--noflag",
    is_flag=True,
    default=False,
    help="Disable birdie and rfiflag options in atlod and avoid target flagging.",
)
@click.option(
    "-I",
    "--interactive",
    is_flag=True,
    default=False,
    help="Run calibration pipeline interactively with manual flagging.",
)
@click.option(
    "-o",
    "--out-dir",
    type=Path,
    default=Path("."),
    help="Path to store calibrated MeasurementSet.",
)
@click.option(
    "-S",
    "--skip-pipeline",
    is_flag=True,
    default=False,
    help="Skip execution of flagging/calibration pipeline.",
)
@click.option(
    "-L",
    "--savelogs",
    is_flag=True,
    default=False,
    help="Store processing logs.",
)
@click.option(
    "-k",
    "--keep-intermediate",
    is_flag=True,
    default=False,
    help="Store intermediate files produced by miriad.",
)
@click.option("-v", "--verbose", is_flag=True, default=False)
@click.argument("data_dir")
@click.argument("project_code")
def main(
    data_dir,
    project_code,
    band,
    primary_cal,
    gain_cal,
    target,
    leakage_cal,
    strong_pol,
    mfinterval,
    bpinterval,
    gpinterval,
    nfbin,
    num_flag_rounds,
    refant,
    shiftra,
    shiftdec,
    noflag,
    interactive,
    out_dir,
    skip_pipeline,
    savelogs,
    keep_intermediate,
    verbose,
):

    os.system(f"mkdir -p {out_dir}")

    logfile = out_dir / "atca-cal.log" if savelogs else None
    setupLogger(verbose=verbose, filename=logfile)

    miriad = MiriadWrapper(
        data_dir=Path(data_dir),
        band=BANDS.get(band),
        project_code=project_code,
        out_dir=out_dir,
        strong_pol=strong_pol,
        mfinterval=mfinterval,
        bpinterval=bpinterval,
        gpinterval=gpinterval,
        nfbin=nfbin,
        refant=refant,
        noflag=noflag,
        verbose=verbose,
    )

    pipeline = CABBContinuumPipeline(
        miriad=miriad,
        shiftra=shiftra,
        shiftdec=shiftdec,
        num_flag_rounds=num_flag_rounds,
        interactive=interactive,
    )

    try:
        pipeline.miriad.set_targets(
            primary_cal=primary_cal,
            leakage_cal=leakage_cal,
            gain_cal=gain_cal,
            target=target,
        )
    except ValueError as e:
        logger.error(e)
        exit(1)

    if not skip_pipeline:
        pipeline.run()


    if not keep_intermediate:
        miriad.cleanup()

    return


if __name__ == "__main__":
    main()
