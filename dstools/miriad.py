import logging
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import astropy.units as u
import click
import dstools
import pandas as pd
from casaconfig import config
from dstools.logger import setupLogger, parse_stdout_stderr
from dstools.utils import prompt

config.logfile = "/dev/null"
from casatasks import importuvfits, listobs

logger = logging.getLogger(__name__)

package_root = dstools.__path__[0]


def calc_total_scantime(header):
    scans = [key for key in header.keys() if "scan" in key]

    time = lambda scan: (scan["EndTime"] - scan["BeginTime"]) * u.day.to(u.minute)

    df = pd.DataFrame({scan: header[scan]["0"] for scan in scans}).T
    df["ScanTime"] = (df.EndTime - df.BeginTime) * u.day.to(u.min)
    scantimes = df.groupby("FieldName").ScanTime.sum().reset_index()

    return scantimes.sort_values("ScanTime")


@dataclass
class MiriadWrapper:
    out_dir: Path
    band: str
    project_code: str
    primary_cal: str = "1934-638"
    gain_cal: str = None
    target: str = None
    leakage_cal: str = None
    mfinterval: float = 1.0
    bpinterval: float = 1.0
    gpinterval: float = 0.1
    refant: int = 3
    noflag: bool = False
    strong_pol: bool = False
    IF: str = None
    data_dir: Path = Path(".")
    keep_intermediate: bool = False
    verbose: bool = False

    def __post_init__(self):
        if self.IF is None:
            self.IF = self.band.IF

        self.opts = {
            "mfinterval": self.mfinterval,
            "bpinterval": self.bpinterval,
            "gpinterval": self.gpinterval,
            "ifsel": self.IF,
            "spec": self.band.spec,
            "refant": self.refant,
        }

    def run_command(self, command, args=None):
        if args is not None:
            args = " ".join([f"{arg}" for arg in args])
        else:
            args = ""

        exports = " ".join(f"export {opt}='{val}';" for opt, val in self.opts.items())

        p = subprocess.Popen(
            f"source {package_root}/functions.sh; {exports} {command} {args}",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )

        parse_stdout_stderr(p, logger)

        return

    def load_data(self, shiftra, shiftdec):
        logger.info(f"Loading RPFITS files from {self.data_dir}")
        self.run_command(
            "load_data",
            args=[
                str(self.out_dir.absolute()),
                str(self.data_dir.absolute()),
                "true" if self.noflag else "false",
                str(self.project_code),
                shiftra,
                shiftdec,
            ],
        )

        self._get_header()

        return

    def path(self, target):
        return Path(self.out_dir / "miriad" / f"{target}.{self.band.freq}")

    def _get_header(self):

        uv = self.out_dir / "miriad" / f"{self.project_code}.uv"
        ms = self.generate_ms(uv)

        header = listobs(vis=ms)
        total_times = calc_total_scantime(header)
        fields = {header[f"field_{i}"]["name"]: i for i in range(header["nfields"])}

        if self.primary_cal not in fields:
            raise ValueError(f"Could not locate {self.primary_cal} field in {uv}")

        # Assume target and gain calibrator are two objects with most observing time
        if self.gain_cal is None:
            gain_cal = total_times.iloc[-2].FieldName
            logger.info(f"No gain calibrator specified, defaulting to {gain_cal}")
        if self.target is None:
            target = total_times.iloc[-1].FieldName
            logger.info(f"No science target specified, defaulting to {target}")

        self.target_paths = {
            "primary_cal": self.path(self.primary_cal),
            "gain_cal": self.path(gain_cal),
            "leakage_cal": self.path(self.leakage_cal) if self.leakage_cal else None,
            "target": self.path(target),
        }

        # Clean up intermediate files
        self.run_command(f"rm -r {uv} {ms}")

        return

    def blflag(self, vis, x, y, options):
        logger.info(f"Manually flagging {vis.name} in {y} vs {x}")
        self.run_command("manflag", args=[vis, x, y, options])

    def autoflag(self, vis):
        logger.info(f"Autoflagging {vis.name}")
        self.run_command("autoflag", args=[vis])

    def bandpass(self, vis):
        logger.info(f"Running bandpass calibration on {vis.name}")
        self.run_command("cal_bandpass", args=[vis])

    def bootstrap(self, vis1, vis2):
        if vis1 != vis2:
            logger.info(f"Bootstrapping flux scale from {vis2.name} to {vis1.name}")
            self.run_command("bootstrap", args=[vis1, vis2])

    def gaincal(self, vis, options):
        logger.info(f"Running gain/leakage calibration on {vis.name}")
        self.run_command("cal_gains", args=[vis, options])

    def copycal(self, vis1, vis2):
        if vis1 != vis2:
            logger.info(f"Copying calibration tables from {vis1.name} to {vis2.name}")
            self.run_command("copy_cal", args=[vis1, vis2])

    def gpaver(self, vis):
        logger.info(f"Averaging calibration solutions for {vis.name}")
        self.run_command("average_gains", args=[vis])

    def uvaver(self, vis):
        logger.info(f"Applying calibration solutions to {vis.name}")
        self.run_command("apply_gains", args=[vis])

    def generate_ms(self, uv):
        fitsfile = f"{uv}.fits"
        ms = f"{uv}.ms"

        self.run_command("uvtofits", args=[uv, fitsfile])

        importuvfits(
            fitsfile=fitsfile,
            vis=ms,
        )

        self.run_command(f"rm -r {fitsfile}")

        return ms

    def cleanup(self):
        self.run_command(f"rm -r {self.out_dir / 'miriad'}")


@dataclass
class CABBContinuumPipeline:
    miriad: MiriadWrapper
    shiftra: float
    shiftdec: float
    interactive: bool

    def __post_init__(self):
        self.miriad.load_data(self.shiftra, self.shiftdec)

        self.primary_cal = self.miriad.target_paths.get("primary_cal")
        self.gain_cal = self.miriad.target_paths.get("gain_cal")
        self.target = self.miriad.target_paths.get("target")
        self.leakage_cal = self.miriad.target_paths.get("leakage_cal")

    def flag_sequence(self, target):
        self.miriad.blflag(
            target,
            x="time",
            y="amp",
            options="nofqav,nobase",
        )
        self.miriad.blflag(
            target,
            x="chan",
            y="amp",
            options="nofqav,nobase",
        )
        self.miriad.blflag(
            target,
            x="chan",
            y="amp",
            options="nofqav",
        )
        self.miriad.autoflag(target)
        self.miriad.blflag(
            target,
            x="real",
            y="imag",
            options="nofqav,nobase",
        )

        return

    def run(self):

        # Primary bandpass / flux calibrator
        # ---------------------------------
        self.miriad.bandpass(self.primary_cal)

        # Flag and solve for bandpass on primary calibrator
        if self.interactive:
            while prompt(f"Do more flagging on {self.primary_cal.name}?"):
                self.flag_sequence(self.primary_cal)
                self.miriad.bandpass(self.primary_cal)
        else:
            self.miriad.autoflag(self.primary_cal)
            self.miriad.bandpass(self.primary_cal)

        # Solve for primary calibrator gains / leakage
        self.miriad.gaincal(self.primary_cal, options="xyvary")

        # Set options to work with strong or weakly polarised calibrator
        if self.miriad.strong_pol:
            gp_options = "xyvary,qusolve,vsolve,xyref,polref"
        else:
            gp_options = "xyvary,qusolve"

        # Leakage calibrator
        # ------------------
        if self.leakage_cal is not None:

            self.miriad.copycal(self.primary_cal, self.leakage_cal)

            # Flag and solve for gains / leakages / xy-phase on leakage calibrator
            if self.interactive:
                while prompt(f"Do more flagging on {self.leakage_cal.name}?"):
                    self.flag_sequence(self.leakage_cal)
                    self.miriad.gaincal(self.leakage_cal, options=gp_options)
            else:
                self.miriad.autoflag(self.leakage_cal)
                self.miriad.gaincal(self.leakage_cal, options=gp_options)

            # To avoid corruption of Stokes V zero-point, we copy solutions
            # back to primary calibrator and repeat the sequence
            self.miriad.copycal(self.leakage_cal, self.primary_cal)
            self.miriad.gaincal(self.primary_cal, options="noxy")
            self.miriad.copycal(self.primary_cal, self.leakage_cal)
            self.miriad.gaincal(self.leakage_cal, options=gp_options)

            # Now we turn off xy-phase / leakage calibration for subsequent gain calibration
            self.miriad.copycal(self.leakage_cal, self.gain_cal)
            gp_options = "qusolve,noxy,nopol"
        else:
            self.miriad.copycal(self.primary_cal, self.gain_cal)

        # Secondary gain calibrator
        # -------------------------
        if self.interactive:
            while prompt(f"Do more flagging on {self.gain_cal.name}?"):
                self.flag_sequence(self.gain_cal)
                self.miriad.gaincal(self.gain_cal, options=gp_options)
        else:
            self.miriad.autoflag(self.gain_cal)
            self.miriad.gaincal(self.gain_cal, options=gp_options)

        self.miriad.bootstrap(self.gain_cal, self.primary_cal)
        self.miriad.copycal(self.gain_cal, self.target)

        # Science target
        # --------------
        if not self.miriad.noflag:
            self.miriad.autoflag(self.target)

        # Average solutions and apply
        self.miriad.gpaver(self.target)
        self.miriad.uvaver(self.gain_cal)
        self.miriad.uvaver(self.target)

        self.miriad.generate_ms(f"{self.target}.cal")
        self.miriad.run_command(
            f"mv {self.target}.cal.ms {self.miriad.out_dir}/{self.target.name}.ms"
        )

        return


if __name__ == "__main__":
    main()
