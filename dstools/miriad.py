import glob
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
from dstools.logger import parse_stdout_stderr, setupLogger
from dstools.utils import prompt

config.logfile = "/dev/null"
from casatasks import importuvfits, listobs

logger = logging.getLogger(__name__)

package_root = dstools.__path__[0]


@dataclass
class Target:
    name: str
    path: Path


@dataclass
class MiriadWrapper:
    data_dir: Path
    band: str
    project_code: str
    mfinterval: float = 1.0
    bpinterval: float = 1.0
    gpinterval: float = 0.1
    nfbin: int = 4
    refant: int = 3
    noflag: bool = False
    strong_pol: bool = False
    IF: str = None
    out_dir: Path = Path(".")
    verbose: bool = False

    def __post_init__(self):
        if self.IF is None:
            self.IF = self.band.IF

        self.uvfile = self.out_dir / "miriad" / f"{self.project_code}.uv"
        self.opts = {
            "mfinterval": self.mfinterval,
            "bpinterval": self.bpinterval,
            "gpinterval": self.gpinterval,
            "nfbin": self.nfbin,
            "ifsel": self.IF,
            "spec": self.band.spec,
            "refant": self.refant,
        }
        logger.debug(f"Miriad options set to:")
        for k, v in self.opts.items():
            logger.debug(f"{k}={v}")

    def run_command(self, command, args=None, print_stdout=False):
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

        parse_stdout_stderr(p, logger, print_stdout)

        return

    def load_data(self, shiftra: str, shiftdec: str):
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

        return

    def path(self, target: str):
        if target is None:
            return Target(name=None, path=None)

        path = Path(self.out_dir / "miriad" / f"{target}.{self.band.freq}")
        return Target(name=target, path=path)

    def set_targets(
        self,
        primary_cal: str = "1934-638",
        gain_cal: str = None,
        target: str = None,
        leakage_cal: str = None,
    ):

        ms = self.generate_ms(self.uvfile)

        # Read and store scan summary
        header = listobs(vis=ms)
        scan_keys = [key for key in header.keys() if "scan" in key]
        df = pd.DataFrame({scan: header[scan]["0"] for scan in scan_keys}).T
        df["ScanTime"] = (df.EndTime - df.BeginTime) * u.day.to(u.min)
        self.scans = df.sort_values("scanId").reset_index()

        fields = [
            str(Path(p).name).replace(f".{self.band.freq}", "")
            for p in glob.glob(f"{self.out_dir}/miriad/*.{self.band.freq}")
        ]

        # Clean up intermediate files
        os.system(f"rm -r {ms}")

        if primary_cal not in fields:
            raise ValueError(f"Could not locate {primary_cal} field in {self.uvfile}")

        # Assume target and gain calibrator are two objects with most observing time
        total_times = (
            df.groupby("FieldName").ScanTime.sum().reset_index().sort_values("ScanTime")
        )

        if gain_cal is None:
            gain_cal = total_times.iloc[-2].FieldName
            logger.info(f"No gain calibrator specified, defaulting to {gain_cal}")
        if target is None:
            target = total_times.iloc[-1].FieldName
            logger.info(f"No science target specified, defaulting to {target}")

        self.target_paths = {
            "primary_cal": self.path(primary_cal),
            "gain_cal": self.path(gain_cal),
            "leakage_cal": self.path(leakage_cal),
            "target": self.path(target),
        }

        return

    def blflag(self, vis, x, y, options):
        per_bl = "" if "nobase" in options else "per baseline"
        logger.info(f"Manually flagging {vis} in {y} vs {x} {per_bl}")
        self.run_command(
            "manflag",
            args=[vis, x, y, options],
            print_stdout=True,
        )

    def autoflag(self, vis):
        logger.info(f"Autoflagging {vis}")
        self.run_command("autoflag", args=[vis])

    def flag_times(self, vis, start_time, end_time):
        logger.info(f"Flagging {vis} between {start_time}-{end_time}")
        self.run_command("flag_timerange", args=[vis, start_time, end_time])

    def bandpass(self, vis):
        logger.info(f"Running bandpass calibration on {vis}")
        interpolate = "true" if self.noflag else "false"
        self.run_command("cal_bandpass", args=[vis, interpolate])

    def bootstrap(self, vis1, vis2):
        if vis1 != vis2:
            logger.info(f"Bootstrapping flux scale from {vis2.name} to {vis1.name}")
            self.run_command("bootstrap", args=[vis1, vis2])

    def gaincal(self, vis, options):
        logger.info(f"Running gain/leakage calibration on {vis}")
        self.run_command("cal_gains", args=[vis, options])

    def copycal(self, vis1, vis2):
        logger.info(f"Copying calibration tables from {vis1.name} to {vis2.name}")
        self.run_command("copy_cal", args=[vis1, vis2])

    def gpaver(self, vis):
        logger.info(f"Averaging calibration solutions for {vis}")
        self.run_command("average_gains", args=[vis])

    def uvaver(self, vis):
        logger.info(f"Applying calibration solutions to {vis}")
        self.run_command("apply_gains", args=[vis])

    def generate_ms(self, uv):
        fitsfile = f"{uv}.fits"
        ms = f"{uv}.ms"

        self.run_command("uvtofits", args=[uv, fitsfile])

        importuvfits(
            fitsfile=fitsfile,
            vis=ms,
        )

        os.system(f"rm -r {fitsfile}")

        return ms

    def cleanup(self):
        os.system(f"rm -r {self.out_dir / 'miriad'}")


@dataclass
class CABBContinuumPipeline:
    miriad: MiriadWrapper
    shiftra: float
    shiftdec: float
    num_flag_rounds: int
    interactive: bool

    def __post_init__(self):
        if os.path.exists(self.miriad.uvfile):
            logger.warning(f"{self.miriad.uvfile} already exists, will not overwrite")
            return

        self.miriad.load_data(self.shiftra, self.shiftdec)

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

    def make_diagnostics(self):
        logger.info("Generating calibration diagnostic plots")

        cwd = Path(".").absolute()
        os.system(f"mkdir -p {self.miriad.out_dir}/diagnostics")
        os.chdir(self.miriad.out_dir)

        pcal = self.miriad.target_paths.get("primary_cal").path.relative_to(
            self.miriad.out_dir
        )
        scal = self.miriad.target_paths.get("gain_cal").path.relative_to(
            self.miriad.out_dir
        )

        calibrators = set([pcal, scal])

        # Primary / secondary plots
        for vis in calibrators:
            self.miriad.run_command(
                f"uvfmeas vis={vis} \
                stokes=i \
                device=diagnostics/{vis.name}_spectrum.png/PNG"
            )
            self.miriad.run_command(
                f"uvplt vis={vis} \
                stokes=i  \
                axis=re,im \
                options=nofqav,nobase \
                device=diagnostics/{vis.name}_real_imag.png/PNG"
            )
            self.miriad.run_command(
                f"uvplt vis={vis} \
                stokes=i \
                axis=freq,phase \
                nxy=4 \
                options=nofqav \
                device=diagnostics/{vis.name}_phase_spectrum.png/PNG"
            )
            self.miriad.run_command(
                f"gpplt vis={vis} \
                options=gains \
                yaxis=amp,phase \
                log=diagnostics/{vis.name}_gains.txt"
            )
            self.miriad.run_command(
                f"gpplt vis={vis} \
                options=polarization \
                yaxis=amp,phase \
                log=diagnostics/{vis.name}_leakages.txt"
            )

        # Primary plots
        self.miriad.run_command(
            f"gpplt vis={pcal} \
            options=dots,bandpass \
            log=diagnostics/{vis.name}_bandpass.txt \
            device=diagnostics/{pcal.name}_bandpass.png/PNG"
        )

        # Secondary plots
        self.miriad.run_command(
            f"uvplt vis={scal} \
            stokes=i \
            axis=time,amp \
            options=nofqav,nobase \
            device=diagnostics/{scal.name}_amp_time.png/PNG"
        )
        self.miriad.run_command(
            f"uvplt vis={scal} \
            stokes=xx,yy \
            axis=time,phase \
            options=nofqav,nobase \
            device=diagnostics/{scal.name}_phase_time.png/PNG"
        )
        self.miriad.run_command(
            f"uvplt vis={scal} \
            stokes=i \
            axis=uc,vc \
            options=nofqav,nobase \
            device=diagnostics/{scal.name}_uv_coverage.png/PNG"
        )
        self.miriad.run_command(
            f"uvplt vis={scal} \
            stokes=i \
            axis=time,parang \
            options=nofqav,nobase \
            device=diagnostics/{scal.name}_parang_coverage.png/PNG"
        )
        self.miriad.run_command(
            f"gpplt vis={scal} \
            options=dots,wrap \
            yaxis=phase \
            yrange=-180,180 \
            device=diagnostics/phase_solutions.png/PNG"
        )
        self.miriad.run_command(
            f"gpplt vis={scal} \
            options=dots,xygains,wrap \
            yaxis=phase \
            yrange=-10,10 \
            device=diagnostics/xy-phase_solutions.png/PNG"
        )

        # Rename bandpass plots for XX and YY polarisations
        self.miriad.run_command(
            f"mv \
            diagnostics/{pcal.name}_bandpass.png \
            diagnostics/{pcal.name}_bandpass_xx.png"
        )
        self.miriad.run_command(
            f"mv \
            diagnostics/{pcal.name}_bandpass.png_2 \
            diagnostics/{pcal.name}_bandpass_yy.png"
        )

        # Remove all sub-band phase/xy-phase solution plots
        for plot in glob.glob(f"diagnostics/*png_*"):
            self.miriad.run_command(f"rm {plot}")

        os.chdir(cwd)

        return

    def run(self):

        primary_cal = self.miriad.target_paths.get("primary_cal").path
        gain_cal = self.miriad.target_paths.get("gain_cal").path
        target = self.miriad.target_paths.get("target").path
        leakage_cal = self.miriad.target_paths.get("leakage_cal").path

        if os.path.exists(f"{self.miriad.out_dir}/{target.name}.ms"):
            logger.warning(
                f"{self.miriad.out_dir}/{target.name}.ms already exists, will not overwrite"
            )
            return

        # Primary bandpass / flux calibrator
        # ---------------------------------
        self.miriad.bandpass(primary_cal)

        # Flag and solve for bandpass on primary calibrator
        if self.interactive:
            while prompt(f"Start interactive flagging of {primary_cal.name}?"):
                self.flag_sequence(primary_cal)
                self.miriad.bandpass(primary_cal)
        else:
            for _ in range(self.num_flag_rounds):
                self.miriad.autoflag(primary_cal)
                self.miriad.bandpass(primary_cal)

        # Solve for primary calibrator gains / leakage
        self.miriad.gaincal(primary_cal, options="xyvary")

        # Set options to work with strong or weakly polarised calibrator
        if self.miriad.strong_pol:
            gp_options = "xyvary,qusolve,vsolve,xyref,polref"

        else:
            gp_options = "xyvary,qusolve"

        # Leakage calibrator
        # ------------------
        if leakage_cal is not None:

            self.miriad.copycal(primary_cal, leakage_cal)

            # Flag and solve for gains / leakages / xy-phase on leakage calibrator
            if self.interactive:
                while prompt(f"Start interactive flagging of {leakage_cal.name}?"):
                    self.flag_sequence(leakage_cal)
                    self.miriad.gaincal(leakage_cal, options=gp_options)
            else:
                for _ in range(self.num_flag_rounds):
                    self.miriad.autoflag(leakage_cal)
                    self.miriad.gaincal(leakage_cal, options=gp_options)

            # To avoid corruption of Stokes V zero-point, we copy solutions
            # back to primary calibrator and repeat the sequence
            self.miriad.copycal(leakage_cal, primary_cal)
            self.miriad.gaincal(primary_cal, options="noxy")
            self.miriad.copycal(primary_cal, leakage_cal)
            self.miriad.gaincal(leakage_cal, options=gp_options)

            # Now we turn off xy-phase / leakage calibration for subsequent gain calibration
            self.miriad.copycal(leakage_cal, gain_cal)
            gp_options = "qusolve,noxy,nopol"
        else:
            if primary_cal != gain_cal:
                self.miriad.copycal(primary_cal, gain_cal)

        # Secondary gain calibrator
        # -------------------------
        if primary_cal != gain_cal:
            if self.interactive:
                while prompt(f"Start interactive flagging of {gain_cal.name}?"):
                    self.flag_sequence(gain_cal)
                    self.miriad.gaincal(gain_cal, options=gp_options)
            else:
                for _ in range(self.num_flag_rounds):
                    self.miriad.autoflag(gain_cal)
                    self.miriad.gaincal(gain_cal, options=gp_options)

            self.miriad.bootstrap(gain_cal, primary_cal)

        # Science target
        # --------------
        self.miriad.copycal(gain_cal, target)

        if not self.miriad.noflag:
            self.miriad.autoflag(target)

        # Average solutions and apply
        # --------------------------
        self.miriad.gpaver(target)
        self.miriad.uvaver(gain_cal)
        self.miriad.uvaver(target)

        self.miriad.generate_ms(f"{target}.cal")
        self.miriad.generate_ms(f"{gain_cal}.cal")
        self.miriad.run_command(
            f"mv {target}.cal.ms {self.miriad.out_dir}/{target.name}.ms"
        )

        return


if __name__ == "__main__":
    main()
