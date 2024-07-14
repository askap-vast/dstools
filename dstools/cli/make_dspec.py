import click
import h5py
import numpy as np
from casacore.tables import table
import itertools as it


@click.command()
@click.option(
    "-F",
    "--noflag",
    is_flag=True,
    default=False,
    help="Remove flagging mask.",
)
@click.option(
    "-d",
    "--datacolumn",
    type=click.Choice(["data", "corrected"]),
    default="data",
    help="Selection of DATA or CORRECTED_DATA column.",
)
@click.argument("ms")
@click.argument("outfile")
def main(noflag, datacolumn, ms, outfile):

    with h5py.File(outfile, "w") as f:

        t = table(ms)

        # Throw away autocorrelations
        t = t.query("ANTENNA1 != ANTENNA2")

        # Get time, frequency, and baseline axes
        tf = table(f"{ms}/SPECTRAL_WINDOW")
        times = np.unique(t.getcol("TIME"))
        freqs = tf[0]["CHAN_FREQ"]
        antennas = np.unique(
            np.append(
                t.getcol("ANTENNA1"),
                t.getcol("ANTENNA2"),
            ),
        )
        nbaselines = len(antennas) * (len(antennas) - 1) // 2

        # Calculate UVrange for each baseline
        uvw = t.getcol("UVW")
        uvdist = np.sqrt(np.sum(np.square(uvw), axis=1))
        uvdist = np.mean(uvdist.reshape(-1, nbaselines), axis=0)

        # Write time and frequency arrays, discarding duplicate timesamples
        datacolumn = "CORRECTED_DATA" if datacolumn == "corrected" else "DATA"
        bl_waterfall = t.getcol(datacolumn)

        # Apply flags
        if not noflag:
            vis_flag = t.getcol("FLAG")
            bl_waterfall[vis_flag] = np.nan

        # Reshape data into <NPOL,NFREQ,NTIME,NBL> shape array
        waterfall = bl_waterfall.T.reshape((4, len(freqs), -1, nbaselines)).T

        # Trim array to start/end time of data
        timeaxis = np.nansum(waterfall, axis=(0, 2, 3))
        timeaxis[timeaxis == 0 + 0j] = np.nan
        nonnan_vals = np.isfinite(timeaxis)

        mintime = np.argmax(nonnan_vals)
        maxtime = nonnan_vals.size - np.argmax(nonnan_vals[::-1])

        waterfall = waterfall[:, mintime:maxtime, :, :]
        times = times[mintime:maxtime]

        # Write all data to file
        f.create_dataset("time", data=times)
        f.create_dataset("frequency", data=freqs)
        f.create_dataset("uvdist", data=uvdist)
        f.create_dataset("flux", data=waterfall)

        t.close()


if __name__ == "__main__":
    main()
