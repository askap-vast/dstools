import click
import numpy as np


@click.command()
@click.option(
    "-u",
    "--uvrange",
    type=str,
    default="0m",
    help="Minimum uv distance to select from data column (e.g. -u 200m).",
)
@click.option(
    "-c",
    "--datacolumn",
    type=click.Choice(["data", "corrected"]),
    default="data",
    help="Selection of DATA or CORRECTED_DATA column.",
)
@click.argument("msname")
def main(msname, datacolumn, uvrange):
    uvrange = ">{}".format(uvrange)
    savename = msname.replace(".ms", ".baseavg.ms")

    intab = tbtool()
    intab.open(msname, nomodify=False)

    ant1 = intab.getcol("ANTENNA1")
    ant2 = intab.getcol("ANTENNA2")

    # Set all antenna pairs equal for baseline averaging
    nrows = intab.nrows()
    intab.putcol("ANTENNA1", np.zeros(nrows))
    intab.putcol("ANTENNA2", np.ones(nrows))

    # Average over baselines by setting timeaverage interval to less than one scan cycle
    interval = intab.getcol("INTERVAL")
    timebin = "{}s".format(min(interval) * 1e-2)

    mstransform(
        vis=msname,
        outputvis=savename,
        datacolumn=datacolumn,
        uvrange=uvrange,
        timeaverage=True,
        timebin=timebin,
        keepflags=False,
    )

    intab.putcol("ANTENNA1", ant1)
    intab.putcol("ANTENNA2", ant2)

    intab.unlock()
    intab.close()


if __name__ == "__main__":
    main()
