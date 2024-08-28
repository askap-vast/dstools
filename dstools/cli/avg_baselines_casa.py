import click
import numpy as np


@click.command()
@click.option(
    "-u",
    "--minuvdist",
    type=float,
    default=0,
    help="Minimum UV distance in meters to retain if averaging over baseline axis.",
)
@click.argument("msname")
def main(msname, minuvdist):
    savename = msname.replace(".ms", ".dstools-temp.baseavg.ms")

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
        datacolumn="all",
        uvrange=f">{minuvdist}m",
        timeaverage=True,
        timebin=timebin,
        keepflags=False,
    )

    # Replace original antenna names
    intab.putcol("ANTENNA1", ant1)
    intab.putcol("ANTENNA2", ant2)

    intab.unlock()
    intab.close()

    return


if __name__ == "__main__":
    main()
