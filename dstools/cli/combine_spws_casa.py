import click


@click.command()
@click.argument("ms")
@click.argument("outputvis")
def main(ms, outputvis):

    # Determine number of spectral windows
    tab = tbtool()
    tab.open(f"{ms}::SPECTRAL_WINDOW")

    nspws = len(tab.getcol("FREQ_GROUP_NAME"))

    tab.unlock()
    tab.close()

    # Combine spectral windows if more than 1
    combine = nspws > 1
    mstransform(
        vis=ms,
        combinespws=combine,
        datacolumn="all",
        outputvis=outputvis,
    )

    return


if __name__ == "__main__":
    main()
