import pytest

pytest.importorskip("casacore.tables")
pytestmark = pytest.mark.fullstack

from dstools.ms import CalTable  # noqa: E402


@pytest.mark.parametrize(("prop", "val"), [("nspws", 16), ("npols", 1)])
def test_caltable_dimensions(prop, val, copy_paths_into_workspace, caltable_sources):
    caltable_path = copy_paths_into_workspace(caltable_sources, ("fred",))["fred"]
    caltable = CalTable(path=caltable_path)

    assert getattr(caltable, prop) == val


@pytest.mark.parametrize("calmode", ["a", "p"])
def test_caltable_plot_solutions(calmode, copy_paths_into_workspace, caltable_sources):
    caltable_path = copy_paths_into_workspace(caltable_sources, ("fred",))["fred"]
    caltable = CalTable(path=caltable_path)

    fig = caltable.plot_solutions(calmode=calmode)

    assert len(fig.axes) == 6


def test_caltable_plot_solutions_invalid_calmode_raises_error(
    copy_paths_into_workspace,
    caltable_sources,
):
    caltable_path = copy_paths_into_workspace(caltable_sources, ("fred",))["fred"]
    caltable = CalTable(path=caltable_path)

    with pytest.raises(ValueError):
        caltable.plot_solutions(calmode="x")
