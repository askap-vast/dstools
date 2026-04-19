import pytest

pytest.importorskip("casacore.tables")
pytestmark = pytest.mark.fullstack

from dstools.casa import flagdata  # noqa: E402


def test_casacore_casatasks_bindings_dont_clash(ms, mocker):
    flagstats = flagdata(vis=ms.path.as_posix(), mode="summary")

    assert "antenna" in flagstats.keys()
