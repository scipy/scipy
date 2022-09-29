from scipy.datasets._registry import registry
from scipy.datasets._fetchers import fetch_data, data_fetcher
from scipy.datasets import ascent, face, electrocardiogram
from numpy.testing import assert_equal, assert_almost_equal, suppress_warnings
import os
import pytest

try:
    # https://github.com/scipy/scipy/pull/15607#issuecomment-1176457275
    # TODO: Remove warning filter after next certifi release
    with suppress_warnings() as sup:
        sup.filter(category=DeprecationWarning)
        import pooch
except ImportError:
    raise ImportError("Missing optional dependency 'pooch' required "
                      "for scipy.datasets module. Please use pip or "
                      "conda to install 'pooch'.")


data_dir = data_fetcher.path  # type: ignore


def _has_hash(path, expected_hash):
    """Check if the provided path has the expected hash."""
    if not os.path.exists(path):
        return False
    return pooch.file_hash(path) == expected_hash


class TestDatasets:

    @pytest.fixture(scope='module', autouse=True)
    def test_download_all(self):
        # This fixture requires INTERNET CONNECTION

        # test_setup phase
        for dataset in registry:
            fetch_data(dataset)

        yield

        # test_teardown phase

        # TODO: Clear all cache from test_download_all
        # and then test downloading all the datasets
        # individually through test_<dataset-name> methods instead

        # See case 2 for more details in the following gh-comment
        # https://github.com/scipy/scipy/pull/15607#discussion_r943557161

    def test_existence_all(self):
        assert len(os.listdir(data_dir)) >= len(registry)

    def test_ascent(self):
        assert_equal(ascent().shape, (512, 512))

        # hash check
        assert _has_hash(os.path.join(data_dir, "ascent.dat"),
                         registry["ascent.dat"])

    def test_face(self):
        assert_equal(face().shape, (768, 1024, 3))

        # hash check
        assert _has_hash(os.path.join(data_dir, "face.dat"),
                         registry["face.dat"])

    def test_electrocardiogram(self):
        # Test shape, dtype and stats of signal
        ecg = electrocardiogram()
        assert_equal(ecg.dtype, float)
        assert_equal(ecg.shape, (108000,))
        assert_almost_equal(ecg.mean(), -0.16510875)
        assert_almost_equal(ecg.std(), 0.5992473991177294)

        # hash check
        assert _has_hash(os.path.join(data_dir, "ecg.dat"),
                         registry["ecg.dat"])
