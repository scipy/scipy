import pooch
from scipy.datasets._registry import registry
from scipy.datasets._fetchers import data
from scipy.datasets import ascent, face, electrocardiogram
from numpy.testing import assert_equal, assert_almost_equal
import os


def _has_hash(path, expected_hash):
    """Check if the provided path has the expected hash."""
    if not os.path.exists(path):
        return False
    return pooch.file_hash(path) == expected_hash


def test_download_all():
    # This test requires INTERNET CONNECTION
    data_dir = data.path
    for dataset in registry:
        data.fetch(dataset)

    assert len(os.listdir(data_dir)) >= len(registry)


def test_ascent():
    assert_equal(ascent().shape, (512, 512))

    # hash check
    assert _has_hash(os.path.join(data.path, "ascent.dat") , registry["ascent.dat"])


def test_face():
    assert_equal(face().shape, (768, 1024, 3))

    # hash check
    assert _has_hash(os.path.join(data.path, "face.dat") , registry["face.dat"])


def test_electrocardiogram():
    # Test shape, dtype and stats of signal
    ecg = electrocardiogram()
    assert_equal(ecg.dtype, float)
    assert_equal(ecg.shape, (108000,))
    assert_almost_equal(ecg.mean(), -0.16510875)
    assert_almost_equal(ecg.std(), 0.5992473991177294)

    # hash check
    assert _has_hash(os.path.join(data.path, "ecg.dat") , registry["ecg.dat"])
