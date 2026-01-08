""" Test tmpdirs module """
import os
import pytest
from scipy._lib._tmpdirs import in_tempdir


@pytest.mark.thread_unsafe(reason="in_tempdir is not thread-safe")
def test_in_tempdir():
    my_cwd = os.getcwd()
    with in_tempdir() as tmpdir:
        with open('test.txt', "w") as f:
            f.write('some text')
        assert os.path.isfile('test.txt')
        assert os.path.isfile(os.path.join(tmpdir, 'test.txt'))
    assert not os.path.exists(tmpdir)
    assert os.getcwd() == my_cwd
