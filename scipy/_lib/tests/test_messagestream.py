import pytest

from scipy._lib import messagestream
from scipy._lib.messagestream import MessageStream

# Number of streams kept for reuse (`_MAX_SPARE` in messagestream.pyx)
MAX_SPARE = 16

pytestmark = pytest.mark.thread_unsafe(
    reason="changes the global state of MessageStream")


@pytest.fixture(autouse=True)
def tempfile_streams():
    # Use temporary files, as on platforms without open_memstream, and start
    # with no streams kept for reuse.
    messagestream._close_spares()
    messagestream._force_tempfile(True)
    try:
        yield
    finally:
        messagestream._force_tempfile(False)
        messagestream._close_spares()


def test_stream_is_reused():
    stream = MessageStream()
    assert messagestream._spare_count() == 0
    stream.close()
    assert messagestream._spare_count() == 1

    stream = MessageStream()
    assert messagestream._spare_count() == 0
    stream.close()
    assert messagestream._spare_count() == 1


def test_dealloc_keeps_stream():
    stream = MessageStream()
    del stream
    assert messagestream._spare_count() == 1


def test_close_twice():
    stream = MessageStream()
    stream.close()
    stream.close()
    assert messagestream._spare_count() == 1


def test_number_of_kept_streams_is_limited():
    n = MAX_SPARE + 4
    for _ in range(2):
        streams = [MessageStream() for _ in range(n)]
        assert messagestream._spare_count() == 0
        for stream in streams:
            stream.close()
        assert messagestream._spare_count() == MAX_SPARE


def test_empty_stream():
    stream = MessageStream()
    assert stream.get() == ""
    stream.clear()
    assert stream.get() == ""
    stream.close()

    # also after reuse
    stream = MessageStream()
    assert stream.get() == ""
    stream.close()
