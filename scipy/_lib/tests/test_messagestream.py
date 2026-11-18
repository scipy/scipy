import pytest

from scipy._lib import messagestream
from scipy._lib.messagestream import MessageStream

MAX_SPARE = messagestream._fhcache_info()["max_spare"]


@pytest.fixture
def tempfile_streams():
    # Use temporary files, as on platforms without open_memstream, and start
    # with no streams kept for reuse.
    messagestream._fhcache_clear()
    messagestream._force_tempfile(True)
    try:
        yield
    finally:
        messagestream._force_tempfile(False)
        messagestream._fhcache_clear()


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
def test_stream_is_reused():
    stream = MessageStream()
    assert messagestream._fhcache_size() == 0
    stream.close()
    assert messagestream._fhcache_size() == 1

    stream = MessageStream()
    assert messagestream._fhcache_size() == 0
    stream.close()
    assert messagestream._fhcache_size() == 1


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
def test_dealloc_keeps_stream():
    stream = MessageStream()
    del stream
    assert messagestream._fhcache_size() == 1


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
def test_close_twice():
    stream = MessageStream()
    stream.close()
    stream.close()
    assert messagestream._fhcache_size() == 1


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
@pytest.mark.parametrize("tempfile", [True, False])
def test_closed_stream(tempfile):
    messagestream._force_tempfile(tempfile)
    stream = MessageStream()
    stream.close()
    with pytest.raises(ValueError, match="closed"):
        stream.get()
    with pytest.raises(ValueError, match="closed"):
        stream.clear()


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
def test_number_of_kept_streams_is_limited():
    n = MAX_SPARE + 4
    for _ in range(2):
        streams = [MessageStream() for _ in range(n)]
        assert messagestream._fhcache_size() == 0
        for stream in streams:
            stream.close()
        assert messagestream._fhcache_size() == MAX_SPARE


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
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


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
def test_fhcache_info():
    start = messagestream._fhcache_info()
    assert start["n_spare"] == 0

    stream = MessageStream()  # no spare: a new file
    info = messagestream._fhcache_info()
    assert info["misses"] == start["misses"] + 1
    assert info["hits"] == start["hits"]
    assert info["n_messagestream"] == start["n_messagestream"] + 1
    assert info["n_tempfile"] == start["n_tempfile"] + 1

    stream.close()  # kept as a spare, still open
    info = messagestream._fhcache_info()
    assert info["n_spare"] == 1
    assert info["n_messagestream"] == start["n_messagestream"]
    assert info["n_tempfile"] == start["n_tempfile"] + 1

    stream = MessageStream()  # reuses the spare
    info = messagestream._fhcache_info()
    assert info["hits"] == start["hits"] + 1
    assert info["misses"] == start["misses"] + 1
    assert info["n_spare"] == 0
    stream.close()

    messagestream._fhcache_clear()
    info = messagestream._fhcache_info()
    assert info["n_spare"] == 0
    assert info["n_messagestream"] == start["n_messagestream"]
    assert info["n_tempfile"] == start["n_tempfile"]
    assert info["close_errors"] == start["close_errors"]


@pytest.mark.thread_unsafe(reason="changes the global state of MessageStream")
@pytest.mark.usefixtures("tempfile_streams")
def test_fhcache_info_no_lost_files():
    # with more coexisting streams than are kept, after all are closed only
    # the spares remain open
    start = messagestream._fhcache_info()
    for _ in range(3):
        streams = [MessageStream() for _ in range(MAX_SPARE + 4)]
        for stream in streams:
            stream.close()
    info = messagestream._fhcache_info()
    assert info["n_spare"] == MAX_SPARE
    assert info["n_messagestream"] == start["n_messagestream"]
    assert (info["n_tempfile"] - info["n_spare"]
            == start["n_tempfile"] - start["n_spare"])


@pytest.mark.thread_unsafe(reason="checks the global state of MessageStream")
def test_fhcache_info_counts_all_streams():
    # in memory where open_memstream exists, or a temporary file
    start = messagestream._fhcache_info()["n_messagestream"]
    stream = MessageStream()
    assert messagestream._fhcache_info()["n_messagestream"] == start + 1
    stream.close()
    assert messagestream._fhcache_info()["n_messagestream"] == start
