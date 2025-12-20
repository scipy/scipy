import io
from scipy.sparse import csr_matrix
from scipy.io import hb_write, hb_read


def test_hb_write_empty_matrix_roundtrip():
    m = csr_matrix([[0.0, 0.0], [0.0, 0.0]])   # shape (2,2) and nnz == 0
    buf = io.StringIO()
    hb_write(buf, m)          # should NOT raise
    buf.seek(0)
    out = hb_read(buf)
    assert out.shape == (2, 2)
    assert out.nnz == 0


test_hb_write_empty_matrix_roundtrip()
