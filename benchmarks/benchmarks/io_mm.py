from .common import set_mem_rlimit, run_monitored, get_mem_info

from io import BytesIO, StringIO
import os
import tempfile

import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    import scipy.sparse
    import scipy.io._mmio
    import scipy.io._fast_matrix_market
    from scipy.io._fast_matrix_market import mmwrite


def generate_coo(size):
    nnz = int(size / (4 + 4 + 8))
    rows = np.arange(nnz, dtype=np.int32)
    cols = np.arange(nnz, dtype=np.int32)
    data = np.random.default_rng().uniform(low=0, high=1.0, size=nnz)
    return scipy.sparse.coo_matrix((data, (rows, cols)), shape=(nnz, nnz))


def generate_csr(size):
    nrows = 1000
    nnz = int((size - (nrows + 1) * 4) / (4 + 8))
    indptr = (np.arange(nrows + 1, dtype=np.float32) / nrows * nnz).astype(np.int32)
    indptr[-1] = nnz
    indices = np.arange(nnz, dtype=np.int32)
    data = np.random.default_rng().uniform(low=0, high=1.0, size=nnz)
    return scipy.sparse.csr_matrix((data, indices, indptr), shape=(nrows, nnz))


def generate_dense(size):
    nnz = size // 8
    return np.random.default_rng().uniform(low=0, high=1.0, size=(1, nnz))


class MemUsage(Benchmark):
    param_names = ['size', 'implementation', 'matrix_type']
    timeout = 4*60
    unit = "actual/optimal memory usage ratio"

    @property
    def params(self):
        return [
            list(self._get_size().keys()),
            ['scipy.io', 'scipy.io._mmio', 'scipy.io._fast_matrix_market'],
            ['dense', 'coo']  # + ['csr']
        ]

    def _get_size(self):
        size = {
            '1M': int(1e6),
            '10M': int(10e6),
            '100M': int(100e6),
            '300M': int(300e6),
            # '500M': int(500e6),
            # '1000M': int(1000e6),
        }
        return size

    def setup(self, size, implementation, matrix_type):
        set_mem_rlimit()
        self.size = self._get_size()
        size = self.size[size]

        mem_info = get_mem_info()
        try:
            mem_available = mem_info['memavailable']
        except KeyError:
            mem_available = mem_info['memtotal']

        max_size = int(mem_available * 0.7)//4

        if size > max_size:
            raise NotImplementedError()

        # Setup temp file
        f = tempfile.NamedTemporaryFile(delete=False, suffix='.mtx')
        f.close()
        self.filename = f.name

    def teardown(self, size, implementation, matrix_type):
        os.unlink(self.filename)

    def track_mmread(self, size, implementation, matrix_type):
        size = self.size[size]

        if matrix_type == 'coo':
            a = generate_coo(size)
        elif matrix_type == 'dense':
            a = generate_dense(size)
        elif matrix_type == 'csr':
            # cannot read directly into csr, only coo
            return 0
        else:
            raise NotImplementedError

        mmwrite(self.filename, a, symmetry='general')
        del a

        code = f"""
        from {implementation} import mmread
        mmread('{self.filename}')
        """
        time, peak_mem = run_monitored(code)
        return peak_mem / size

    def track_mmwrite(self, size, implementation, matrix_type):
        size = self.size[size]

        code = f"""
        import numpy as np
        import scipy.sparse
        from {implementation} import mmwrite
        
        def generate_coo(size):
            nnz = int(size / (4 + 4 + 8))
            rows = np.arange(nnz, dtype=np.int32)
            cols = np.arange(nnz, dtype=np.int32)
            data = np.random.default_rng().uniform(low=0, high=1.0, size=nnz)
            return scipy.sparse.coo_matrix((data, (rows, cols)), shape=(nnz, nnz))

        def generate_csr(size):
            nrows = 1000
            nnz = int((size - (nrows + 1) * 4) / (4 + 8))
            indptr = (np.arange(nrows + 1, dtype=np.float32) / nrows * nnz).astype(np.int32)
            indptr[-1] = nnz
            indices = np.arange(nnz, dtype=np.int32)
            data = np.random.default_rng().uniform(low=0, high=1.0, size=nnz)
            return scipy.sparse.csr_matrix((data, indices, indptr), shape=(nrows, nnz))
        
        def generate_dense(size):
            nnz = size // 8
            return np.random.default_rng().uniform(low=0, high=1.0, size=(1, nnz))


        a = generate_{matrix_type}({size})
        mmwrite('{self.filename}', a, symmetry='general')
        """
        time, peak_mem = run_monitored(code)
        return peak_mem / size


class IOSpeed(Benchmark):
    """
    Basic speed test. Does not show full potential as
    1) a relatively small matrix is used to keep test duration reasonable
    2) StringIO/BytesIO are noticeably slower than native C++ I/O to an SSD.
    """
    param_names = ['implementation', 'matrix_type']
    params = [
        ['scipy.io', 'scipy.io._mmio', 'scipy.io._fast_matrix_market'],
        ['dense', 'coo']  # + ['csr']
    ]

    def setup(self, implementation, matrix_type):
        # Use a 10MB matrix size to keep the runtimes somewhat short
        self.size = int(10e6)

        if matrix_type == 'coo':
            self.a = generate_coo(self.size)
        elif matrix_type == 'dense':
            self.a = generate_dense(self.size)
        elif matrix_type == 'csr':
            self.a = generate_csr(self.size)
        else:
            raise NotImplementedError

        bio = BytesIO()
        mmwrite(bio, self.a, symmetry='general')
        self.a_str = bio.getvalue().decode()

    def time_mmread(self, implementation, matrix_type):
        if matrix_type == 'csr':
            # cannot read directly into csr, only coo
            return

        if implementation == 'scipy.io':
            impl_module = scipy.io
        elif implementation == 'scipy.io._mmio':
            impl_module = scipy.io._mmio
        elif implementation == 'scipy.io._fast_matrix_market':
            impl_module = scipy.io._fast_matrix_market
        else:
            raise NotImplementedError

        impl_module.mmread(StringIO(self.a_str))

    def time_mmwrite(self, implementation, matrix_type):
        if implementation == 'scipy.io':
            impl_module = scipy.io
        elif implementation == 'scipy.io._mmio':
            impl_module = scipy.io._mmio
        elif implementation == 'scipy.io._fast_matrix_market':
            impl_module = scipy.io._fast_matrix_market
        else:
            raise NotImplementedError

        impl_module.mmwrite(BytesIO(), self.a, symmetry='general')
