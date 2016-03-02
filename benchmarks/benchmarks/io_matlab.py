from __future__ import division, absolute_import, print_function
from .common import set_mem_rlimit, run_monitored, get_mem_info

import os
import tempfile
import collections
from io import BytesIO

import numpy as np

try:
    from scipy.io import savemat, loadmat
except ImportError:
    pass

from .common import Benchmark


class MemUsage(Benchmark):
    param_names = ['size', 'compressed']
    timeout = 4*60
    unit = "actual/optimal memory usage ratio"

    @property
    def params(self):
        return [self._get_sizes().keys(), [True, False]]

    def _get_sizes(self):
        sizes = collections.OrderedDict([
            ('1M', 1e6),
            ('10M', 10e6),
            ('100M', 100e6),
            ('300M', 300e6),
            # ('500M', 500e6),
            # ('1000M', 1000e6),
        ])
        return sizes

    def setup(self, size, compressed):
        set_mem_rlimit()
        self.sizes = self._get_sizes()
        size = int(self.sizes[size])

        mem_info = get_mem_info()
        try:
            mem_available = mem_info['memavailable']
        except KeyError:
            mem_available = mem_info['memtotal']

        max_size = int(mem_available * 0.7)//4

        if size > max_size:
            raise NotImplementedError()

        # Setup temp file
        f = tempfile.NamedTemporaryFile(delete=False, suffix='.mat')
        f.close()
        self.filename = f.name

    def teardown(self, size, compressed):
        os.unlink(self.filename)

    def track_loadmat(self, size, compressed):
        size = int(self.sizes[size])

        x = np.random.rand(size//8).view(dtype=np.uint8)
        savemat(self.filename, dict(x=x), do_compression=compressed, oned_as='row')
        del x

        code = """
        from scipy.io import loadmat
        loadmat('%s')
        """ % (self.filename,)
        time, peak_mem = run_monitored(code)

        return peak_mem / size

    def track_savemat(self, size, compressed):
        size = int(self.sizes[size])

        code = """
        import numpy as np
        from scipy.io import savemat
        x = np.random.rand(%d//8).view(dtype=np.uint8)
        savemat('%s', dict(x=x), do_compression=%r, oned_as='row')
        """ % (size, self.filename, compressed)
        time, peak_mem = run_monitored(code)
        return peak_mem / size


class StructArr(Benchmark):
    params = [
        [(10, 10, 20), (20, 20, 40), (30, 30, 50)],
        [False, True]
    ]
    param_names = ['(vars, fields, structs)', 'compression']

    @staticmethod
    def make_structarr(n_vars, n_fields, n_structs):
        var_dict = {}
        for vno in range(n_vars):
            vname = 'var%00d' % vno
            end_dtype = [('f%d' % d, 'i4', 10) for d in range(n_fields)]
            s_arrs = np.zeros((n_structs,), dtype=end_dtype)
            var_dict[vname] = s_arrs
        return var_dict

    def setup(self, nvfs, compression):
        n_vars, n_fields, n_structs = nvfs

        self.var_dict = StructArr.make_structarr(n_vars, n_fields, n_structs)
        self.str_io = BytesIO()

        savemat(self.str_io, self.var_dict, do_compression=compression)

    def time_savemat(self, nvfs, compression):
        savemat(self.str_io, self.var_dict, do_compression=compression)

    def time_loadmat(self, nvfs, compression):
        loadmat(self.str_io)
