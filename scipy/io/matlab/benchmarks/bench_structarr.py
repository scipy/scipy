from __future__ import division
from numpy.testing import *

from cStringIO import StringIO

import numpy as np
import scipy.io as sio


def make_structarr(n_vars, n_fields, n_structs):
    var_dict = {}
    for vno in range(n_vars):
        vname = 'var%00d' % vno
        end_dtype = [('f%d' % d, 'i4', 10) for d in range(n_fields)]
        s_arrs = np.zeros((n_structs,), dtype=end_dtype)
        var_dict[vname] = s_arrs
    return var_dict


def bench_run(): 
    str_io = StringIO()
    print
    print 'Read / writing matlab structs'
    print '='*60
    print ' write |  read |   vars | fields | structs '
    print '-'*60
    print
    for n_vars, n_fields, n_structs in (
        (10, 10, 20),):
        var_dict = make_structarr(n_vars, n_fields, n_structs)
        str_io = StringIO()
        write_time = measure('sio.savemat(str_io, var_dict)')
        read_time = measure('sio.loadmat(str_io)')
        print '%.5f | %.5f | %5d | %5d | %5d ' % (
            write_time,
            read_time,
            n_vars,
            n_fields,
            n_structs)


if __name__ == '__main__' :
    bench_run()
