
from numpy.testing import *

from cStringIO import StringIO
from zlib import compress, decompress

from numpy import linspace, inf
from scipy.io.matlab.zlibstreams import ZlibInputStream, \
     TwoShotZlibInputStream

def print_header(klass):
    print
    print '    %s reading gzip streams' % klass
    print '='*40
    print 'time(s) |  nbytes'
    print '-'*40
    print

def bench_run():
    test_sizes = [1e6, 5e6]
    for klass in (ZlibInputStream, TwoShotZlibInputStream):
        print_header(klass)
        for size in test_sizes:
            data = linspace(1, 100, num=size).astype('float32')
            zdata = compress(data.tostring())
            zstream = klass(StringIO(zdata), len(zdata))
            read_time = measure('zstream.read()')
            raw_read_time = measure('decompress(zdata)')
            if raw_read_time > 0:
                ratio = read_time / raw_read_time
            else:
                ratio = inf
            print '%.3f | %.3f | %d' % (read_time,
                                        read_time / raw_read_time,
                                        data.nbytes)

if __name__ == '__main__' :
    run_module_suite()
