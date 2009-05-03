from __future__ import division
from numpy.testing import *

from math import ceil
from StringIO import StringIO
from zlib import compress, decompress

from numpy import linspace, inf
from scipy.io.matlab.zlibstreams import ZlibInputStream, \
     StubbyZlibInputStream

def print_header(klass):
    print
    print '    %s reading gzip streams' % klass
    print '='*40
    print 'time(s) | ratio | nbytes | block size'
    print '-'*40
    print


def raw_read(block_size, zstr, n):
    sstream = StringIO(decompress(zstr))
    reads = int(ceil(n / block_size))
    for i in range(reads):
        sstream.read(block_size)


def zstream_read(block_size, zstr, n, klass):
    zstream = klass(StringIO(zstr), len(zstr))
    reads = int(ceil(n / block_size))
    for i in range(reads):
        zstream.read(block_size)
    

def bench_run():
    test_sizes = [1e5]
    for klass in (ZlibInputStream, StubbyZlibInputStream):
        print_header(klass)
        for size in test_sizes:
            data = linspace(1, 100, num=size).tostring()
            n = len(data)
            zdata = compress(data)
            for block_size in (1, 10, 100, 1000):
                read_time = measure(
                    'zstream_read(block_size, zdata, n, klass)')
                raw_read_time = measure(
                    'raw_read(block_size, zdata, n)')
                if raw_read_time > 0:
                    ratio = read_time / raw_read_time
                else:
                    ratio = inf
                print '%.3f | %4.2f | %d | %d ' % (
                    read_time,
                    ratio,
                    n,
                    block_size)

if __name__ == '__main__' :
    run_module_suite()
