#!/usr/bin/env python 

import numpy as np
np.seterr(all='raise')

from numpy import random 

from scipy.linalg.expokit import expm as expm_expokit

#from scipy.linalg._expokit import dgpadm, zgpadm, dspadm, zhpadm, \
#                                  dgexpv, zgexpv

from scipy.linalg import expm as pyexpm
from scipy.sparse.base import isspmatrix

import logging
import unittest
from timeit import Timer


REPEAT = 1000   #: Number of times to repeat TimeExpokit test

MINARRSIZE  = 5
MAXARRSIZE = 50
STRIDE = 5
PRECISION = 10

class CompareAccuracy(unittest.TestCase):

    def compare_precision(self, array, precision):
        """Call :func:`scipy.linalg.expm` and :func:`scipy.linalg.expokit.expm`.

        Compare results to specified precision, and raise
        :class:`AssertionError` if there are differences.
        """
        padm = expm_expokit(array)
        pypd = pyexpm(array)
        #test = (padm == pypd).all()
        test = np.array_equal(padm, pypd)
        if test is False:
            logging.warn("np.array_equal says matrices differ!")
            #print('EXPOKIT result:\n {0}'.format(padm))
            #print('SciPy   result:\n {0}'.format(pypd))
            assert padm.shape == pypd.shape, \
                    'Shapes differ: EXPOKIT: {0}, SciPy: {1}'\
                    .format( padm.shape, pypd.shape )
            assert np.allclose(padm, pypd), "Matrices are not even close!"
            for xi, x in enumerate(pypd):
                for yi, xy in enumerate(x):
                    if '{0:.{1}f}'.format(padm[xi][yi], precision) != \
                            '{0:.{1}f}'.format(pypd[xi][yi], precision):
                        print('At {xy}:\n' \
                              '    expokit:     {expoexpm},\n'
                              '    linalg.expm: {pyexpm}'.format( \
                                prec=precision,
                                xy=(xi, yi),
                                expoexpm=padm[xi][yi],
                                pyexpm=pypd[xi][yi]) )
        logging.info('arrays equal to {0} decimal places :)'.format(precision))

    def test_small(self):
        array, carray = make_small_arrays(MINARRSIZE)
        logging.info("DGPADM small test on dtype: {0}".format(array.dtype))
        self.compare_precision(array, precision=PRECISION)
        logging.info("ZGPADM small test on dtype: {0}".format(carray.dtype))
        self.compare_precision(carray, precision=PRECISION)

    def test_dense(self):
        array, carray = make_arrays(MINARRSIZE)
        logging.info("DGPADM test on dtype: {0}".format(array.dtype))
        self.compare_precision(array, precision=PRECISION)
        logging.info("ZGPADM test on dtype: {0}".format(carray.dtype))
        self.compare_precision(carray, precision=PRECISION)


    def test_sparse(self):
        array, carray = make_sparse_arrays()
        logging.info("DGEXPV sparse test on dtype: {0}".format(array.dtype))
        self.compare_precision(array, precision=PRECISION)
        logging.info("ZGEXPV sparse test on dtype: {0}".format(carray.dtype))
        self.compare_precision(carray, precision=PRECISION)


TEST_ARRAY = np.array(())
TEST_CARRAY = np.array(())

class TimeExpokit(unittest.TestCase):

    def test_times(self):
        global TEST_ARRAY, TEST_CARRAY
        print("Times averaged over {0} iterations".format(REPEAT))
        print("Array size".ljust(25) +
              "EXPOKIT times".ljust(25) +
              "SciPy times".ljust(25) )
        for arrsize in xrange(MINARRSIZE, MAXARRSIZE + 1, STRIDE):
            TEST_ARRAY, TEST_CARRAY = make_arrays(arrsize)
            print('{0}'.format(arrsize).ljust(25)),
            self.time_expokit()
            self.time_scipy()

    def time_expokit(self):
        timer = Timer("expm_expokit(TEST_ARRAY)",
                setup="from __main__ import expm_expokit, TEST_ARRAY")
        times = timer.repeat(number=REPEAT)
        print('{0:.03f} += {1:.03f}'.format(np.mean(times), np.std(times)).ljust(25)),

    def time_scipy(self):
        timer = Timer("pyexpm(TEST_ARRAY)",
                setup="from __main__ import pyexpm, TEST_ARRAY")
        times = timer.repeat(number=REPEAT)
        print('{0:.03f} += {1:.03f}'.format(np.mean(times), np.std(times) ).ljust(25))

def make_arrays(arrsize):
    """Return a real and complex array, using host's default data-types."""
    array = random.rand(arrsize,arrsize)
    carray = np.dot(array, np.eye(arrsize)*np.complex(1+1j))
    return array, carray

def make_small_arrays(arrsize):
    """Make a real and complex array, using data-types used by default 
    on 32bit machines."""
    array = np.array( random.rand(arrsize,arrsize), dtype = np.float32 )
    carray = np.array( np.dot( array, np.eye(arrsize)*np.complex(1+1j) ), dtype=np.complex64 )
    return array, carray

def make_sparse_arrays(arrsize=50):
    from scipy.sparse import diags
    # real sparse matrix
    ddiags = [ random.rand( arrsize ),
               random.rand(arrsize-1),
               random.rand(arrsize-1) ]
    offsets = np.array([0, 1, -1])
    array = diags( ddiags, offsets)

    # complex sparse matrix
    cdata  = [ np.dot( ddiags[0],  np.eye(arrsize)   * np.complex(1+1j) ),
               np.dot( ddiags[1],  np.eye(arrsize-1) * np.complex(1+1j) ),
               np.dot( ddiags[1],  np.eye(arrsize-1) * np.complex(1+1j) ) ]
    carray = diags( cdata, offsets )
    return array, carray


if __name__ == '__main__':
    unittest.main()
    exit()

