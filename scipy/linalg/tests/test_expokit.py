#!/usr/bin/env python 
#
# Created by: Alex Leach, November 2012
#

import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_array_almost_equal, \
    assert_array_almost_equal_nulp

from scipy.linalg.expowrap import expm

from scipy.linalg.expokit import dgchbv, dgexpv, dgpadm, dgphiv, dmexpv, dnchbv, \
    dschbv, dsexpv, dspadm, dsphiv, zgchbv, zgexpv, zgpadm, zgphiv, zhexpv,  \
    zhpadm, zhphiv, znchbv


class TestPADM(TestCase):

    def test_zero(self):
        a = array([[0.,0],[0,0]])
        assert_array_almost_equal(expm(a),[[1,0],[0,1]])

    def test_simple(self):
        a = array([[0.,1.],[0,0]])
        assert_array_almost_equal(expm(a),[[1,1],[0,1]] )


def time_dgpadm():
    return expmwrap( array )

def time_pyexpm():
    return pyexpm( array )

def sanity_test(array):
    padm = expmwrap( array )
    pypd = pyexpm( array )
    test = np.all(padm == pypd)
    if test is False:
        print "Matrices differ!\n\nexpokit:\n{0}\n\nexpm:\n{1}\n\n{2}"\
        .format( padm, pypd, test )
        assert padm.shape == pypd.shape
        for xi,x in enumerate(pypd):
            for yi,y in enumerate(x):
                assert '{0:.10f}'.format(padm[xi][yi]) == '{0:.10f}'.format(y), \
                        'expokit:{0}, expm: {1} at ({2},{3})'.format(padm[xi][yi], y, xi, yi)

def make_arrays(arrsize):
    array = random.rand(arrsize,arrsize)
    carray = np.dot( array, np.eye(arrsize)*np.complex(1+1j) )
    return array, carray

if __name__ == '__main__':
    print '{0} {1}'.format( ''.ljust(25), '----- Mean time +- S.D. -----'.rjust(37) )
    print '{0} {1} {2}'.format( 'Array size'.ljust(25), \
            'Expokit'.ljust(25), 'scipy..expm'.ljust(25) )
    array, carray = make_arrays(minarrsize)
    sanity_test(array)
    sanity_test(carray)
    for arrsize in xrange( minarrsize, maxarrsize+1, stride ):
        array, carray = make_arrays(arrsize)
        print '{0}'.format(arrsize).ljust(25),

        timer = timeit.Timer("time_dgpadm()", setup="from __main__ import time_dgpadm, np, array")
        times = timer.repeat(number=1000)
        print '{0:.03f} += {1:.03f}'.format( np.mean(times), np.std(times) ).ljust(25),

        timer = timeit.Timer("time_pyexpm()", setup="from __main__ import time_pyexpm, np, array")
        times = timer.repeat(number = 1000)
        print '{0:.03f} += {1:.03f}'.format( np.mean(times), np.std(times) ).ljust(25)

