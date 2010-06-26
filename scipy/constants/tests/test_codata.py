
import warnings

import codata
import constants

from scipy.constants import find
from numpy.testing import assert_equal, run_module_suite

def test_find():

    warnings.simplefilter('ignore', DeprecationWarning)

    keys = find('weak mixing', disp=False)
    assert_equal(keys, ['weak mixing angle'])

    keys = find('qwertyuiop', disp=False)
    assert_equal(keys, [])

    keys = find('natural unit', disp=False)
    assert_equal(keys, sorted(['natural unit of velocity',
                                'natural unit of action',
                                'natural unit of action in eV s',
                                'natural unit of mass',
                                'natural unit of energy',
                                'natural unit of energy in MeV',
                                'natural unit of momentum',
                                'natural unit of momentum in MeV/c',
                                'natural unit of length',
                                'natural unit of time']))


def test_basic_table_parse():
    c = 'speed of light in vacuum'
    assert_equal(codata.value(c), constants.c)
    assert_equal(codata.value(c), constants.speed_of_light)

def test_basic_lookup():
    assert_equal('%d %s' % (codata.c, codata.unit('speed of light in vacuum')),
                 '299792458 m s^-1')

if __name__ == "__main__":
    run_module_suite()
