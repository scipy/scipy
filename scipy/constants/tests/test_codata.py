
import warnings

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

if __name__ == "__main__":
    run_module_suite()
