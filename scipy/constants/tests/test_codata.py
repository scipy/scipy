from numpy.testing import (assert_equal, assert_, assert_almost_equal)
from scipy.constants import _codata


def test_find():
    keys = _codata.find('weak mixing', disp=False)
    assert_equal(keys, ['weak mixing angle'])

    keys = _codata.find('qwertyuiop', disp=False)
    assert_equal(keys, [])

    keys = _codata.find('natural unit', disp=False)
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


def test_basic_lookup():
    assert_equal('%d %s' % (_codata.value('speed of light in vacuum'),
                            _codata.unit('speed of light in vacuum')),
                 '299792458 m s^-1')


def test_find_all():
    assert_(len(_codata.find(disp=False)) > 300)


def test_find_single():
    assert_equal(_codata.find('Wien freq', disp=False)[0],
                 'Wien frequency displacement law constant')


def test_2002_vs_2006():
    assert_almost_equal(_codata.value('magn. flux quantum'),
                        _codata.value('mag. flux quantum'))


def test_exact_values():
    # Check that updating stored values with exact ones worked.
    exact = dict((k, v[0]) for k, v in _codata._physical_constants_2018.items())
    replace = _codata.exact2018(exact)
    for key, val in _codata.replace.items():
        assert_equal(val, _codata.value(key))
        assert_(precision(key) == 0)
