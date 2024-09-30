import pytest
from scipy.constants import find, value, c, speed_of_light, precision
from numpy.testing import assert_equal, assert_, assert_almost_equal
import scipy.constants._codata as _cd


def test_find():
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
    c_s = 'speed of light in vacuum'
    assert_equal(value(c_s), c)
    assert_equal(value(c_s), speed_of_light)


def test_basic_lookup():
    assert_equal('%d %s' % (_cd.value('speed of light in vacuum'),
                            _cd.unit('speed of light in vacuum')),
                 '299792458 m s^-1')

def test_find_all():
    assert_(len(find(disp=False)) > 300)


def test_find_single():
    assert_equal(find('Wien freq', disp=False)[0],
                 'Wien frequency displacement law constant')


def test_2002_vs_2006():
    assert_almost_equal(value('magn. flux quantum'),
                        value('mag. flux quantum'))


def test_exact_values():
    # Check that updating stored values with exact ones worked.
    exact = dict((k, v[0]) for k, v in _cd._physical_constants_2018.items())
    replace = _cd.exact2018(exact)
    for key, val in replace.items():
        assert_equal(val, value(key))
        assert precision(key) == 0


def test_trunc_not_marked_exact_value_2002to2014():
    exact = dict(
        (k, v[0]) for k, v in _cd._physical_constants_2002.items()
    )

    with pytest.raises(Warning):
        assert _cd.parse_constants_2002to2014(
            _cd.txt2002.replace("(exact)", "0000000"), exact
        )


def test_trunc_not_marked_exact_value_2018to2022():
    exact = dict(
        (k, v[0]) for k, v in _cd._physical_constants_2018.items()
    )

    with pytest.raises(Warning):
        assert _cd.parse_constants_2018toXXXX(
            _cd.txt2018.replace("(exact)", "0000000"), exact
        )


def test_not_listed_as_exact():
    exact = dict(
        (k, v[0]) for k, v in _cd._physical_constants_2018.items()
    )

    with pytest.raises(Warning):
        assert _cd.replace_exact(
            _cd.txt2018, {"fictitious constant"}, exact
        )


def test_not_correctly_calculated_constant():
    with pytest.raises(Warning):
        assert _cd.replace_exact(
            _cd._physical_constants_2002,
            {"magn. constant"},
            {"magn. constant": 0},
        )


def test_unmatched_exact_constants():
    exact = dict(
        (k, v[0]) for k, v in _cd._physical_constants_2018.items()
    )

    with pytest.raises(Warning):
        assert _cd.replace_exact(
            _cd._physical_constants_2018,
            {"empty set"},
            exact,
        )
