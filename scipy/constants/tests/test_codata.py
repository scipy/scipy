from scipy.constants import codata, find, value
from numpy.testing import (assert_equal, assert_,
                           assert_almost_equal)


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


def test_defining_constants():
    assert_equal('%d %s' %
                 (codata.value('hyperfine transition frequency of Cs-133'),
                  codata.unit('hyperfine transition frequency of Cs-133')),
                 '9192631770 Hz')
    assert_equal('%d %s' % (codata.value('speed of light in vacuum'),
                            codata.unit('speed of light in vacuum')),
                 '299792458 m s^-1')
    assert_equal('%f %s' % (codata.value('Planck constant'),
                            codata.unit('Planck constant')),
                 '6.62607015e-34 J Hz^-1')
    assert_equal('%f %s' % (codata.value('elementary charge'),
                            codata.unit('elementary charge')),
                 '1.602176634e-19 C')
    assert_equal('%f %s' % (codata.value('Boltzmann constant'),
                            codata.unit('Boltzmann constant')),
                 '1.380649e-23 J K^-1')
    assert_equal('%d %s' % (codata.value('Avogadro constant'),
                            codata.unit('Avogadro constant')),
                 '6.02214076e23 mol^-1')
    assert_equal('%d %s' % (codata.value('luminous efficacy'),
                            codata.unit('luminous efficacy')),
                 '683 lm W^-1')


def test_find_all():
    assert_(len(codata.find(disp=False)) > 300)


def test_find_single():
    assert_equal(codata.find('Wien freq', disp=False)[0],
                 'Wien frequency displacement law constant')


def test_2002_vs_2006():
    assert_almost_equal(codata.value('magn. flux quantum'),
                        codata.value('mag. flux quantum'))


def test_exact_values():
    # Check that updating stored values with exact ones worked.
    for key in codata.exact_values:
        assert_((codata.exact_values[key][0] - value(key)) / value(key) == 0)
