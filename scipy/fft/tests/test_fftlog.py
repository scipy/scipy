import numpy as np
from numpy.testing import assert_allclose
import pytest

from scipy.fft._fftlog import fht, ifht, fhtoffset


def test_fht_agrees_with_fftlog():
    def f(r, mu):
        return r**(mu+1)*np.exp(-r**2/2)

    # test 1: do not change to optimal offset
    r = np.logspace(-4, 4, 16)
    dln = np.log(r[1]/r[0])
    mu = 0.3
    offset = 0.0
    bias = 0.0
    a = f(r, mu)
    ours = fht(a, dln, mu, offset=offset, bias=bias)
    theirs = [ -0.1159922613593045E-02,  0.1625822618458832E-02,
               -0.1949518286432330E-02,  0.3789220182554077E-02,
                0.5093959119952945E-03,  0.2785387803618774E-01,
                0.9944952700848897E-01,  0.4599202164586588    ,
                0.3157462160881342    , -0.8201236844404755E-03,
               -0.7834031308271878E-03,  0.3931444945110708E-03,
               -0.2697710625194777E-03,  0.3568398050238820E-03,
               -0.5554454827797206E-03,  0.8286331026468585E-03 ]
    assert_allclose(ours, theirs)

    # test 2: change to optimal offset
    offset = fhtoffset(dln, mu, bias=bias)
    ours = fht(a, dln, mu, offset=offset, bias=bias)
    theirs = [  0.4353768523152057E-04, -0.9197045663594285E-05,
                0.3150140927838524E-03,  0.9149121960963704E-03,
                0.5808089753959363E-02,  0.2548065256377240E-01,
                0.1339477692089897    ,  0.4821530509479356    ,
                0.2659899781579785    , -0.1116475278448113E-01,
                0.1791441617592385E-02, -0.4181810476548056E-03,
                0.1314963536765343E-03, -0.5422057743066297E-04,
                0.3208681804170443E-04, -0.2696849476008234E-04 ]
    assert_allclose(ours, theirs)

    # test 3: positive bias
    bias = 0.8
    offset = fhtoffset(dln, mu, bias=bias)
    ours = fht(a, dln, mu, offset=offset, bias=bias)
    theirs = [ -0.5269564404680044E-03,  0.6621809022379771E-03,
               -0.8650417270041665E-03,  0.1228803626821765E-02,
               -0.1585545936440683E-02,  0.6459647456325081E-02,
                0.3827405357707621E-01,  0.4517325873081965    ,
                0.9078447643770889E-01, -0.1492998636505558    ,
               -0.7951038767491089E-02, -0.5637860690939526E-02,
               -0.2350071850573260E-03, -0.7638491233662864E-04,
               -0.3149905120494838E-03,  0.4071038541859640E-03 ]
    assert_allclose(ours, theirs)

    # test 4: negative bias
    bias = -0.8
    offset = fhtoffset(dln, mu, bias=bias)
    ours = fht(a, dln, mu, offset=offset, bias=bias)
    theirs = [  0.1032111432323560E-01,  0.1861355637905322E-01,
                0.3501096038827679E-01,  0.6380101165463653E-01,
                0.1193350913910257    ,  0.2175062758501763    ,
                0.4041134221771817    ,  0.6219659277167918    ,
                0.3624727167380778    ,  0.5711362574323833E-01,
                0.1115827451790637E-01,  0.2803499258222323E-02,
                0.1353303187235311E-02,  0.1639742462916552E-02,
                0.3090587168643993E-02,  0.5393830473972794E-02 ]
    assert_allclose(ours, theirs)


@pytest.mark.parametrize('optimal', [True, False])
@pytest.mark.parametrize('offset', [0.0, 1.0, -1.0])
@pytest.mark.parametrize('bias', [0, 0.5, -0.5])
@pytest.mark.parametrize('n', [64, 63])
def test_fht_inverse(n, bias, offset, optimal):
    a = np.random.randn(n)
    dln = np.random.uniform(-1, 1)
    mu = np.random.uniform(-2, 2)

    if optimal:
        offset = fhtoffset(dln, mu, initial=offset, bias=bias)

    A = fht(a, dln, mu, offset=offset, bias=bias)
    a_ = ifht(A, dln, mu, offset=offset, bias=bias)

    assert_allclose(a, a_)


def test_fht_special_cases():
    a = np.random.randn(64)
    dln = np.random.uniform(-1, 1)

    # let xp = (mu+1+q)/2, xm = (mu+1-q)/2, M = {0, -1, -2, ...}

    # case 1: xp in M, xm in M => well-defined transform
    mu, bias = -4.0, 1.0
    with pytest.warns(None) as record:
        fht(a, dln, mu, bias=bias)
        assert not record, 'fhtq warned about a well-defined transform'

    # case 2: xp not in M, xm in M => well-defined transform
    mu, bias = -2.5, 0.5
    with pytest.warns(None) as record:
        fht(a, dln, mu, bias=bias)
        assert not record, 'fhtq warned about a well-defined transform'

    # case 3: xp in M, xm not in M => singular transform
    mu, bias = -3.5, 0.5
    with pytest.warns(Warning) as record:
        fht(a, dln, mu, bias=bias)
        assert record, 'fhtq did not warn about a singular transform'
