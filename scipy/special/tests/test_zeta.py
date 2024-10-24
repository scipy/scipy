import scipy.special as sc
import numpy as np
import pytest

from numpy.testing import assert_equal, assert_allclose


def test_zeta():
    assert_allclose(sc.zeta(2,2), np.pi**2/6 - 1, rtol=1e-12)


def test_zetac():
    # Expected values in the following were computed using Wolfram
    # Alpha's `Zeta[x] - 1`
    x = [-2.1, 0.8, 0.9999, 9, 50, 75]
    desired = [
        -0.9972705002153750,
        -5.437538415895550,
        -10000.42279161673,
        0.002008392826082214,
        8.881784210930816e-16,
        2.646977960169853e-23,
    ]
    assert_allclose(sc.zetac(x), desired, rtol=1e-12)


def test_zetac_special_cases():
    assert sc.zetac(np.inf) == 0
    assert np.isnan(sc.zetac(-np.inf))
    assert sc.zetac(0) == -1.5
    assert sc.zetac(1.0) == np.inf

    assert_equal(sc.zetac([-2, -50, -100]), -1)


def test_riemann_zeta_special_cases():
    assert np.isnan(sc.zeta(np.nan))
    assert sc.zeta(np.inf) == 1
    assert sc.zeta(0) == -0.5

    # Riemann zeta is zero add negative even integers.
    assert_equal(sc.zeta([-2, -4, -6, -8, -10]), 0)

    assert_allclose(sc.zeta(2), np.pi**2/6, rtol=1e-12)
    assert_allclose(sc.zeta(4), np.pi**4/90, rtol=1e-12)


def test_riemann_zeta_avoid_overflow():
    s = -260.00000000001
    desired = -5.6966307844402683127e+297  # Computed with Mpmath
    assert_allclose(sc.zeta(s), desired, atol=0, rtol=5e-14)


@pytest.mark.parametrize(
    "z, desired",
    # # Reference cases taken from mpmath with the script:
    # import numpy as np
    # import scipy.stats as stats

    # from mpmath import mp

    # # seed = np.random.SeedSequence().entropy
    # seed = 154689806791763421822480125722191067828
    # rng = np.random.default_rng(seed)

    # # A small point in each quadrant outside of the critical strip
    # cases = []
    # for x_sign, y_sign in [1, 1], [1, -1], [-1, 1], [-1, -1]:
    #     x = x_sign * rng.uniform(2, 8)
    #     y = y_sign * rng.uniform(2, 8)
    #     z = x + y*1j
    #     reference = complex(mp.zeta(z))
    #     cases.append((z, reference))

    # # Moderately large imaginary part in each quadrant outside of critical strip
    # for x_sign, y_sign in [1, 1], [1, -1], [-1, 1], [-1, -1]:
    #     x = x_sign * rng.uniform(2, 8)
    #     y = y_sign * rng.uniform(50, 80)
    #     z = x + y*1j
    #     reference = complex(mp.zeta(z))
    #     cases.append((z, reference))

    # # points in critical strip
    # x = rng.uniform(-1.0, 1.0, size=5)
    # y = np.exp(rng.uniform(0, 5, size=5))
    # z = x + y*1j
    # for t in z:
    #     reference = complex(mp.zeta(t))
    #     cases.append((complex(t), reference))
    # z = x - y*1j
    # for t in z:
    #     reference = complex(mp.zeta(t))
    #     cases.append((complex(t), reference))
    [
        # small values in each quadrant outside critical strip
        ((3.12838509346655+7.111085974836645j),
         (1.0192654793474945+0.08795174413289127j)),
        ((7.06791362314716-7.219497492626728j),
         (1.0020740683598117-0.006752725913243711j)),
        ((-6.806227077655519+2.724411451005281j),
         (0.06312488213559667-0.061641496333765956j)),
        ((-3.0170751511621026-6.3686522550665945j),
         (-0.10330747857150148-1.214541994832571j)),
        # Moderately large imaginary part, outside critical strip
        ((6.133994402212294+60.03091448000761j),
         (0.9885701843417336+0.009636925981078128j)),
        ((6.17268142822657-64.74883149743795j),
         (1.0080474225840865+0.012032804974965354j)),
        ((-3.462191939791879+76.16258975567534j),
         (18672.072070850158+2908.5104826247184j)),
        ((-6.955735216531752-74.75791554155748j),
         (-77672258.72276545+71625206.0401107j)),
        # Inside critical strip, varying imaginary part.
        ((-0.1823923420352156+1.4596830498094384j),
         (0.15641024400724224-0.29104726162518146j)),
        ((0.9346987902419266+4.918968547259143j),
         (0.7452215845376634+0.17593462255710995j)),
        ((0.7384965359955509+66.6142398421354j),
         (0.5365014352325163-0.35081725113825424j)),
        ((-0.14456304559992472+21.783747851715468j),
         (0.1373304916093494+1.6091267223996006j)),
        ((-0.5904101064314209+33.17656449538932j),
         (-4.409223214592868-1.6987333585451423j)),
         ((-0.1823923420352156-1.4596830498094384j),
          (0.15641024400724224+0.29104726162518146j)),
        ((0.9346987902419266-4.918968547259143j),
         (0.7452215845376634-0.17593462255710995j)),
        ((0.7384965359955509-66.6142398421354j),
         (0.5365014352325163+0.35081725113825424j)),
        ((-0.14456304559992472-21.783747851715468j),
         (0.1373304916093494-1.6091267223996006j)),
        ((-0.5904101064314209-33.17656449538932j),
         (-4.409223214592868+1.6987333585451423j)),
    ]
)
def test_riemann_zeta_complex(z, desired):
    assert_allclose(sc.zeta(z), desired, rtol=1e-13)
