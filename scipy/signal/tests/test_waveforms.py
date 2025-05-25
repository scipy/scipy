import math
from typing import Literal

import numpy as np
import pytest
from pytest import raises as assert_raises
from scipy._lib._array_api import (
    array_namespace, assert_almost_equal, xp_assert_equal, xp_assert_close
)

import scipy.signal._waveforms as waveforms
from scipy.fft import irfft


# These chirp_* functions are the instantaneous frequencies of the signals
# returned by chirp().

def chirp_linear(t, f0, f1, t1):
    f = f0 + (f1 - f0) * t / t1
    return f


def chirp_quadratic(t, f0, f1, t1, vertex_zero=True):
    if vertex_zero:
        f = f0 + (f1 - f0) * t**2 / t1**2
    else:
        f = f1 - (f1 - f0) * (t1 - t)**2 / t1**2
    return f


def chirp_geometric(t, f0, f1, t1):
    f = f0 * (f1/f0)**(t/t1)
    return f


def chirp_hyperbolic(t, f0, f1, t1):
    f = f0*f1*t1 / ((f0 - f1)*t + f1*t1)
    return f


def compute_frequency(t, theta):
    """
    Compute theta'(t)/(2*pi), where theta'(t) is the derivative of theta(t).
    """
    # Assume theta and t are 1-D NumPy arrays.
    # Assume that t is uniformly spaced.
    dt = t[1] - t[0]
    f = np.diff(theta)/(2*np.pi) / dt
    tf = 0.5*(t[1:] + t[:-1])
    return tf, f


class TestChirp:

    def test_exceptions(self):
        """Raise all exceptions in used functions. """
        t = np.arange(3)
        with pytest.raises(ValueError, match="^For a logarithmic chirp, f0 and f1 "):
            waveforms._chirp_phase(t, f0=-1, t1=2, f1=1, method='log')
        for f0, f1 in [(0, 1), (1, 0)]:
            with pytest.raises(ValueError, match="^For a hyperbolic chirp, f0 and f1"):
                waveforms._chirp_phase(t, f0=f0, t1=2, f1=f1, method='hyp')
        with pytest.raises(ValueError, match="^method must be 'linear', 'quadratic',"):
            waveforms._chirp_phase(t, f0=f0, t1=2, f1=f1, method='INVALID')

    def test_linear_at_zero(self):
        w = waveforms.chirp(t=0, f0=1.0, f1=2.0, t1=1.0, method='linear')
        assert_almost_equal(w, 1.0)

    def test_linear_freq_01(self):
        method = 'linear'
        f0 = 1.0
        f1 = 2.0
        t1 = 1.0
        t = np.linspace(0, t1, 100)
        phase = waveforms._chirp_phase(t, f0, t1, f1, method)
        tf, f = compute_frequency(t, phase)
        abserr = np.max(np.abs(f - chirp_linear(tf, f0, f1, t1)))
        assert abserr < 1e-6

    def test_linear_freq_02(self):
        method = 'linear'
        f0 = 200.0
        f1 = 100.0
        t1 = 10.0
        t = np.linspace(0, t1, 100)
        phase = waveforms._chirp_phase(t, f0, t1, f1, method)
        tf, f = compute_frequency(t, phase)
        abserr = np.max(np.abs(f - chirp_linear(tf, f0, f1, t1)))
        assert abserr < 1e-6

    def test_linear_complex_power(self):
        method = 'linear'
        f0 = 1.0
        f1 = 2.0
        t1 = 1.0
        t = np.linspace(0, t1, 100)
        w_real = waveforms.chirp(t, f0, t1, f1, method, complex=False)
        w_complex = waveforms.chirp(t, f0, t1, f1, method, complex=True)
        w_pwr_r = np.var(w_real)
        w_pwr_c = np.var(w_complex)

        # Making sure that power of the real part is not affected with
        # complex conversion operation
        err = w_pwr_r - np.real(w_pwr_c)

        assert(err < 1e-6)

    def test_linear_complex_at_zero(self):
        w = waveforms.chirp(t=0, f0=-10.0, f1=1.0, t1=1.0, method='linear',
                            complex=True)
        xp_assert_close(w, 1.0+0.0j)  # dtype must match

    def test_quadratic_at_zero(self):
        w = waveforms.chirp(t=0, f0=1.0, f1=2.0, t1=1.0, method='quadratic')
        assert_almost_equal(w, 1.0)

    def test_quadratic_at_zero2(self):
        w = waveforms.chirp(t=0, f0=1.0, f1=2.0, t1=1.0, method='quadratic',
                            vertex_zero=False)
        assert_almost_equal(w, 1.0)

    def test_quadratic_complex_at_zero(self):
        w = waveforms.chirp(t=0, f0=-1.0, f1=2.0, t1=1.0, method='quadratic',
                            complex=True)
        xp_assert_close(w, 1.0+0j)

    def test_quadratic_freq_01(self):
        method = 'quadratic'
        f0 = 1.0
        f1 = 2.0
        t1 = 1.0
        t = np.linspace(0, t1, 2000)
        phase = waveforms._chirp_phase(t, f0, t1, f1, method)
        tf, f = compute_frequency(t, phase)
        abserr = np.max(np.abs(f - chirp_quadratic(tf, f0, f1, t1)))
        assert abserr < 1e-6

    def test_quadratic_freq_02(self):
        method = 'quadratic'
        f0 = 20.0
        f1 = 10.0
        t1 = 10.0
        t = np.linspace(0, t1, 2000)
        phase = waveforms._chirp_phase(t, f0, t1, f1, method)
        tf, f = compute_frequency(t, phase)
        abserr = np.max(np.abs(f - chirp_quadratic(tf, f0, f1, t1)))
        assert abserr < 1e-6

    def test_logarithmic_at_zero(self):
        w = waveforms.chirp(t=0, f0=1.0, f1=2.0, t1=1.0, method='logarithmic')
        assert_almost_equal(w, 1.0)

    def test_logarithmic_freq_01(self):
        method = 'logarithmic'
        f0 = 1.0
        f1 = 2.0
        t1 = 1.0
        t = np.linspace(0, t1, 10000)
        phase = waveforms._chirp_phase(t, f0, t1, f1, method)
        tf, f = compute_frequency(t, phase)
        abserr = np.max(np.abs(f - chirp_geometric(tf, f0, f1, t1)))
        assert abserr < 1e-6

    def test_logarithmic_freq_02(self):
        method = 'logarithmic'
        f0 = 200.0
        f1 = 100.0
        t1 = 10.0
        t = np.linspace(0, t1, 10000)
        phase = waveforms._chirp_phase(t, f0, t1, f1, method)
        tf, f = compute_frequency(t, phase)
        abserr = np.max(np.abs(f - chirp_geometric(tf, f0, f1, t1)))
        assert abserr < 1e-6

    def test_logarithmic_freq_03(self):
        method = 'logarithmic'
        f0 = 100.0
        f1 = 100.0
        t1 = 10.0
        t = np.linspace(0, t1, 10000)
        phase = waveforms._chirp_phase(t, f0, t1, f1, method)
        tf, f = compute_frequency(t, phase)
        abserr = np.max(np.abs(f - chirp_geometric(tf, f0, f1, t1)))
        assert abserr < 1e-6

    def test_hyperbolic_at_zero(self):
        w = waveforms.chirp(t=0, f0=10.0, f1=1.0, t1=1.0, method='hyperbolic')
        assert_almost_equal(w, 1.0)

    def test_hyperbolic_freq_01(self):
        method = 'hyperbolic'
        t1 = 1.0
        t = np.linspace(0, t1, 10000)
        #           f0     f1
        cases = [[10.0, 1.0],
                 [1.0, 10.0],
                 [-10.0, -1.0],
                 [-1.0, -10.0]]
        for f0, f1 in cases:
            phase = waveforms._chirp_phase(t, f0, t1, f1, method)
            tf, f = compute_frequency(t, phase)
            expected = chirp_hyperbolic(tf, f0, f1, t1)
            xp_assert_close(f, expected, atol=1e-7)

    def test_hyperbolic_zero_freq(self):
        # f0=0 or f1=0 must raise a ValueError.
        method = 'hyperbolic'
        t1 = 1.0
        t = np.linspace(0, t1, 5)
        assert_raises(ValueError, waveforms.chirp, t, 0, t1, 1, method)
        assert_raises(ValueError, waveforms.chirp, t, 1, t1, 0, method)

    def test_hyperbolic_const_freq(self):
        """Test case parameter f1 == f2.

        This test is required to achieve 100% coverage.
        """
        t, x_ref = np.arange(3), np.array([1., 1., 1.])

        x = waveforms.chirp(t, f0=2, t1=2, f1=2, method='hyp')
        xp_assert_close(x, x_ref)

    def test_unknown_method(self):
        method = "foo"
        f0 = 10.0
        f1 = 20.0
        t1 = 1.0
        t = np.linspace(0, t1, 10)
        assert_raises(ValueError, waveforms.chirp, t, f0, t1, f1, method)

    def test_integer_t1(self):
        f0 = 10.0
        f1 = 20.0
        t = np.linspace(-1, 1, 11)
        t1 = 3.0
        float_result = waveforms.chirp(t, f0, t1, f1)
        t1 = 3
        int_result = waveforms.chirp(t, f0, t1, f1)
        err_msg = "Integer input 't1=3' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)

    def test_integer_f0(self):
        f1 = 20.0
        t1 = 3.0
        t = np.linspace(-1, 1, 11)
        f0 = 10.0
        float_result = waveforms.chirp(t, f0, t1, f1)
        f0 = 10
        int_result = waveforms.chirp(t, f0, t1, f1)
        err_msg = "Integer input 'f0=10' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)

    def test_integer_f1(self):
        f0 = 10.0
        t1 = 3.0
        t = np.linspace(-1, 1, 11)
        f1 = 20.0
        float_result = waveforms.chirp(t, f0, t1, f1)
        f1 = 20
        int_result = waveforms.chirp(t, f0, t1, f1)
        err_msg = "Integer input 'f1=20' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)

    def test_integer_all(self):
        f0 = 10
        t1 = 3
        f1 = 20
        t = np.linspace(-1, 1, 11)
        float_result = waveforms.chirp(t, float(f0), float(t1), float(f1))
        int_result = waveforms.chirp(t, f0, t1, f1)
        err_msg = "Integer input 'f0=10, t1=3, f1=20' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)


class TestSweepPoly:

    def test_sweep_poly(self):
        """This test ensures that function `sweep_poly` has 100% coverage. """
        t, p = np.arange(3), np.ones(2)
        x_ref = np.array([1., -1., 1.])

        x = waveforms.sweep_poly(t, p)
        xp_assert_close(x, x_ref)


    def test_sweep_poly_quad1(self):
        p = np.poly1d([1.0, 0.0, 1.0])
        t = np.linspace(0, 3.0, 10000)
        phase = waveforms._sweep_poly_phase(t, p)
        tf, f = compute_frequency(t, phase)
        expected = p(tf)
        abserr = np.max(np.abs(f - expected))
        assert abserr < 1e-6

    def test_sweep_poly_const(self):
        p = np.poly1d(2.0)
        t = np.linspace(0, 3.0, 10000)
        phase = waveforms._sweep_poly_phase(t, p)
        tf, f = compute_frequency(t, phase)
        expected = p(tf)
        abserr = np.max(np.abs(f - expected))
        assert abserr < 1e-6

    def test_sweep_poly_linear(self):
        p = np.poly1d([-1.0, 10.0])
        t = np.linspace(0, 3.0, 10000)
        phase = waveforms._sweep_poly_phase(t, p)
        tf, f = compute_frequency(t, phase)
        expected = p(tf)
        abserr = np.max(np.abs(f - expected))
        assert abserr < 1e-6

    def test_sweep_poly_quad2(self):
        p = np.poly1d([1.0, 0.0, -2.0])
        t = np.linspace(0, 3.0, 10000)
        phase = waveforms._sweep_poly_phase(t, p)
        tf, f = compute_frequency(t, phase)
        expected = p(tf)
        abserr = np.max(np.abs(f - expected))
        assert abserr < 1e-6

    def test_sweep_poly_cubic(self):
        p = np.poly1d([2.0, 1.0, 0.0, -2.0])
        t = np.linspace(0, 2.0, 10000)
        phase = waveforms._sweep_poly_phase(t, p)
        tf, f = compute_frequency(t, phase)
        expected = p(tf)
        abserr = np.max(np.abs(f - expected))
        assert abserr < 1e-6

    def test_sweep_poly_cubic2(self):
        """Use an array of coefficients instead of a poly1d."""
        p = np.array([2.0, 1.0, 0.0, -2.0])
        t = np.linspace(0, 2.0, 10000)
        phase = waveforms._sweep_poly_phase(t, p)
        tf, f = compute_frequency(t, phase)
        expected = np.poly1d(p)(tf)
        abserr = np.max(np.abs(f - expected))
        assert abserr < 1e-6

    def test_sweep_poly_cubic3(self):
        """Use a list of coefficients instead of a poly1d."""
        p = [2.0, 1.0, 0.0, -2.0]
        t = np.linspace(0, 2.0, 10000)
        phase = waveforms._sweep_poly_phase(t, p)
        tf, f = compute_frequency(t, phase)
        expected = np.poly1d(p)(tf)
        abserr = np.max(np.abs(f - expected))
        assert abserr < 1e-6


class TestGaussPulse:

    def test_exceptions(self):
        """Raise all exceptions in function. """
        with pytest.raises(ValueError, match="^Center frequency "):
             waveforms.gausspulse('cutoff', fc=-1.0)
        with pytest.raises(ValueError, match="^Fractional bandwidth "):
            waveforms.gausspulse('cutoff', fc=1000.0, bw=-0.5)
        with pytest.raises(ValueError, match="^Reference level for bandwidth "):
            waveforms.gausspulse('cutoff', bwr=1.0)
        with pytest.raises(ValueError, match="^Reference level for time cutoff must"):
            waveforms.gausspulse('cutoff', tpr=1.0)
        with pytest.raises(ValueError, match="^If `t` is a string, it must be "):
            waveforms.gausspulse('INVALID')

    def test_integer_fc(self):
        float_result = waveforms.gausspulse('cutoff', fc=1000.0)
        int_result = waveforms.gausspulse('cutoff', fc=1000)
        err_msg = "Integer input 'fc=1000' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)

    def test_integer_bw(self):
        float_result = waveforms.gausspulse('cutoff', bw=1.0)
        int_result = waveforms.gausspulse('cutoff', bw=1)
        err_msg = "Integer input 'bw=1' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)

    def test_integer_bwr(self):
        float_result = waveforms.gausspulse('cutoff', bwr=-6.0)
        int_result = waveforms.gausspulse('cutoff', bwr=-6)
        err_msg = "Integer input 'bwr=-6' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)

    def test_integer_tpr(self):
        float_result = waveforms.gausspulse('cutoff', tpr=-60.0)
        int_result = waveforms.gausspulse('cutoff', tpr=-60)
        err_msg = "Integer input 'tpr=-60' gives wrong result"
        xp_assert_equal(int_result, float_result, err_msg=err_msg)

    def test_parameter_t_values(self):
        """Pass numeric array instead of string "cutoff". """
        t = np.array([-0.5, 0, 0.5])
        x0_ref = np.array([-0.7999183981317266, 1, -0.7999183981317266])
        e0_ref = np.array([0.7999183981317266, 1, 0.7999183981317266])
        y0_ref = np.array([-9.796175058510971e-17, 0, 9.796175058510971e-17])


        x = waveforms.gausspulse(t, fc=1, retquad=False, retenv=False)
        print(f"{x[0]=:0.16f}")
        xp_assert_close(x, x0_ref, rtol=1e-8,
                        err_msg="Invalid result x for retquad=False, retenv=False)")

        x, e = waveforms.gausspulse(t, fc=1, retquad=False, retenv=True)
        xp_assert_close(x, x0_ref, rtol=1e-8,
                        err_msg="Invalid result x for retquad=False, retenv=True)")
        xp_assert_close(e, e0_ref, rtol=1e-8,
                        err_msg="Invalid result e for retquad=False, retenv=True)")

        x, y = waveforms.gausspulse(t, fc=1, retquad=True, retenv=False)
        print(f"{y[0]=:.16g}")
        xp_assert_close(x, x0_ref, rtol=1e-8,
                        err_msg="Invalid result x for retquad=True, retenv=False)")
        xp_assert_close(y, y0_ref, rtol=1e-8,
                        err_msg="Invalid result y for retquad=True, retenv=False)")

        x, y, e = waveforms.gausspulse(t, fc=1, retquad=True, retenv=True)
        print(f"{y[0]=:.16g}")
        xp_assert_close(x, x0_ref, rtol=1e-8,
                        err_msg="Invalid result x for retquad=True, retenv=True)")
        xp_assert_close(y, y0_ref, rtol=1e-8,
                        err_msg="Invalid result y for retquad=True, retenv=True)")
        xp_assert_close(e, e0_ref, rtol=1e-8,
                        err_msg="Invalid result e for retquad=True, retenv=True)")

class TestUnitImpulse:

    def test_no_index(self):
        xp_assert_equal(waveforms.unit_impulse(7),
                        np.asarray([1.0, 0, 0, 0, 0, 0, 0]))
        xp_assert_equal(waveforms.unit_impulse((3, 3)),
                        np.asarray([[1.0, 0, 0], [0, 0, 0], [0, 0, 0]]))

    def test_index(self):
        xp_assert_equal(waveforms.unit_impulse(10, 3),
                        np.asarray([0.0, 0, 0, 1, 0, 0, 0, 0, 0, 0]))
        xp_assert_equal(waveforms.unit_impulse((3, 3), (1, 1)),
                        np.asarray([[0.0, 0, 0], [0, 1, 0], [0, 0, 0]]))

        # Broadcasting
        imp = waveforms.unit_impulse((4, 4), 2)
        xp_assert_equal(imp, np.asarray([[0.0, 0, 0, 0],
                                         [0.0, 0, 0, 0],
                                         [0.0, 0, 1, 0],
                                         [0.0, 0, 0, 0]]))

    def test_mid(self):
        xp_assert_equal(waveforms.unit_impulse((3, 3), 'mid'),
                        np.asarray([[0.0, 0, 0], [0, 1, 0], [0, 0, 0]]))
        xp_assert_equal(waveforms.unit_impulse(9, 'mid'),
                        np.asarray([0.0, 0, 0, 0, 1, 0, 0, 0, 0]))

    def test_dtype(self):
        imp = waveforms.unit_impulse(7)
        assert np.issubdtype(imp.dtype, np.floating)

        imp = waveforms.unit_impulse(5, 3, dtype=int)
        assert np.issubdtype(imp.dtype, np.integer)

        imp = waveforms.unit_impulse((5, 2), (3, 1), dtype=complex)
        assert np.issubdtype(imp.dtype, np.complexfloating)


class TestSawtoothWaveform:
    def test_dtype(self):
        waveform = waveforms.sawtooth(
            np.array(1, dtype=np.float32), width=np.float32(1)
        )
        assert waveform.dtype == np.float64

        waveform = waveforms.sawtooth(1)
        assert waveform.dtype == np.float64

    @pytest.mark.parametrize("width, x_ref", ([
                              (   0, [ 1,  1/2,   0, -1/2]),
                              (0.25, [-1,    1, 1/3, -1/3]),
                              ( 0.5, [-1,    0,   1,    0]),
                              (0.75, [-1, -1/3, 1/3,    1]),
                              (   1, [-1, -1/2,   0,  1/2])]))
    def test_values(self, width, x_ref):
        t = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        x = waveforms.sawtooth(t, width)
        xp_assert_close(x, np.asarray(x_ref, dtype=float), atol=1e-12)


class TestSawtoothRFFT:
    """Unit tests for function `signal.sawtooth_rfft`. """

    def setup_class(self):
        """Reference values are calculated here. """
        xp = array_namespace()
        # For self.test_values():
        vv0 = -8 / 3 * xp.array([1 + 1j, 2., 1 - 1j, 0, 1 + 1j])
        self.Xa_ref = xp.zeros((6,), dtype=vv0.dtype)
        self.Xa_ref[1:] = vv0 / (xp.arange(1, 6) * xp.pi) ** 2
        self.Xb_ref = xp.zeros_like(self.Xa_ref)
        self.Xb_ref[2::2] = self.Xa_ref[1:3]

        # For self.test_norm_parameter():
        self.x0_ref = irfft(self.Xa_ref, norm='forward')

        # For self.test_duty01():
        self.X01a_ref = xp.zeros(6, dtype=complex)
        self.X01a_ref[1:] = 1j / xp.pi / xp.arange(1, 6)
        self.X01b_ref = xp.zeros_like(self.X01a_ref)
        self.X01b_ref[2::2] = self.X01a_ref[1:3]

    def test_exceptions(self):
        """Raise all exceptions in function. """
        for n in (10.5, -10):
            with pytest.raises(ValueError, match=f"^Parameter {n=} is not a positive"):
                waveforms.sawtooth_rfft(n, 1, 0.5)
        for m_cyc in (5.5, 0, 10):
            with pytest.raises(ValueError, match=f"^Parameter {m_cyc=} is not a "):
                waveforms.sawtooth_rfft(10, m_cyc, 0.5)
        for duty in (-0.1, 1.1):
            with pytest.raises(ValueError, match="^0 <= duty <= 1 does not hold for "):
                waveforms.sawtooth_rfft(10, 1, duty)
        norm = 'INVALID'
        with pytest.raises(ValueError, match=f"^Parameter {norm=} not in "):
            # noinspection PyTypeChecker
            waveforms.sawtooth_rfft(10, 1, 0.5, norm=norm)

    @pytest.mark.parametrize('n', (10, 11))
    def test_values(self, n):
        """Test output against reference values. """
        Xa = waveforms.sawtooth_rfft(n, 1, duty=0.25)
        xp_assert_close(Xa, n * self.Xa_ref, atol=1e-12)

        Xb = waveforms.sawtooth_rfft(n, 2, duty=0.25)
        xp_assert_close(Xb, n * self.Xb_ref, atol=1e-12)

    @pytest.mark.parametrize('norm', ('backward', 'ortho', 'forward'))
    def test_norm_parameter(self, norm: Literal['backward', 'ortho', 'forward']):
        """Verify that parameter `norm` is compatible with `scipy.fft.irfft`. """
        X0 = waveforms.sawtooth_rfft(10, 1, duty=0.25, norm=norm)
        x0 = irfft(X0, norm=norm)
        xp_assert_close(x0, self.x0_ref, atol=1e-12)

    @pytest.mark.parametrize('norm', ('backward', 'ortho', 'forward'))
    @pytest.mark.parametrize('n', (10, 11))
    def test_duty01(self, n, norm: Literal['backward', 'ortho', 'forward']):
        """Test for parameter ``duty=0`` and ``duty=1``. """
        s = {'backward': n, 'ortho': math.sqrt(n), 'forward': 1}[norm]

        X0a = waveforms.sawtooth_rfft(n, 1, duty=0, norm=norm)
        xp_assert_close(X0a, -s * self.X01a_ref, atol=1e-12)

        X0b = waveforms.sawtooth_rfft(n, 2, duty=0, norm=norm)
        xp_assert_close(X0b, -s * self.X01b_ref, atol=1e-12)

        X1a = waveforms.sawtooth_rfft(n, 1, duty=1, norm=norm)
        xp_assert_close(X1a, s * self.X01a_ref, atol=1e-12)

        X1b = waveforms.sawtooth_rfft(n, 2, duty=1, norm=norm)
        xp_assert_close(X1b, s * self.X01b_ref, atol=1e-12)



class TestSquareWaveform:
    def test_dtype(self):
        waveform = waveforms.square(np.array(1, dtype=np.float32),
                                    duty=np.float32(0.5))
        assert waveform.dtype == np.float64

        waveform = waveforms.square(1)
        assert waveform.dtype == np.float64

    @pytest.mark.parametrize('duty, x_ref', ([
                             (   0, [-1, -1, -1, -1]),
                             (0.25, [ 1, -1, -1, -1]),
                             ( 0.5, [ 1,  1, -1, -1]),
                             (0.75, [ 1,  1,  1, -1]),
                             (   1, [ 1,  1,  1,  1])]))
    def test_values(self, duty, x_ref):
        t = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        x = waveforms.square(t, duty)
        xp_assert_equal(x, np.asarray(x_ref, dtype=np.float64))


class TestSquareRFFT:
    """Unit tests for function `signal.square_rfft`. """

    def setup_class(self):
        """Reference values are calculated here. """
        # For self.test_values():
        xp = array_namespace()
        self.X0 = xp.array([0, 1 - 1j, -1j, -(1 + 1j) / 3, 0., (1 - 1j) / 5]) / xp.pi
        self.X0[0] = -0.5

        # For self.test_norm_parameter():
        X_ref = xp.asarray([4., 0., -(8 + 8j) / xp.pi, 0., -8j / xp.pi])
        self.x_ref = irfft(X_ref, norm='backward')

    def test_exceptions(self):
        """Raise all exceptions in function. """
        for n in (10.5, -10):
            with pytest.raises(ValueError, match=f"^Parameter {n=} is not a positive"):
                waveforms.square_rfft(n, 1, 0.5)
        for m_cyc in (5.5, 0, 10):
            with pytest.raises(ValueError, match=f"^Parameter {m_cyc=} is not a "):
                waveforms.square_rfft(10, m_cyc, 0.5)
        for duty in (0, 1):
            with pytest.raises(ValueError, match="^0 < duty < 1 does not hold for "):
                waveforms.square_rfft(10, 1, duty)
        norm = 'INVALID'
        with pytest.raises(ValueError, match=f"^Parameter {norm=} not in "):
            # noinspection PyTypeChecker
            waveforms.square_rfft(10, 1, 0.5, norm=norm)

    @pytest.mark.parametrize('n', (10, 11))
    def test_values(self, n):
        """Test output against reference values. """
        X = waveforms.square_rfft(n, 1, 0.25)
        xp_assert_close(X, self.X0 * n, atol=1e-12)

    @pytest.mark.parametrize('norm', ('backward', 'ortho', 'forward'))
    def test_norm_parameter(self, norm: Literal['backward', 'ortho', 'forward']):
        """Test that parameter `norm` is compatible with `scipy.fft.irfft`. """
        X = waveforms.square_rfft(8, 2, 0.75, norm=norm)
        x = irfft(X, norm=norm)
        xp_assert_close(x, self.x_ref, atol=1e-12)
