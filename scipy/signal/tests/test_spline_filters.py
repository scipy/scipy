# pylint: disable=missing-docstring
import numpy as np
from numpy import array
from numpy.testing import (assert_allclose, assert_array_equal,
                           assert_almost_equal)
import pytest
from pytest import raises

import scipy.signal._spline_filters as bsp
from scipy import signal

from scipy.signal._spline import (
    symiirorder1_ic, symiirorder2_ic_fwd, symiirorder2_ic_bwd)
from scipy.signal._spline_filters import symiirorder1, symiirorder2


class TestBSplines:
    """Test behaviors of B-splines. Some of the values tested against were
    returned as of SciPy 1.1.0 and are included for regression testing
    purposes. Others (at integer points) are compared to theoretical
    expressions (cf. Unser, Aldroubi, Eden, IEEE TSP 1993, Table 1)."""

    def test_spline_filter(self):
        np.random.seed(12457)
        # Test the type-error branch
        raises(TypeError, bsp.spline_filter, array([0]), 0)
        # Test the real branch
        np.random.seed(12457)
        data_array_real = np.random.rand(12, 12)
        # make the magnitude exceed 1, and make some negative
        data_array_real = 10*(1-2*data_array_real)
        result_array_real = array(
            [[-.463312621, 8.33391222, .697290949, 5.28390836,
              5.92066474, 6.59452137, 9.84406950, -8.78324188,
              7.20675750, -8.17222994, -4.38633345, 9.89917069],
             [2.67755154, 6.24192170, -3.15730578, 9.87658581,
              -9.96930425, 3.17194115, -4.50919947, 5.75423446,
              9.65979824, -8.29066885, .971416087, -2.38331897],
             [-7.08868346, 4.89887705, -1.37062289, 7.70705838,
              2.51526461, 3.65885497, 5.16786604, -8.77715342e-03,
              4.10533325, 9.04761993, -.577960351, 9.86382519],
             [-4.71444301, -1.68038985, 2.84695116, 1.14315938,
              -3.17127091, 1.91830461, 7.13779687, -5.35737482,
              -9.66586425, -9.87717456, 9.93160672, 4.71948144],
             [9.49551194, -1.92958436, 6.25427993, -9.05582911,
              3.97562282, 7.68232426, -1.04514824, -5.86021443,
              -8.43007451, 5.47528997, 2.06330736, -8.65968112],
             [-8.91720100, 8.87065356, 3.76879937, 2.56222894,
              -.828387146, 8.72288903, 6.42474741, -6.84576083,
              9.94724115, 6.90665380, -6.61084494, -9.44907391],
             [9.25196790, -.774032030, 7.05371046, -2.73505725,
              2.53953305, -1.82889155, 2.95454824, -1.66362046,
              5.72478916, -3.10287679, 1.54017123, -7.87759020],
             [-3.98464539, -2.44316992, -1.12708657, 1.01725672,
              -8.89294671, -5.42145629, -6.16370321, 2.91775492,
              9.64132208, .702499998, -2.02622392, 1.56308431],
             [-2.22050773, 7.89951554, 5.98970713, -7.35861835,
              5.45459283, -7.76427957, 3.67280490, -4.05521315,
              4.51967507, -3.22738749, -3.65080177, 3.05630155],
             [-6.21240584, -.296796126, -8.34800163, 9.21564563,
              -3.61958784, -4.77120006, -3.99454057, 1.05021988e-03,
              -6.95982829, 6.04380797, 8.43181250, -2.71653339],
             [1.19638037, 6.99718842e-02, 6.72020394, -2.13963198,
              3.75309875, -5.70076744, 5.92143551, -7.22150575,
              -3.77114594, -1.11903194, -5.39151466, 3.06620093],
             [9.86326886, 1.05134482, -7.75950607, -3.64429655,
              7.81848957, -9.02270373, 3.73399754, -4.71962549,
              -7.71144306, 3.78263161, 6.46034818, -4.43444731]])
        assert_allclose(bsp.spline_filter(data_array_real, 0),
                        result_array_real)

    def test_spline_filter_complex(self):
        # Test the complex branch : restore the test which was removed
        # in https://github.com/scipy/scipy/pull/17455 but should not
        # have been. cf. https://github.com/scipy/scipy/pull/9208/
        np.random.seed(12457)
        data_array_complex = np.random.rand(7, 7) + np.random.rand(7, 7)*1j
        # make the magnitude exceed 1, and make some negative
        data_array_complex = 10*(1+1j-2*data_array_complex)
        result_array_complex = array(
            [[-4.61489230e-01-1.92994022j, 8.33332443+6.25519943j,
              6.96300745e-01-9.05576038j, 5.28294849+3.97541356j,
              5.92165565+7.68240595j, 6.59493160-1.04542804j,
              9.84503460-5.85946894j],
             [-8.78262329-8.4295969j, 7.20675516+5.47528982j,
              -8.17223072+2.06330729j, -4.38633347-8.65968037j,
              9.89916801-8.91720295j, 2.67755103+8.8706522j,
              6.24192142+3.76879835j],
             [-3.15627527+2.56303072j, 9.87658501-0.82838702j,
              -9.96930313+8.72288895j, 3.17193985+6.42474651j,
              -4.50919819-6.84576082j, 5.75423431+9.94723988j,
              9.65979767+6.90665293j],
             [-8.28993416-6.61064005j, 9.71416473e-01-9.44907284j,
              -2.38331890+9.25196648j, -7.08868170-0.77403212j,
              4.89887714+7.05371094j, -1.37062311-2.73505688j,
              7.70705748+2.5395329j],
             [2.51528406-1.82964492j, 3.65885472+2.95454836j,
              5.16786575-1.66362023j, -8.77737999e-03+5.72478867j,
              4.10533333-3.10287571j, 9.04761887+1.54017115j,
              -5.77960968e-01-7.87758923j],
             [9.86398506-3.98528528j, -4.71444130-2.44316983j,
              -1.68038976-1.12708664j, 2.84695053+1.01725709j,
              1.14315915-8.89294529j, -3.17127085-5.42145538j,
              1.91830420-6.16370344j],
             [7.13875294+2.91851187j, -5.35737514+9.64132309j,
              -9.66586399+0.70250005j, -9.87717438-2.0262239j,
              9.93160629+1.5630846j, 4.71948051-2.22050714j,
              9.49550819+7.8995142j]])
        # FIXME: for complex types, the computations are done in
        # single precision (reason unclear). When this is changed,
        # this test needs updating.
        assert_allclose(bsp.spline_filter(data_array_complex, 0),
                        result_array_complex, rtol=1e-6)

    def test_gauss_spline(self):
        np.random.seed(12459)
        assert_almost_equal(bsp.gauss_spline(0, 0), 1.381976597885342)
        assert_allclose(bsp.gauss_spline(array([1.]), 1), array([0.04865217]))

    def test_gauss_spline_list(self):
        # regression test for gh-12152 (accept array_like)
        knots = [-1.0, 0.0, -1.0]
        assert_almost_equal(bsp.gauss_spline(knots, 3),
                            array([0.15418033, 0.6909883, 0.15418033]))

    def test_cspline1d(self):
        np.random.seed(12462)
        assert_array_equal(bsp.cspline1d(array([0])), [0.])
        c1d = array([1.21037185, 1.86293902, 2.98834059, 4.11660378,
                     4.78893826])
        # test lamda != 0
        assert_allclose(bsp.cspline1d(array([1., 2, 3, 4, 5]), 1), c1d)
        c1d0 = array([0.78683946, 2.05333735, 2.99981113, 3.94741812,
                      5.21051638])
        assert_allclose(bsp.cspline1d(array([1., 2, 3, 4, 5])), c1d0)

    def test_qspline1d(self):
        np.random.seed(12463)
        assert_array_equal(bsp.qspline1d(array([0])), [0.])
        # test lamda != 0
        raises(ValueError, bsp.qspline1d, array([1., 2, 3, 4, 5]), 1.)
        raises(ValueError, bsp.qspline1d, array([1., 2, 3, 4, 5]), -1.)
        q1d0 = array([0.85350007, 2.02441743, 2.99999534, 3.97561055,
                      5.14634135])
        assert_allclose(bsp.qspline1d(array([1., 2, 3, 4, 5])), q1d0)

    def test_cspline1d_eval(self):
        np.random.seed(12464)
        assert_allclose(bsp.cspline1d_eval(array([0., 0]), [0.]), array([0.]))
        assert_array_equal(bsp.cspline1d_eval(array([1., 0, 1]), []),
                           array([]))
        x = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
        dx = x[1]-x[0]
        newx = [-6., -5.5, -5., -4.5, -4., -3.5, -3., -2.5, -2., -1.5, -1.,
                -0.5, 0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.,
                6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12.,
                12.5]
        y = array([4.216, 6.864, 3.514, 6.203, 6.759, 7.433, 7.874, 5.879,
                   1.396, 4.094])
        cj = bsp.cspline1d(y)
        newy = array([6.203, 4.41570658, 3.514, 5.16924703, 6.864, 6.04643068,
                      4.21600281, 6.04643068, 6.864, 5.16924703, 3.514,
                      4.41570658, 6.203, 6.80717667, 6.759, 6.98971173, 7.433,
                      7.79560142, 7.874, 7.41525761, 5.879, 3.18686814, 1.396,
                      2.24889482, 4.094, 2.24889482, 1.396, 3.18686814, 5.879,
                      7.41525761, 7.874, 7.79560142, 7.433, 6.98971173, 6.759,
                      6.80717667, 6.203, 4.41570658])
        assert_allclose(bsp.cspline1d_eval(cj, newx, dx=dx, x0=x[0]), newy)

    def test_qspline1d_eval(self):
        np.random.seed(12465)
        assert_allclose(bsp.qspline1d_eval(array([0., 0]), [0.]), array([0.]))
        assert_array_equal(bsp.qspline1d_eval(array([1., 0, 1]), []),
                           array([]))
        x = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6]
        dx = x[1]-x[0]
        newx = [-6., -5.5, -5., -4.5, -4., -3.5, -3., -2.5, -2., -1.5, -1.,
                -0.5, 0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.,
                6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12.,
                12.5]
        y = array([4.216, 6.864, 3.514, 6.203, 6.759, 7.433, 7.874, 5.879,
                   1.396, 4.094])
        cj = bsp.qspline1d(y)
        newy = array([6.203, 4.49418159, 3.514, 5.18390821, 6.864, 5.91436915,
                      4.21600002, 5.91436915, 6.864, 5.18390821, 3.514,
                      4.49418159, 6.203, 6.71900226, 6.759, 7.03980488, 7.433,
                      7.81016848, 7.874, 7.32718426, 5.879, 3.23872593, 1.396,
                      2.34046013, 4.094, 2.34046013, 1.396, 3.23872593, 5.879,
                      7.32718426, 7.874, 7.81016848, 7.433, 7.03980488, 6.759,
                      6.71900226, 6.203, 4.49418159])
        assert_allclose(bsp.qspline1d_eval(cj, newx, dx=dx, x0=x[0]), newy)


# i/o dtypes with scipy 1.9.1, likely fixed by backwards compat
sepfir_dtype_map = {np.uint8: np.float32, int: np.float64,
                    np.float32: np.float32, float: float,
                    np.complex64: np.complex64, complex: complex}

class TestSepfir2d:
    def test_sepfir2d_invalid_filter(self):
        filt = np.array([1.0, 2.0, 4.0, 2.0, 1.0])
        image = np.random.rand(7, 9)
        # No error for odd lengths
        signal.sepfir2d(image, filt, filt[2:])

        # Row or column filter must be odd
        with pytest.raises(ValueError, match="odd length"):
            signal.sepfir2d(image, filt, filt[1:])
        with pytest.raises(ValueError, match="odd length"):
            signal.sepfir2d(image, filt[1:], filt)

        # Filters must be 1-dimensional
        with pytest.raises(ValueError, match="object too deep"):
            signal.sepfir2d(image, filt.reshape(1, -1), filt)
        with pytest.raises(ValueError, match="object too deep"):
            signal.sepfir2d(image, filt, filt.reshape(1, -1))

    def test_sepfir2d_invalid_image(self):
        filt = np.array([1.0, 2.0, 4.0, 2.0, 1.0])
        image = np.random.rand(8, 8)

        # Image must be 2 dimensional
        with pytest.raises(ValueError, match="object too deep"):
            signal.sepfir2d(image.reshape(4, 4, 4), filt, filt)

        with pytest.raises(ValueError, match="object of too small depth"):
            signal.sepfir2d(image[0], filt, filt)

    @pytest.mark.parametrize('dtyp',
        [np.uint8, int, np.float32, float, np.complex64, complex]
    )
    def test_simple(self, dtyp):
        # test values on a paper-and-pencil example
        a = np.array([[1, 2, 3, 3, 2, 1],
                      [1, 2, 3, 3, 2, 1],
                      [1, 2, 3, 3, 2, 1],
                      [1, 2, 3, 3, 2, 1]], dtype=dtyp)
        h1 = [0.5, 1, 0.5]
        h2 = [1]
        result = signal.sepfir2d(a, h1, h2)
        expected = array([[2.5, 4. , 5.5, 5.5, 4. , 2.5],
                          [2.5, 4. , 5.5, 5.5, 4. , 2.5],
                          [2.5, 4. , 5.5, 5.5, 4. , 2.5],
                          [2.5, 4. , 5.5, 5.5, 4. , 2.5]])

        assert_allclose(result, expected, atol=1e-16)
        assert result.dtype == sepfir_dtype_map[dtyp]

        result = signal.sepfir2d(a, h2, h1)
        expected = array([[2., 4., 6., 6., 4., 2.],
                          [2., 4., 6., 6., 4., 2.],
                          [2., 4., 6., 6., 4., 2.],
                          [2., 4., 6., 6., 4., 2.]])
        assert_allclose(result, expected, atol=1e-16)
        assert result.dtype == sepfir_dtype_map[dtyp]

    @pytest.mark.parametrize('dtyp',
        [np.uint8, int, np.float32, float, np.complex64, complex]
    )
    def test_strided(self, dtyp):
        a = np.array([[1, 2, 3, 3, 2, 1, 1, 2, 3],
                     [1, 2, 3, 3, 2, 1, 1, 2, 3],
                     [1, 2, 3, 3, 2, 1, 1, 2, 3],
                     [1, 2, 3, 3, 2, 1, 1, 2, 3]])
        h1, h2 = [0.5, 1, 0.5], [1]
        result_strided = signal.sepfir2d(a[:, ::2], h1, h2)
        result_contig = signal.sepfir2d(a[:, ::2].copy(), h1, h2)
        assert_allclose(result_strided, result_contig, atol=1e-15)
        assert result_strided.dtype == result_contig.dtype

    @pytest.mark.xfail(reason="XXX: filt.size > image.shape: flaky")
    def test_sepfir2d_strided_2(self):
        # XXX: this test is flaky: fails on some reruns, with
        # result[0, 1] and result[1, 1] being ~1e+224.
        np.random.seed(1234)
        filt = np.array([1.0, 2.0, 4.0, 2.0, 1.0, 3.0, 2.0])
        image = np.random.rand(4, 4)

        expected = np.array([[36.018162, 30.239061, 38.71187 , 43.878183],
                             [38.180999, 35.824583, 43.525247, 43.874945],
                             [43.269533, 40.834018, 46.757772, 44.276423],
                             [49.120928, 39.681844, 43.596067, 45.085854]])
        assert_allclose(signal.sepfir2d(image, filt, filt[::3]), expected)

    @pytest.mark.parametrize('dtyp',
        [np.uint8, int, np.float32, float, np.complex64, complex]
    )
    def test_sepfir2d_strided_3(self, dtyp):
        # NB: 'image' and 'filt' dtypes match here. Otherwise we can run into
        # unsafe casting errors for many combinations. Historically, dtype handling
        # in `sepfir2d` is a tad baroque; fixing it is an enhancement.
        filt = np.array([1, 2, 4, 2, 1, 3, 2], dtype=dtyp)
        image = np.array([[0, 3, 0, 1, 2],
                          [2, 2, 3, 3, 3],
                          [0, 1, 3, 0, 3],
                          [2, 3, 0, 1, 3],
                          [3, 3, 2, 1, 2]], dtype=dtyp)

        expected = array([[123., 101.,  91., 136., 127.],
             [133., 125., 126., 152., 160.],
             [136., 137., 150., 162., 177.],
             [133., 124., 132., 148., 147.],
             [173., 158., 152., 164., 141.]])
        result = signal.sepfir2d(image, filt, filt[::3])
        assert_allclose(result, expected, atol=1e-15)
        assert result.dtype == sepfir_dtype_map[dtyp]

        expected = array([[22., 35., 41., 31., 47.],
             [27., 39., 48., 47., 55.],
             [33., 42., 49., 53., 59.],
             [39., 44., 41., 36., 48.],
             [67., 62., 47., 34., 46.]])
        result = signal.sepfir2d(image, filt[::3], filt[::3])
        assert_allclose(result, expected, atol=1e-15)
        assert result.dtype == sepfir_dtype_map[dtyp]


def test_cspline2d():
    np.random.seed(181819142)
    image = np.random.rand(71, 73)
    signal.cspline2d(image, 8.0)


def test_qspline2d():
    np.random.seed(181819143)
    image = np.random.rand(71, 73)
    signal.qspline2d(image)


# ### Test SymIIR-related functionality

def _compute_symiirorder2_bwd_hs(k, cs, rsq, omega):
    cssq = cs * cs
    k = np.abs(k)
    rsupk = np.power(rsq, k / 2.0)

    c0 = (cssq * (1.0 + rsq) / (1.0 - rsq) /
          (1 - 2 * rsq * np.cos(2 * omega) + rsq * rsq))
    gamma = (1.0 - rsq) / (1.0 + rsq) / np.tan(omega)
    return c0 * rsupk * (np.cos(omega * k) + gamma * np.sin(omega * k))


class TestSymIIR:
    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64, np.complex64, np.complex128])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir1_ic(self, dtype, precision):
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Symmetrical initial conditions for a IIR filter of order 1 are:
        # x[0] + z1 * \sum{k = 0}^{n - 1} x[k] * z1^k

        # Check the initial condition for a low-pass filter
        # with coefficient b = 0.85 on a step signal. The initial condition is
        # a geometric series: 1 + b * \sum_{k = 0}^{n - 1} u[k] b^k.

        # Finding the initial condition corresponds to
        # 1. Computing the index n such that b**n < precision, which
        # corresponds to ceil(log(precision) / log(b))
        # 2. Computing the geometric series until n, this can be computed
        # using the partial sum formula: (1 - b**n) / (1 - b)
        # This holds due to the input being a step signal.
        b = 0.85
        n_exp = int(np.ceil(np.log(c_precision) / np.log(b)))
        expected = np.asarray([[(1 - b ** n_exp) / (1 - b)]], dtype=dtype)
        expected = 1 + b * expected

        # Create a step signal of size n + 1
        x = np.ones(n_exp + 1, dtype=dtype)
        assert_allclose(symiirorder1_ic(x, b, precision), expected,
                        atol=2e-6, rtol=2e-7)

        # Check the conditions for a exponential decreasing signal with base 2.
        # Same conditions hold, as the product of 0.5^n * 0.85^n is
        # still a geometric series
        b_d = array(b, dtype=dtype)
        expected = np.asarray(
            [[(1 - (0.5 * b_d) ** n_exp) / (1 - (0.5 * b_d))]], dtype=dtype)
        expected = 1 + b_d * expected

        # Create an exponential decreasing signal of size n + 1
        x = 2 ** -np.arange(n_exp + 1, dtype=dtype)
        assert_allclose(symiirorder1_ic(x, b, precision), expected,
                        atol=2e-6, rtol=2e-7)

    def test_symiir1_ic_fails(self):
        # Test that symiirorder1_ic fails whenever \sum_{n = 1}^{n} b^n > eps
        b = 0.85
        # Create a step signal of size 100
        x = np.ones(100, dtype=np.float64)

        # Compute the closed form for the geometrical series
        precision = 1 / (1 - b)
        pytest.raises(ValueError, symiirorder1_ic, x, b, precision)

        # Test that symiirorder1_ic fails when |z1| >= 1
        pytest.raises(ValueError, symiirorder1_ic, x, 1.0, -1)
        pytest.raises(ValueError, symiirorder1_ic, x, 2.0, -1)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64, np.complex64, np.complex128])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir1(self, dtype, precision):
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Test for a low-pass filter with c0 = 0.15 and z1 = 0.85
        # using an unit step over 200 samples.
        c0 = 0.15
        z1 = 0.85
        n = 200
        signal = np.ones(n, dtype=dtype)

        # Find the initial condition. See test_symiir1_ic for a detailed
        # explanation
        n_exp = int(np.ceil(np.log(c_precision) / np.log(z1)))
        initial = np.asarray((1 - z1 ** n_exp) / (1 - z1), dtype=dtype)
        initial = 1 + z1 * initial

        # Forward pass
        # The transfer function for the system 1 / (1 - z1 * z^-1) when
        # applied to an unit step with initial conditions y0 is
        # 1 / (1 - z1 * z^-1) * (z^-1 / (1 - z^-1) + y0)

        # Solving the inverse Z-transform for the given expression yields:
        # y[n] = y0 * z1**n * u[n] +
        #        -z1 / (1 - z1) * z1**(k - 1) * u[k - 1] +
        #        1 / (1 - z1) * u[k - 1]
        # d is the Kronecker delta function, and u is the unit step

        # y0 * z1**n * u[n]
        pos = np.arange(n, dtype=dtype)
        comp1 = initial * z1**pos

        # -z1 / (1 - z1) * z1**(k - 1) * u[k - 1]
        comp2 = np.zeros(n, dtype=dtype)
        comp2[1:] = -z1 / (1 - z1) * z1**pos[:-1]

        # 1 / (1 - z1) * u[k - 1]
        comp3 = np.zeros(n, dtype=dtype)
        comp3[1:] = 1 / (1 - z1)

        expected_fwd = comp1 + comp2 + comp3

        # Reverse condition
        sym_cond = -c0 / (z1 - 1.0) * expected_fwd[-1]

        # Backward pass
        # The transfer function for the forward result is equivalent to
        # the forward system times c0 / (1 - z1 * z).

        # Computing a closed form for the complete expression is difficult
        # The result will be computed iteratively from the difference equation
        exp_out = np.zeros(n, dtype=dtype)
        exp_out[0] = sym_cond

        for i in range(1, n):
            exp_out[i] = c0 * expected_fwd[n - 1 - i] + z1 * exp_out[i - 1]

        exp_out = exp_out[::-1]

        out = symiirorder1(signal, c0, z1, precision)
        assert_allclose(out, exp_out, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize('dtyp', [np.float32, np.float64])
    def test_symiir1_values(self, dtyp):
        rng = np.random.RandomState(1234)
        s = rng.uniform(size=16).astype(dtyp)
        res = symiirorder1(s, 0.5, 0.1)

        # values from scipy 1.9.1
        exp_res = np.array([0.14387447, 0.35166047, 0.29735238, 0.46295986, 0.45174927,
                            0.19982875, 0.20355805, 0.47378628, 0.57232247, 0.51597393,
                           0.25935107, 0.31438554, 0.41096728, 0.4190693 , 0.25812255,
                           0.33671467])
        assert res.dtype == dtyp
        atol = {np.float64: 1e-15, np.float32: 1e-7}[dtyp]
        assert_allclose(res, exp_res, atol=atol)

        # complex now
        s = s + 1j*s
        res = symiirorder1(s, 0.5, 0.1)
        assert res.dtype == np.complex64 if dtyp == np.float32 else np.complex128
        assert_allclose(res, exp_res + 1j*exp_res, atol=atol)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir2_initial_fwd(self, dtype, precision):
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        # Compute the initial conditions for a order-two symmetrical low-pass
        # filter with r = 0.5 and omega = pi / 3 for an unit step input.
        r = np.asarray(0.5, dtype=dtype)
        omega = np.asarray(np.pi / 3.0, dtype=dtype)
        cs = 1 - 2 * r * np.cos(omega) + r**2

        # The index n for the initial condition is bound from 0 to the
        # first position where sin(omega * (n + 2)) = 0 => omega * (n + 2) = pi
        # For omega = pi / 3, the maximum initial condition occurs when
        # sqrt(3) / 2 * r**n < precision.
        # => n = log(2 * sqrt(3) / 3 * precision) / log(r)
        ub = np.ceil(np.log(c_precision / np.sin(omega)) / np.log(c_precision))
        lb = np.ceil(np.pi / omega) - 2
        n_exp = min(ub, lb)

        # The forward initial condition for a filter of order two is:
        # \frac{cs}{\sin(\omega)} \sum_{n = 0}^{N - 1} {
        #    r^(n + 1) \sin{\omega(n + 2)}} + cs
        # The closed expression for this sum is:
        # s[n] = 2 * r * np.cos(omega) -
        #        r**2 - r**(n + 2) * np.sin(omega * (n + 3)) / np.sin(omega) +
        #        r**(n + 3) * np.sin(omega * (n + 2)) / np.sin(omega) + cs
        fwd_initial_1 = (
            cs +
            2 * r * np.cos(omega) -
            r**2 -
            r**(n_exp + 2) * np.sin(omega * (n_exp + 3)) / np.sin(omega) +
            r**(n_exp + 3) * np.sin(omega * (n_exp + 2)) / np.sin(omega))

        # The second initial condition is given by
        # s[n] = 1 / np.sin(omega) * (
        #        r**2 * np.sin(3 * omega) -
        #        r**3 * np.sin(2 * omega) -
        #        r**(n + 3) * np.sin(omega * (n + 4)) +
        #        r**(n + 4) * np.sin(omega * (n + 3)))
        ub = np.ceil(np.log(c_precision / np.sin(omega)) / np.log(c_precision))
        lb = np.ceil(np.pi / omega) - 3
        n_exp = min(ub, lb)

        fwd_initial_2 = (
            cs + cs * 2 * r * np.cos(omega) +
            (r**2 * np.sin(3 * omega) -
             r**3 * np.sin(2 * omega) -
             r**(n_exp + 3) * np.sin(omega * (n_exp + 4)) +
             r**(n_exp + 4) * np.sin(omega * (n_exp + 3))) / np.sin(omega))

        expected = np.r_[fwd_initial_1, fwd_initial_2][None, :]

        n = 100
        signal = np.ones(n, dtype=dtype)

        out = symiirorder2_ic_fwd(signal, r, omega, precision)
        assert_allclose(out, expected, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir2_initial_bwd(self, dtype, precision):
        c_precision = precision
        if precision <= 0.0 or precision > 1.0:
            if dtype in {np.float32, np.complex64}:
                c_precision = 1e-6
            else:
                c_precision = 1e-11

        r = np.asarray(0.5, dtype=dtype)
        omega = np.asarray(np.pi / 3.0, dtype=dtype)
        cs = 1 - 2 * r * np.cos(omega) + r * r
        a2 = 2 * r * np.cos(omega)
        a3 = -r * r

        n = 100
        signal = np.ones(n, dtype=dtype)

        # Compute initial forward conditions
        ic = symiirorder2_ic_fwd(signal, r, omega, precision)
        out = np.zeros(n + 2, dtype=dtype)
        out[:2] = ic[0]

        # Apply the forward system cs / (1 - a2 * z^-1 - a3 * z^-2))
        for i in range(2, n + 2):
            out[i] = cs * signal[i - 2] + a2 * out[i - 1] + a3 * out[i - 2]

        # Find the backward initial conditions
        ic2 = np.zeros(2, dtype=dtype)
        idx = np.arange(n)

        diff = (_compute_symiirorder2_bwd_hs(idx, cs, r * r, omega) +
                _compute_symiirorder2_bwd_hs(idx + 1, cs, r * r, omega))
        ic2_0_all = np.cumsum(diff * out[:1:-1])
        pos = np.where(diff ** 2 < c_precision)[0]
        ic2[0] = ic2_0_all[pos[0]]

        diff = (_compute_symiirorder2_bwd_hs(idx - 1, cs, r * r, omega) +
                _compute_symiirorder2_bwd_hs(idx + 2, cs, r * r, omega))
        ic2_1_all = np.cumsum(diff * out[:1:-1])
        pos = np.where(diff ** 2 < c_precision)[0]
        ic2[1] = ic2_1_all[pos[0]]

        out_ic = symiirorder2_ic_bwd(out, r, omega, precision)[0]
        assert_allclose(out_ic, ic2, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize(
        'dtype', [np.float32, np.float64])
    @pytest.mark.parametrize('precision', [-1.0, 0.7, 0.5, 0.25, 0.0075])
    def test_symiir2(self, dtype, precision):
        r = np.asarray(0.5, dtype=dtype)
        omega = np.asarray(np.pi / 3.0, dtype=dtype)
        cs = 1 - 2 * r * np.cos(omega) + r * r
        a2 = 2 * r * np.cos(omega)
        a3 = -r * r

        n = 100
        signal = np.ones(n, dtype=dtype)

        # Compute initial forward conditions
        ic = symiirorder2_ic_fwd(signal, r, omega, precision)
        out1 = np.zeros(n + 2, dtype=dtype)
        out1[:2] = ic[0]

        # Apply the forward system cs / (1 - a2 * z^-1 - a3 * z^-2))
        for i in range(2, n + 2):
            out1[i] = cs * signal[i - 2] + a2 * out1[i - 1] + a3 * out1[i - 2]

        # Find the backward initial conditions
        ic2 = symiirorder2_ic_bwd(out1, r, omega, precision)[0]

        # Apply the system cs / (1 - a2 * z - a3 * z^2)) in backwards
        exp = np.empty(n, dtype=dtype)
        exp[-2:] = ic2[::-1]

        for i in range(n - 3, -1, -1):
            exp[i] = cs * out1[i] + a2 * exp[i + 1] + a3 * exp[i + 2]

        out = symiirorder2(signal, r, omega, precision)
        assert_allclose(out, exp, atol=4e-6, rtol=6e-7)

    @pytest.mark.parametrize('dtyp', [np.float32, np.float64])
    def test_symiir2_values(self, dtyp):
        rng = np.random.RandomState(1234)
        s = rng.uniform(size=16).astype(dtyp)
        res = symiirorder2(s, 0.1, 0.1, precision=1e-10)

        # values from scipy 1.9.1
        exp_res = np.array([0.26572609, 0.53408018, 0.51032696, 0.72115829, 0.69486885,
           0.3649055 , 0.37349478, 0.74165032, 0.89718521, 0.80582483,
           0.46758053, 0.51898709, 0.65025605, 0.65394321, 0.45273595,
           0.53539183])

        assert res.dtype == dtyp
        # The values in SciPy 1.14 agree with those in SciPy 1.9.1 to this
        # accuracy only. Implementation differences are twofold:
        # 1. boundary conditions are computed differently
        # 2. the filter itself uses sosfilt instead of a hardcoded iteration
        # The boundary conditions seem are tested separately (see
        # test_symiir2_initial_{fwd,bwd} above, so the difference is likely
        # due to a different way roundoff errors accumulate in the filter.
        # In that respect, sosfilt is likely doing a better job.
        atol = {np.float64: 2e-6, np.float32: 2e-6}[dtyp]
        assert_allclose(res, exp_res, atol=atol)

        # complex now
        s = s + 1j*s
        with pytest.raises(TypeError):
            res = symiirorder2(s, 0.5, 0.1)

    def test_symiir1_integer_input(self):
        s = np.where(np.arange(100) % 2, -1, 1)
        expected = symiirorder1(s.astype(float), 0.5, 0.5)
        out = symiirorder1(s, 0.5, 0.5)
        assert_allclose(out, expected)

    def test_symiir2_integer_input(self):
        s = np.where(np.arange(100) % 2, -1, 1)
        expected = symiirorder2(s.astype(float), 0.5, np.pi / 3.0)
        out = symiirorder2(s, 0.5, np.pi / 3.0)
        assert_allclose(out, expected)
