"""
Test Cython optimize zeros API functions: ``bisect``, ``ridder``, ``brenth``,
and ``brentq``, in `scipy.optimize.cython_optimize`, by finding the roots of a
3rd order polynomial given a sequence of constant terms, ``a0``, and fixed 1st,
2nd, and 3rd order terms in ``args``.

.. math::

    f(x, a0, args) =  ((args[2]*x + args[1])*x + args[0])*x + a0

The 3rd order polynomial function is written in Cython and called in a Python
wrapper named after the zero function. See the private ``_zeros`` Cython module
in `scipy.optimize.cython_optimze` for more information. 
"""

from __future__ import division, print_function, absolute_import
import numpy.testing as npt
from scipy.optimize.cython_optimize import _zeros

# constants
A0 = tuple(-2.0 - x/10.0 for x in range(10))  # constant term
ARGS = (0.0, 0.0, 1.0)  # 1st, 2nd, and 3rd order terms
XLO, XHI = 0.0, 2.0  # first and second bounds of zeros functions
# absolute and relative tolerances and max iterations for zeros functions
XTOL, RTOL, MITR = 0.001, 0.001, 10


EXPECTED_BISECT = [
    1.259765625,
    1.279296875,
    1.298828125,
    1.318359375,
    1.337890625,
    1.357421875,
    1.376953125,
    1.392578125,
    1.408203125,
    1.427734375]


# test bisect
def test_bisect():
    npt.assert_allclose(
        EXPECTED_BISECT,
        list(
            _zeros.loop_example('bisect', A0, ARGS, XLO, XHI, XTOL, RTOL, MITR)
        )
    )


EXPECTED_RIDDER = [
    1.2588478785767947,
    1.2795040615075954,
    1.299514441316524,
    1.318927124420269,
    1.3377847289304623,
    1.356125211864335,
    1.373982543571637,
    1.3913872624129802,
    1.408366934614972,
    1.4249465383291897]


# test ridder
def test_ridder():
    npt.assert_allclose(
        EXPECTED_RIDDER,
        list(
            _zeros.loop_example('ridder', A0, ARGS, XLO, XHI, XTOL, RTOL, MITR)
        )
    )


EXPECTED_BRENT = [
    1.259872799904563,
    1.28042866862737,
    1.3003083443276644,
    1.3195775152132223,
    1.3382912269178606,
    1.3576742112124918,
    1.375418396530433,
    1.3927268676453197,
    1.4096289546536032,
    1.4261502005552444]


# test brenth
def test_brenth():
    npt.assert_allclose(
        EXPECTED_BRENT,
        list(
            _zeros.loop_example('brenth', A0, ARGS, XLO, XHI, XTOL, RTOL, MITR)
        ),
        rtol=RTOL, atol=XTOL
    )


# test brentq
def test_brentq():
    npt.assert_allclose(
        EXPECTED_BRENT,
        list(
            _zeros.loop_example('brentq', A0, ARGS, XLO, XHI, XTOL, RTOL, MITR)
        ),
        rtol=RTOL, atol=XTOL
    )
                       

# test brentq with full output
def test_brentq_full_output():
    output = _zeros.full_output_example(
        (A0[0],) + ARGS, XLO, XHI, XTOL, RTOL, MITR)
    npt.assert_allclose(EXPECTED_BRENT[0], output['root'])
    npt.assert_equal(6, output['iterations'])
    npt.assert_equal(7, output['funcalls'])
    npt.assert_equal(0, output['error_num'])
