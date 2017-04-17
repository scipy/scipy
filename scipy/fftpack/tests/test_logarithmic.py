# Created by Dieter Werthmüller, January 2017

from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.fftpack import fftlog, fftlogargs
from numpy.testing import assert_allclose, assert_raises, run_module_suite

__usage__ = """
Build fftpack:
  python setup.py install
Run tests if scipy is installed:
  python -c 'from scipy import fftpack; fftpack.test()'
Run tests if fftpack is not installed:
  python tests/test_logarithmic.py
"""


class TestFFTLog():

    def test_fct1_sine(self):

        # FFTLog parameters
        w, t, kr, rk = fftlogargs(512, 0.01, 0, 'sine', -0.2, 1, 1)

        # Analytical solution
        fw = np.sqrt(np.pi/2)/2*np.exp(-w)*np.sin(w)  # Frequency domain
        ft = t/(4 + t**4)                             # Time domain

        # FFTLog
        fftl = fftlog(fw, 0.01, 'sine', -0.2, kr, rk)

        # Check
        assert_allclose(fftl[120:-150], ft[120:-150], rtol=1e-3)

    def test_fct2_sine(self):

        # FFTLog parameters
        w, t, kr, rk = fftlogargs(512, 0.01, 0, 'sine', 0.5, 1, 1)

        # Analytical solution
        fw = np.sqrt(np.pi/(2*np.ones(w.size)))  # Frequency domain
        ft = 1/t                                 # Time domain

        # FFTLog
        fftl = fftlog(fw, 0.01, 'sine', 0.5, kr, rk)

        # Check
        assert_allclose(fftl, ft)

    def test_fct3_cosine(self):

        # FFTLog parameters
        w, t, kr, rk = fftlogargs(512, 0.01, 0, 'cosine', 0, 1, 1)

        # Analytical solution
        fw = np.sqrt(1/w)  # Frequency domain
        ft = 1/np.sqrt(t)  # Time domain

        # FFTLog
        fftl = fftlog(fw, 0.01, 'cosine', 0, kr, rk)

        # Check
        assert_allclose(fftl, ft)

    def test_empty(self):
        with assert_raises(ValueError):
            fftlog([], 0.01)


class TestFFTLogargs():
    # Default values in fftlogargs:
    # - dlogr : 0.01
    # - logrc : 0.0
    # - mu : 'sine'
    # - q : 0
    # - kr : 1
    # - kropt : 0

    def test_kr2(self):  # Test 1: even n, kr = 2
        inppts = [0.96605088, 0.98855309, 1.01157945, 1.03514217]
        outpts = [1.93210176, 1.97710619, 2.02315891, 2.07028433]
        kr = 2.0
        rk = 0.5
        out = fftlogargs(n=4, kr=2)
        assert_allclose(inppts, out[0])
        assert_allclose(outpts, out[1])
        assert_allclose(kr, out[2])
        assert_allclose(rk, out[3])

    def test_kropt1(self):  # Test 2: even n, kropt = 1
        inppts = [0.96605088, 0.98855309, 1.01157945, 1.03514217]
        outpts = [0.97306236, 0.9957279, 1.01892138, 1.04265511]
        kr = 1.0072578812188107
        rk = 0.99279441605358465
        out = fftlogargs(n=4, kropt=1)
        assert_allclose(inppts, out[0])
        assert_allclose(outpts, out[1])
        assert_allclose(kr, out[2])
        assert_allclose(rk, out[3])

    def test_cosine(self):  # Test 3: odd n, kr = pi, mu = 'cosine'
        inppts = [0.95499259, 0.97723722, 1., 1.02329299, 1.04712855]
        outpts = [3.00019769, 3.07008127, 3.14159265, 3.21476975, 3.28965135]
        kr = 3.141592653589793
        rk = 0.31830988618379069
        out = fftlogargs(5, mu='cosine', kr=np.pi)
        assert_allclose(inppts, out[0])
        assert_allclose(outpts, out[1])
        assert_allclose(kr, out[2])
        assert_allclose(rk, out[3])

    def test_logrc1(self):  # Test 4: odd n, logrc = 1, kropt=1
        inppts = [9.54992586, 9.77237221, 10., 10.23292992, 10.47128548]
        outpts = [0.09619238, 0.09843299, 0.10072579, 0.10307199, 0.10547285]
        kr = 1.0072578812188107
        rk = 99.279441605358485
        out = fftlogargs(n=5, logrc=1, kropt=1)
        assert_allclose(inppts, out[0])
        assert_allclose(outpts, out[1])
        assert_allclose(kr, out[2])
        assert_allclose(rk, out[3])

    def test_falsemu(self):  # Test 5: wrong mu
        with assert_raises(ValueError):
            fftlogargs(5, mu=1)

if __name__ == "__main__":
    run_module_suite()
