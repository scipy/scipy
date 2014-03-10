import numpy as np

from numpy.testing import TestCase, run_module_suite, assert_allclose

from scipy.signal.hankel import HankelTransform

def run_test(transform, func, true, B):

    HB = transform(func, B)
    TH = true(B)

    assert_allclose(HB, TH)

class TestHankelTransforms(TestCase):

    def test_all(self, name=None):
        # With no arguments, test the 801-point filter due to Anderson (1979)
        HT = HankelTransform(name)

        # First, generate the vector of B values.
        B = np.r_[0.0001, 0.001, 0.01:0.05:0.01,  0.05:5.01:0.05 ]

        #  Test 1, int(g*exp(-g*g)*J0(g*b),g=0..infinity)=exp(-b^2/4)/2.
        run_test(HT.hankel0, lambda g: g * np.exp(-g**2),
                          lambda b: np.exp(-b**2/4)/2, B)

        #  Test 2, int(g*g*exp(-g*g)*J1(g*b),g=0..infinity)=b*exp(-b*b/4)/4
        run_test(HT.hankel1, lambda g: g**2 * np.exp(-g**2),
                          lambda b: b * np.exp(-b**2/4)/4, B)

        #  Test 3, int(exp(-g)*J1(g*B),g=0..infinity)=(sqrt(1+B^2)-1)/(B*sqrt(1+B^2))
        run_test(HT.hankel1, lambda g: np.exp(-g),
                          lambda b: (np.sqrt(1+b**2)-1)/(b*np.sqrt(1+b**2)), B);

        #  Test 4, int(g*exp(-2*g)*J1(g*B),g=0..infinity)=
        run_test(HT.hankel1, lambda g: g * np.exp(-2*g),
                          lambda b: b / (4+b**2)**(3./2.), B)

        #  Test 5, int(exp(-2*g)*J0(g*B),g=0..infinity)=1/(sqrt(4+B^2))
        run_test(HT.hankel0, lambda g: np.exp(-2*g),
                          lambda b: 1. / np.sqrt(4+b**2), B)



if __name__ == "__main__":
    run_module_suite()
