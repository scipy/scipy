from pytest import raises

from scipy._lib._finite_differences import _derivative

class TestDerivative:
    def test_not_enough_points(self):
        with raises(ValueError, match="derivative order \\+ 1."):
            _derivative(lambda x: x, 0, order=5, n=3)
    
    def test_odd_n(self):
        with raises(ValueError, match="must be odd."):
            _derivative(lambda x: x,0, n=2)