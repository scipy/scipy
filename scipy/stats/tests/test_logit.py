import numpy as np
from numpy.testing import *
from scipy.stats import logit, expit

class TestLogit(TestCase):
    
    def check_logit_out(self, dtype, expected_out):
        a = np.linspace(0,1,10)
        a = np.array(a, dtype=dtype)
        actual_out = logit(a)
        assert_almost_equal(actual_out, expected_out)
        assert_equal(actual_out.dtype, np.dtype(dtype))

    def test_float16(self):
        expected = np.array([-np.inf, -2.08007812, 
                            -1.25292969, -0.69335938, 
                            -0.22363281, 0.22363281,  
                            0.69238281,  1.25292969,  
                            2.078125  , np.inf], dtype='float16')
        self.check_logit_out('f2', expected)

    def test_float32(self):
        expected = np.array([-np.inf, -2.07944155, 
                            -1.25276291, -0.69314718, 
                            -0.22314353, 0.22314365,  
                            0.6931473 ,  1.25276303,  
                            2.07944155, np.inf], dtype='float32')
        self.check_logit_out('f4', expected)

    def test_float64(self):
        expected = np.array([-np.inf, -2.07944154, 
                            -1.25276297, -0.69314718, 
                            -0.22314355, 0.22314355,  
                            0.69314718,  1.25276297,  
                            2.07944154,         np.inf])
        self.check_logit_out('f8', expected)

class TestExpit(TestCase):
    def check_expit_out(self, dtype, expected_out):
        a = np.linspace(-4,4,10)
        a = np.array(a, dtype=dtype)
        actual_out = expit(a)
        assert_almost_equal(actual_out, expected_out)
        assert_equal(actual_out.dtype, np.dtype(dtype))

    def test_float16(self):
        expected = np.array([ 0.01799011,  0.04263306,  
                            0.09771729,  0.20861816,  
                            0.390625  , 0.609375  ,  
                            0.79150391,  0.90234375,  
                            0.95751953,  0.98193359], dtype='float16')
        self.check_expit_out('f2', expected)

    def test_float32(self):
        expected = np.array([ 0.01798621,  0.04265125,  
                            0.09777259,  0.20860852,  
                            0.39068246, 0.60931754,  
                            0.79139149,  0.9022274 ,  
                            0.95734876,  0.98201376], dtype='float32')
        self.check_expit_out('f4',expected)

    def test_float64(self):
        expected = np.array([ 0.01798621,  0.04265125,  
                            0.0977726 ,  0.20860853,  
                            0.39068246, 0.60931754,  
                            0.79139147,  0.9022274 ,  
                            0.95734875,  0.98201379])
        self.check_expit_out('f8', expected)


