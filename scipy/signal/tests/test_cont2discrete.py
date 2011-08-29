import numpy as np
from numpy.testing import TestCase, run_module_suite, \
                          assert_array_almost_equal, assert_almost_equal

from scipy.signal import cont2discrete as c2d

# Author: Jeffrey Armstrong <jeff@approximatrix.com>
# March 29, 2011

class TestC2D(TestCase):

    def test_zoh(self):

        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0,], [0.0,], [-0.33,]])
        
        ad_truth = 1.648721270700128 * np.eye(2)
        bd_truth = 0.324360635350064 * np.ones((2, 1))
        # c and d in discrete should be equal to their continuous counterparts
        dt_requested = 0.5

        ad, bd, cd, dd, dt = c2d((ac, bc, cc, dc), dt_requested, method='zoh')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cc, cd)
        assert_array_almost_equal(dc, dd)
        
        assert_almost_equal(dt_requested, dt)

    def test_gbt(self):

        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0,], [0.0,], [-0.33,]])
        
        dt_requested = 0.5
        alpha = 1.0/3.0

        ad_truth = 1.6 * np.eye(2)
        bd_truth = 0.3 * np.ones((2, 1))
        cd_truth = np.array([[0.9, 1.2],
                             [1.2, 1.2],
                             [1.2, 0.3]])
        dd_truth = np.array([[0.175,],
                             [0.2,],
                             [-0.205,]])

        ad, bd, cd, dd, dt = c2d((ac, bc, cc, dc), dt_requested,
                                 method='gbt',alpha=alpha)

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)

    def test_euler(self):
    
        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0,], [0.0,], [-0.33,]])

        dt_requested = 0.5

        ad_truth = 1.5 * np.eye(2)
        bd_truth = 0.25 * np.ones((2, 1))
        cd_truth = np.array([[0.75, 1.0],
                             [1.0, 1.0],
                             [1.0, 0.25]])
        dd_truth = dc

        ad, bd, cd, dd, dt = c2d((ac, bc, cc, dc), dt_requested, 
                                 method='euler')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)
        
        assert_almost_equal(dt_requested, dt)


    def test_backward_diff(self):

        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0,], [0.0,], [-0.33,]])

        dt_requested = 0.5

        ad_truth = 2.0 * np.eye(2)
        bd_truth = 0.5 * np.ones((2, 1))
        cd_truth = np.array([[1.5, 2.0],
                             [2.0, 2.0],
                             [2.0, 0.5]])
        dd_truth = np.array([[0.875,],
                             [1.0,],
                             [0.295,]])

        ad, bd, cd, dd, dt = c2d((ac, bc, cc, dc), dt_requested,
                                 method='backward_diff')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)

    def test_bilinear(self):

        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0,], [0.0,], [-0.33,]])

        dt_requested = 0.5

        ad_truth = (5.0 / 3.0) * np.eye(2)
        bd_truth = (1.0 / 3.0) * np.ones((2, 1))
        cd_truth = np.array([[1.0, 4.0 / 3.0],
                             [4.0 / 3.0, 4.0 / 3.0],
                             [4.0 / 3.0, 1.0 / 3.0]])
        dd_truth = np.array([[0.291666666666667,],
                             [1.0 / 3.0,],
                             [-0.121666666666667,]])

        ad, bd, cd, dd, dt = c2d((ac, bc, cc, dc), dt_requested, 
                                 method='bilinear')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)
        
        assert_almost_equal(dt_requested, dt)

        # Same continuous system again, but change sampling rate
        
        ad_truth = 1.4 * np.eye(2)
        bd_truth = 0.2 * np.ones((2, 1))
        cd_truth = np.array([[0.9, 1.2], [1.2, 1.2], [1.2, 0.3]])
        dd_truth = np.array([[0.175,], [0.2,], [-0.205,]])

        dt_requested = 1.0 / 3.0

        ad, bd, cd, dd, dt = c2d((ac, bc, cc, dc), dt_requested, 
                                 method='bilinear')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)
        
        assert_almost_equal(dt_requested, dt)

    def test_transferfunction(self):
    
        numc = np.array([0.25, 0.25, 0.5])
        denc = np.array([0.75, 0.75, 1.0])

        numd = np.array([[1.0 / 3.0, -0.427419169438754, 0.221654141101125],])
        dend = np.array([1.0, -1.351394049721225, 0.606530659712634])

        dt_requested = 0.5
    
        num, den, dt = c2d((numc, denc), dt_requested, method='zoh')

        assert_array_almost_equal(numd, num)
        assert_array_almost_equal(dend, den)
        
        assert_almost_equal(dt_requested, dt)
        
    def test_zerospolesgain(self):
        
        zeros_c = np.array([0.5, -0.5])
        poles_c = np.array([1.j / np.sqrt(2), -1.j / np.sqrt(2)])
        k_c = 1.0
        
        zeros_d = [1.23371727305860, 0.735356894461267]
        polls_d = [0.938148335039729 + 0.346233593780536j,
                   0.938148335039729 - 0.346233593780536j]
        k_d = 1.0
        
        dt_requested = 0.5
        
        zeros, poles, k, dt = c2d((zeros_c, poles_c, k_c), dt_requested, 
                                  method='zoh')
                                
        assert_array_almost_equal(zeros_d, zeros)
        assert_array_almost_equal(polls_d, poles)
        assert_almost_equal(k_d, k)
        assert_almost_equal(dt_requested, dt)

if __name__ == "__main__":
    run_module_suite()
