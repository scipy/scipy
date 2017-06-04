######### Unit test section #########
from numpy import *
from numpy.testing import *

from scipy.optimize._lp import lp_solve

def test_lp(prt=False):
    m1 = 20
    m2 = 50
    probs = [
        {
            'A': array([
                [2.,  5., 3., -1.,  0.,  0.],
                [3., 2.5, 8.,  0., -1.,  0.],
                [8.,10.,  4.,  0.,  0., -1.]]),
            'b': array([185., 155., 600.]),
            'c': array([4., 8., 3., 0., 0., 0.]),
            'result': [
                    array([ 66.25, 0., 17.5, 0., 183.75, 0.]),
                    317.5,
                    True,
                    True,
                    array([2, 0, 4])            
                ]
        },
        {        
            'A': array([
                [-1., -1., -1.,  0.,  0.,  0.],
                [ 0.,  0.,  0.,  1.,  1.,  1.],
                [ 1.,  0.,  0.,  1.,  0.,  0.],
                [ 0.,  1.,  0.,  0.,  1.,  0.],
                [ 0.,  0.,  1.,  0.,  0.,  1.]]),
            'b': array([-0.5, 0.4, 0.3, 0.3, 0.3]),
            'c': array([2.8, 6.3, 10.8, -2.8, -6.3, -10.8]),
            'result': [
                    array([0.3, 0.2, 0.0, 0.0, 0.1, 0.3]),
                    -1.77,
                    True,
                    True,
                    array([1, 7, 0, 4, 5])            
                ]
        },
        {   # with degeneracy
            'A': array([[cos(2*pi*i/(m1+1))-1., sin(2*pi*i/(m1+1))]
                        for i in xrange(1,m1+1)]).T,
            'b': zeros(2).T,
            'c': -ones(m1).T,
            'result': [
                    zeros(m1),
                    0.,
                    True,
                    True,
                    array([0,19])
                ]
            
        },
        {   # with unboundedness (0 is a member of the convex hull
            # of these vectors)
            'A': array([[cos(2*pi*i/(m2+1))-1., sin(2*pi*i/(m2+1))]
                        for i in xrange(0,m2)]).T,
            'b': zeros(2).T,
            'c': -ones(m2).T,
            'result': [
                    None,   # unchecked when unbounded
                    -Inf,   # unchecked when unbounded
                    False,
                    True,
                    array([2, 49])
                ]
            
        }, 
        {   # Unsolvable
            'A': array([[cos(2*pi*i/(m2+1))-1., sin(2*pi*i/(m2+1))]
                        for i in xrange(0,m2)]).T,
            'b': ones(2).T,
            'c': -ones(m2).T,
            'result': [
                    None,   # unchecked when unsolvable
                    None,   # unchecked when unsolvable
                    None,   # unchecked when unsolvable
                    False,
                    array([50, 1])
                ]
            
        }, # add other test cases here...
    ]


    for prob in probs:
        lpsol = lp_solve(prob['c'],prob['A'],prob['b'])
        optx = lpsol.x
        zmin = lpsol.fun
        bounded = lpsol.is_bounded
        solvable = lpsol.is_solvable
        basis = lpsol.basis
        if prt:
            print "A:\n",prob['A']
            print "b:",prob['b']
            print "c:",prob['c']
            print " ---->"
            print "optx:",optx
            print "zmin:",zmin
            print "bounded:",bounded
            print "solvable:",solvable
            print "basis:",basis
            print "-------------------------------------------"
        else:
            expected_res = prob['result']
            assert_equal(solvable, expected_res[3], err_msg=repr(prob))
            assert_equal(basis, expected_res[4], err_msg=repr(prob))
            if solvable:
                assert_equal(bounded, expected_res[2], err_msg=repr(prob))
                if bounded:
                    assert_almost_equal(optx, expected_res[0],
                                        err_msg=repr(prob))
                assert_almost_equal(zmin, expected_res[1], err_msg=repr(prob))
                # when unbounded zmin == -Inf, but -Inf != -Inf
                # so we won't check it...
