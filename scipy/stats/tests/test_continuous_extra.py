# contains additional tests for continuous distributions
#
# NOTE: one test, _est_cont_skip, that is renamed so that nose doesn't
#       run it,
#       6 distributions return nan for entropy
#       truncnorm fails by design for private method _ppf test


import numpy.testing as npt
import numpy as np
import nose

from scipy import stats

from test_continuous_basic import distcont

DECIMAL = 5

@npt.dec.slow
def test_cont_extra():
    for distname, arg in distcont[:]:
        distfn = getattr(stats, distname)

        yield check_ppf_limits, distfn, arg, distname + \
              ' ppf limit test'
        yield check_isf_limits, distfn, arg, distname + \
              ' isf limit test'
        yield check_loc_scale, distfn, arg, distname + \
              ' loc, scale test'

@npt.dec.slow
def _est_cont_skip():
    for distname, arg in distcont:
        distfn = getattr(stats, distname)        
        #entropy test checks only for isnan, currently 6 isnan left
        yield check_entropy, distfn, arg, distname + \
              ' entropy nan test'
        # _ppf test has 1 failure be design
        yield check_ppf_private, distfn, arg, distname + \
              ' _ppf private test'

def test_540_567():
    # test for nan returned in tickets 540, 567
    npt.assert_almost_equal(stats.norm.cdf(-1.7624320982),0.03899815971089126,
                            decimal=10, err_msg = 'test_540_567')
    npt.assert_almost_equal(stats.norm.cdf(-1.7624320983),0.038998159702449846,
                            decimal=10, err_msg = 'test_540_567')
    npt.assert_almost_equal(stats.norm.cdf(1.38629436112, loc=0.950273420309,
                            scale=0.204423758009),0.98353464004309321,
                            decimal=10, err_msg = 'test_540_567')


def check_ppf_limits(distfn,arg,msg):
    below,low,upp,above = distfn.ppf([-1,0,1,2], *arg)
    #print distfn.name, distfn.a, low, distfn.b, upp
    #print distfn.name,below,low,upp,above
    assert_equal_inf_nan(distfn.a,low, msg + 'ppf lower bound')
    assert_equal_inf_nan(distfn.b,upp, msg + 'ppf upper bound')
    assert np.isnan(below), msg + 'ppf out of bounds - below'
    assert np.isnan(above), msg + 'ppf out of bounds - above'

def check_ppf_private(distfn,arg,msg):
    #fails by design for trunk norm self.nb not defined
    ppfs = distfn._ppf(np.array([0.1,0.5,0.9]), *arg)
    assert not np.any(np.isnan(ppfs)), msg + 'ppf private is nan'


def check_isf_limits(distfn,arg,msg):
    below,low,upp,above = distfn.isf([-1,0,1,2], *arg)
    #print distfn.name, distfn.a, low, distfn.b, upp
    #print distfn.name,below,low,upp,above
    assert_equal_inf_nan(distfn.a,upp, msg + 'isf lower bound')
    assert_equal_inf_nan(distfn.b,low, msg + 'isf upper bound')
    assert np.isnan(below), msg + 'isf out of bounds - below'
    assert np.isnan(above), msg + 'isf out of bounds - above'


def check_loc_scale(distfn,arg,msg):
    m,v = distfn.stats(*arg)
    loc, scale = 10.0, 10.0
    mt,vt = distfn.stats(loc=loc, scale=scale, *arg)
    assert_equal_inf_nan(m*scale+loc,mt,msg + 'mean')
    assert_equal_inf_nan(v*scale*scale,vt,msg + 'var')

def check_entropy(distfn,arg,msg):
    ent = distfn.entropy(*arg)
    #print 'Entropy =', ent
    assert not np.isnan(ent), msg + 'test Entropy is nan'\

def assert_equal_inf_nan(v1,v2,msg):
    assert not np.isnan(v1)
    if not np.isinf(v1):
        npt.assert_almost_equal(v1, v2, decimal=DECIMAL, err_msg = msg + \
                                   ' - finite')
    else:
        assert np.isinf(v2) or np.isnan(v2), \
               msg + ' - infinite, v2=%s' % str(v2)

if __name__ == "__main__":
    import nose
    #nose.run(argv=['', __file__])
    nose.runmodule(argv=[__file__,'-s'], exit=False)

