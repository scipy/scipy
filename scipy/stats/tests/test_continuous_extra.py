import numpy.testing as npt
import numpy as np
import nose

from scipy import stats

from test_continuous_basic import distcont

DECIMAL = 5

def test_cont_extra():
    for distname, arg in distcont[:]:
        distfn = getattr(stats, distname)
##        rvs = distfn.rvs(size=1000,*arg)
##        sm = rvs.mean()
##        sv = rvs.var()
##        skurt = stats.kurtosis(rvs)
##        sskew = stats.skew(rvs)
        
        yield check_ppf_limits, distfn, arg, distname + \
              ' ppf limit test'
        yield check_isf_limits, distfn, arg, distname + \
              ' isf limit test'
        yield check_loc_scale, distfn, arg, distname + \
              ' loc, scale test'
        #entropy test checks only for isnan, currently 6 isnan
##        yield check_entropy, distfn, arg, distname + \
##              ' entropy nan test'


def check_ppf_limits(distfn,arg,msg):
    below,low,upp,above = distfn.ppf([-1,0,1,2], *arg)
    #print distfn.name, distfn.a, low, distfn.b, upp
    #print distfn.name,below,low,upp,above
    assert_equal_inf_nan(distfn.a,low, msg + 'ppf lower bound')
    assert_equal_inf_nan(distfn.b,upp, msg + 'ppf upper bound')
    assert np.isnan(below), msg + 'ppf out of bounds - below'
    assert np.isnan(above), msg + 'ppf out of bounds - above'

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
        npt.assert_almost_equal(v1, v2, decimal=DECIMAL, err_msg= msg + \
                                   ' - finite')
    else:
        assert np.isinf(v2) or np.isnan(v2), \
               msg + ' - infinite, v2=%s' % str(v2)

if __name__ == "__main__":
    import nose
    #nose.run(argv=['', __file__])
    #print distcont[:5]
    #test_cont_extra()
    nose.runmodule(argv=[__file__,'-s'], exit=False)

