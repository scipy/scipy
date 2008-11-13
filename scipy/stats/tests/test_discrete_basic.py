import numpy.testing as npt
import numpy as np
import nose

from scipy import stats

distdiscrete = [
    ['bernoulli',(0.3,)],
    ['binom',    (5, 0.4)],
    ['boltzmann',(1.4, 19)],
    ['dlaplace', (0.8,)], #0.5
    ['geom',     (0.5,)],
    ['hypergeom',(30, 12, 6)],
    ['logser',   (0.6,)],
    ['nbinom',   (5, 0.5)],
    ['planck',   (0.51,)],   #4.1
    ['poisson',  (0.6,)],
    ['randint',  (7, 31)],
    ['zipf',     (4,)] ]   # arg=4 is ok,
                           # Zipf broken for arg = 2, e.g. weird .stats

def test_discrete_basic():
    for distname, arg in distdiscrete:
        distfn = getattr(stats,distname)
        #assert stats.dlaplace.rvs(0.8) != None
        rvs = distfn.rvs(size=10000,*arg)
        m,v = distfn.stats(*arg)
        #yield npt.assert_almost_equal(rvs.mean(), m, decimal=4,err_msg='mean')
        #yield npt.assert_almost_equal, rvs.mean(), m, 2, 'mean' # does not work
        yield check_sample_mean, rvs.mean(), m, distname + 'sample mean test'
        yield check_sample_mean, rvs.var(), v, distname + 'sample var test'
        yield check_cdf_ppf, distfn, arg
        yield check_pmf_cdf, distfn, arg
        yield check_oth, distfn, arg

def test_discrete_private():
    #testing private methods mostly for debugging
    #   some tests might fail by design,
    #   e.g. incorrect definition of distfn.a and distfn.b
    for distname, arg in distdiscrete:
        distfn = getattr(stats,distname)
        rvs = distfn.rvs(size=10000,*arg)
        m,v = distfn.stats(*arg)

        yield check_ppf_ppf, distfn, arg
        yield check_cdf_ppf_private, distfn, arg
        yield check_generic_moment, distfn, arg, m, 1, 3   # last is decimal
        yield check_generic_moment, distfn, arg, v+m*m, 2, 3 # last is decimal
        yield check_moment_frozen, distfn, arg, m, 1, 3   # last is decimal
        yield check_moment_frozen, distfn, arg, v+m*m, 2, 3 # last is decimal


def check_sample_mean(sm,m,msg):
    if m < np.inf:
        npt.assert_almost_equal(sm, m, decimal=0, err_msg=msg + \
                                ' - finite moment')
    else:
        assert sm > 10000, 'infinite moment, sm = ' + str(sm)

def check_sample_var(sm,m):
    npt.assert_almost_equal(sm, m, decimal=0, err_msg= 'var')

def check_cdf_ppf(distfn,arg):
    ppf05 = distfn.ppf(0.5,*arg)
    cdf05 = distfn.cdf(ppf05,*arg)
    npt.assert_almost_equal(distfn.ppf(cdf05-1e-6,*arg),ppf05,
                            err_msg=str(distfn) + 'ppf-cdf-median')
    assert (distfn.ppf(cdf05+1e-4,*arg)>ppf05), str(distfn) + 'ppf-cdf-next'

def check_cdf_ppf_private(distfn,arg):
    ppf05 = distfn._ppf(0.5,*arg)
    cdf05 = distfn.cdf(ppf05,*arg)
    npt.assert_almost_equal(distfn._ppf(cdf05-1e-6,*arg),ppf05,
                            err_msg=str(distfn) + 'ppf-cdf-median')
    assert (distfn._ppf(cdf05+1e-4,*arg)>ppf05), str(distfn) + 'ppf-cdf-next'

def check_ppf_ppf(distfn, arg):
    assert distfn.ppf(0.5,*arg) < np.inf
    ppfs = distfn.ppf([0.5,0.9],*arg)
    ppf_s = [distfn._ppf(0.5,*arg), distfn._ppf(0.9,*arg)]
    assert np.all(ppfs < np.inf)
    assert ppf_s[0] == distfn.ppf(0.5,*arg)
    assert ppf_s[1] == distfn.ppf(0.9,*arg)
    assert ppf_s[0] == ppfs[0]
    assert ppf_s[1] == ppfs[1]

def check_pmf_cdf(distfn, arg):
    startind = np.int(distfn._ppf(0.01,*arg)-1)
    index = range(startind,startind+10)
    cdfs = distfn.cdf(index,*arg)
    npt.assert_almost_equal(cdfs, distfn.pmf(index, *arg).cumsum() + \
                            cdfs[0] - distfn.pmf(index[0],*arg),
                            decimal=4, err_msg='pmf-cdf')

def check_generic_moment(distfn, arg, m, k, decim):
    npt.assert_almost_equal(distfn.generic_moment(k,*arg), m, decimal=decim,
                            err_msg= str(distfn) + ' generic moment test')

def check_moment_frozen(distfn, arg, m, k, decim):
    npt.assert_almost_equal(distfn(*arg).moment(k), m, decimal=decim,
                            err_msg= str(distfn) + ' frozen moment test')

def check_oth(distfn, arg):
    #checking other methods of distfn
    meanint = round(distfn.stats(*arg)[0]) # closest integer to mean
    npt.assert_almost_equal(distfn.sf(meanint, *arg), 1 - \
                            distfn.cdf(meanint, *arg), decimal=8)
    median_sf = distfn.isf(0.5, *arg)

    assert distfn.sf(median_sf - 1, *arg) > 0.5
    assert distfn.cdf(median_sf + 1, *arg) > 0.5
    npt.assert_equal(distfn.isf(0.5, *arg), distfn.ppf(0.5, *arg))


if __name__ == "__main__":
    nose.run(argv=['', __file__])
