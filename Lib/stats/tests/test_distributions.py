""" Test functions for stats module

"""


import unittest
from scipy_test.testing import assert_array_equal, assert_equal
from scipy_test.testing import assert_almost_equal, rand
from scipy_test.testing import assert_array_almost_equal
from scipy_test.testing import ScipyTestCase as TestCase
import scipy
from scipy import stats

def kolmogorov_test(diststr,args=(),N=20,significance=0.01):
    qtest = stats.ksoneisf(significance,N)
    cdf = eval('stats.'+diststr+'cdf')
    dist = eval('stats.'+diststr)
    # Get random numbers
    kwds = {'size':N}
    vals = scipy.sort(dist(*args,**kwds))
    cdfvals = cdf(vals,*args)
    q = max(abs(cdfvals-arange(1.0,N+1)/N))
    assert (q < qtest), "Failed q=%f, bound=%f, alpha=%f" % (q, qtest, significance)
    return


# generate test cases to test cdf and distribution consistency

dists = ['uniform','stnorm','norm','lognorm','expon','beta',
         'power','bradford','burr','fisk','cauchy','halfcauchy',
         'foldcauchy','gamma','gengamma','loggamma',
         'alpha','anglit','arcsine','betaprime','erlang',
         'dgamma','extreme3','exponweib','exponpow','frechet',
         'gilbrat','f','ncf','chi2','chi','nakagami','genpareto',
         'genextreme','genhalflogistic','pareto','lomax','halfnorm',
         'halflogistic','fatiguelife','foldnorm','ncx2','t','nct',
         'weibull','dweibull','maxwell','rayleigh','genlogistic',
         'logistic','gumbel','gompertz','hypsecant','laplace',
         'reciprocal','triang','tukeylambda']

for dist in dists:
    distfunc = eval('stats.'+dist)
    nargs = distfunc.func_code.co_argcount - 1  # size argument is last
    alpha = 0.01
    if dist == 'fatiguelife':
        alpha = 0.001        
    if dist == 'erlang':
        args = str((4,)+tuple(rand(2)))
    elif dist == 'frechet':
        args = str(tuple(2*rand(1))+(0,)+tuple(2*rand(2)))
    elif dist == 'triang':
        args = str(tuple(rand(nargs)))
    elif dist == 'reciprocal':
        vals = rand(nargs)
        vals[1] = vals[0] + 1.0
        args = str(tuple(vals))
    else:
        args = str(tuple(1.0+rand(nargs)))
    exstr = r"""
class test_%s(TestCase):
    def check_cdf(self):
        print "Testing %s"
        D,pval = stats.kstest('%s','',args=%s,N=30)
        if (pval < %f):
            D,pval = stats.kstest('%s','',args=%s,N=30)
            #if (pval < %f):
            #    D,pval = stats.kstest('%s','',args=%s,N=30)
        assert (pval > %f), "D = " + str(D) + "; pval = " + str(pval) + "; alpha = " + str(alpha) + "\nargs = " + str(%s)
""" % (dist,dist,dist,args,alpha,dist,args,alpha,dist,args,alpha,args)
    exec exstr

def test_suite(level=1):
    suites = []
    if level > 0:
        for dist in dists:
            thistest = eval("test_%s"%dist)
            suites.append( unittest.makeSuite(thistest,'check_') )
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level=level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
