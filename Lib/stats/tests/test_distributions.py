""" Test functions for stats module

"""


import unittest,sys
from scipy_test.testing import *
from scipy_test.testing import ScipyTestCase as TestCase
set_package_path()
import scipy
from scipy import stats
del sys.path[0]

import types

def kolmogorov_test(diststr,args=(),N=20,significance=0.01):
    qtest = stats.ksoneisf(significance,N)
    cdf = eval('stats.'+diststr+'.cdf')
    dist = eval('stats.'+diststr)
    # Get random numbers
    kwds = {'size':N}
    vals = scipy.sort(dist.rvs(*args,**kwds))
    cdfvals = cdf(vals,*args)
    q = max(abs(cdfvals-arange(1.0,N+1)/N))
    assert (q < qtest), "Failed q=%f, bound=%f, alpha=%f" % (q, qtest, significance)
    return


# generate test cases to test cdf and distribution consistency

dists = ['uniform','norm','lognorm','expon','beta',
         'powerlaw','bradford','burr','fisk','cauchy','halfcauchy',
         'foldcauchy','gamma','gengamma','loggamma',
         'alpha','anglit','arcsine','betaprime','erlang',
         'dgamma','exponweib','exponpow','frechet_l','frechet_r',
         'gilbrat','f','ncf','chi2','chi','nakagami','genpareto',
         'genextreme','genhalflogistic','pareto','lomax','halfnorm',
         'halflogistic','fatiguelife','foldnorm','ncx2','t','nct',
         'weibull_min','weibull_max','dweibull','maxwell','rayleigh',
         'genlogistic', 'logistic','gumbel_l','gumbel_r','gompertz',
         'hypsecant', 'laplace', 'reciprocal','triang','tukeylambda']

for dist in dists:
    distfunc = eval('stats.'+dist)
    nargs = distfunc.numargs
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


class test_randint(TestCase):
    def check_rvs(self):
        vals = stats.randint.rvs(5,30,size=100)
        assert(scipy.all(vals < 30) & scipy.all(vals >= 5))
        assert(len(vals) == 100)
        vals = stats.randint.rvs(5,30,size=(2,50))
        assert(scipy.shape(vals) == (2,50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.randint.rvs(15,46)
        assert((val >= 15) & (val < 46))
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])

    def check_pdf(self):
        k = scipy.r_[0:36]
        out = scipy.where((k >= 5) & (k < 30), 1.0/(30-5), 0)
        vals = stats.randint.pmf(k,5,30)
        assert(scipy.all(vals == out))

    def check_cdf(self):
        x = scipy.r_[0:36:100j]
        k = scipy.floor(x)
        out = scipy.select([k>=30,k>=5],[1.0,(k-5.0+1)/(30-5.0)],0)
        vals = stats.randint.cdf(x,5,30)
        assert_array_almost_equal(vals, out, decimal=12)

class test_binom(TestCase):
    def check_rvs(self):
        vals = stats.binom.rvs(10, 0.75, size=(2, 50))
        assert(scipy.all(vals >= 0) & scipy.all(vals <= 10))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.binom.rvs(10, 0.75)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
        
class test_bernoulli(TestCase):
    def check_rvs(self):
        vals = stats.bernoulli.rvs(0.75, size=(2, 50))
        assert(scipy.all(vals >= 0) & scipy.all(vals <= 1))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.bernoulli.rvs(0.75)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
class test_nbinom(TestCase):
    def check_rvs(self):
        vals = stats.nbinom.rvs(10, 0.75, size=(2, 50))
        assert(scipy.all(vals >= 0))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.nbinom.rvs(10, 0.75)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
class test_geom(TestCase):
    def check_rvs(self):
        vals = stats.geom.rvs(0.75, size=(2, 50))
        assert(scipy.all(vals >= 0))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.geom.rvs(0.75)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
class test_hypergeom(TestCase):
    def check_rvs(self):
        vals = stats.hypergeom.rvs(20, 10, 3, size=(2, 50))
        assert(scipy.all(vals >= 0) &
               scipy.all(vals <= 3))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.hypergeom.rvs(20, 3, 10)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
class test_logser(TestCase):
    def check_rvs(self):
        vals = stats.logser.rvs(0.75, size=(2, 50))
        assert(scipy.all(vals >= 1))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.logser.rvs(0.75)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
class test_poisson(TestCase):
    def check_rvs(self):
        vals = stats.poisson.rvs(0.5, size=(2, 50))
        assert(scipy.all(vals >= 0))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.poisson.rvs(0.5)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
class test_zipf(TestCase):
    def check_rvs(self):
        vals = stats.zipf.rvs(1.5, size=(2, 50))
        assert(scipy.all(vals >= 1))
        assert(scipy.shape(vals) == (2, 50))
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.zipf.rvs(1.5)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
class test_dlaplace(TestCase):
    def check_rvs(self):
        vals = stats.dlaplace.rvs(1.5 , size=(2, 50))
        assert(scipy.shape(vals) == (2, 50))        
        assert(vals.typecode() in scipy.typecodes['AllInteger'])
        val = stats.dlaplace.rvs(1.5)
        assert(isinstance(val, scipy.ArrayType))
        assert(val.typecode() in scipy.typecodes['AllInteger'])
        
if __name__ == "__main__":
    ScipyTest('stats.distributions').run()

