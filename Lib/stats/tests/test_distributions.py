""" Test functions for stats module

"""


from numpy.testing import *

set_package_path()
import numpy.oldnumeric as numpy
from numpy import typecodes
import stats
restore_path()

import types

def kolmogorov_test(diststr,args=(),N=20,significance=0.01):
    qtest = stats.ksoneisf(significance,N)
    cdf = eval('stats.'+diststr+'.cdf')
    dist = eval('stats.'+diststr)
    # Get random numbers
    kwds = {'size':N}
    vals = numpy.sort(dist.rvs(*args,**kwds))
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
class test_%s(ScipyTestCase):
    def check_cdf(self):
        D,pval = stats.kstest('%s','',args=%s,N=30)
        if (pval < %f):
            D,pval = stats.kstest('%s','',args=%s,N=30)
            #if (pval < %f):
            #    D,pval = stats.kstest('%s','',args=%s,N=30)
        assert (pval > %f), "D = " + str(D) + "; pval = " + str(pval) + "; alpha = " + str(alpha) + "\nargs = " + str(%s)
""" % (dist,dist,args,alpha,dist,args,alpha,dist,args,alpha,args)
    exec exstr


class test_randint(ScipyTestCase):
    def check_rvs(self):
        vals = stats.randint.rvs(5,30,size=100)
        assert(numpy.all(vals < 30) & numpy.all(vals >= 5))
        assert(len(vals) == 100)
        vals = stats.randint.rvs(5,30,size=(2,50))
        assert(numpy.shape(vals) == (2,50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.randint.rvs(15,46)
        assert((val >= 15) & (val < 46))
        assert isinstance(val, numpy.ScalarType),`type(val)`
        assert(val.dtype.char in typecodes['AllInteger'])

    def check_pdf(self):
        k = numpy.r_[0:36]
        out = numpy.where((k >= 5) & (k < 30), 1.0/(30-5), 0)
        vals = stats.randint.pmf(k,5,30)
        assert_array_almost_equal(vals,out)

    def check_cdf(self):
        x = numpy.r_[0:36:100j]
        k = numpy.floor(x)
        out = numpy.select([k>=30,k>=5],[1.0,(k-5.0+1)/(30-5.0)],0)
        vals = stats.randint.cdf(x,5,30)
        assert_array_almost_equal(vals, out, decimal=12)

class test_binom(ScipyTestCase):
    def check_rvs(self):
        vals = stats.binom.rvs(10, 0.75, size=(2, 50))
        assert(numpy.all(vals >= 0) & numpy.all(vals <= 10))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.binom.rvs(10, 0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])


class test_bernoulli(ScipyTestCase):
    def check_rvs(self):
        vals = stats.bernoulli.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 0) & numpy.all(vals <= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.bernoulli.rvs(0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class test_nbinom(ScipyTestCase):
    def check_rvs(self):
        vals = stats.nbinom.rvs(10, 0.75, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.nbinom.rvs(10, 0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class test_geom(ScipyTestCase):
    def check_rvs(self):
        vals = stats.geom.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.geom.rvs(0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class test_hypergeom(ScipyTestCase):
    def check_rvs(self):
        vals = stats.hypergeom.rvs(20, 10, 3, size=(2, 50))
        assert(numpy.all(vals >= 0) &
               numpy.all(vals <= 3))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.hypergeom.rvs(20, 3, 10)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class test_logser(ScipyTestCase):
    def check_rvs(self):
        vals = stats.logser.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.logser.rvs(0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class test_poisson(ScipyTestCase):
    def check_rvs(self):
        vals = stats.poisson.rvs(0.5, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.poisson.rvs(0.5)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class test_zipf(ScipyTestCase):
    def check_rvs(self):
        vals = stats.zipf.rvs(1.5, size=(2, 50))
        assert(numpy.all(vals >= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.zipf.rvs(1.5)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class test_dlaplace(ScipyTestCase):
    def check_rvs(self):
        vals = stats.dlaplace.rvs(1.5 , size=(2, 50))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.dlaplace.rvs(1.5)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

if __name__ == "__main__":
    ScipyTest('stats.distributions').run()
