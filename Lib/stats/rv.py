from __future__ import nested_scopes

import rand
import Numeric
import sys
import math
from types import *
Num = Numeric
import scipy.special
special = scipy.special
from scipy_base.fastumath import *

ArgumentError = "ArgumentError"

_scalerr = "Scale must be positive."
_parmerr = "Parameters must be positive."

def _chkpos(parm):
    if (parm <=0):
        raise ValueError, _parmerr
    return

def seed(x=0,y=0):
    """seed(x, y), set the seed using the integers x, y; 
    Set a random one from clock if  y == 0
    """
    if type (x) != IntType or type (y) != IntType :
        raise ArgumentError, "seed requires integer arguments."
    if y == 0:
        import time
        t = time.time()
        ndigits = int(math.log10(t))
        base = 10**(ndigits/2)
        x = int(t/base)
        y = 1 + int(t%base)
    rand.set_seeds(x,y)

seed()

def get_seed():
    "Return the current seed pair"
    return rand.get_seeds()

def _build_random_array(fun, args, size=None):
# Build an array by applying function fun to
# the arguments in args, creating an array with
# the specified shape.
# Allows an integer shape n as a shorthand for (n,).
    if isinstance(size, IntType): 
        size = [size]
    if size is not None and len(size) != 0:
        n = Num.multiply.reduce(size)
        s = apply(fun, args + (n,))
        s.shape = size
        return s
    else:
        n = 1
        s = apply(fun, args + (n,))
        return s[0]    

def random(size=None):
    "Returns array of random numbers between 0 and 1"
    return _build_random_array(rand.sample, (), size)

def uniform(loc=0.0, scale=1.0, size=None):
    """Returns array of given shape of uniformly distributed reals
    in the range loc to scale+loc (exclusive of endpoints).
    """
    return _build_random_array(rand.uniform, (loc, scale+loc), size)

def randint(min, max=None, size=None):
    """random integers >=min, < max.  If max not given, random integers >= 0, <min"""
    if max is None:
        max = min
        min = 0
    a = Num.floor(uniform(min, max, size))
    if isinstance(a, Num.ArrayType):
        return a.astype(Num.Int)
    else:
        return int(a)
     
def random_integers(max, min=1, size=None):
    """random_integers(max, min=1, size=None) = random integers in range min-max inclusive"""
    return randint(min, max+1, size) 
     
def permutation(arg):
    """If arg is an integer, a permutation of indices arange(n), otherwise
    a permuation of the sequence"""
    if isinstance(arg,IntType):
        arg = Num.arange(arg)
    return rand.permutation(arg)

def stnorm(size=None):
    """returns array of random numbers normally distributed with mean 0
    and standard deviation 1"""
    return _build_random_array(rand.standard_normal, (), size)

def norm(loc=0.0, scale=1.0, size=None):
        """returns array of random numbers randomly distributed with
        specified mean (location) and standard deviation (scale)"""
        return _build_random_array(rand.normal, (loc, scale), size)

def lognorm(std=0.0, loc=0.0, scale=1.0, size=None):
    """returns array of random numbers lognormally distributed where
    the underlying normal has the given mean and standard-deviation.
    """
    if (std <= 0.0) or (scale <= 0.0):
        raise ValueError, 'std and scale must be positive'
    else:
        return exp(std * norm(size=size))*scale + loc
    return

def multivariate_normal(mean, cov, size=None):
       """returns an array containing multivariate normally distributed random
          numbers with specified mean and covariance.

          mean must be a 1 dimensional array. cov must be a square two dimensional
          array with the same number of rows and columns as mean has elements.

          The first form returns a single 1-D array containing a multivariate
          normal.

          The second form returns an array of shape (m, n, ..., cov.shape[0]).
          In this case, output[i,j,...,:] is a 1-D array containing a multivariate
          normal."""
       if isinstance(size, IntType): 
           size = [size]
       if size is None:
           n = 1
       else:
           n = Num.product(size)
       output = rand.multivariate_normal(mean, cov, n)
       if size is not None:
           final_shape = list(size[:])
           final_shape.append(mean.shape[0])
           output.shape = final_shape
       return output       

def stexpon(size=None):
    """returns array of random numbers with a standard exponential distribution (mean=1)"""
    return _build_random_array(rand.standard_exp, (), size)

def expon(loc=0.0, scale=1.0, size=None):
    """returns array of random numbers exponentially distributed scale is equivalent to 1.0 / lambda"""
    return _build_random_array(rand.exponential, (scale,), size) + loc

def beta(a, b, loc=0.0, scale=1.0, size=None):
    """returns array of beta distributed random numbers."""
    return _build_random_array(rand.beta, (a, b), size)*scale + loc

def power(a, loc=0.0, scale=1.0, size=None):
    """returns array of power-function distributed random numbers

    Special case of beta with second shape parameter = 1.0
    """
    if (a <= 0) or (scale <=0):
        raise ValueError, _parmerr
    return pow(random(size=size),1.0/a)*scale + loc

def bradford(c, loc=0.0, scale=1.0, size=None):
    """returns array of bradford distributed random numbers."""
    if (c<=0) or (scale<=0):
        raise ValueError, "The shape (first) and scale (third) parameters must be larger than 0." 
    U = random(size=size) 
    return 1.0/c*((1.0+c)**U - 1)*scale + loc


# Generalization of Fisk (d=1)
#
def burr(c, d, loc=0.0, scale=1.0, size=None):
    """returns array of Burr distributed random numbers."""
    if (c<=0) or (d<=0) or (scale<=0):
        raise ValueError, "The shape parameters (first and second arguments) and the scale must all be positive."
    U = random(size=size)
    return (U**(-1.0/d)-1)**(-1.0/c)*scale + loc

def fisk(c, loc=0.0, scale=1.0, size=None):
    """Fisk random numbers (special case of Burr for d=1.0)
    """
    return burr(c, 1.0, loc, scale, size)


## Sometimes called the Lorentz
def cauchy(loc=0.0, scale=1.0, size=None):
    if (scale <=0):
        raise ValueError, _scalerr
    U = random(size=size)
    return tan(pi*(U-0.5))*scale+loc

def halfcauchy(loc=0.0, scale=1.0, size=None):
    if (scale <=0):
        raise ValueError, _scalerr
    return abs(cauchy(size=size))*scale + loc

def foldcauchy(c, loc=0.0, scale=1.0, size=None):
    if (scale <=0) or (c<0):
        raise ValueError, _parmerr
    return abs(cauchy(loc=c,size=size))*scale + loc

def stgamma(a, size=None):
    """returns array of random numbers with a standard gamma distribution with mean a"""
    return _build_random_array(rand.standard_gamma, (a,), size)

def gamma(a, loc=0.0, scale=1.0, size=None):
    """returns array of gamma distributed random numbers with
    shape (alpha) and scale (beta) and location parameter as given."""
    # Note the different parameter definitions as the ones used in
    #  ranlib.
    if (a <=0) or (scale <=0):
        raise ValueError, _parmerr
    return _build_random_array(rand.gamma, (1.0/scale, a), size) + loc

def gengamma(a, c, loc=0.0, scale=1.0, size=None):
    """Generalized Gamma random numbers.
    """
    if (a<=0) or (c==0) or (scale<=0):
        raise ValueError, _parmerr
    u = random(size=size)
    vals = special.gammaincinv(a, special.gamma(a)*q)**(1.0/c)
    return vals*scale + loc

def loggamma(c, loc=0.0, scale=1.0, size=None):
    """Log Gamma distribution.
    """
    if (c <=0) or (scale<=0):
        return ValueErro, _parmerr
    u = random(size=size)
    vals = log(special.gammaincinv(c,u*special.gamma(c)))
    return vals*scale + loc

######
#####

def alpha(a, loc=0.0, scale=1.0, size=None):
    """Alpha distributed random variables.
    """
    u = random(size=size)
    return (1.0/(a-special.ndtri(u*special.ndtr(a))))*scale + loc

def anglit(loc=0.0, scale=1.0, size=None):
    """Anglit random variables
    """
    u = random(size=size)
    return (arcsin(sqrt(u))-pi/4)*scale + loc

def arcsine(loc=0.0, scale=1.0, size=None):
    u = random(size=size)
    return sin(pi*u/2.0)**2 * scale + loc


def betaprime(a, b, loc=0.0, scale=1.0, size=None):
    u1 = gamma(a,size=size)
    u2 = gamma(b,size=size)
    return (u1 / u2)*scale + loc

def erlang(n, loc=0.0, scale=1.0, size=None):
    """Erlang distributed random variables."""
    if (n != floor(n)):
        raise TypeError, "Erlang distribution only defined for integer n."
    return gamma(n, loc=loc, scale=scale, size=size)

def dgamma(a, loc=0.0, scale=1.0, size=None):
    u = random(size=size)
    return (gamma(a, 0.0, 1.0, size=size)*(Num.where(u>=0.5,1,-1)))*scale+loc

def extreme3(c, loc=0.0, scale=1.0, size=None):
    """Extreme Value Type III
    """
    if (c <=0) or (scale <=0):
        raise ValueError, _parmerr
    u = random(size=size)
    return -pow(-log(u),1.0/c)*scale + loc

def exponweib(a, c, loc=0.0, scale=1.0, size=None):
    """Exponentiated Weibull
    """
    if (c <=0) or (scale <=0) or (a <=0):
        raise ValueError, _parmerr
    u = random(size=size)
    return -pow(-log(1-pow(u,1.0/a)),1.0/c)*scale + loc

def exponpow(b, loc=0.0, scale=1.0, size=None):
    """Exponentiated Weibull
    """
    if (b <=0) or (scale <=0):
        raise ValueError, _parmerr
    u = random(size=size)
    return pow(log(1.0-log(1-u)),1.0/b)*scale + loc

def frechet(c, loc=0.0, scale=1.0, size=None):
    if (c <=0) or (scale <=0):
        raise ValueError, _parmerr
    u = random(size=size)
    return pow(-log(u),-1.0/c)*scale + loc


def gilbrat(size=None):
    """returns array of gilbert distributed random numbers (special case of
    lognormal distribution with std=1.0, loc=0.0, and scale=1.0)
    """
    return lognorm(std=1.0, loc=0.0, scale=1.0, size=size)

def f(dfn, dfd, size=None):
    """returns array of F distributed random numbers with dfn degrees of freedom in the numerator and dfd degrees of freedom in the denominator."""
    return _build_random_array(rand.f, (dfn, dfd), size)

def ncf(dfn, dfd, nc, size=None):
    """returns array of noncentral F distributed random numbers with dfn degrees of freedom in the numerator and dfd degrees of freedom in the denominator, and noncentrality parameter nonc."""
    return _build_random_array(rand.noncentral_f, (dfn, dfd, nc), size)

def chi2(df, loc=0.0, scale=1.0, size=None):
    """returns array of chi squared distributed random numbers with df degrees of freedom."""
    return _build_random_array(rand.chi2, (df,), size)*scale + loc

def chi(df, loc=0.0, scale=1.0, size=None):
    """returns array of chi distributed random numbers with df degrees of freedom."""
    return sqrt(chi2(df,0,1.0,size=size))*scale + loc

def nakagami(df, loc=0.0, scale=1.0, size=None):
    """returns array of Nakagami distributed random numbers with df degrees of freedom."""
    if (df <=0):
        raise ValueError, _parmerr
    U = random(size=size)
    return loc + scale*sqrt(1.0/df*special.gammainccinv(df,1-U))

def genpareto(c, loc=0.0, scale=1.0, size=None):
    """Generalized Pareto
    """
    if (scale <=0) or (c==0):
        raise ValueError, "Scale must be positive and parameter non-zero."
    u = random(size=size)
    return 1.0/c*(pow(1.0/(1-q),c)-1)

def genextreme(c, loc=0.0, scale=1.0, size=None):
    if (scale <=0):
        raise ValueError, _scalerr
    u = random(size=size)
    if c == 0:
        vals= -log(-log(u))
    else:
        vals = 1.0/c*(1-(-log(u))^c)
    return vals*scale + loc

def genhalflogistic(c, loc=0.0, scale=1.0, size=None):
    if (c <=0) or (scale <=0):
        raise ValueError, _parmerr
    u = random(size=size)
    return 1.0/c * (1-((1-u)/(1+u))**c)

def pareto(c, loc=0.0, scale=1.0, size=None):
    """Pareto First Kind.
    """
    if (c<=0) or (scale <= 0):
        raise ValueError, _parmerr
    u = random(size=size)
    return (1-u)**(-1.0/c)

# special case with no location parameter.  Use pareto if you need
#   location parameter.
def lomax(c, scale=1.0, size=None):
    """Pareto of the second kind
    """
    return pareto(c, -1, scale, size)

def halfnorm(loc=0.0, scale=1.0, size=None):
    """HalfNormal Random variates."""
    # return chi(1, loc, scale, size)
    return abs(stnorm(size=size))*scale + loc

def halflogistic(loc=0.0, scale=1.0, size=None):
    """HalfLogistic Random variates."""
    u = random(size=size)
    return log((1+u)/(1-u))*scale + loc

def fatiguelife(c, loc=0.0, scale=1.0, size=None):
    z = stnorm(size=size)
    U = random(size=size)
    fac = 2 + c*c*z*z
    det = sqrt(fac*fac - 4)
    t1 = fac + det
    t2 = fac - det
    vals = t1*(U>0.5) + t2*(U<0.5)
    return vals*scale + loc

def foldnorm(c=0.0, loc=0.0, scale=1.0, size=None):
    """Folded Normal Random Variates."""
    return abs(norm(c*scale, scale)) + loc

def ncx2(df, nc, size=None):
    """returns array of noncentral chi squared distributed random numbers
    with df degrees of freedom and noncentrality parameter."""
    return _build_random_array(rand.noncentral_chi2, (df, nc), size)

def t(df, loc=0.0, scale=1.0, size=None):
    """returns array of student_t distributed random numbers
    with df degrees of freedom."""
    Y = f(df, df, size=size)
    sY = sqrt(Y)
    return 0.5*sqrt(df)*(sY - 1.0/sY)*scale + loc

def nct(df, nc, size=None):
    """returns array of noncentral student_t distributed random numbers
    with df degrees of freedom."""    
    return norm(nc)*Num.sqrt(df) / Num.sqrt(noncentral_chi2(df,nc,size))

def weibull(c, loc=0.0, scale=1.0, size=None):
    """Weibull distributed random variates."""
    u = random(size=size)
    return pow(log(1.0/(1.0-u)), 1.0/c)*scale + loc

def dweibull(c, loc=0.0, scale=1.0, size=None):
    """Weibull distributed random variates."""
    u = random(size=size)
    return (weibull(c, 0.0, 1.0, size=size)*(Num.where(u>=0.5,1,-1)))*scale+loc

def maxwell(scale=1.0, size=None):
    """return array of Maxwell distributed random numbers"""
    return chi(3, loc=0.0, scale=scale, size=size)
    
def rayleigh(scale=1.0, size=None):
    """returns array of Rayleigh distributed random numbers with mode given
    by the scale parameter."""
    return chi(2,loc=0.0,scale=scale,size=size)

def genlogistic(shape, loc=0.0, scale=1.0, size=None):
    """Generalized Logistic Random Variates."""
    if (shape<=0) or (scale <= 0.0):
        raise ValueError, "shape and scale must be positive."
    u = random(size=size)
    return -log(pow(u,-1.0/shape)-1.0)*scale + loc


def logistic(loc=0.0, scale=1.0, size=None):
    """Logistic Random Numbers."""
    return genlogistic(1.0, loc=loc, scale=scale, size=size)

def gumbel(loc=0.0, scale=1.0, size=None):
    """Gumbel (Log-Weibull, Fisher-Tippett) RN"""
    u = random(size=size)
    return loc - scale*log(-log(u))

def gompertz(c, loc=0.0, scale=1.0, size=None):
    if (c<=0) or (scale <=0):
        raise ValueError, _parmerr
    u = random(size=size)
    return log(1.0-1.0/c*log(1.0-u))*scale + loc


def hypsecant(loc=0.0, scale=1.0, size=None):
    """Hyperbolic secant random numbers."""
    u = random(size=size)
    return log(tan(pi*u/2))*scale + loc

def laplace(loc=0.0, scale=1.0, size=None):
    """Laplace (Double Exponential, Bilateral Exponential) random numbers."""
    u = random(size=size)
    return (Num.where(u<=0.5,log(2*u),-log(2*(1-u))))*scale + loc

def reciprocal(a, b, loc=0.0, scale=1.0, size=None):
    if (a<=0) or (b<=0) or (scale<=0):
        raise ValueError, _parmerr
    u = random(size=size)
    return a*pow(b/a,u)*scale + loc

def triang(c, loc=0.0, scale=1.0, size=None):
    if (c<0) or (c>1):
        raise ValueError, "C must be in [0,1]."
    u = random(size=size)
    return Num.where(u<c, sqrt(c*u), 1-sqrt((1-c)*(1-u)))*scale + loc

def tukeylambda(lam, loc=0.0, scale=1.0, size=None):
    """Tukey-Lambda random numbers.
    """
    _chkpos(scale)
    u = random(size=size)
    if (lam == 0):
        return scale * log(u/(1.0-u)) + loc
    else:
        return scale*1.0/lam*(u**lam - (1-u)**lam) + loc

    
#####################################
# General purpose continuous
######################################

def randwppf(ppf, args=(), size=None):
    """returns an array of randomly distributed integers of a distribution
    whose percent point function (inverse of the CDF) is given.

    args is a tuple of extra arguments to the ppf function (i.e. shape,
    location, scale), and size is the size of the output.  Note the ppf
    function must accept an array of q values to compute over.
    """
    U = random(size=size)
    return apply(ppf, (U,)+args)

def randwcdf(cdf, mean=1.0, args=(), size=None):
    """returns an array of randomly distributed integers of a distribution
    whose cumulative distribution function (CDF) is given.

    mean is the mean of the distribution (helps the solver).
    args is a tuple of extra arguments to the cdf function (i.e. shape,
    location, scale), and size is the size of the output.  Note the
    cdf function needs to accept a single value to compute over.
    """
    import scipy.optimize as optimize
    def _ppfopt(x, q, *nargs):
        newargs = (x,)+nargs
        return cdf(*newargs) - q

    def _ppf(q, *nargs):
        return optimize.fsolve(_ppfopt, mean, args=(q,)+nargs)

    _vppf = scipy.special.general_function(_ppf)
    U = random(size=size)
    return apply(_vppf,(U,)+args)


#################################################
## DISCRETE
##################################################

def bernoulli(pr=0.5, size=None):
    return binom(1, pr, size)

def binom(trials, pr=0.5, size=None):
    """returns array of binomially distributed random integers.

           trials is the number of trials in the binomial distribution.
           p is the probability of an event in each trial of the binomial distribution."""
    return _build_random_array(rand.binomial, (trials, pr), size)

def nbinom(trials, pr=0.5, size=None):
    """returns array of negative binomially distributed random integers.
    
           trials is the number of trials in the negative binomial
                  distribution.
           p is the probability of an event in each trial of the
                  negative binomial distribution.
    """
    return _build_random_array(rand.negative_binomial, (trials, pr), size)


def multinom(trials, probs, size=None):
    """returns array of multinomial distributed integer vectors.

           trials is the number of trials in each multinomial distribution.
           probs is a one dimensional array. There are len(prob)+1 events. 
           prob[i] is the probability of the i-th event, 0<=i<len(prob).
           The probability of event len(prob) is 1.-Numeric.sum(prob).

       The first form returns a single 1-D array containing one multinomially
           distributed vector.

           The second form returns an array of size (m, n, ..., len(probs)).
           In this case, output[i,j,...,:] is a 1-D array containing a multinomially
           distributed integer 1-D array."""
        # Check preconditions on arguments
    probs = Num.array(probs)
    if len(probs.shape) != 1:
        raise ArgumentError, "probs must be 1 dimensional."
        # Compute shape of output
    if type(size) == type(0): size = [size]
    final_shape = size[:]
    final_shape.append(probs.shape[0]+1)
    x = rand.multinomial(trials, probs.astype(Num.Float32), Num.multiply.reduce(size))
        # Change its shape to the desire one
    x.shape = final_shape
    return x

def poisson(mu, size=None):
    """returns array of poisson distributed random integers with specifed mean."""
    if (mu < 0):
        raise ValueError, "mu must be > 0."
    return _build_random_array(rand.poisson, (mu,), size)

def mean_var_test(x, type, mean, var, skew=[]):
    n = len(x) * 1.0
    x_mean = Num.sum(x)/n
    x_minus_mean = x - x_mean
    x_var = Num.sum(x_minus_mean*x_minus_mean)/(n-1.0)
    print "\nAverage of ", len(x), type
    print "(should be about ", mean, "):", x_mean
    print "Variance of those random numbers (should be about ", var, "):", x_var
    if skew != []:
       x_skew = (Num.sum(x_minus_mean*x_minus_mean*x_minus_mean)/9998.)/x_var**(3./2.)
       print "Skewness of those random numbers (should be about ", skew, "):", x_skew

def test():
    x, y = get_seed()
    print "Initial seed", x, y
    seed(x, y)
    x1, y1 = get_seed()
    if x1 != x or y1 != y:
        raise SystemExit, "Failed seed test."
    print "First random number is", random()
    print "Average of 10000 random numbers is", Num.sum(random(10000))/10000.
    x = random([10,1000])
    if len(x.shape) != 2 or x.shape[0] != 10 or x.shape[1] != 1000:
        raise SystemExit, "random returned wrong shape"
    x.shape = (10000,)
    print "Average of 100 by 100 random numbers is", Num.sum(x)/10000.
    y = uniform(0.5,0.6, (1000,10))
    if len(y.shape) !=2 or y.shape[0] != 1000 or y.shape[1] != 10:
        raise SystemExit, "uniform returned wrong shape"
    y.shape = (10000,)
    if Num.minimum.reduce(y) <= 0.5 or Num.maximum.reduce(y) >= 0.6:
        raise SystemExit, "uniform returned out of desired range"
    print "randint(1, 10, size=[50])"
    print randint(1, 10, size=[50])
    print "permutation(10)", permutation(10)
    print "randint(3,9)", randint(3,9)
    print "random_integers(10, size=[20])"
    print random_integers(10, size=[20])
    s = 3.0
    x = norm(2.0, s, [10, 1000])
    if len(x.shape) != 2 or x.shape[0] != 10 or x.shape[1] != 1000:
        raise SystemExit, "standard_normal returned wrong shape"
    x.shape = (10000,)
    mean_var_test(x, "normally distributed numbers with mean 2 and variance %f"%(s**2,), 2, s**2, 0)
    x = exponential(3, 10000)
    mean_var_test(x, "random numbers exponentially distributed with mean %f"%(s,), s, s**2, 2)
    x = multivariate_normal(Num.array([10,20]), Num.array(([1,2],[2,4])))
    print "\nA multivariate normal", x
    if x.shape != (2,): raise SystemExit, "multivariate_normal returned wrong shape"
    x = multivariate_normal(Num.array([10,20]), Num.array([[1,2],[2,4]]), [4,3])
    print "A 4x3x2 array containing multivariate normals"
    print x
    if x.shape != (4,3,2): raise SystemExit, "multivariate_normal returned wrong shape"
    x = multivariate_normal(Num.array([-100,0,100]), Num.array([[3,2,1],[2,2,1],[1,1,1]]), 10000)
    x_mean = Num.sum(x)/10000.
    print "Average of 10000 multivariate normals with mean [-100,0,100]"
    print x_mean
    x_minus_mean = x - x_mean
    print "Estimated covariance of 10000 multivariate normals with covariance [[3,2,1],[2,2,1],[1,1,1]]"
    print Num.matrixmultiply(Num.transpose(x_minus_mean),x_minus_mean)/9999.
    x = beta(5.0, 10.0, 10000)
    mean_var_test(x, "beta(5.,10.) random numbers", 0.333, 0.014)
    x = gamma(.01, 2., 10000)
    mean_var_test(x, "gamma(.01,2.) random numbers", 2*100, 2*100*100)
    x = chi_square(11., 10000)
    mean_var_test(x, "chi squared random numbers with 11 degrees of freedom", 11, 22, 2*Num.sqrt(2./11.))
    x = F(5., 10., 10000)
    mean_var_test(x, "F random numbers with 5 and 10 degrees of freedom", 1.25, 1.35)
    x = poisson(50., 10000)
    mean_var_test(x, "poisson random numbers with mean 50", 50, 50, 0.14)
    print "\nEach element is the result of 16 binomial trials with probability 0.5:"
    print binomial(16, 0.5, 16)
    print "\nEach element is the result of 16 negative binomial trials with probability 0.5:"
    print negative_binomial(16, 0.5, [16,])
    print "\nEach row is the result of 16 multinomial trials with probabilities [0.1, 0.5, 0.1 0.3]:"
    x = multinomial(16, [0.1, 0.5, 0.1], 8)
    print x
    print "Mean = ", Num.sum(x)/8.

if __name__ == '__main__': 
    test()
