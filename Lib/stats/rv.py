import rand
import Numeric
import sys
import math
from types import *

ArgumentError = "ArgumentError"

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
        n = Numeric.multiply.reduce(size)
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

def uniform(minimum, maximum, size=None):
    """Returns array of given shape of random reals in given range (exclusive of enpoints)
    """
    return _build_random_array(rand.uniform, (minimum, maximum), size)

def randint(minimum, maximum=None, size=None):
    """random integers >=min, < max.  If max not given, random integers >= 0, <min"""
    if maximum is None:
        maximum = minimum
        minimum = 0
    a = Numeric.floor(uniform(minimum, maximum, size))
    if isinstance(a, Numeric.ArrayType):
        return a.astype(Numeric.Int)
    else:
        return int(a)
     
def random_integers(maximum, minimum=1, size=None):
    """random_integers(max, min=1, size=None) = random integers in range min-max inclusive"""
    return randint(minimum, maximum+1, size) 
     
def permutation(arg):
    """If arg is an integer, a permutation of indices arange(n), otherwise
    a permuation of the sequence"""
    if isinstance(arg,IntType):
        arg = Numeric.arange(arg)
    return rand.permutation(arg)

def standard_normal(size=None):
    """standard_normal(n) or standard_normal([n, m, ...]) returns array of
           random numbers normally distributed with mean 0 and standard
           deviation 1"""
    return _build_random_array(rand.standard_normal, (), size)

def normal(mean, std, size=None):
        """normal(mean, std, n) or normal(mean, std, [n, m, ...]) returns
           array of random numbers randomly distributed with specified mean and
           standard deviation"""
        return _build_random_array(rand.normal, (mean, std), size)
    
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
           n = Numeric.product(size)
       output = rand.multivariate_normal(mean, cov, n)
       if size is not None:
           final_shape = list(size[:])
           final_shape.append(mean.shape[0])
           output.shape = final_shape
       return output       

def standard_exp(size=None):
    """returns array of random numbers with a standard exponential distribution (mean=1)"""
    return _build_random_array(rand.standard_exp, (), size)

def exponential(mean, size=None):
    """returns array of random numbers exponentially distributed with specified mean"""
    return _build_random_array(rand.exponential, (av,), size)

def beta(a, b, size=None):
    """returns array of beta distributed random numbers."""
    return _build_random_array(rand.beta, (a, b), size)

def standard_gamma(a, size=None):
    """returns array of random numbers with a standard gamma distribution with mean a"""
    return _build_random_array(rand.standard_gamma, (a,), size)

def gamma(a, r, size=None):
    """returns array of gamma distributed random numbers."""
    return _build_random_array(rand.gamma, (a, r), size)

def f(dfn, dfd, size=None):
    """returns array of F distributed random numbers with dfn degrees of freedom in the numerator and dfd degrees of freedom in the denominator."""
    return _build_random_array(rand.f, (dfn, dfd), size)

def noncentral_f(dfn, dfd, nonc, size=None):
    """returns array of noncentral F distributed random numbers with dfn degrees of freedom in the numerator and dfd degrees of freedom in the denominator, and noncentrality parameter nonc."""
    return _build_random_array(rand.noncentral_f, (dfn, dfd, nonc), size)

def chi2(df, size=None):
    """returns array of chi squared distributed random numbers with df degrees of freedom."""
    return _build_random_array(rand.chi2, (df,), size)

def noncentral_chi2(df, nonc, size=None):
    """returns array of noncentral chi squared distributed random numbers with df degrees of freedom and noncentrality parameter."""
    return _build_random_array(rand.noncentral_chi2, (df, nonc), size)

def t(df, size=None):
    """returns array of noncentral student_t distributed random numbers with df degrees of freedom."""
    return standard_normal(size)*Numeric.sqrt(df) / Numeric.sqrt(chi2(df,size))

def noncentral_t(df, nonc, size=None):
    return normal(nonc)*Numeric.sqrt(df) / Numeric.sqrt(noncentral_chi2(df,nonc,size))

def binomial(trials, p, size=None):
    """returns array of binomially distributed random integers.

           trials is the number of trials in the binomial distribution.
           p is the probability of an event in each trial of the binomial distribution."""
    return _build_random_array(rand.binomial, (trials, p), size)

def negative_binomial(trials, p, size=None):
    """returns array of negative binomially distributed random integers.

           trials is the number of trials in the negative binomial distribution.
           p is the probability of an event in each trial of the negative binomial distribution."""
    return _build_random_array(rand.negative_binomial, (trials, p), size)

def multinomial(trials, probs, size=None):
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
    probs = Numeric.array(probs)
    if len(probs.shape) != 1:
        raise ArgumentError, "probs must be 1 dimensional."
        # Compute shape of output
    if type(size) == type(0): size = [size]
    final_shape = size[:]
    final_shape.append(probs.shape[0]+1)
    x = rand.multinomial(trials, probs.astype(Numeric.Float32), Numeric.multiply.reduce(size))
        # Change its shape to the desire one
    x.shape = final_shape
    return x

def poisson(mean, size=None):
    """returns array of poisson distributed random integers with specifed mean."""
    return _build_random_array(rand.poisson, (mean,), size)


def mean_var_test(x, type, mean, var, skew=[]):
    n = len(x) * 1.0
    x_mean = Numeric.sum(x)/n
    x_minus_mean = x - x_mean
    x_var = Numeric.sum(x_minus_mean*x_minus_mean)/(n-1.0)
    print "\nAverage of ", len(x), type
    print "(should be about ", mean, "):", x_mean
    print "Variance of those random numbers (should be about ", var, "):", x_var
    if skew != []:
       x_skew = (Numeric.sum(x_minus_mean*x_minus_mean*x_minus_mean)/9998.)/x_var**(3./2.)
       print "Skewness of those random numbers (should be about ", skew, "):", x_skew

def test():
    x, y = get_seed()
    print "Initial seed", x, y
    seed(x, y)
    x1, y1 = get_seed()
    if x1 != x or y1 != y:
        raise SystemExit, "Failed seed test."
    print "First random number is", random()
    print "Average of 10000 random numbers is", Numeric.sum(random(10000))/10000.
    x = random([10,1000])
    if len(x.shape) != 2 or x.shape[0] != 10 or x.shape[1] != 1000:
        raise SystemExit, "random returned wrong shape"
    x.shape = (10000,)
    print "Average of 100 by 100 random numbers is", Numeric.sum(x)/10000.
    y = uniform(0.5,0.6, (1000,10))
    if len(y.shape) !=2 or y.shape[0] != 1000 or y.shape[1] != 10:
        raise SystemExit, "uniform returned wrong shape"
    y.shape = (10000,)
    if Numeric.minimum.reduce(y) <= 0.5 or Numeric.maximum.reduce(y) >= 0.6:
        raise SystemExit, "uniform returned out of desired range"
    print "randint(1, 10, size=[50])"
    print randint(1, 10, size=[50])
    print "permutation(10)", permutation(10)
    print "randint(3,9)", randint(3,9)
    print "random_integers(10, size=[20])"
    print random_integers(10, size=[20])
    s = 3.0
    x = normal(2.0, s, [10, 1000])
    if len(x.shape) != 2 or x.shape[0] != 10 or x.shape[1] != 1000:
        raise SystemExit, "standard_normal returned wrong shape"
    x.shape = (10000,)
    mean_var_test(x, "normally distributed numbers with mean 2 and variance %f"%(s**2,), 2, s**2, 0)
    x = exponential(3, 10000)
    mean_var_test(x, "random numbers exponentially distributed with mean %f"%(s,), s, s**2, 2)
    x = multivariate_normal(Numeric.array([10,20]), Numeric.array(([1,2],[2,4])))
    print "\nA multivariate normal", x
    if x.shape != (2,): raise SystemExit, "multivariate_normal returned wrong shape"
    x = multivariate_normal(Numeric.array([10,20]), Numeric.array([[1,2],[2,4]]), [4,3])
    print "A 4x3x2 array containing multivariate normals"
    print x
    if x.shape != (4,3,2): raise SystemExit, "multivariate_normal returned wrong shape"
    x = multivariate_normal(Numeric.array([-100,0,100]), Numeric.array([[3,2,1],[2,2,1],[1,1,1]]), 10000)
    x_mean = Numeric.sum(x)/10000.
    print "Average of 10000 multivariate normals with mean [-100,0,100]"
    print x_mean
    x_minus_mean = x - x_mean
    print "Estimated covariance of 10000 multivariate normals with covariance [[3,2,1],[2,2,1],[1,1,1]]"
    print Numeric.matrixmultiply(Numeric.transpose(x_minus_mean),x_minus_mean)/9999.
    x = beta(5.0, 10.0, 10000)
    mean_var_test(x, "beta(5.,10.) random numbers", 0.333, 0.014)
    x = gamma(.01, 2., 10000)
    mean_var_test(x, "gamma(.01,2.) random numbers", 2*100, 2*100*100)
    x = chi_square(11., 10000)
    mean_var_test(x, "chi squared random numbers with 11 degrees of freedom", 11, 22, 2*Numeric.sqrt(2./11.))
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
    print "Mean = ", Numeric.sum(x)/8.

if __name__ == '__main__': 
    test()
