from __future__ import nested_scopes

import string, time, array
import rand
import Numeric
import sys
import math
from types import *
import types
Num = Numeric
import scipy.special
special = scipy.special
from scipy_base.fastumath import *
from scipy_base import vectorize
acos = arccos

SequenceType = [types.TupleType, types.ListType, array.ArrayType, Num.ArrayType]

def _check_shape(sh):
   if type(sh) not in SequenceType:
      sh = [sh]
   for val in sh:
      if not isinstance(val,types.IntType):
         raise ValueError, "Each element of the shape parameter must be an integer."
   prod = Num.product(sh)
   return tuple(sh), prod


ArgumentError = "ArgumentError"

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

    _vppf = vectorize(_ppf)
    U = random(size=size)
    return apply(_vppf,(U,)+args)


#################################################
## DISCRETE
##################################################

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


##############################
## Functions from old rv2.py
##############################


class _pranv:
   """Encapsulating class for all random number generator methods and data

   All <_pranv> data and method attributes are private.  The single
   instance <_inst> is created and public aliases to external methods are
   provided at the end of the module.  To get the alias, delete the leading
   underscore from the name of the corresponding method in class <_pranv>.
   There are no public data attributes.  This class should have only one
   instance, <_inst>."""

#                    Internal (private) data attributes:
# Line
# Number
#     _algorithm_name   one-line description of uniform(0,1) algorithm used
#     _count            number of delivered random()s; Python long integer
#     _index            integer index to next uniform(0,1) in ranbuf[]
#     _iterator         a tuple of tuples used by some uniform generators
#     _ln_fac           a double array containing ln(n!) for n = 0, ..., 129
#     _ranbuf           a double array to store batched uniform(0,1) randoms
#     _ranbuf_size      integer, length of _ranbuf
#     _fillbuf          a function variable, aliasing a uniform(0,1) generator
#     _seed             initial Python long integer from user or clock
#     _series           array of pseudo-random series values, or single value
#     _second_nrv       Normal pseudo-randoms are produced in pairs; 2nd here.

#                    Internal (private) utility methods:

#     __init__()            Initialize <_pranv> instance _inst; called once.
#     _build_iterator()     Construct _iterator; allocate _series, _ranbuf.
#     _cmrg(size=None)    Fill ranbuf with random #s from algorithm 'cmrg'.
#     _fill_ln_fac()        Calculate and store ln(n!) for n = 0 through 129.
#     _flip(size=None)    Fill ranbuf with randoms from algorithm 'flip'.
#     _ln_factorial(n)      Return ln(n!), from ln_fac[n] or by calculation.
#     _smrg(size=None)    Fill ranbuf with random #s from algorithm 'sMRG'.
#     _twister(size=None) Fill ranbuf with randoms from algorithm 'twister'.

#            External (private, but aliased public) utility methods:

#     _initial_seed()        Return the seed used to start this random series.
#     _initialize(file_string=None, algorithm='cmrg', seed=0L)        Restart.
#     _random_algorithm()    Return a 1-line description of current generator.
#     _random_count()    Return number of randoms generated in current series.
#     _random(size=None)             Return next pseudo-random uniform(0,1).
#     _save_state(file_string='prvstate.dat')      Save <_inst> state to disk.

#  External (private, but aliased public) non-uniform random variate methods:

#     _choice(seq=(0,1), size=None)       <seq> item
#     _geom(pr_failure=0.5,
#                              size=None) integer     non-negative
#     _hypergeom(tot=35, good=25,
#                   sample=10, size=None) integer     >= max(0, sample-good)
#                                                       <= min(sample, bad)
#     _logser(p=0.5, size=None)      integer     positive
#     _von_Mises(mode=0.0, shape=1.0,
#                              size=None) double      > -pi, < +pi
#     _Wald(mean=1.0, scale=1.0,
#                              size=None) double      positive
#     _Zipf(a=4.0, size=None)             integer     positive

#           External (private, but aliased public) geometrical,
#                permutation, and subsampling routines:

#     _in_simplex(mseq=5*[0.0], bound=1.0)            Return point in simplex.
#     _in_sphere(center=5*[0.0], radius=1.0)           Return point in sphere.
#     _on_simplex(mseq=5*[0.0], bound=1.0)            Return point on simplex.
#     _on_sphere(center=5*[0.0], radius=1.0)           Return point on sphere.
#     _sample(inlist=[0,1,2], sample_size=2)      Return simple random sample.
#     _smart_sample(inlist=[0,1,2], sample_size=2)    Return 'no-dups' sample.


   def __init__(self):
      """Initialize the single class <_pranv> instance.

      """
#  Declaration of class _pranv data attributes, all introduced here mostly for
#  documentation purposes.  Most are allocated or set in _initialize().

      self._algorithm_name = '' # one-line aLgorithm description
      self._count = 0L     #  number of pseudorandoms generated in this series
      self._index = 0      #  index to output buffer array ranbuf[]
      self._iterator = ((),)  # tuple of tuples used by 'flip' and 'twister'.
      self._ln_fac = array.array('d', 130*[0.0])    # ln(n!); n = 0 to 129;
                           #  double array used by some discrete generators.
      self._ranbuf = []    #  uniform(0,1) output buffer array; allocated in
                           #  _build_iterator(), called by _initialize().
      self._ranbuf_size = 0    #  This is just len(_ranbuf), after allocation.
      self._fillbuf = None #  function variable pointing to uniform(0,1) gen.
      self._seed = 0L      #  starting seed from user or system clock
      self._series = []    #  rotating array of random integers, longs or
                           #  doubles; allocated in _build_iterator.
      self._second_nrv = None  # 2nd of generated pair of normal random vars.
#            End of class <_pranv> data attribute declarations

      self._fill_ln_fac()  # Calculate values of  _ln_fac[n] = ln(n!)
                           # for n = 0 to 129; only need to do this once.
      self._initialize()   # Do detailed initialization.  <_initialize()>
                           # may also be called by the user (via its global
                           # alias that drops the leading '_') to change
                           # generators.

   def _build_iterator(self):
      """Construct <_iterator[]> if necessary; allocate _ranbuf and _series.

      _build_iterator()

      'flip' and 'twister', pseudo-random integer generators which manipulate
      long series of integers, utilize a tuple of tuples to simplify and speed
      up the logic.  This function constructs the (global) tuple of tuples,
      <_iterator>, if it's not already there.  The array <_ranbuf> is
      allocated if needed and its length, <_ranbuf_size>, is set.  Also, the
      array <_series> is allocated if necessary.  <_build_iterator> is an
      internal (private) <_pranv> method, called whenever the selected
      uniform(0,1) pseudo-random generator changes, by <_initialize()>."""

      if self._fillbuf == self._cmrg:
         if self._ranbuf_size != 660:    # 'cmrg' does not use _iterator[].
            del self._ranbuf
            self._ranbuf = array.array('d', 660*[0.0])     # output buffer
            self._ranbuf_size = 660
         if len(self._series) != 6:
            del self._series
            self._series = array.array('d', 6*[0.0])

      elif self._fillbuf == self._flip:
         if len(self._iterator) != 55:  # May need to recreate _ranbuf ,
            del self._series            # _iterator and _series.
            self._series = array.array('l',56*[0])
            del self._iterator
            self._iterator = 55 * [()]
            i = 1                  # The 'flip' _iterator is:
            j = 32                 # ( (1,32), (2,33), ..., (24,55), (25,1),
            while j < 56:          #    (26,2), ..., (55,31) ).  Note that
               self._iterator[i-1] = (i, j) # _iterator is global.  A list
               i = i + 1           # here, _iterator converts to tuple below.
               j = j + 1

            j = 1
            while j < 32:
               self._iterator[i-1] = (i, j)
               i = i + 1
               j = j + 1

            self._iterator = tuple(self._iterator)

         if self._ranbuf_size != 660:
            del self._ranbuf
            self._ranbuf = array.array('d', 660*[0.0])  # output buffer
            self._ranbuf_size = 660

      elif self._fillbuf == self._smrg:
         if self._ranbuf_size != 660:    # 'SMRG' does not use _iterator[].
            del self._ranbuf
            self._ranbuf = array.array('d', 660*[0.0])  # output buffer
            self._ranbuf_size = 660

         del self._series
         self._series = 0L               # _series is a single Python long.

      elif self._fillbuf == self._twister:
          if self._ranbuf_size != 624:   # May need to reload _ranbuf,
            del self._ranbuf             # _series and _iterator.
            self._ranbuf = array.array('d', 624*[0.0])
            self._ranbuf_size = 624
            del self._series
            self._series = array.array('L', 624*[0L])
            del self._iterator
            self._iterator = 624 * [()]
                                   #  _iterator is:
            k = 0                  # ( (0,1,397),(1,2,398),...,(226,227,623),
            while k < 227:         # (227,228,0),228,229,1),...,(622,623,395),
               self._iterator[k] =  (k, k + 1, k + 397) # (623,0,386) ).A list
               k = k + 1           # here, _iterator converts to tuple below.

            while k < 623:
               self._iterator[k] =  (k, k + 1, k - 227)
               k = k + 1

            self._iterator[623] =   (623, 0, 396)
            self._iterator = tuple(self._iterator)

      else:
         pass                   # Can't get here.

   def _cmrg(self):
      """Produce batch of uniform(0,1) RVs from L'Ecuyer's MRG32k3a generator

      _cmrg()

      See L'Ecuyer, Pierre., "Good parameters and implementations for combined
      multiple recursive random number generators," May, 1998, to appear in
      "Operations Research."  To download a Postscript copy of the article:
          http://www.iro.umontreal.ca/~lecuyer/myftp/papers/combmrg2.ps
      The period is about 2**191, (or 10**57), and it is well-behaved in all
      dimensions through 45. The generator has been tested extensively by
      L'Ecuyer; no statistical faults were found."""

    # This algorithm is much simpler than the implementation appearing below,
    # which has been hand-optimized to minimize array references. The savings
    # in execution time was about 25%.  Here is the underlying algorithm:
    #
    #     s = self.series                         # local alias
    #
    #     p1 = 1403580.0 * s[1] - 810728.0 * s[0]
    #     k =         floor(p1  /  4294967087.0)  # first modulus, 2**32 - 209
    #     p1 =        p1 - k    *  4294967087.0
    #     if p1 < 0.0: p1 = p1  +  4294967087.0
    #
    #     p2 = 527612.0 * s[5] - 1370589.0 * s[3]
    #     k =         floor(p2  /  4294944443.0)  # 2nd modulus, 2**32 - 22853
    #     p2 =          p2 - k  *  4294944443.0
    #     if p2 < 0.0: p2 = p2  +  4294944443.0
    #     s[3] = s[4]
    #     s[4] = s[5]
    #
    #     dif = p1 - p2
    #     if dif <= 0.0: return (dif + 4294967087 * 2.328306549295728e-10
    #     else:          return               dif * 2.328306549295728e-10
    #

      ranlim = self._ranbuf_size                  # local alias
      buf = self._ranbuf                          # local alias
      s = self._series                            # local alias
      (s0, s1, s2, s3, s4, s5) = s                # Unpack to local aliases.
      i = 0
      while i < ranlim:
         p11 = 1403580.0 * s1 - 810728.0 * s0
         p11 = p11  -   floor(p11 / 4294967087.0) * 4294967087.0
         if p11 < 0.0: p11 = p11 +  4294967087.0  # first mod, 2**32-209

         p2 = 527612.0 * s5 - 1370589.0 * s3
         p2 = p2  -  floor(p2 / 4294944443.0) * 4294944443.0
         if p2 < 0.0: p2 = p2 + 4294944443.0      # 2nd mod, 2**32-22853

         s0 = p11
         s3 = p2
         dif = p11 - p2
         if dif < 0: dif = dif + 4294967087.0
         buf[i] = dif * 2.328306549295728e-10
         i = i + 1

         p1 = 1403580.0 * s2 - 810728.0 * s1
         p1 = p1  -  floor(p1 / 4294967087.0) * 4294967087.0
         if p1 < 0.0: p1 = p1 + 4294967087.0      # first mod, 2**32-209

         p2 = 527612.0 * p2 - 1370589.0 * s4
         p2 = p2  -  floor(p2 / 4294944443.0) * 4294944443.0
         if p2 < 0.0: p2 = p2 + 4294944443.0      # 2nd mod, 2**32-22853

         s1 = p1
         s4 = p2
         dif = p1 - p2
         if dif < 0: dif = dif + 4294967087.0
         buf[i] = dif * 2.328306549295728e-10
         i = i + 1

         p1 = 1403580.0 * p11 - 810728.0 * s2
         p1 = p1  -  floor(p1 / 4294967087.0) * 4294967087.0
         if p1 < 0.0: p1 = p1 + 4294967087.0      # first mod, 2**32-209

         p2 = 527612.0 * p2 - 1370589.0 * s5
         p2 = p2  -  floor(p2 / 4294944443.0) * 4294944443.0
         if p2 < 0.0: p2 = p2 + 4294944443.0      # 2nd mod, 2**32-22853

         s2 = p1
         s5 = p2
         dif = p1 - p2
         if dif < 0: dif = dif + 4294967087.0
         buf[i] = dif * 2.328306549295728e-10
         i = i + 1

      s[0] = s0                                # Save the coefficients
      s[1] = s1                                # for the next batch of
      s[2] = s2                                # 660 uniform randoms.
      s[3] = s3                                # s = [s0,s1,s2,s3,s4,s5]
      s[4] = s4                                # doesn't work here because
      s[5] = s5                                # s is an array, not a list.
      self._index = 0
      self._count = self._count + ranlim


   def _fill_ln_fac(self):
      """Calculate and store array values of  <_ln_fac[n]> = ln(n!).

      _fill_ln_fac()

      This internal routine is called only once, by  __init__() in
      class <_pranv> during the instantiation of <_inst>."""
      self._ln_fac[0] = 0.0
      i = 0
      sum = 0.0
      while i < 129:
         i = i + 1
         sum = sum + log( float(i) )
         self._ln_fac[i] = sum


   def _flip(self):
      """Produce a batch of uniform(0,1) RVs with Knuth's 1993 generator.

      _flip()

      See Knuth, D.E., "The Stanford GraphBase: A Platform for Combinatorial
      Computing," ACM Press, Addison-Wesley, 1993, pp 216-221.  The period for
      the underlying pseudo-random integer series is at least 3.6e16, and may
      be nearly 4e25. The low-order bits of the integers are just as random
      as the high-order bits. The requirements are 32-bit integers and two's
      complement arithmetic.  The recurrence is quite fast, but for sensitive
      applications, autocorrelations or other statistical defects may dictate
      the use of a more complex generator, a longer series or both."""
      ranlim = self._ranbuf_size                  # local alias
      buf = self._ranbuf                          # local alias
      s = self._series                            # local alias
            # Here is the algorithm, C-language style:
            #    ii = 1             # Refresh series of 55 random integers.
            #    jj = 32
            #    while jj < 56:     # Calculate (s[ii] - s[ii+31]) % 2**31.
            #       s[ii] = (s[ii] - s[jj]) & 0x7fffffff
            #       ii = ii + 1
            #       jj = jj + 1
            #    jj = 1
            #    while jj < 32:     # Calculate (s[ii] - s[ii-24]) % 2**31.
            #       s[ii] = (s[ii] - s[jj]) & 0x7fffffff
            #       ii = ii + 1
            #       jj = jj + 1
            #
            # And here it is using <_iterator>:

      i = 0
      while i < ranlim:
         for (ii, jj) in self._iterator:
            s[ii] = (s[ii] - s[jj]) & 0x7fffffff
            buf[i] = ( float(s[ii]) + 1.0 ) * 4.6566128709089882e-10
            # The constant is 1/(2**31+1), so stored value is in (0,1).
            i = i + 1

      self._index = 0
      self._count = self._count + ranlim


   def _ln_factorial(self, n=0.0):
      """Return ln(n!)from the array <_ln_fac[n]> or use Stirling's formula.

      _ln_factorial(n=0.0)

      For 0 <= n <= 129, the value pre-calculated and stored in <_ln_fac[]>
      by internal routine <_fill_ln_fac()> is returned.  For n >= 130, the
      Stirling asymptotic approximation with 3 terms is used.
      Approximation error is negligible.  NOTE:  n must be non-negative
      and should be integral.  There is no error-checking in this function
      because it is internal, and we don't make mistakes inside this class."""

      if n < 130:
         return self._ln_fac[int(n)]
      else:
         y = n + 1.0
         yi = 1.0 / y
         yi2 = yi * yi
         return ((((  0.793650793650793650e-3) * yi2  #  1 / (1260 * (n+1)**5)
                    - 0.277777777777777778e-2) * yi2  # -1 / ( 360 * (n+1)**3)
                    + 0.833333333333333333e-1) * yi   # +1 / ( 8   * (n+1)   )
                    + 0.918938533204672742            # +(log(sqrt(2*math.pi))
                    + (y - 0.5) * log(y) - y )

   def _smrg(self):
      """Produce a batch of uniform(0,1) RVs from Wu's mod 2**61 - 1 generator

      _smrg()

      See Wu, Pei-Chi, "Multiplicative, Congruential Random-Number Generators
      with multiplier (+ or -) 2**k1 (+ or -) 2**k2 and Modulus 2**p - 1, ACM
      Transactions on mathematical Software, June, 1997, Vol. 23, No. 2,
      pp 255 - 265.  The generator has modulus 2**61 - 1, which is the
      Mersenne prime immediately following 2**31 - 1, and multiplier
      37**458191 % (2**61 - 1). Because 37 is the minimal primitive root of
      2**61 -1, the generator has full period (2**61 -2, about 2.3e18).
      It was found by Wu to have the best performance on spectral tests
      (through dimension 8) of any multiplier of the type 37**k % (2**61 - 1)
      with k in the range from 1 to 1,000,000.
      Generated pseudo-random numbers range from 4.337e-19 to 1.0 - 2**53
      -- all bits are random, an advantage over 31- or 32-bit generators. """
      buf = self._ranbuf                          # local alias
      s = self._series                            # local alias
      ranlim = self._ranbuf_size                  # local alias
      for i in xrange(ranlim):
         s = s * 2137866620694229420L % 0x1fffffffffffffffL  # The first
                                 # constant is 37**458191 % (2**61 - 1),
                                 # and the second is 2**61-1.
         buf[i] = float(s) * 4.3368086899420168e-19 # For s = 2**61 - 2,
                                                    # this is 1 - 2**-53.
      self._index = 0
      self._series = s           # save the current seed in <_series>
      self._count = self._count + ranlim

   def _twister(self):
      """Produce a batch of uniform(0,1) RVs from the Mersenne Twister MT19937

      _twister()

      See M. matsumoto and T. Nishamura, "Mersenne Twister," ACM Transactions
      on Modeling and Computer Simulation, Jan. 1998, vol. 8, no. 1, pp 3-30.
      The period is 2**19937 - 1,(> 1e6000); the generator has a 623 dim-
      ensional equi-distributional property.  This means that every sequence
      of less than 624 pseudo-random 32-bit integers occurs equally often. It
      has passed a number of statistical tests, including those in Marsaglia's
      "Diehard" package. """
      buf = self._ranbuf                          # local alias
      s = self._series                            # local alias
    # Here is the first part of the algorithm, C-language style:
    # if i == 624:             # self._index points to long ints in _series[].
    #    k = 0                 # Generate 624 new pseudo-random long integers.
    #    while k < 227:        # 227 = 624 - 397
    #       y = (s[ k ] & 0x80000000L) | (s[k + 1] & 0x7fffffffL)
    #       s[k] = s[k + 397] ^ (y >> 1) ^ (y & 0x1) * 0x9908b0dfL
    #       k = k + 1
    #
    #    while k < 623:
    #       y = (s[ k ] & 0x80000000L) | (s[k + 1] & 0x7fffffffL)
    #       s[k] = s[k - 227] ^ (y >> 1) ^ (y & 0x1) * 0x9908b0dfL
    #       k = k + 1
    #
    #    y =    (s[623] & 0x80000000L) | (s[  0  ] & 0x7fffffffL)
    #     s[623] =  s[  396  ] ^ (y >> 1) ^ (y & 0x1) * 0x9908b0dfL
    #
    # And here is is in Python, using the global iterator, just three lines:
      for (k, kp1, koff) in self._iterator:
         y = (s[k] & 0x80000000L) | (s[kp1] & 0x7fffffffL)
         y = s[koff] ^ (y >>1) ^ (y & 0x1) * 0x9908b0dfL
         s[k] = y
         y = y ^ (y >> 11)
         y = y ^ (y <<  7) & 0x9d2c5680L
         y = y ^ (y << 15) & 0xefc60000L
         y = y ^ (y >> 18)
         buf[k] = (float(y) + 1.0) * 2.3283064359965952e-10
                                 # The constant is 1.0 / (2.0**32 +1.0).
      self._index = 0
      self._count = self._count + self._ranbuf_size

#            External (private; aliased public) utility methods:

   def _initialize(self, file_string=None, algorithm='cmrg', seed=0L):
      """Set uniform algorithm and seed, or restore state from disk.

      Called by <__init__()> or by user, via the alias <initialize()>.
      if <file_string> is specified, it can be a null string, indicating
      that class _pranv state is not to be retrieved from a file; a single
      blank or the string 'default', indicating that _pranv state is to be
      restored from the file with the default file name ('prvstate.dat'); or
      the path/file name of a file containing a _pranv state.  <algorithm> may
      be 'cmrg', 'flip', 'smrg', or 'twister'; see the documentation for
      guidance.  <seed> is by default 0, indicating the system clock is to
      be used to determine algorithm starting value(s), or a Python long
      integer.  For some of the algorithms, a negative seed has a special
      meaning; see the code below.  <algorithm> and <seed are ignored if
      <file_string is non-null."""
      if file_string:                 # User wants to restore generator state.
         if (file_string == ' ' or    # a space (not the null string)
          string.strip( string.lower(file_string) ) == 'default'):
            file_string = 'prvstate.dat'
         f = open(file_string, 'r')   # Get saved state of uniform generator.
         inlist = f.readlines()
         f.close()
         self._algorithm_name = inlist[1]    # Skip the header (inlist[0]).
         self._algorithm_name = self._algorithm_name[:-1]  # Drop the \n.
         self._seed = eval(inlist[2])
         self._count = eval(inlist[3])
         self._second_nrv = eval(inlist[4])
         self._index = eval(inlist[5])

         pos = string.find(self._algorithm_name, ":")
         if algorithm == 'cmrg':
            self._fillbuf = self._cmrg       # Set generator function pointer.

         elif algorithm == 'flip':
            self._fillbuf = self._flip

         elif algorithm == 'smrg':
            self._fillbuf = self._smrg

         elif algorithm == 'twister':
            self._fillbuf = self._twister

         else:
            pass                     # can't get here

         self._build_iterator()      # Construct  _iterator, and allocate
                                     # _series and _ranbuf if necessary.
         if type(self._series) == array.ArrayType:
            offset = 6
            for i in xrange( len(self._series) ):         # 'CMRG', 'Flip', or
               self._series[i] = eval(inlist[offset + i]) # 'Twister'

            offset = len(self._series) + 6
            for i in xrange(self._ranbuf_size):
               self._ranbuf[i] = eval(inlist[offset + i])

         else:
            self._series = eval( inlist[6] )              # 'SMRG'
            offset = 7
            for i in xrange(self._ranbuf_size):
               self._ranbuf[i] = eval(inlist[offset + i])

      else:                            # User wants a new generator.
         self._second_nrv = None       # Initialize storage for 2nd normal RV.
         self._count = 0L              # Initialize random() delivery count.
         algorithm = string.strip( string.lower(algorithm) )
         if algorithm == 'cmrg':
            self._algorithm_name = \
        "'CMRG': Combined Multiplicative Recursion  MRG32k3a (L'Ecuyer, 1998)"
            self._fillbuf = self._cmrg # Set generator function pointer.
            self._build_iterator()     # Allocate _series and _ranbuf.
            if seed < 0:
               self._seed = long(seed)                # User wants all initial
               seed = float((-seed) & 0x7fffffffL)    # values equal.
               for i in range(6): self._series[i] = seed

            else:
               if seed == 0:          # Generate initial seed from epoch time.
                  t = (long(time.time()*128) * 1048576L) & 0x7fffffffffffffffL
                  seed = ((t & 0x7fffffffL)^(t >> 32)) & 0x7fffffffL # 31 bits
                  seed = max(seed, 1L) # <seed> must be positive in CMRG
                  seed = min(seed, 4294944442L)  # 2**32 - 22853 - 1

               self._seed = long(seed) # Save for user inquiry.
               self._series[0] = float(seed)
               for j in range(1,6):    # Use a standard linear multiplicative
                  k = seed / 127773    # congruential generator to get initial
                  seed = 16807 * (seed - k * 127773) - 2836 * k  # values for
                  if seed < 0: seed = seed + 0x7fffffff # other 5 multipliers.
                  self._series[j] = float(seed)

            self._index = self._ranbuf_size-1  # Request new batch of randoms.

         elif algorithm == 'flip':
            self._algorithm_name = \
               "'Flip': Subtractive Series Mod 2**31 (Knuth,1993)"
            self._fillbuf = self._flip  # Set generator function pointer.
            self._build_iterator()      # define_iterator;_ranbuf and _series
            if seed == 0L:              # Get seed from epoch time.
               t = (long(time.time() * 128) * 1048576L) & 0x7fffffffffffffffL
               seed = int( (t & 0x7fffffffL)^(t >> 32) ) & 0x7fffffffL # 31 b.
               seed = max(seed, 1)      # Insure seed is positive.

            self._seed = long(seed)  # Save for possible user inquiry.
            seed = seed & 0x7fffffff # Need positive 4-byte integer here.
            s = self._series         # local alias
            s[0] = -1                # <_series[0]> is a sentinel unused here.
            prev = seed              # <_series[]> contains signed integers.
            next = 1
            s[55] = prev
            i = 21
            while i:
               s[i] = next
               next = (prev - next) & 0x7fffffff     # difference mod 2**31
               if seed & 0x1: seed = 0x40000000 + (seed >> 1)
               else: seed = seed >> 1                # cyclic shift right
               next = (next - seed) & 0x7fffffff     # difference mod 2**31
               prev = s[i]
               i = (i + 21) % 55                     # Jump arround in array.

            self._index = self._ranbuf_size-1 # Request new batch of uniforms.
            self._random()                    # Exercize the generator by
            self._index = self._index + 219   # "using" 4 sets of 55 integers.

         elif algorithm == 'smrg':
            self._algorithm_name = \
               "'SMRG': Single Multiplicative Recursion m61-p3019 (Wu, 1997)"
            self._fillbuf = self._smrg       # Set generator function pointer.
            self._build_iterator()           # Allocate _ranbuf and _series
            if seed == 0L:
               t = (long( time.time() * 128 ) * 544510892281524561475242L)
               t = t & 0x3ffffffffffffffffffffffffffffffL  # Keep 122 bits.
               seed = ( t ^ (t >> 61) ) & 0x1fffffffffffffffL # 61 bits
                                 # The time() multiplier above is arbitrary
                                 # so long as t >= (but close to) 2**122.
            self._seed = seed    # Save for user inquiry.
            seed = abs(seed)     # Insure seed is positive.
            seed = min(seed, 2L**61 - 2)  # 2**61 - 1 is the modulus
            self._series = seed  #_series is not an array here; a Python long.
            self._index = self._ranbuf_size-1 # Request new batch of uniforms.
            self._random()       # Exercize the generator by "using" 16
            self._index = self._index + 15    # pseudo-random numbers.

         elif algorithm == 'twister':
            self._algorithm_name = \
         "'Twister': Mersenne Twister MT19937 (matsumoto and Nishamura, 1998)"
            self._fillbuf = self._twister # Set generator function pointer.
            self._build_iterator()        # Define _iterator;_ranbuf, _series.
            if seed == 0L:
                                          #  Construct a seed from epoch time.
                                          #  32-bit unsigned integers needed.
               t = (long(time.time() * 128) * 2097152L) & 0xffffffffffffffffL
               seed = ( (t & 0xffffffffL)^(t >> 32) ) & 0xffffffffL # 32 bits

            self._series[0] = seed   # <_series[]> contains unsigned integers.
            seed = abs(seed)
            self._seed = long(seed)  # Save the seed for later user inquiry.
            seed = seed & 0xffffffffL

            # Use an algorithm from Knuth(1981) to continue filling <_series>
            # with 623 more unsigned 32-bit pseudo-random integers.

            for k in xrange(1, 624):
               seed = (69069L * seed) & 0xffffffffL    # Retain lower 32 bits.
               self._series[k] = seed

            self._index = self._ranbuf_size - 1
                                # Force generation of first batch of 624
                                # pseudo-randoms.  <_index> points to the
                                # pseudo-random uniform just delivered from
                                # <_ranbuf>; none remain when <_index> is 623.

         else: raise ValueError, \
           "algorithm? --must be 'CMRG', 'Flip', 'SMRG', or 'Twister'"


   def _initial_seed(self):
      """Return seed for this series; was either user-supplied or from clock.

      initial_seed()

      The result is a Python long integer."""
      return self._seed

   def _random_algorithm(self):
      """Return 1-line description of uniform random number algorithm used.

      random_algorithm()

      The result is a string."""
      return self._algorithm_name

   def _random_count(self):
      """Return number of uniform(0,1) RVs used since the seed was last reset.

      random_count()

      A Python long integer is returned.  This the number of uniform(0,1)
      random numbers delivered to the user or consumed by _pranv methods."""
      return self._count - (self._ranbuf_size - 1 - self._index)


   def _save_state(self, file_string='prvstate.dat'):
      """Save <_pranv> state to disk for later recall and continuation.

      save_state(file_string='prvstate.dat')

      A backup file ('filename.bak' or 'prvstate.bak') is created if
      the file <file_string> can be opened, where <filename*> is
      the leading portion of <file_string> before the '.', if present."""

      try:
         f = open(file_string, 'r')     # If file exists, create backup file.
         temp = f.read()
         dot_pos = string.rfind(file_string,'.')
         if dot_pos == -1: bak_file_string = file_string + '.bak'
         else: bak_file_string = file_string[:dot_pos] + '.bak'
         f.close()
         f = open(bak_file_string, 'w') # Write backup file.
         f.write(temp)
         f.close()
      except IOError:                   # No file was found.
         pass
      f = open(file_string, 'w')        # Open file for current state.
                                        # Save state.  Arrays can't be
      outlist = []                      # pickled, so just write out the data.
      outlist.append('Python module rv 1.1 random numbers save_state: '
         + time.ctime( time.time() ) + '\n')  # save_state file header
      outlist.append(self._algorithm_name + '\n')
      outlist.append(`self._seed` + '\n')
      outlist.append(`self._count` + '\n')
      outlist.append(`self._second_nrv` + '\n')
      outlist.append(`self._index` + '\n')
      if type(self._series) is array.ArrayType:
         for i in xrange( len(self._series) ): # 'CMRG', 'Flip', and 'Twister'
            outlist.append(`self._series[i]` + '\n')

      else:                                          # In 'SMRG', _series is a
         outlist.append(`self._series` + '\n')       # single Python long.

      for i in xrange (self._ranbuf_size):     # buffer of uniform(0,1) RVs
            outlist.append(`self._ranbuf[i]` + '\n')

      f.writelines(outlist)
      f.close()

#             Uniform(0,1) Pseudo-random Variate Generator

   def _random(self, size=None):
      """Return a single random uniform(0,1), or <None> and fill buffer.

      """
      i = self._index                        # local alias
      if size is not None:
         size, Ns = _check_shape(size)
         buffer = Num.zeros(Ns,Num.Float)
         buf = self._ranbuf                  # local alias
         ranlim = self._ranbuf_size          # local alias
         for j in xrange( Ns ):
            i = i + 1
            if i == ranlim: self._fillbuf(); i = 0
            buffer[j] = buf[i]
         return Num.reshape(buffer, size)

      else:
         i = i + 1
         if i == self._ranbuf_size: self._fillbuf(); i = 0
         self._index = i
         return self._ranbuf[i]

      self._index = i       # Restore _ranbuf index and return (buffer filled).

   def _choice(self, seq=(0,1), size=None):
      """Return element k with probability 1/len(seq) from non-empty sequence.

      choice(seq=(0,1), size=None)

      Default is a 0, 1 coin flip.  <seq> must not be empty.  If <buffer> is
      specified, it must be a mutable sequence.  It is filled with random
      <seq> elements and <None> is returned.  Otherwise, a single <seq>
      element is returned."""
      lenseq = len(seq)
      if lenseq == 0:
         raise ValueError, '<seq> must not be empty'

      else:
         i = self._index                          # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            buflim = self._ranbuf_size            # local alias
            buf = self._ranbuf                    # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               x = buf[i]
               buffer[j] = seq[int( floor(x * lenseq) )]
            return Num.reshape(buffer, size)

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            x = self._ranbuf[i]
            self._index = i          # Restore updated value of _ranbuf index.
            return seq[int( floor(x * lenseq) )]

      self._index = i  # Update _ranbuf index and return (buffer[] is filled).

   def _geom(self, pr=0.5, size=None):
      """Return non-negative random integers from a geometric distribution.

      geom(pr=0.5, size=None)

      Z has the geometric distribution if it is the number of successes
      before the first failure in a series of independent Bernouli trials
      with probability of success 1 - <pr>.  0 <= pr < 1.
      The result is a non-negative integer less than or equal to 2**31 -1.
      If <buffer> is specified, it must be a mutable sequence.
      <buffer> is filled with geometric(pr) pseudo-randoms.
      Otherwise, a single geometric pseudo-random variate is returned."""

      pr = 1.0-pr    # Added to be consistent with distributions.py 
      if not 0.0 <= pr < 1.0:
         raise ValueError, '0.0 <=  <pr>  < 1.0'
      else:
         i = self._index                  # local alias
         buf = self._ranbuf               # local alias
         buflim = self._ranbuf_size       # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns            
            j = 0

         else:
            out_len = 0
            j = -1

         if pr <= 0.9:            # <pr> is small or moderate;
            while j < out_len:            # use inverse transformation.
               pk = sum = 1.0 - pr
               successes = 0
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               u = buf[i]
               while sum < u:
                  successes = successes + 1
                  pk = pk * pr
                  sum = sum + pk

               if out_len:
                  buffer[j] = successes
                  j = j + 1

               else:
                  self._index = i    # Restore updated value of _ranbuf index.
                  return successes

         else:                       # <pr> larger than 0.9.
            while j < out_len:
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               u = buf[i]
               successes = floor( log(u) / log(pr) )
               successes = int( min(successes, 2147483647.0) )     # 2**31 - 1
               if out_len:
                  buffer[j] = successes
                  j = j + 1

               else:
                  self._index = i    # Restore updated value of _ranbuf index.
                  return successes

         self._index = i # Restore _ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)         

   def _hypergeom(self, tot=35, good=25, N=10, size=None):
      """Return hypergeometric pseudorandom variates: #"bad" in <sample>

      hypergeom(tot=35, good=25, N=10)

      Z has a hypergeometric distribution if it is the number of "bad"
      (tot-good) items in a sample of size <N> from a population of
      size tot items.  <tot>, <good> and <N> must be positive integers.  Also,
      <tot> >= <sample> >= 1.  The result  Z satisfies:
      max(0, <sample> - <good>) <= Z <= min(<sample>, <tot>-<good>).
      if size is not None is
      supplied, it must be a mutable sequence (list or array).  It is filled
      with hypergeometric pseudo-random integers and <None> is returned.
      Otherwise, a single hypergeometric pseudo-random integer is returned.
      See Fishman, "Monte Carlo," pp 218-221.  Algorithms HYP and HRUA* were
      used."""
      bad = tot - good
      alpha = int(bad); beta=int(good); n = int(N)
      if (alpha < 1) or (beta < 1) or (n < 1) or (alpha + beta) < n:
         raise ValueError, '<bad>, <good>, or <sample> out of range'

      else:
         i = self._index                  # local alias
         buf = self._ranbuf               # local alias
         buflim = self._ranbuf_size       # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns            
            j = 0

         else:
            out_len = 0
            j = -1

         if n <= 10:                      # 10 is an empirical estimate.
            d1 = alpha + beta - n         # Use Fishman's HYP algorithm.
            d2 = float( min(alpha, beta) )
            while j < out_len:
               y = d2
               k = n
               while y > 0.0:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u = buf[i]
                  y = y - floor( u + y / (d1 + k) )
                  k = k - 1
                  if k == 0:
                     break

               z = int(d2 - y)
               if alpha > beta: z = n - z
               if out_len:
                  buffer[j] = z
                  j = j + 1

               else:
                  self._index = i
                  return z

         else:                            # Use Fishman's HRUA* algorithm.
            minalbe = min(alpha, beta)
            popsize = alpha + beta
            maxalbe = popsize - minalbe
            lnfac = self._ln_factorial    # function alias
            m = min(n, popsize - n)
            d1 = 1.715527770              # 2 * sqrt(2 / math.e)
            d2 = 0.898916162              # 3 - 2 * sqrt(3 / math.e)
            d4 = float(minalbe) / popsize
            d5 = 1.0 - d4
            d6 = m * d4 + 0.5
            d7 = sqrt((popsize - m) * n * d4 * d5 / (popsize - 1) + 0.5)
            d8 = d1 * d7 + d2
            d9 = floor( float((m + 1) * (minalbe + 1)) / (popsize + 2) )
            d10 = ( lnfac(d9) + lnfac(minalbe - d9) +
                    lnfac(m - d9) + lnfac(maxalbe - m + d9) )
            d11 = min( min(m, minalbe) + 1.0, floor(d6 + 7 * d7) )  # 7 for 9-
            while j < out_len:                               # digit precision
               while 1:                                      # in d1 and d2
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  x = buf[i]
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  y = buf[i]
                  w = d6 + d8 * (y - 0.5) / x
                  if w < 0.0 or w >= d11:
                     continue             # fast rejection; try another x, y

                  else:
                     z = floor(w)
                     t = d10 - ( lnfac(z) + lnfac(minalbe - z) +
                                 lnfac(m - z) + lnfac(maxalbe - m + z) )
                     if x * (4.0 - x) - 3.0 <= t:
                        break             # fast acceptance

                     else:
                        if x * (x - t) >= 1:
                           continue       # fast rejection, try another x, y

                        else:
                           if 2.0 * log(x) <= t:
                              break       # acceptance

               if alpha > beta:           # Error in HRUA*, this is correct.
                  z = m - z

               if m < n:                  # This fix allows n to exceed
                  z = alpha - z           # popsize / 2.

               if out_len:
                  buffer[j] = z
                  j = j + 1

               else:
                  self._index = i
                  return z

         self._index = i # Restore ranbuf index and return (buffer is filled).
         return Num.reshape(buffer, size)
      

   def _logser(self, pr=0.5, size=None):
      """Return logarithmic (log-series) positive pseudo-random integers;
      0<p<1.

      logser(pr=0.5, size=None)

      Z has the logseries(logarithmic)distribution if Pr{Z=i} = (a/i) * p ** i,
      for i = 1, 2, ... and a = -1.0 / log(1.0 - p). If <buffer> is specified,
      it must be a mutable sequence (list or array).  It is filled with
      logseries (positive) pseudo-random integers and <None> is returned.
      Otherwise, a single logarithmic series pseudo-random is returned.
      The algorithm is by A.W. Kemp; see Devroye, 1986, pp 547-548."""
      p = pr
      if not(0.0 < p < 1.0):
         raise ValueError, '<p> must be in (0.0, 1.0)'

      else:
         r = log(1.0 - p)
         i = self._index                  # local alias
         buf = self._ranbuf               # local alias
         buflim = self._ranbuf_size       # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns            
            j = 0

         else:
            out_len = 0
            j = -1

         if p < 0.95:
            pr_one = -p / r
            while j < out_len:
               sum = pr_one
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               u = buf[i]
               k = 1
               while u > sum:
                  u = u - sum
                  k = k + 1
                  sum = sum * p * (k - 1) / k

               if out_len:
                  buffer[j] = k
                  j = j + 1

               else:
                  self._index = i
                  return k

         else:                                                 # p is >= 0.95.
            while j < out_len:
               k = 1
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               v = buf[i]
               if v < p:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u = buf[i]
                  q = 1.0 - exp(r * u)
                  q2 = q * q
                  if v <= q2 :
                     k = int( floor(1.0 + log(v) / log(q)) )

                  elif (q2 < v <= q):
                     k = 1

                  else:
                     k = 2

               if out_len:
                  buffer[j] = k
                  j = j + 1

               else:
                  self._index = i
                  return k

         self._index = i  # Restore ranbuf index and return (<buffer> filled).
         return Num.reshape(buffer, size)         

   def _von_Mises(self, b, loc=0.0, size=None):
      """Return von Mises distribution pseudo-random variates on [-pi, +pi].

      von_Mises(b, loc=0.0, size=None)

      <mean> must be in the open interval (-math.pi, +math.pi).  If <buffer>
      is specified, it must be a mutable sequence (list or array).  It is
      filled with von Mises(mean, shape) pseudo-random variates and <None> is
      returned.  Otherwise, a single von Mises RV is returned. The method is
      an algorithm of Best and Fisher, 1979; see Fisher, N. I., "Statistical
      Analysis of Circular Data," Cambridge University Press, 1995, p. 49."""

      shape, mean = b, loc
      z = exp(1j*mean)
      mean = arctan2(z.imag, z.real)
      if not (-3.1415926535897931 < mean < +3.1415926535897931):
         raise ValueError, \
            '<loc> must be in the open interval (-math.pi, math.pi)'

      else:
         a = 1.0 + sqrt( 1.0 + 4.0 * shape*shape )
         b = ( a - sqrt(a + a) ) / (shape + shape)
         r = (1.0 + b * b) / (b + b)
         i = self._index                          # local alias
         buf = self._ranbuf                       # local alias
         buflim = self._ranbuf_size               # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns            
            j = 0

         else:
            out_len = 0
            j = - 1

         while j < out_len:
            while (1):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               z = cos(3.1415926535897931 * buf[i])
               f = (1.0 + r * z) / (r + z)
               c = shape * (r - f)
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               u = buf[i]
               if (c * (2.0 - c) - u > 0.0):
                  break                           # quick acceptance

               if log(c / u) + 1.0 - c < 0.0:
                  continue                        # quick rejection

               else:
                  break                           # acceptance

            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            if out_len:
               if buf[i] > 0.5:
                  buffer[j] = fmod(acos(f) + mean, 6.2831853071795862)

               else:
                  buffer[j] = -fmod(acos(f) + mean, 6.2831853071795862)

               j = j + 1

            else:
               self._index = i # Restore updated value of _ranbuf index.
               if buf[i] > 0.5: return  fmod(acos(f)+mean, 6.2831853071795862)
               else:            return -fmod(acos(f)+mean, 6.2831853071795862)

         self._index = i # Restore _ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)
               
   def _Wald(self, mu, loc=0.0, scale=1.0, size=None):
      """Return Inverse gaussian pseudo-random variates.

      invnorm(mu, loc=0.0, scale=1.0, size=None)

      Both <mean> and <scale> must be positive.  If <buffer> is specified, it
      must be a mutable sequence (list or array).  It is filled with inverse
      gaussian (wald if mu=1.0) pseudo-random variates, and <None> is returned.
      Otherwise, a single Wald pseudo-random variate is returned."""
      mean = mu
      oscale = scale
      scale = 1.0
      if (mean <= 0.0) or (scale <= 0.0):
         raise ValueError, 'both <mean> and <scale> must be positive'

      else:
         i = self._index                          # local alias
         buf = self._ranbuf                       # local alias
         buflim = self._ranbuf_size               # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns            
            j = 0

         else:
            out_len = 0
            j = -1

         normal = self._normal                    # local alias
         r = 0.5 * mean / scale
         con = 4.0 * scale
         meansq = mean * mean
         while j < out_len:
            n = normal()
            muy = mean * n * n
            x1 = mean + r * ( muy -  sqrt(muy * (con + muy)) )
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            if out_len:
               if buf[i] <= mean / (mean + x1):
                  buffer[j] = x1

               else:
                  buffer[j] =  meansq / x1

               j = j + 1

            else:
               self._index = i       # Restore updated value of _ranbuf index.
               if buf[i] <= mean / (mean + x1):
                  return x1

               else:
                  return meansq / x1

      self._index = i    # Restore _ranbuf index and return (buffer[] filled).
      return Num.reshape(buffer*oscale+loc, size)      
      
   def _Zipf(self, a=4.0, size=None):
      """Return positive pseudo-random integers from the Zipf distribution

      Zipf(a=4.0, size=None)

      Z has the Zipf ditribution with parameter a if
         Pr{Z=i} = 1.0 / (zeta(a) * i ** a) for i = 1, 2, ... and a > 1.0.
      Here zeta(a) is the Riemann zeta function, the sum from 1 to infinity of
      1.0 / i **a.  If <buffer> is specified, it must be a mutable sequence
      (list or array).  It is filled with Zipf pseudo-random integers and
      <None> is returned.  Otherwise, a single (positive) Zipf pseudo-random
      integer is returned; see Devroye, 1986, pp 550-551 for the algorithm."""
      if a <= 1.0:
         raise ValueError, '<a> must be larger than 1.0'

      else:
         i = self._index                  # local alias
         buf = self._ranbuf               # local alias
         buflim = self._ranbuf_size       # local alias
         am1 = a - 1.0
         b = 2.0 ** am1
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns
            j = 0

         else:
            out_len = 0
            j = -1

         while j < out_len:
            while 1:
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               x = floor( buf[i] ** (-1 / am1) )
               t = (1.0 + 1.0 / x) ** am1
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               if buf[i] * x * (t - 1.0) / (b - 1.0)  <= t / b:
                  break

            if out_len:
               buffer[j] = int(x)
               j = j + 1

            else:
               self._index = i
               return int(x)

         self._index = i  # Restore ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)         

#               End of Non-uniform(0,1) Generators

#               Geometric and Subsampling Routines

   def _in_simplex(self, mseq=5*[0.0], bound=1.0):
      """Return a pseudorandom point in a simplex.

      in_simplex(mseq=2*[0.0], bound=1.0)

      The simplex of dimension d and bound b > 0 is defined by all points
      z = (z1, z2, ..., zd) with z1 + z2 + ... + zd <= b and all zi >= 0.
      <bound> must be positive.
      A two-dimensional simplex with bound b is the region enclosed by the
      horizontal (x) and vertical (y) axes and the line y = b - x.  A
      one dimensional simplex with bound b is the single point b.
      <mseq> must be a mutable sequence (list or array); if it is an array, it
      must be typed double or single float.  A sequence of the same type and
      length as <mseq> containing the coordinates of a pseudo-random point in
      the simplex is returned.
      See Fishman, George, "Monte Carlo," Springer-Verlag, 1996, pp. 232-233
      and 226-229.  Algorithms EUNIFORM and OSIMP were used."""
      d = len(mseq)
      point = mseq[:]
      if bound <= 0.0:
         raise ValueError, '<bound> must be positive.'

      elif d == 1:
         point[0] = bound

      else:
         i = self._index                          # local alias
         buf = self._ranbuf                       # local alias
         buflim = self._ranbuf_size               # local alias
         range_d = range(d)
         i = i + 1
         if i == buflim: self._fillbuf(); i = 0
         point[0] = -log(buf[i])
         for j in range_d[1:]: # Accumulate running sum of exponential(1) RVs.
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            point[j] = point[j-1] - log(buf[i])

         i = i + 1
         if i == buflim: self._fillbuf(); i = 0
         total = point[d-1] - log(buf[i])
         for j in range_d:
            point[j] = point[j] / total

         point[0] = point[0] * bound
         j = d - 1
         while j > 0:
            point[j] = bound * (point[j] - point[j-1])
            j = j - 1

      self._index = i     # Restore updated value of _ranbuf index and return.
      return point


   def _in_sphere(self, center=5*[0.0], radius=1.0):
      """Return a pseudo-random point in a sphere centered at <center>.

      in_sphere(center=5*[0.0], radius=1.0)

      A one-dimensional sphere centered at 0 is the two points <+radius> and
      <-radius>.  A two-dimensional sphere centered at 0 is the circle with
      radius <radius>.  <center> must be a sequence containing double
      coordinates for the sphere's center. <radius> must be positive.
      pseudo-random point in the "sphere" with center 0 and radius <radius>.
      See Fishman, George, "Monte Carlo," Springer-Verlag, 1996, pp 234-235.
      Algorithm OSPHERE was used."""
      d = len(center)
      if radius <= 0.0:
         raise ValueError, '<radius> must be positive.'

      else:
         return self._on_sphere( center, radius * self._beta(float(d), 1.0) )

   def _on_simplex(self, mseq=5*[0.0], bound=1.0):
      """Return a pseudo_random point on the boundary of a simplex.

      on_simplex(mseq=5*[0.0], bound=1.0)

      The simplex of dimension d and bound b > 0 is defined by all vectors
      z = (z1, z2, ..., zd) with z1 + z2 + ... + zd <= b and all zi >= 0.
      <mseq> must be a mutable sequence (list or array) and <bound> must be
      positive.  A boundary point on the simplex has the sum of its
      coefficients equal to b.
      The coordinates of a pseudo-random point on the boundary of the
      simplex of dimension len(seq) and bound <bound> are returned as a
      sequence of the same type and length as <mseq>.  The input values in
      <mseq> are ignored, and undisturbed.
      See Fishman, George, "Monte Carlo," Springer-Verlag, 1996,
      pp 232-233 and 226-229.  Algorithms EUNIFORM and OSIMP were used."""
      d = len(mseq)
      point = mseq[:]
      if bound <= 0.0:
         raise ValueError, '<bound> must be positive.'

      elif d == 1:                         # The point is just <bound>.
         point[0] = bound
      else:
         i = self._index                           # local alias
         buf = self._ranbuf                        # local alias
         buflim = self._ranbuf_size                # local alias
         range_d = range(d)                        # [0,1,...,d-1]
         range_dm1 = range_d[:(d - 1)]             # [0,1,...,d-2]
         range_1_dm1 = range_dm1[1:]               # [1,2,...,d-2]
         i = i + 1
         if i == buflim: self._fillbuf(); i = 0
         point[0] = -log(buf[i])
         for j in range_1_dm1:             # range_1_dm1 is empty for d = 2.
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            point[j] = point[j - 1] - log(buf[i])

         i = i + 1
         if i == buflim: self._fillbuf(); i = 0
         total = point[d - 2] - log(buf[i])
         for j in range_dm1:
            point[j] = point[j] / total

         last = point[d - 2]
         point[0] = bound * point[0]
         j = d - 2
         while j > 0:                      # not executed for d = 2.
            point[j] = bound * (point[j] - point[j - 1])
            j = j - 1

         point[d - 1] = bound - last

      self._index = i     # Restore updated value of _ranbuf index and return.
      return point

   def _on_sphere(self, center=5*[0.0], radius=1.0):
      """Return a pseudo-random point on the sphere centered at point[].

      on_sphere(center=5*[0.0], radius=1.0)

      The sphere of dimension d centered at 0 is all points (x1,x2,...,xd)
      with x1**2 + x2**2 + ... + xd**2 <= <radius>**2.  The surface of this
      sphere is any point with the sum of its components squared equal to
      the radius squared.  A 2-dimensional sphere is a circle, and a
      1-dimensional sphere is the line from <-radius> to <+radius>.  <center>
      must be a sequence and <radius> must be positive.  A sequence of the
      same type as <center[]> containing the coordinates of a pseudo-random
      point on the "sphere" centered at <center[]> is returned.
      See Fishman, George, "Monte Carlo," Springer-Verlag, 1996, pp 234-235.
      Algorithm OSPHERE was used.  For dimensions 2, 3, and 4, rejection
      algorithms were used. See Marsaglia, G., "Choosing a point from the
      surface of a sphere," Annals of mathematical Statistics, vol. 43,
      pp. 645-646, 1972."""
      d = len(center)
      buflim = self._ranbuf_size        # local alias
      buf = self._ranbuf                # local alias
      i = self._index                   # local alias
      point = center[:]
      if  radius <= 0.0:
         raise ValueError, '<radius> must be positive.'

      elif d == 0:                      # Just return center[:].
         pass

      elif d == 1:
         i = i + 1
         if i == buflim: self._fillbuf(); i = 0
         if buf[i] > 0.5:
            point[0] = radius + center[0]

         else:
            point[0] = -radius + center[0]

      elif d == 2:                      # rejection from circumscribed square.
         ss = 2.0
         while ss > 1.0:
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            x = buf[i]
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            y = buf[i]
            x = x + x - 1.0                # (x, y) is in the square circum-
            y = y + y - 1.0                # scribed on the unit circle
            ss = x * x + y * y             # centered at (0, 0).

         ss = radius / sqrt(ss)            # (x, y) now restricted to unit
         point[0] = x * ss + center[0]     # circle.
         point[1] = y * ss + center[1]
         self._index = i                   # Restore ranbuf index and return.

      elif d == 3:                         # algorithm by Marsaglia, 1972
         ss = 2.0
         while ss > 1.0:
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            u = buf[i]
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            v = buf[i]
            u = u + u - 1.0
            v = v + v - 1.0
            ss = u * u + v * v

         con = sqrt(1.0 - ss) * radius
         point[0] = 2.0 * u * con + center[0]
         point[1] = 2.0 * v * con + center[1]
         point[2] = (1.0 - ss - ss) * radius + center[2]
         self._index = i                   # Restore ranbuf index and return.

      elif d == 4:
         ssuv = 2.0                        # algorithm by Marsaglia, 1972
         while ssuv > 1.0:
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            u = buf[i]
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            v = buf[i]
            u = u + u - 1.0
            v = v + v - 1.0
            ssuv = u * u + v * v

         ssrs = 2.0
         while ssrs > 1.0:
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            r = buf[i]
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            s = buf[i]
            r = r + r - 1.0
            s = s + s - 1.0
            ssrs = r * r + s * s

         con = sqrt((1.0 - ssuv) / ssrs) * radius
         point[0] = u * radius + center[0]
         point[1] = v * radius + center[1]
         point[2] = r * con + center[2]
         point[3] = s * con + center[3]
         self._index = i                   # Restore ranbuf index and return.

      else:                                # Use radially symmetric normals
         range_d = range(d)                # (Fishmans algorithm OSPHERE).
         sum = 0.0
         normal = self._normal
         for j in range_d:
            z = normal()
            sum = sum + z * z
            point[j] = z

         constant = radius / sqrt(sum)
         for j in range_d:
           point[j] = constant * point[j] + center[j]

      return point

   def _sample(self, inlist=[0,1,2], sample_size=2):
      """Return simple random sample of <sample_size> items from <inlist>.

      sample(inlist=[0,1,2], sample_size=2)

      <inlist> must be a mutable sequence (list or array).  The
      <sample_size> <= len(inlist) elements of the returned mutable sequence
      (which has the same type as <inlist>) are a simple random sample from
      all the elements in <inlist>.
      Sampling is with replacement; duplicates may occur in the returned
      sequence.  The length of the returned list or array is <sample_size>,
      which may be any positive integer."""
      pop_size = len(inlist)
      n = int(sample_size)
      if pop_size == 0:
         outlist = inlist[:]

      elif n <= 0:
         outlist = 0 * inlist[:1]

      else:
         i = self._index                          # local alias
         buf = self._ranbuf                       # local alias
         buflim = self._ranbuf_size               # local alias
         outlist = n*inlist[:1]                   # <sample_size> >= 1
         for j in xrange(n):
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            outlist[j] = inlist[int( floor(buf[i] * pop_size) )]

      self._index = i                     # Restore _ranbuf index arnd return.
      return outlist

   def _smart_sample(self, inlist=[0,1,2], sample_size=2):
      """Return a random sample without replacement of elements in <inlist>.

      smart_sample(inlist=[0,1,2], sample_size=2)

      <inlist> must be a mutable sequence (list or array).  The
      <sample_size> <= len(inlist) elements of the returned mutable sequence
      (which has the same type as <inlist>) are a sample without replacement
      of all the elements in <inlist>.  Since sampling is without replacement,
      there is no possibility of duplicates in the returned list or array
      unless there are duplicates in <inlist>.  The length of the
      returned mutable sequence is the smaller of <sample_size> and
      len(inlist), and its type is the type of <inlist>."""
      pop_size = len(inlist)
      n = int( min(sample_size, pop_size) )  # can't have sample > population.
      if pop_size <= 1:
         outlist = inlist[:]

      elif n < 1:
         outlist = 0*inlist[:1]

      else:
         outlist = inlist[:]            # <outlist> has same type as <inlist>.
         i = self._index                # local alias
         buf = self._ranbuf             # local alias
         buflim = self._ranbuf_size     # local alias
         for j in xrange(n):
            temp = outlist[j]           # Save item in ith position.
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            k = j + int( floor(buf[i] * (pop_size - j)) )
            outlist[j] = outlist[k]     # Move sampled item to ith position.
            outlist[k] = temp           # copy to sampled item's position.

         if n < pop_size:
            del outlist[n:pop_size]     # Remove any unsampled elements.

      self._index = i                   # Restore _ranbuf index and return.
      return outlist

_inst = _pranv()   # Initialize uniform(0,1) generator to default; clock seed.

initial_seed      =  _inst._initial_seed
initialize        =  _inst._initialize
random_algorithm  =  _inst._random_algorithm
random_count      =  _inst._random_count
save_state        =  _inst._save_state

choice            =  _inst._choice

#     Geometrical Point Generators, Permutations and Subsampling Routines
in_simplex        =  _inst._in_simplex
in_sphere         =  _inst._in_sphere
on_simplex        =  _inst._on_simplex
on_sphere         =  _inst._on_sphere
sample            =  _inst._sample
smart_sample      =  _inst._smart_sample


##def mean_var_test(x, type, mean, var, skew=[]):
##    n = len(x) * 1.0
##    x_mean = Num.sum(x)/n
##    x_minus_mean = x - x_mean
##    x_var = Num.sum(x_minus_mean*x_minus_mean)/(n-1.0)
##    print "\nAverage of ", len(x), type
##    print "(should be about ", mean, "):", x_mean
##    print "Variance of those random numbers (should be about ", var, "):", x_var
##    if skew != []:
##       x_skew = (Num.sum(x_minus_mean*x_minus_mean*x_minus_mean)/9998.)/x_var**(3./2.)
##       print "Skewness of those random numbers (should be about ", skew, "):", x_skew

##def test():
##    x, y = get_seed()
##    print "Initial seed", x, y
##    seed(x, y)
##    x1, y1 = get_seed()
##    if x1 != x or y1 != y:
##        raise SystemExit, "Failed seed test."
##    print "First random number is", random()
##    print "Average of 10000 random numbers is", Num.sum(random(10000))/10000.
##    x = random([10,1000])
##    if len(x.shape) != 2 or x.shape[0] != 10 or x.shape[1] != 1000:
##        raise SystemExit, "random returned wrong shape"
##    x.shape = (10000,)
##    print "Average of 100 by 100 random numbers is", Num.sum(x)/10000.
##    y = uniform(0.5,0.6, (1000,10))
##    if len(y.shape) !=2 or y.shape[0] != 1000 or y.shape[1] != 10:
##        raise SystemExit, "uniform returned wrong shape"
##    y.shape = (10000,)
##    if Num.minimum.reduce(y) <= 0.5 or Num.maximum.reduce(y) >= 0.6:
##        raise SystemExit, "uniform returned out of desired range"
##    print "randint(1, 10, size=[50])"
##    print randint(1, 10, size=[50])
##    print "permutation(10)", permutation(10)
##    print "randint(3,9)", randint(3,9)
##    print "random_integers(10, size=[20])"
##    print random_integers(10, size=[20])
##    s = 3.0
##    x = norm(2.0, s, [10, 1000])
##    if len(x.shape) != 2 or x.shape[0] != 10 or x.shape[1] != 1000:
##        raise SystemExit, "standard_normal returned wrong shape"
##    x.shape = (10000,)
##    mean_var_test(x, "normally distributed numbers with mean 2 and variance %f"%(s**2,), 2, s**2, 0)
##    x = exponential(3, 10000)
##    mean_var_test(x, "random numbers exponentially distributed with mean %f"%(s,), s, s**2, 2)
##    x = multivariate_normal(Num.array([10,20]), Num.array(([1,2],[2,4])))
##    print "\nA multivariate normal", x
##    if x.shape != (2,): raise SystemExit, "multivariate_normal returned wrong shape"
##    x = multivariate_normal(Num.array([10,20]), Num.array([[1,2],[2,4]]), [4,3])
##    print "A 4x3x2 array containing multivariate normals"
##    print x
##    if x.shape != (4,3,2): raise SystemExit, "multivariate_normal returned wrong shape"
##    x = multivariate_normal(Num.array([-100,0,100]), Num.array([[3,2,1],[2,2,1],[1,1,1]]), 10000)
##    x_mean = Num.sum(x)/10000.
##    print "Average of 10000 multivariate normals with mean [-100,0,100]"
##    print x_mean
##    x_minus_mean = x - x_mean
##    print "Estimated covariance of 10000 multivariate normals with covariance [[3,2,1],[2,2,1],[1,1,1]]"
##    print Num.matrixmultiply(Num.transpose(x_minus_mean),x_minus_mean)/9999.
##    x = beta(5.0, 10.0, 10000)
##    mean_var_test(x, "beta(5.,10.) random numbers", 0.333, 0.014)
##    x = gamma(.01, 2., 10000)
##    mean_var_test(x, "gamma(.01,2.) random numbers", 2*100, 2*100*100)
##    x = chi_square(11., 10000)
##    mean_var_test(x, "chi squared random numbers with 11 degrees of freedom", 11, 22, 2*Num.sqrt(2./11.))
##    x = F(5., 10., 10000)
##    mean_var_test(x, "F random numbers with 5 and 10 degrees of freedom", 1.25, 1.35)
##    x = poisson(50., 10000)
##    mean_var_test(x, "poisson random numbers with mean 50", 50, 50, 0.14)
##    print "\nEach element is the result of 16 binomial trials with probability 0.5:"
##    print binomial(16, 0.5, 16)
##    print "\nEach element is the result of 16 negative binomial trials with probability 0.5:"
##    print negative_binomial(16, 0.5, [16,])
##    print "\nEach row is the result of 16 multinomial trials with probabilities [0.1, 0.5, 0.1 0.3]:"
##    x = multinomial(16, [0.1, 0.5, 0.1], 8)
##    print x
##    print "Mean = ", Num.sum(x)/8.

##if __name__ == '__main__': 
##    test()
