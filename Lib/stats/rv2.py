"""Module <rv>:  Random Numbers From Uniform and Other Distributions

                          Release 1.1  27-Nov-98

External (public) globals, which are aliases for public method attributes in
class <_pranv>, are listed below in the order in which their parent
methods appear (alphabetically within groups) in the class <_pranv> source:

                            U t i l i t i e s
------------------------------------------------------------------------------
initialize(file_string=None, algorithm='cmrg', seed=0L).....Restart generator.
initial_seed()........Return the seed used to start the current random series.
random_algorithm()...........Return a 1-line description of current generator.
random_count().......Return the number of randoms generated in current series.
save_state(file_string='prvstate.dat')......Save class <_pranv> state to disk.
------------------------------------------------------------------------------

       U n i f o r m   R a n d o m   N u m b e r   G e n e r a t o r
------------------------------------------------------------------------------
random2(size=None).....Return float from the chosen generator('CMRG', 'Flip',
    'SMRG' or 'Twister'), or return <None> and fill <buffer> if it was
    supplied.  The expected mean and standard deviation from any of the
    uniform(0,1) generators in the test program below are 0.500 and 0.289,
    respectively.  See the HTML documentation for more information on
    individual generator characteristics.
------------------------------------------------------------------------------

   N o n  -  U n i f o r m   R a n d o m   N u m b e r   G e n e r a t o r s
------------------------------------------------------------------------------
beta(mu=4.0, nu=2.0, size=None)..................default mean 0.67, SD 0.18.
binomial(trials=25, pr_success=0.5, size=None)...default mean 12.5, SD 2.50.
cauchy(median=0.0, scale=1.0, size=None).................(no moments exist).
chi2(df=10.0, size=None)...................default mean 10.0, SD 4.50.
choice(seq=(0,1), size=None).....................default mean 0.50, SD 0.50.
exponential(scale=1.0, size=None)................default mean 1.00, SD 1.00.
fratio(numdf=2.0, denomdf=10.0, size=None).....default mean 1.25, SD 1.60.
gamma(mu=2.0, size=None).........................default mean 2.00, SD 1.40.
geom(pr=0.5, size=None).............default mean 1.00, SD 1.40.
gumbel(mode=0.0, scale=1.0, size=None)...........default mean 0.58, SD 1.30.
hypergeom(tot=35, good=25, sample=10, size=None)deflt mn 2.86, SD 1.22
laplace(mu=0.0, scale=1.0, size=None)............default mean 0.00, SD 1.40.
logseries(p=0.5, size=None)....................default mean 1.44, SD 0.90
logistic(mu=0.0, scale=1.0, size=None)...........default mean 0.00, SD 1.80.
lognormal(mu=0.0, sigma=1.0, size=None)..........default mean 1.65, SD 2.20.
negative_binomial(r=1.0, pr=0.5, size=None)  deflt mn 1.00, SD 1.40.
normal(mu=0.0, sigma=1.0, size=None).............default mean 0.00, SD 1.00.
pareto(mode=1.0, shape=4.0, size=None)...........default mean 1.33, SD 0.47.
poisson(rate=1.0, size=None).....................default mean 1.00, SD 1.00.
randint(lowint=0, upint=1, size=None)............default mean 0.50, SD 0.50.
rayleigh(scale=1.0, size=None)....................default mean 1.25, SD 0.65.
student_t(df=100.0, size=None)...................default mean 0.00, SD 1.01.
triangular(left=0.0, mode=0.5, right=1.0, size=None) deflt mn 0.50, SD 0.20.
uniform(lower=-0.5, upper=0.5, size=None)........default mean 0.00, SD 0.29.
von_mises(mean=0.0, shape=1.0, size=None)........default mean 0.00, SD 1.30.
invnorm(mu, loc=0.0, scale=1.0, size=None).......
wald(loc=0.0, scale=1.0, size=None).............default mean 1.00, SD 1.00.
weibull(scale=1.0, shape=0.5, size=None).........default mean 2.00, SD 4.50.
zipf(a=4.0, size=None)...........................default mean 1.11, SD 0.54
------------------------------------------------------------------------------

 Geometrical Point Generators, and Routines for Sampling and Permutations
------------------------------------------------------------------------------

in_simplex(mseq=5*[0.0], bound=1.0)............Return random point in simplex.
in_sphere(center=5*[0.0], radius=1.0).........Return random point in "sphere".
on_simplex(mseq=5*[0.0], bound=1.0)............Return random point on simplex.
on_sphere(center=5*[0.0], radius=1.0).........Return random point on "sphere".
permutation(elements=[0,1,2,3,4]).....Return random permutation of <elements>.
sample(inlist=[0,1,2], sample_size=2).............Return simple random sample.
smart_sample(inlist=[0,1,2], sample_size=2)....Return "no dups" random sample.
------------------------------------------------------------------------------


           Copyright 1998 by Ivan Frohne; Wasilla, Alaska, U.S.A.
                          All Rights Reserved

Permission to use, copy, modify and distribute this software and its
documentation for any purpose, free of charge, is granted subject to the
following conditions:
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the software.

   THE SOFTWARE AND DOCUMENTATION IS PROVIDED WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO MERCHANTABILITY, FITNESS
   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHOR
   OR COPYRIGHT HOLDER BE LIABLE FOR ANY CLAIM OR DAMAGES IN A CONTRACT
   ACTION, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
   SOFTWARE OR ITS DOCUMENTATION.
------------------------------------------------------------------------------
"""

# Changed Feb-2002  by Travis E. Oliphant  <oliphant.travis@ieee.org>
#  Altered interface so that functions take size argument and return an
#  array of that size.


#  Changed  27-Nov-98       begun 19-Sep-98  (Release 1.1 changes)
#  Author:  Ivan Frohne
#           HC31 Box 5236E
#           Wasilla, AK 99654-9802
#           frohne@alaska.net

import Numeric as Num
import types
import string, time, array
from math import sqrt, log, floor, acos, fmod, cos, exp, tan

SequenceType = [types.TupleType, types.ListType, array.ArrayType, Num.ArrayType]

def _check_shape(sh):
   if type(sh) not in SequenceType:
      sh = [sh]
   for val in sh:
      if not isinstance(val,types.IntType):
         raise ValueError, "Each element of the shape parameter must be an integer."
   prod = Num.product(sh)
   return tuple(sh), prod

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

#                                           returns...  range...
#     ------------------------------------  ----------  ----------------
#     _beta(mu=4.0, nu=2.0, size=None)    double      > 0.0, < 1.0
#     _binomial(trials=25, pr_success=0.5,
#                              size=None) integer     >= 0, <= trials
#     _Cauchy(median=0.0, scale=1.0,
#                              size=None) double      unbounded
#     _chi_square(df=10.0, size=None)     double      positive
#     _choice(seq=(0,1), size=None)       <seq> item
#     _exponential(scale=1.0, size=None)  double      positive
#     _Fisher_F(numdf=2.0, denomdf=10.0,
#                              size=None) double      positive
#     _gamma(mu=2.0, size=None)           double      positive
#     _geom(pr_failure=0.5,
#                              size=None) integer     non-negative
#     _Gumbel(mode=1.0, scale=1.0,
#                              size=None) double      unbounded
#     _hypergeom(tot=35, good=25,
#                   sample=10, size=None) integer     >= max(0, sample-good)
#                                                       <= min(sample, bad)
#     _Laplace(mu=0.0, scale=1.0,
#                              size=None) double      unbounded
#     _logser(p=0.5, size=None)      integer     positive
#     _logistic(mu=0.0, scale=1.0,
#                              size=None) double      unbounded
#     _lognormal(mu=0.0, sigma=1.0,
#                              size=None) double      positive
#     _negative_binomial(r=1.0,
#              pr_failure=0.5, size=None) integer     non-negative
#     _normal(mu=0.0, sigma=1.0,
#                              size=None) double      unbounded
#     _Pareto(mode=1.0, shape=4.0,
#                              size=None) double      > mode
#     _Poisson(rate=1.0, size=None)       integer     non-negative
#     _randint(lowint=0, upint=1,
#                              size=None) integer     >= lowint, <= upint
#     _Rayleigh(mode=1.0, size=None)      double      positive
#     _Student_t(df=100.0, size=None)     double      unbounded
#     _triangular(left=0.0, mode=0.5,
#                   right=1.0, size=None) double      > left, < right
#     _uniform(lower=-0.5, upper=+0.5,
#                              size=None) double      > lower, < upper
#     _von_Mises(mode=0.0, shape=1.0,
#                              size=None) double      > -pi, < +pi
#     _Wald(mean=1.0, scale=1.0,
#                              size=None) double      positive
#     _Weibull(scale=1.0, shape=0.5,
#                              size=None) double      positive
#     _Zipf(a=4.0, size=None)             integer     positive

#           External (private, but aliased public) geometrical,
#                permutation, and subsampling routines:

#     _in_simplex(mseq=5*[0.0], bound=1.0)            Return point in simplex.
#     _in_sphere(center=5*[0.0], radius=1.0)           Return point in sphere.
#     _on_simplex(mseq=5*[0.0], bound=1.0)            Return point on simplex.
#     _on_sphere(center=5*[0.0], radius=1.0)           Return point on sphere.
#     _permutation(elements=[0,1,2,3,4])     Return permutation of <elements>.
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
      Transactions on Mathematical Software, June, 1997, Vol. 23, No. 2,
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

      See M. Matsumoto and T. Nishamura, "Mersenne Twister," ACM Transactions
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
         "'Twister': Mersenne Twister MT19937 (Matsumoto and Nishamura, 1998)"
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

#             Non-Uniform(0,1) Pseudo-random Variate Generators

   def _beta(self, mu=4.0, nu=2.0, size=None):
      """Return incomplete beta function random variates; mean mu/(mu+nu).

      beta(mu=4.0, nu=2.0, size=None)

      Both <mu> and <nu> must be positive.  <buffer>, if specified, must be a
      mutable sequence (list or array).  It is filled with beta(mu,nu)
      pseudo-random variates and <None> is returned.  If <buffer> is not
      specified, a single beta(mu,nu) variate is returned.  See Fishman,
      "Monte Carlo," pp 200-296, algorithms AW and BB*; and Devroye, "Non-
      uniform Random Variate Generation," pp 428-443, for the second algorithm
      of Atkinson and Whittaker (1976, 1979) on p. 443. """
      if (mu <= 0.0) or (nu <= 0.0):
         raise ValueError, 'both <mu> and <nu> must be positive'

      else:
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns
            j = 0

         else:
            out_len = 0
            j = -1

         i = self._index                 # local alias
         buflim = self._ranbuf_size      # local alias
         buf = self._ranbuf              # local alias
         if nu > mu:
            maxmunu, minmunu = nu, mu

         else:
            maxmunu, minmunu = mu, nu

         sum = mu + nu

         if maxmunu < 1.0:           # Use Fishman's algorithm AW
            one_less_mu = 1.0 - mu
            one_less_nu = 1.0 - nu
            t = 1.0 / ( 1.0 + sqrt(nu * one_less_nu / (mu * one_less_mu)) )
            one_less_t = 1.0 - t
            nut = nu * t
            p = nut / ( nut + mu * one_less_t )
            invmu = 1.0 / mu
            invnu = 1.0 / nu
            while j < out_len:
               while 1:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u = buf[i]
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  y = -log(buf[i])       # exponential(1) RV
                  if u <= p:
                     z = t * (u / p) ** invmu
                     if y >= one_less_nu * (t - z) / one_less_t:
                        break            # fast acceptance

                     if y >= one_less_nu * log((1.0 - z) / one_less_t):
                        break            # acceptance

                  else:
                     z = 1.0 - one_less_t * ( (1.0 - u) / (1.0 - p) ) ** invnu
                     if y >= one_less_mu * (z / t - 1.0):
                        break            # fast acceptance

                     if y >= one_less_mu * log(z / t):
                        break            # acceptance

               if out_len:
                  buffer[j] = z
                  j = j + 1

               else:
                  self._index = i
                  return z

         elif minmunu >= 1.0:        # Use Fishman's algorithm BB*
            d4 = sqrt( (sum - 2.0) / (2.0 * mu * nu - sum) )
            d5 = minmunu + 1.0 / d4
            while j < out_len:
               while 1:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u1 = buf[i]
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u2 = buf[i]
                  v = d4 * log( u1 / (1.0 - u1) )
                  w = minmunu * exp(v)
                  z1 = u1 * u1 * u2
                  r = d5 * v - 1.38629436
                  s = minmunu + r - w
                  if s + 2.60943791 > 5.0 * z1:
                     break           # fast acceptance 1, ln(4) = 1.386...
                                     # and 1 + ln(5) = 2.609...
                  else:
                     t = log(z1)
                     if s >= t:
                        break        # fast acceptance

                     else:
                        if r + sum * log( sum / (maxmunu + w) ) >= t:
                           break     # acceptance

               if minmunu == mu:
                  z = w / (maxmunu + w)

               else:
                  z = maxmunu / (maxmunu + w)

               if out_len:
                  buffer[j] = z
                  j = j + 1

               else:
                  self._index = i
                  return z

         else:                          # mu and nu straddle 1.0; use the 2nd.
            inv_minmunu = 1.0 / minmunu # algorithm of Atkinson and Whittaker.
            inv_maxmunu = 1.0 / maxmunu
            c = (1.0 - minmunu)
            t = c / (c + maxmunu)
            onemt = 1.0 - t
            c = maxmunu * t
            p = c / (c + minmunu * onemt ** maxmunu)
            while j < out_len:
               while 1:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u = buf[i]
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  y = -log(buf[i])
                  if u <= p:
                     x = t * (u / p) ** inv_minmunu
                     if (1.0 - maxmunu) * log(1.0 - x) <= y:
                        break

                  else:
                     x = 1.0 - onemt * ((1.0 - u) / (1.0 - p)) ** inv_maxmunu
                     if (1.0 - minmunu) * log(x / t) <= y:
                        break

               if minmunu == nu:
                  x = 1.0 - x

               if out_len:
                  buffer[j] = x
                  j = j + 1

               else:
                  self._index = i
                  return x

         self._index = i    # Restore ranbuf index and return (buffer filled).
         return Num.reshape(buffer,size)

   def _binomial(self, trials=25, pr_success=0.5, size=None):
      """Return integer pseudo-random variates from a binomial distribution.

      binomial(trials=25, pr_success=0.5, size=None)

      <trials> must be an integer >= 1, and <pr_success> must be in the closed
      interval [0.0, 1.0].  <buffer>, if specified, must be a mutable sequence
      (list or array).  It is filled with binomial(trials, pr_success) pseudo-
      random variates (non-negative integers) and <None> is returned. If
      <buffer> is not specified, a single binomial pseudo-random integer
      is returned.
      See Fishman, George S., "Monte Carlo Concepts, Algorithms and
      Applications," Springer-Verlag, New York, 1996.  Algorithms ITR
      (pp. 154-155) and BRUA* (pp. 215-218) are used.  Mean computing time
      is bounded."""
      n = int(trials)                         # n is a (positive) integer.
      np1 = n + 1.0
      d3 = min(pr_success, 1.0 - pr_success)  # 0.0 <= d3 <= 0.5
      d4 = n * d3 + 0.5
      d5 = 1.0 - d3
      d8 = d3 / d5
      buflim  = self._ranbuf_size             # local alias
      i = self._index                         # local alias
      buf = self._ranbuf                      # local alias
      if (n <= 0.0) or not (0.0 <= pr_success <= 1.0):
         raise ValueError, \
          '<trials> must be a positive integer; 0.0 <= <pr_success> <= 1.0'

      elif d4 <= 10.5:      # n*min(p, 1-p) small; use inverse transformation.
         initial_pr = (1.0 - d3) ** n
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)            
            out_len = Ns
            j = 0

         else:
            out_len = 0
            j = -1

         while j < out_len:
            pk = sum = initial_pr
            successes = 0
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            u = buf[i]
            while sum < u:
               successes = successes + 1
               pk = d8 * pk * (np1 - successes) /  successes
               sum = sum + pk

            if pr_success > 0.5: successes = n - successes
            if size is not None:
               buffer[j] = successes
               j = j + 1

            else:
               self._index = i    # Restore updated value of _ranbuf index.
               return successes

      else:                                  # n * min(p, 1-p) moderate or
                                             # large; use acceptance/rejection
         d6 = sqrt(n * d3 + 0.5)             # method augmented by ratio-of-
                                             # uniforms fast tests.
         d7 = (d6 * 1.715527770              # 2*sqrt(2/math.e) to 9 places
                  + 0.909016162)             # 3-2*sqrt(3/math.e) to 9 places
         d9 = log(d8)
         d10 = floor(np1 * d3)
         d11 = min( np1, floor(d4 + 7 * d6) )     # 7 is a precision constant.
         lf = self._ln_factorial                  # local alias
         d12 = lf(d10) + lf(n - d10)
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns
         else:
            out_len = 1
         for j in xrange(out_len):
            while 1:
               i = i + 1                          # Generate two pseudo-random
               if i == buflim: self._fillbuf(); i = 0          # uniforms.
               x = buf[i]
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               y = buf[i]
               w = d4 + d7 * (y - 0.5) / x        # ratio of uniforms
               if w < 0.0 or w > d11:
                  continue                        # fast rejection; try
                                                  # another x and y
               else:
                  successes = floor(w)
                  t = (successes-d10)*d9+d12 - lf(successes) - lf(n-successes)
                  if x * (4.0 - x) - 3.0 < t:
                     break                        # fast acceptance

                  elif x * (x - t) >= 1.0:
                     continue                     # fast rejection, try
                                                  # another x and y
                  elif t >= 2 * log(x):
                     break                        # acceptance

            if size is not None:
               buffer[j] = int(successes)

            else:
               self._index = i       # Restore updated value of _ranbuf index.
               return int(successes)

      self._index = i             # Restore updated _ranbuf index and return.
      return Num.reshape(buffer, size)

   def _Cauchy(self, median=0.0, scale=1.0, size=None):
      """Return Cauchy random variates; parameters <median> and <scale>.

      Cauchy(median=0.0, scale=1.0, size=None)

      <scale> must be positive.  If <buffer> is specified, it must be a
      mutable sequence (list or array).  It is filled with
      Cauchy(median, scale) pseudo_random variates and <None> is returned.
      Otherwise, a single Cauchy pseudo-random variate is returned."""
      if scale <= 0.0:
         raise ValueError, '<scale> must be positive'

      else:
         i = self._index                             # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            buflim = self._ranbuf_size               # local alias
            buf = self._ranbuf                       # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               x = buf[i]
               buffer[j] = median + scale / tan(3.1415926535897931 * x)

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            x = self._ranbuf[i]
            self._index = i          # Restore updated value of _ranbuf index.
            return median + scale / tan(3.1415926535897931 * x)

      self._index = i  # Restore _ranbuf index and return; buffer[] is filled.
      return Num.reshape(buffer, size) 

   def _chi_square(self, df=10.0, size=None):
      """Return chi-square pseudo-random variates; <df> degrees of freedom.

      chi_square(df=10.0, size=None)

      <df> must be positive.  If <buffer> is specified, it must be a mutable
      sequence (list or array).  It is filled with chi-square(df) pseudo-
      random variates and <None> is returned.  Otherwise, a single chi-square
      pseudo-random variate is returned."""
      if df <= 0.0:
         raise ValueError, '<df> (degrees of freedom) must be positive'

      else:
         gamma = self._gamma                   # local alias
         return 2.0*gamma(0.5*df,size=size)

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

   def _exponential(self, scale=1.0, size=None):
      """Return exponential pseudo-random variates; scale-factor <scale>.

      exponential(scale=1.0, size=None)

      <scale> is positive;  the density function is (1/scale)*exp(-x/scale).
      If <buffer is specified, it must be a mutable sequence (list or array).
      It is filled with exponential RVs and <None> is returned.  Otherwise,
      a single exponential pseudo-random variate is returned."""
      scal = abs(scale)             # Fix scale if the user guessed wrong.
      if scal == 0.0:
         raise ValueError, '<scale> must be positive'

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
               buffer[j] = -scal * log(x)

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            x = self._ranbuf[i]
            self._index = i          # Restore updated value of _ranbuf index.
            return -scal * log(x)

      self._index = i  # Restore _ranbuf index and return (buffer is filled).
      return Num.reshape(buffer, size)
      

   def _Fisher_F(self, numdf=2.0, denomdf=10.0, size=None):
      """Return F(numdf, denomdf) distribution pseudo-random variates.

      Fisher_F(numdf=2.0, denomdf=10.0, size=None)

      Both <numdf> and <denomdf> must be positive.  if <buffer> is specified,
      it must be a mutable sequence (list or array).  <buffer> is filled with
      F(numdf, denomdf) pseudo-randoms and <None> is returned.  Otherwise,
      a single F pseudo-random variate is returned."""
      if (numdf <= 0.0) or (denomdf <= 0.0):
         raise ValueError, 'both numdf and denomdf must be positive'

      else:
         beta = self._beta           # local alias
         x = beta(0.5 * numdf, 0.5 * denomdf, size=size)
         return denomdf * x / (numdf * (1.0-x))

   def _gamma(self, mu=2.0, size=None):
      """Return pseudo-random variates from gamma distribution; mean <mu>.

       gamma(mu=2.0, size=None)

       The gamma distribution is also known as the incomplete gamma function.
       <mu> must be positive.  If <buffer> is specified, it must be a mutable
       sequence (list or array).  <buffer> is filled with gamma(mu) RVs and
       <None> is returned.  Otherwise, a single gamma(mu) variate is returned.
       """
      if mu <= 0.0:
         raise ValueError, 'mu must be positive'

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

         if mu <= 1.0:
            t = 0.07 + 0.75 * sqrt(1.0 - mu) # Use Best's (1983)  RGS alg.:
            b = 1.0 + exp(-t) * mu / t       # Devroye, Luc, "Non-Uniform
            c = 1.0 / mu                     # Random Variate Generation,"
            while j < out_len:               # Sprinter-Verlag, 1985, p. 426.
               while 1:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u = buf[i]
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  w = buf[i]
                  v = b * u
                  if v <= 1.0:
                     x = t * v ** c
                     if w <= (2.0 - x) / (2.0 + x) or w <= exp(-x): break

                  else:
                     x = -log(c * t * (b - v))
                     y = x / t
                     if w * (mu+y-mu*y) <= 1.0 or w <= y**(mu-1.0): break

               if out_len:
                  buffer[j] = x
                  j = j + 1

               else:
                  self._index = i    # Restore updated value of _ranbuf index.
                  return x

         else:                  # Use Best's rejection algorithm ( Best(1978);
            b = mu - 1.0        # Devoyre, ibid. 1985, p. 410).
            c = 3.0 * mu - 0.75
            while j < out_len:
               while 1:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u = buf[i]
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  v = buf[i]
                  w = u * (1.0 - u)
                  y = sqrt(c/w) * (u - 0.5)
                  x = b + y
                  if x > 0.0:
                     z = 64.0 * w * w * w * v * v
                     if ( (z <= (1.0 - 2.0 * y * y / x)) or
                        (log(z) <= 2.0 * (b * log(x/b) - y)) ): break

               if out_len:
                  buffer[j] = x
                  j = j + 1

               else:
                  self._index = i    # Restore updated value of _ranbuf index.
                  return x

         self._index = i # Restore _ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)        


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

   def _Gumbel(self, mode=0.0, scale=1.0, size=None):
      """Return Gumbel (extreme value) distribution pseudo-random variates.

      Gumbel(mode=0.0, scale=1.0, size=None)

      <scale> must be positive.  If <buffer> is specified, it must be a
      mutable sequence.  <buffer> is filled with Gumbel(mode, scale) RVs and
      <None>is returned.  Otherwise, a single Gumbel RV is returned."""
      if scale <= 0.0:
         raise ValueError, '<scale> must be positive'

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
               u = buf[i]
               buffer[j] =  mode - scale * log( -log(u) )

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            u = self._ranbuf[i]
            self._index = i          # Restore updated value of _ranbuf index.
            return mode - scale * log( -log(u) )

      self._index = i # Restore _ranbuf index and return (buffer[] is filled).
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
         

   def _Laplace(self, loc=0.0, scale=1.0, size=None):
      """Return Laplace (double exponential) pseudo-random variates.

      laplace(loc=0.0, scale=1.0, size=None)

      The mean is <mu> and the (positive) scale factor is <scale>.  If
      <buffer> is specified, it must be a mutable sequence (list or array).
      <buffer> is filled with Laplace(mu, scale) pseudo-random variates and
      <None> is returned.  Otherwise, a single Laplace pseudo-random variate
      is returned."""
      mu = loc
      if scale <= 0.0:
         raise ValueError, 'Scale must be positive.'

      else:
         i = self._index                          # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)            
            buf = self._ranbuf                    # local alias
            buflim = self._ranbuf_size            # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               x = buf[i]
               if x < 0.5:
                  x =  mu + scale * log(x + x)

               else:
                  x =  mu - scale * log(2.0 - x - x)

               buffer[j] = x

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            x = self._ranbuf[i]
            if x < 0.5: x = mu + scale * log(x + x)
            else:       x = mu - scale * log(2.0 - x - x)
            self._index = i          # Restore updated value of _ranbuf index.
            return x

      self._index = i # Restore _ranbuf index and return (buffer[] is filled).
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

   def _logistic(self, mu=0.0, scale=1.0, size=None):
      """Return logistic pseudo-random variates; mean <mu>.

      logistic(mu=0.0, scale=1.0, size=None)

      <scale> must be positive.  If <buffer> is specified, it must be a
      mutable sequence (list or array).  <buffer> is filled with
      logistic(mu, scale) pseudo-random variates and <None> is returned.
      Otherwise, a single logistic pseudo-random variate is returned."""
      if scale <= 0.0:
         raise ValueError, 'scale must be positive'

      else:
         i = self._index                          # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            buf = self._ranbuf                    # local alias
            buflim = self._ranbuf_size            # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               x = buf[i]
               buffer[j] =  mu + scale * log( x / (1.0 - x) )

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            x = self._ranbuf[i]
            self._index = i          # Restore updated value of _ranbuf index.
            return mu + scale * log( x / (1.0 - x) )

      self._index = i # Restore _ranbuf index and return (buffer[] is filled).
      return Num.reshape(buffer, size)
      


   def _lognormal(self, mean=0.0, sigma=1.0, size=None):
      """Return lognormal distribution pseudo-random variates.

      lognormal(mean=0.0, sigma=1.0, size=None)

      The parameters <mean> and <sigma> are the mean and standard deviation
      of the underlying normal distribution.  <sigma> must be positive. If
      <buffer> is specified, it must be a mutable sequence (list or array).
      <buffer> is filled with lognormal(mean, sigma) pseudo-random variates
      and <None> is returned.  Otherwise, a single lognormal variate is
      returned."""
      if sigma <= 0.0:
         raise ValueError, '<sigma> must be positive'

      else:
         normal = self._normal
         return Num.exp( mean + sigma * normal(size=size))

   def _negative_binomial(self, r=1.0, pr=0.5, size=None):
      """Return negative binomial random variates (non-negative integers).

      negative_binomial(r=1.0, pr=0.5, size=None)

      For <r> integral, the negative binomial distribution is also called the
      Pascal distribution.  If <r> = 1, it is the geometric distribution.
      <r> must be positive, and 0 <= pr < 1.  If <buffer> is speci-
      fied, it must be a mutable sequence (list or array).  <buffer> is filled
      with negative binomial(r, pr) pseudo-random variates, and <None>
      is returned.  Otherwise, a single negative binomial pseudo-random
      variate is returned."""
      if ( not (0.0 < pr < 1.0) ) or r <= 0.0:
         raise ValueError, '<r> must be positive and 0.0 < pr < 1.0'

      else:
         pr_success = 1.0 - pr
         ratio = pr_success / pr
         Poisson = self._Poisson
         gamma = self._gamma
         return Poisson(gamma(r) / ratio, size=size)

   def _normal(self, mu=0.0, sigma=1.0, size=None):
      """Return normal pseudo-random variates; mean <mu> & std. dev. <sigma>.

      normal(mu=0.0, sigma=1.0, size=None)

      Two random normal variates are generated on odd-numbered calls;
      the second is returned on even-numbered calls for single normals.
      <sigma> must be positive.  If <buffer> is specified, it must be a
      mutable sequence (list or array).  <buffer> is filled with
      normal(mu, sigma) pseudo-random variates and <None> is returned.
      Otherwise, a single normal pseudo-random variate is returned."""

      if sigma <= 0.0:
         raise ValueError, '<sigma> must be positive'

      else:
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns
            j = 0

         else:
            out_len = 0
            j = - 1

         i = self._index                          # local alias
         buf = self._ranbuf                       # local alias
         buflim = self._ranbuf_size               # local alias
         second = self._second_nrv                # local alias
         while j < out_len:
            x = second
            if x:                   # _second_nrv set to None in _initialize()
               second = None

            else:
               z = 2.0
               while z >= 1.0:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  x = buf[i]
                  x = x + x - 1.0      # -1 < x < 1
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  y = buf[i]
                  y = y + y - 1.0      # -1 < y < 1
                  z = x * x + y * y    # Insure (x,y) is in the unit circle.

               e = -log(z)
               e = sqrt((e + e) / z)
               second = y * e
               x = x * e

            if out_len:
               buffer[j] = x * sigma + mu
               j = j + 1

            else:
               self._index = i               # Restore updated  _ranbuf index.
               self._second_nrv = second     # Restore normal buffer.
               return x * sigma + mu

         self._second_nrv = second           # Restore normal buffer.
         self._index = i # Restore _ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)         

   def _Pareto(self, mode=1.0, shape=4.0, size=None):
      """Return Pareto distribution pseudo-random variates.

      Pareto(mode=1.0, shape=4.0, size=None)

      Both <mode> and <shape> must be positive.  <buffer>, if specified, must
      be a mutable sequence (list or array).  It is filled with
      Pareto(mode, shape) pseudo-random variates and <None> is returned.
      Otherwise, a single Pareto pseudo-random variate is returned."""
      if (mode <= 0.0) or (shape <= 0.0):
         raise ValueError, 'both <mode> and <shape> must be positive'

      else:
         i = self._index                          # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            buflim = self._ranbuf_size            # local alias
            buf = self._ranbuf                    # local alias
            exponent = -1.0 / shape
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               buffer[j] = mode * (1.0 - buf[i]) ** exponent

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            self._index = i          # Restore updated value of _ranbuf index.
            return mode * (1.0 - self._ranbuf[i]) ** (-1.0 / shape)

      self._index = i # Restore _ranbuf index and return (buffer[] is filled).
      return Num.reshape(buffer, size)
      


   def _Poisson(self, rate=1.0, size=None):
      """Return integer pseudo-random variates from a Poisson distribution.

      Poisson(rate=1.0, size=None)

      See Fishman, George S., "Monte Carlo Concepts, Algorithms and
      Applications," Springer-Verlag, New York, 1996.  Algorithms ITR
      (pp. 154-155) and PRUA* (pp. 211-215) are used.  Mean computing time
      is bounded.  <rate> must be positive.  If <buffer> is specified, it
      must be a mutable sequence.  It is filled with Poisson(rate)
      pseudo-random variates (integers) and <None> is returned.  Otherwise,
      a single Poisson pseudo-random variate (an integer) is returned."""

      if rate <= 0.0:
         raise ValueError, '<rate> must be positive'

      else:
         i = self._index                          # local alias
         buflim = self._ranbuf_size               # local alias
         buf = self._ranbuf                       # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            out_len = Ns            
            j = 0

         else:
            out_len = 0
            j = -1

         init = exp(-rate)
         if rate < 10.0:     # small or moderate rate:  Use inverse transform.
            while j < out_len:
               pk = sum = init
               successes = 0
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               u = buf[i]
               while sum < u:
                  successes = successes + 1
                  pk = pk * rate / successes
                  sum = sum + pk

               if out_len:
                  buffer[j] = successes
                  j = j + 1

               else:
                  self._index = i    # Restore updated value of _ranbuf index.
                  return successes

         else:                    # large rate: Use acceptance/rejection
            while j < out_len:    # augmented by the ratio of uniforms method.
               d3 = rate + 0.5
               d4 = sqrt(d3)
               d5 = ( 1.715527770 * d4  # 2 * sqrt(2 / math.e) to 9 digits
                  + 0.898916162 )       # 3 - 2 * sqrt(3 / math.e) to 9 digits
               d6 = log(rate)
               d7 = floor(rate)
               d8 = floor( d3 + 7.0 * (d4 + 1.0) )# 7 is a precision constant.
               lf = self._ln_factorial            # local alias
               d9 = d6 * d7 - lf(d7)
               while 1:
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  x = buf[i]
                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  y = buf[i]
                  z = d3 + d5 * (y - 0.5) / x     # ratio of uniforms
                  if z < 0.0 or z >= d8:
                     continue                     # fast rejection: new x,y

                  else:
                     successes = floor(z)
                     t = successes * d6 - lf(successes) - d9
                     if t >= x * (4.0 - x) - 3.0:
                        break                     # fast acceptance; return k

                     elif x * (x - t) >= 1.0:
                        continue                  # fast rejection; new x, y

                     elif t >= 2.0 * log(x):
                        break                     # acceptance; ret. successes

               if out_len:
                  buffer[j] = int(successes)
                  j = j + 1

               else:
                  self._index = i    # Restore updated value of _ranbuf index.
                  return int(successes)

         self._index = i # Restore _ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)         

   def _randint(self, lowint=0, upint=1, size=None):
      """Return integers k with equal probability; lowint <= k <= upint.

      randint(lowint=0, upint=1, size=None)

      <lowint> must be <= <upint>.  If <buffer> is specified, it must be a
      mutable sequence (list or array).  It is filled with pseudo-random
      integers from the closed interval [lowint, upint].
      Note that the default is a Bernouli (0,1) random variate.  If
      <lowint> == <upint>, the common value is returned."""
      if not (lowint <= upint):
         raise ValueError, '<lowint> must not exceed <upint>'

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
               buffer[j] = lowint + int(floor(buf[i] * (upint - lowint + 1)))

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            self._index = i          # Restore updated value of _ranbuf index.
            return lowint + int(floor(self._ranbuf[i] * (upint - lowint + 1)))

         self._index = i # Restore _ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)         


   def _Rayleigh(self, scale=1.0, size=None):
      """Return Rayleigh distribution pseudo-random variates.

      rayleigh(scale=1.0, size=None)

      <mode> must be positive.  if shape is not None is specified, it must be a mutable
      sequence (list or array).  It is filled with Rayleigh(mode) pseudo-
      random variates and <None> is returned.  Otherwise, a single Rayleigh
      pseudo-random variate is returned."""
      mode = scale
      if mode <= 0.0:
         raise ValueError, '<mode> must be positive'

      else:
         i = self._index
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            buf = self._ranbuf                 # local alias
            buflim = self._ranbuf_size         # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               buffer[j] = mode * sqrt( -2.0 * log(buf[i]) )

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            self._index = i          # Restore updated value of _ranbuf index.
            return mode * sqrt( -2.0 * log(self._ranbuf[i]) )

         self._index = i # Restore _ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)
         

   def _Student_t(self, df=100.0, size=None):
      """Return Student's t pseudo-random variate; <df> degrees of freedom.

      Student_t(df=100.0, size=None)

      <df> must be positive.  If <buffer> is specified, it must be a mutable
      sequence (list or array).  It is filled with Student-t(df) pseudo-random
      variates and <None> is returned.  Otherwise, a single Student-t pseudo-
      random variate is returned.  For df > 3.0, see Fishman, "Monte Carlo,"
      pp 207-209, algorithm T3T*."""
      if df <= 0.0:
         raise ValueError, '<df> must be positive'

      else:
         if df <= 3.0:
            cdf = 0.5 * df
            con = 2.0 / df
            if size is not None:
               size, Ns = _check_shape(size)
               buffer = Num.zeros(Ns,Num.Float)
               normal = self._normal
               gamma = self._gamma
               for j in xrange( Ns ):
                  buffer[j] = normal()/sqrt( con * gamma(cdf) )
               return Num.reshape(buffer, size)
            else:
               return self._normal() / sqrt( con * self._gamma(cdf) )
            
         else:
            if size is not None:
               size, Ns = _check_shape(size)
               buffer = Num.zeros(Ns,Num.Float)
               out_len = Ns
               j = 0

            else:
               out_len = 0
               j = -1

            i = self._index              # local alias
            buflim = self._ranbuf_size   # local alias
            buf = self._ranbuf           # local alias
            while j < out_len:
               while 1:
                  while 1:
                     i = i + 1
                     if i == buflim: self._fillbuf(); i = 0
                     x = buf[i]
                     i = i + 1
                     if i == buflim: self._fillbuf(); i = 0
                     y = buf[i]
                     y = y + y - 1.0     # uniform(-1,+1)
                     z = 0.866025404 * y / x
                     t = z * z
                     if x * (3.0 + t) <= 3.0:
                        break

                  i = i + 1
                  if i == buflim: self._fillbuf(); i = 0
                  u = buf[i]
                  w = u * u
                  if 0.860081359786 - t / 8.59111293933 >= w:
                     break               # 0.86008... = (9*sqrt(math.e)/16)**2
                                         # 8.591... = 256*math.e/81
                  else:
                     v = 1.0 + 0.33333333333333333 * t
                     s = 2.0 * log(v * v * 0.927405714769 / u)  # 9*sqrt(e)/16
                     if s >= t:
                        break

                     else:
                        dfpone = df + 1.0
                        if s >= 1.0 + dfpone * log((t + df)/ dfpone):
                           break

               if out_len:
                  buffer[j] = z
                  j = j + 1

               else:
                  self._index = i
                  return z

            self._index = i # Restore ranbuf index and return (buffer filled).
            return Num.reshape(buffer, size)
            

   def _triangular(self, left=0.0, mode=0.5, right=1.0, size=None):
      """Return triangular pseudo-random variates; left <= mode <= right.

      triangular(left=0.0, mode=0.5, right=1.0, size=None)

      The triangle is ABC; A and B ( <left> and <right> ) are on the
      horizontal axis, and a perpendicular dropped from C intersects the
      horizontal axis at <mode>.  Neither angle with the horizontal axis is
      allowed to exceed 90 degrees.  Thus, left <= mode <= right.  If <buffer>
      is specified, it must be a mutable sequence (list or array).  It is
      filled with triangular pseudo-random varuiates and <None> is returned.
      Otherwise, a single triangular pseudo-random variate is returned."""
      if not(left <= mode <= right):
         raise ValueError, '<left> <= <mode> <= <right>'

      else:
         i = self._index                          # local alias
         base = right - left
         leftbase = mode - left
         ratio = leftbase / base
         leftprod = leftbase * base
         rightprod = (right - mode) * base
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)            
            buf = self._ranbuf                    # local alias
            buflim = self._ranbuf_size            # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               u = buf[i]
               if u <= ratio:
                  buffer[j] = left + sqrt(u * leftprod)

               else:
                  buffer[j] = right - sqrt( (1.0 - u) * rightprod )

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            u = self._ranbuf[i]
            self._index = i          # Restore updated value of _ranbuf index.
            if u <= ratio:
               return left + sqrt(u * leftprod)

            else:
               return right - sqrt( (1.0 - u) * rightprod )

         self._index = i  # Restore ranbuf index and return (buffer[] filled).
         return Num.reshape(buffer, size)
         
   def _uniform(self, lower=-0.5, upper=+0.5, size=None):
      """Return pseudo-randoms from uniform distribution on [lower, upper].

      uniform(lower=-0.5, upper=+0.5, size=None)

      If <buffer> is specified, it must be a mutable sequence.  It is filled
      with uniform(lower, upper) pseudo-randoms and <None> is returned.
      Otherwise, a single uniform pseudo-random variate is returned."""
      length = upper - lower
      if length < 0.0:
         raise ValueError, '<upper> must be at least as large as <lower>'

      else:
         i = self._index                          # local alias
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)
            buf = self._ranbuf                    # local alias
            buflim = self._ranbuf_size            # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               buffer[j] = lower + length * buf[i]
         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            self._index = i          # Restore updated value of _ranbuf index.
            return lower + length * self._ranbuf[i]

      self._index = i # Restore _ranbuf index and return (buffer[] is filled).
      return Num.reshape(buffer, size)      

   def _von_Mises(self, mean=0.0, shape=1.0, size=None):
      """Return von Mises distribution pseudo-random variates on [-pi, +pi].

      von_Mises(mean=0.0, shape=1.0, size=None)

      <mean> must be in the open interval (-math.pi, +math.pi).  If <buffer>
      is specified, it must be a mutable sequence (list or array).  It is
      filled with von Mises(mean, shape) pseudo-random variates and <None> is
      returned.  Otherwise, a single von Mises RV is returned. The method is
      an algorithm of Best and Fisher, 1979; see Fisher, N. I., "Statistical
      Analysis of Circular Data," Cambridge University Press, 1995, p. 49."""
      if not (-3.1415926535897931 < mean < +3.1415926535897931):
         raise ValueError, \
            '<mean> must be in the open interval (-math.pi, math.pi)'

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
      
   def _Weibull(self, scale=1.0, shape=0.5, size=None):
      """Return Weibull distribution pseudo-random variates.

      Weibull(scale=1.0, shape=0.5, size=None)

      Both <scale> and <shape> must be positive.  If <buffer> is specified, it
      must be a mutable sequence.  It is filled with Weibull(scale, shape)
      pseudo-random variates and <None> is returned.  Otherwise, a single
      Weibull variate is returned."""
      if ( (scale <= 0.0) or (shape <= 0.0) ):
         raise ValueError, 'both <scale> and <shape> must be positive'

      else:
         i = self._index                          # local alias
         exponent = 1.0 / shape
         if size is not None:
            size, Ns = _check_shape(size)
            buffer = Num.zeros(Ns,Num.Float)            
            buf = self._ranbuf                    # local alias
            buflim = self._ranbuf_size            # local alias
            for j in xrange( Ns ):
               i = i + 1
               if i == buflim: self._fillbuf(); i = 0
               buffer[j] = scale * pow(-log(buf[i]), exponent)

         else:
            i = i + 1
            if i == self._ranbuf_size: self._fillbuf(); i = 0
            self._index = i          # Restore udpated value of _ranbuf index.
            return scale * pow(-log(self._ranbuf[i]), exponent)

      self._index = i # Restore _ranbuf index and return (buffer[] is filled).
      return Num.reshape(buffer, size)
      
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
      surface of a sphere," Annals of Mathematical Statistics, vol. 43,
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


   def _permutation(self, elements=[0, 1, 2, 3, 4]):
      """Return a pseudo-random permutation of the items in elements[].

      permutation(elements=[0, 1, 2, 3, 4])

      <elements> must be a mutable sequence containing the items to be
      permuted.  A sequence of the same type as <elements[]> containing the
      randomly permuted items is returned.

         template = array.array('l', range(5))
         perm = rv.smart_sample(template, 5)
      #  is equivalent to:
         template = array.array('l', range(5))
         perm = rv.permutation(template)
      """
      lng = len(elements)
      outlist = elements[:]
      if lng == 0:
         pass

      else:
         i = self._index                          # local alias
         buf = self._ranbuf                       # local alias
         buflim = self._ranbuf_size               # local alias
         for j in xrange(lng):
            i = i + 1
            if i == buflim: self._fillbuf(); i = 0
            k = int( floor(buf[i] * (lng - j)) )
            temp = outlist[j]
            outlist[j] = outlist[k]
            outlist[k] = temp

         self._index = i              # Restore updated value of _buflim index
         return outlist               # and return.

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


#               End of class <_pranv> definition

#                ---------------GLOBALS---------------
#
#             N.B:  Globals are set here and nowhere else!

_inst = _pranv()   # Initialize uniform(0,1) generator to default; clock seed.

#                  Shortcuts to <_inst> methods

#     Uniform(0,1) Random Number Generator
random2            =  _inst._random

#            Utilities
initial_seed      =  _inst._initial_seed
initialize        =  _inst._initialize
random_algorithm  =  _inst._random_algorithm
random_count      =  _inst._random_count
save_state        =  _inst._save_state

#     Non-Uniform(0,1) Random Number Generators
#beta              =  _inst._beta
#binomial          =  _inst._binomial
#cauchy            =  _inst._Cauchy
#chi2              =  _inst._chi_square
choice            =  _inst._choice
#exponential       =  _inst._exponential
#fratio            =  _inst._Fisher_F
#gamma             =  _inst._gamma
geom         =  _inst._geom
#gumbel            =  _inst._Gumbel
hypergeom    =  _inst._hypergeom
#laplace           =  _inst._Laplace
logser       =  _inst._logser
#logistic      =  _inst._logistic
#lognormal         =  _inst._lognormal
#negative_binomial =  _inst._negative_binomial
#normal            =  _inst._normal
#pareto            =  _inst._Pareto
#poisson           =  _inst._Poisson
#randint           =  _inst._randint
#rayleigh          =  _inst._Rayleigh
#student_t         =  _inst._Student_t
triang            =  _inst._triangular
#uniform           =  _inst._uniform
von_mises         =  _inst._von_Mises
invnorm             =  _inst._Wald
def wald(loc=0.0, scale=1.0, size=None):
   return invnorm(1.0, loc, scale, size=size)
#weibull           =  _inst._Weibull
zipf              =  _inst._Zipf

#     Geometrical Point Generators, Permutations and Subsampling Routines
in_simplex        =  _inst._in_simplex
in_sphere         =  _inst._in_sphere
on_simplex        =  _inst._on_simplex
on_sphere         =  _inst._on_sphere
#permutation       =  _inst._permutation
#sample            =  _inst._sample
smart_sample      =  _inst._smart_sample

#                    End of Module rv Proper

#                    Module rv self-testing code
class _gen_timer:
   """Execution times for pseudo-random generators

   """

   def __init__(self):
      """Set sample size and determine slack (dead) time for timing loop.

      """
#      data attributes                definition
#  self._buffer_deadtime  execution time for buffered null loop, microseconds
#  self._buffer_size      length of buffer specified in function call
#  self._deadtime         execution time for unbuffered null loop, mcs
#  self._repetitions      Number of pseudorandoms generated
#
      self._repetitions = 1000   # Should be >= 100,000 for accurate timing.
      self._buffer_size = 100
      expression = 'None'
      nullcode = compile(expression, expression, 'eval')
      clock = time.clock                    # local alias
      x = 1.0
      sum = sumsq = smallest = largest = 0.0
      t0 = clock()
      for i in xrange(self._repetitions):
         eval(nullcode)
         sum = sum + x
         sumsq = sumsq + x * x
         smallest = min(x, smallest)
         largest = max(x, largest)

      t1 = clock()
      self._deadtime = (t1 - t0) * 1e6 / self._repetitions  # slack time, mcs.

      dummy = array.array('d', [0.0])
      groups = int(self._repetitions / self._buffer_size)
      sum = sumsq = smallest = largest = 0.0
      t0 = clock()
      for i in xrange(groups):
         eval(nullcode)
         for j in xrange(self._buffer_size):
            sum = sum + x
            sumsq = sumsq + x * x
            smallest = min(x, smallest)
            largest = max(x, largest)

      t1 = clock()
      self._buffer_deadtime = (t1 - t0) * 1e6 / self._repetitions  # buffer
                                                            # slack time, mcs.

   def _bench(self, generator, expected, var, buf=None):
      """Print time required and stats for n pseudo-randoms from <generator>.

      """
      gen_call_code = compile(generator, generator, 'eval')
      clock = time.clock                    # local alias
      n = self._repetitions                 # local read-only alias
      sum = sumsq = 0.0
      smallest = 1e300
      largest = -1e300
      if buf is not None:
         buffer_size = self._buffer_size
         groups = n / buffer_size
         t0 = clock()
         for i in xrange(groups):
            eval(gen_call_code)
            for j in xrange(buffer_size):
               x = buf[j]
               sum = sum + x
               sumsq = sumsq + x * x
               smallest = min(smallest, x)
               largest = max(largest, x)

         t1 = clock()
         t = (t1 - t0) * 1e6 / n - self._buffer_deadtime
         mean = sum / n
         stdev = sqrt( (sumsq - sum*sum/n) / (n - 1.0) )
         print "rv.%-50s   theoretical mean = %9.3g; SD = %9.3g" \
           %(generator,expected,sqrt(var))
         print "time: %3.0f mcs; min = %#11.5g;  max = %#14.10g;" \
            "  sample average = %#9.3g; SD = %9.3g" \
            %(t, smallest, largest, mean, stdev)
         limit = 5.0 * sqrt(var / n)
         if not(mean - limit <= expected <= mean + limit):
            print 'Above mean is out of range.  Possible algorithm error?'

         print


      else:
         t0 = clock()
         for i in xrange(n):
            x = eval(gen_call_code)
            sum = sum + x
            sumsq = sumsq + x * x
            smallest = min(smallest, x)
            largest = max(largest, x)

         t1 = clock()
         t = (t1 - t0) * 1.e6 / n - self._deadtime
         mean = sum / n
         stdev = sqrt( (sumsq - sum*sum/n) / (n - 1.0) )
         print "rv.%-50s   theoretical mean = %9.3g; SD = %9.3g" \
           %(generator,expected,sqrt(var))
         print "time: %3.0f mcs; min = %#11.5g;  max = %#14.10g;" \
            "  sample average = %#9.3g; SD = %9.3g" \
            %(t, smallest, largest, mean, stdev)
         limit = 5.0 * sqrt(var / n)
         if not(mean - limit <= expected <= mean + limit):
            print 'Above mean is out of range.  Possible algorithm error?'

         print


def _vergen(gen_abbr, seed, u100):
   """Verify uniform(0,1) pseudo-random <gen_abbr>, initial seed <seed>.

   <u100> is the correct value of the 100th pseudorandom.
   """
   initialize('', gen_abbr, seed)
   x = 105*[0.0]
   random2(x)
   temp = random_algorithm()
   print temp[0:string.find(temp, ':')],
   print ' 100th result is %12.10g' % (x[99],),
   if abs(round(x[99],8) - u100) > 1e-8:
      print "-- Algorithm error!, expected %12.10f" % u100
      for i in range(96, 104):
         print '%3d, %12.10f' %(i, x[i])

   else:
      print '-- verified.'

   assert initial_seed() == seed, "supplied seed stored incorrectly"
   del x

def _versave(gen_abbr):
   """Verify proper operation of _save_state() with generator <gen_abbr>.

   """
   print random_algorithm(),
   initialize(file_string = None, algorithm = gen_abbr)
   for i in xrange(1000): x1000 = random2()
   initialize(file_string = None, algorithm=gen_abbr, seed=initial_seed())
   for i in xrange(500): y1000 = random2()
   save_state()
   initialize(file_string = 'default')
   for i in xrange(500): y1000 = random2()
   assert abs(round(x1000,8) - round(y1000,8)) < 1e-7,\
     "Uniform generator state not saved correctly: " + `x1000` + ' ' + `y1000`
   print "save_state() verified."


def _rvtest():
   """Test and benchmark Module rv pseudorandom generators

   """
   # Test default uniform generator.
   print "The default pseudorandom generator is:"
   print '   ', random_algorithm()
   print "The first uniform(0,1) pseudorandom delivered is: %12.6g" \
          %(random2(),)
   print "The (epoch-time-generated) seed was: %19s" %(initial_seed(),)
   print "%3d uniform(s) generated from this seed." %(random_count(),)
   print

   _vergen('CMRG', -12345L,  0.75923860)      # Verify generator 'CMRG'.
   _vergen('Flip', -314159L, 0.50010253)      # Verify generator 'Flip'.
   _vergen('SMRG', -12345L,  0.62457420)      # Verify generator 'SMRG'.
   _vergen('Twister', 4357L, 0.24882571)      # Verify generator 'Twister'.
   print

   _versave('CMRG')       # Verify save_state('CMRG').
   _versave('Flip')       # Verify save_state('Flip').
   _versave('SMRG')       # Verify save_state('SMRG').
   _versave('Twister')    # Verify save_state('Twister').
   print

   pt = 2*[0.0]
   qt = 2*[0.0]
   print "pt = ", pt
   for j in xrange(700):
      generator='qt = rv.in_simplex(seq=pt,bound=1.0)'; pt=in_simplex(qt,1.0)
      assert (1.0 - pt[0]) >= pt[1], 'point outside triangle'
      assert (pt[0] >= 0.0) and (pt[1] >= 0.0), 'negative coordinate'
      if j < 2: print '%35s, qt is %6.4f, %6.4f' %(generator, pt[0], pt[1])
      generator='qt = rv.in_sphere(center=pt,radius=1.0)';pt=in_sphere(qt,1.0)
      assert (pt[0]**2 + pt[1]**2) <= 1.0, 'point outside circle'
      if j < 2: print '%35s, qt is %6.4f, %6.4f' %(generator,pt[0], pt[1])
      generator='qt = rv.on_simplex(seq=pt,bound=1.0)'; pt=on_simplex(qt, 1.0)
      assert (pt[0] + pt[1] - 1.0)<1.e-12, 'point not on triangle hypotenuse'
      if j < 2: print '%35s, qt is %6.4f, %6.4f' %(generator, pt[0], pt[1])
      generator='qt = rv.on_sphere(center=pt,radius=1.0)';pt=on_sphere(qt,1.0)
      assert (pt[0]**2+pt[1]**2-1.0) < 1.e-12, 'point not on circle boundary'
      if j < 2: print '%35s, qt is %6.4f, %6.4f' %(generator, pt[0], pt[1])
      if j < 2: print

   sum = pt[:]
   sum[0] = sum[1] = 0.0
   for j in xrange(1000):
      pt=on_sphere(qt, 1.0)
      for k in range(2): sum[k] = sum[k] + pt[k]

   for k in range(2): sum[k] = sum[k] / 1000.0
   print 'on_sphere() Average of 1000 points: %6.4f, %6.4f' %(sum[0], sum[1])
   del sum
   print

   pt = array.array('d', 3*[0.0])
   qt = array.array('d', 3*[0.0])
   print 'pt = ', pt
   for j in xrange(700):
      generator='qt = rv.in_simplex(seq=pt,bound=1.0)';pt=in_simplex(qt,1.0)
      assert (pt[0]+pt[1]+pt[2]) <= 1.0, 'point outside simplex'
      assert (pt[0]>=0.0)and(pt[1]>=0.0)and(pt[2]>=0.0), 'negative coordinate'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2])
      generator='at = rv.in_sphere(center=pt,radius=1.0)';pt=in_sphere(qt,1.0)
      assert (pt[0]**2+pt[1]**2+pt[2]**2) <= 1.0, 'point outside sphere'
      if j < 2: print '%35s, pt is %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2])
      generator='qt = rv.on_simplex(seq=pt,bound=1.0)';pt=on_simplex(qt, 1.0)
      assert (pt[0]+pt[1]+pt[2]-1.0)<1.e-12, 'point not on simplex boundary'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2])
      generator='qt = rv.on_sphere(center=pt,radius=1.0)';pt=on_sphere(qt,1.0)
      assert (pt[0]**2+pt[1]**2+pt[2]**2-1.0)<1.e-12, 'point not on surface'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2])
      if j < 2: print

   sum = pt[:]
   sum[0] = sum[1] = sum[2] = 0.0
   for j in xrange(1000):
      pt=on_sphere(qt, 1.0)
      for k in range(3): sum[k] = sum[k] + pt[k]

   for k in range(3): sum[k] = sum[k] / 1000.0
   print 'on_sphere() Average of 1000 points: %6.4f, %6.4f, %6.4f' \
          %(sum[0], sum[1], sum[2])
   del sum
   print

   pt = array.array('d', 4*[0.0])
   qt = array.array('d', 4*[0.0])
   print 'pt = ', pt
   for j in xrange(700):
      generator='qt = rv.in_simplex(seq=pt,bound=1.0)';pt=in_simplex(qt, 1.0)
      assert (pt[0]+pt[1]+pt[2]+pt[3]) <= 1.0, 'point outside simplex'
      assert (pt[0]>=0.0)and(pt[1]>=0.0)and(pt[2]>=0.0)and(pt[3]>=0.0), \
             'negative coordinate'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3])
      generator='qt = rv.in_sphere(center=pt,radius=1.0)';pt=in_sphere(qt,1.0)
      assert (pt[0]**2+pt[1]**2+pt[2]**2+pt[3]**2) <= 1.0, \
             'point outside sphere'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3])
      generator='qt = rv.on_simplex(seq=pt,bound=1.0)';pt=on_simplex(qt, 1.0)
      assert (pt[0]+pt[1]+pt[2]+pt[3]-1.0)<1.e-12, \
             'point not on simplex boundary'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3])
      generator='qt = rv.on_sphere(center=pt,radius=1.0)';pt=on_sphere(qt,1.0)
      assert (pt[0]**2+pt[1]**2+pt[2]**2+pt[3]**2-1.0)<1.e-12, \
             'point not on surface'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3])
      if j < 2: print

   sum = pt[:]
   sum[0] = sum[1] = sum[2] = sum[3] = 0.0
   for j in xrange(1000):
      pt=on_sphere(qt, 1.0)
      for k in range(4): sum[k] = sum[k] + pt[k]

   for k in range(4): sum[k] = sum[k] / 1000.0
   print 'on_sphere() Average of 1000 points: %6.4f, %6.4f, %6.4f, %6.4f' \
           %(sum[0], sum[1], sum[2], sum[3])
   del sum
   print

   pt = array.array('d', 5*[0.0])
   qt = array.array('d', 5*[0.0])
   print 'pt = ', pt
   for j in xrange(700):
      generator='qt = rv.in_simplex(seq=pt,bound=1.0)';pt=in_simplex(qt, 1.0)
      assert (pt[0]+pt[1]+pt[2]+pt[3]+pt[4]) <= 1.0, 'point outside simplex'
      assert (pt[0]>=0.0)and(pt[1]>=0.0)and(pt[2]>=0.0)and(pt[3]>=0.0), \
             'negative coordinate'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3], pt[4])
      generator='qt = rv.in_sphere(center=pt,radius=1.0)';pt=in_sphere(qt,1.0)
      assert (pt[0]**2+pt[1]**2+pt[2]**2+pt[3]**2+pt[4]**2) <= 1.0, \
             'point outside sphere'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3], pt[4])
      generator='qt = rv.on_simplex(seq=pt,bound=1.0)';pt=on_simplex(qt, 1.0)
      assert (pt[0]+pt[1]+pt[2]+pt[3]+pt[4]-1.0)<1.e-12, \
             'point not on simplex boundary'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3], pt[4])
      generator='qt = rv.on_sphere(center=pt,radius=1.0)';pt=on_sphere(qt,1.0)
      assert (pt[0]**2+pt[1]**2+pt[2]**2+pt[3]**2+pt[4]**2-1.0)<1.e-12, \
             'point not on surface'
      if j < 2: print '%35s, qt is %6.4f, %6.4f, %6.4f, %6.4f, %6.4f' \
                      %(generator, pt[0], pt[1], pt[2], pt[3], pt[4])
      if j < 2: print

   sum = pt[:]
   sum[0] = sum[1] = sum[2] = sum[3] = sum[4] = 0.0
   for j in xrange(1000):
      pt=on_sphere(qt, 1.0)
      for k in range(5): sum[k] = sum[k] + pt[k]

   for k in range(5): sum[k] = sum[k] / 1000.0
   print 'on_sphere() Average of 1000 points: %6.4f,%6.4f,%6.4f,%6.4f,%6.4f' \
           %(sum[0], sum[1], sum[2], sum[3], sum[4])
   del sum
   print

   elems = ['this', 'might', 'be', 'a', 'garbled', 'phrase']
   print "elems = ", elems
   generator = 'per = rv.permutation(elements=elems)';  per=permutation(elems)
   print generator, '; per is ', per
   generator = 'p = rv.sample([4,1,15,31,19,10,3,2,8,22], 8)'
   p = sample([4,1,15,31,19,10,3,2,8,22], 8); p.sort()
   print generator, '; p.sort(); p is ', p
   generator = 'p = rv.smart_sample([4,1,15,31,19,10,3,2,8,22], 8)'
   p = smart_sample([4,1,15,31,19,10,3,2,8,22], 8); p.sort()
   print generator, '; p.sort(); p is ', p
   print

   digits = array.array( 'l', range(10) )
   print 'digits = ', digits
   print
   print 'Simple random samples of size 5 from the decimal digits:'
   for i in range(5): print 'p = rv.sample(digits, 5); p is ',sample(digits,5)
   print
   print 'Simple random samples of size 15 from the decimal digits:'
   for i in range(5): print 'p = rv.sample(digits, 15); p is ', \
      sample(digits,15)
   print
   print \
       'Random samples of size 5 from the decimal digits without replacement:'
   for i in range(5):
      print 'p = rv.smart_sample(digits, 5); p is ', smart_sample(digits, 5)
   print
   del digits

   _do = _gen_timer()            # Instance of _gen_timer class
   benchmark = _do._bench        # alias for _bench()
   sz = _do._buffer_size
   buf = array.array('d', sz*[0.0])
   print "        Benchmark Timing and Statistics for Module rv Pseudo-random"
   print "                 Generators, Sample Size: % 10d" \
         % (int(_do._repetitions),)
   print
   if _do._repetitions < 100000:
      print "WARNING:  The sample size should be at least 100,000 for"
      print "two-decimal-digit accuracy execution time estimates."
      print "Change the variable <self._repetitions> near line number 2810."
   print
   print "buf = array.array('d',%5d*[0.0])" %(_do._buffer_size,)
   print
   initialize(file_string=None, algorithm='cmrg')
   print "rv.initialize(file_string=None, algorithm='cmrg')"
   benchmark('random2()', 0.5, 0.083333)
   benchmark('random2(size=sz)', 0.5, 0.083333, buf)

   initialize(file_string=None, algorithm='flip')
   print "rv.initialize(file_string=None, algorithm='flip')"
   benchmark('random2()', 0.5, 0.083333)
   benchmark('random2(size=sz)', 0.5, 0.083333, buf)

   initialize(file_string=None, algorithm='smrg')
   print "rv.initialize(file_string=None, algorithm='smrg')"
   benchmark('random2()', 0.5, 0.083333)
   benchmark('random2(size=sz)', 0.5, 0.083333, buf)

   initialize(file_string=None, algorithm='twister')
   print "rv.initialize(file_string=None, algorithm='twister')"
   benchmark('random2()', 0.5, 0.083333)
   benchmark('random2(size=sz)', 0.5, 0.083333, buf)

   initialize(file_string=None, algorithm='flip')  # Use fastest uniform(0,1)
                                                   # generator.
   benchmark('beta(mu=4.0, nu=2.0)', 0.66667, 0.031746)
   benchmark('beta(mu=4.0, nu=2.0, size=sz)', 0.66667, 0.031746, buf)
   benchmark('beta(mu=0.4, nu=0.2)', 0.66667, 0.138889)
   benchmark('beta(mu=0.4, nu=0.2, size=sz)', 0.66667, 0.138889, buf)
   benchmark('beta(mu=0.4, nu=2.0)', 0.16667, 0.040850)
   benchmark('beta(mu=0.4, nu=2.0, size=sz)', 0.16667, 0.040850, buf)

   benchmark('binomial(trials=25, pr_success=0.5)', 12.5, 6.25)
   benchmark('binomial(trials=25, pr_success=0.5, size=sz)',12.5,6.25,buf)
   benchmark('binomial(trials=25, pr_success=0.2)', 5.0, 4.0)
   benchmark('binomial(trials=25, pr_success=0.2, size=sz)',5.0,4.00,buf)
   benchmark('binomial(trials=500, pr_success=0.5)', 250.0, 125.0)

   benchmark('Cauchy(median=0.0, scale=1.0)', 0.0, 1e100)
   benchmark('Cauchy(median=0.0, scale=1.0, size=sz)', 0.0, 1e100, buf)

   benchmark('chi_square(df=1.0)', 1.0, 2.0)
   benchmark('chi_square(df=1.0,  size=sz)',   1.0,  2.0, buf)
   benchmark('chi_square(df=10.0)', 10.0, 20.0)
   benchmark('chi_square(df=10.0, size=sz)', 10.0, 20.0, buf)

   benchmark('choice( seq=(0,1) )', 0.5, 0.25)
   benchmark('choice(seq=(0,1), size=sz)', 0.5, 0.25,buf)
   benchmark('choice( seq=(2, 4, 6, 8) )', 5.0, 5.0)
   benchmark('choice(seq=(2, 4, 6, 8), size=sz)', 5.0, 5.0, buf)
   benchmark('choice(seq=[1, 2, 3])', 2.0, 0.66667)
   benchmark('choice(seq=[1, 2, 3], size=sz)', 2.0, 0.66667, buf)

   benchmark('exponential(scale=1.0)', 1.0, 1.0)
   benchmark('exponential(scale=1.0, size=sz)', 1.0, 1.0, buf)

   benchmark('Fisher_F(numdf=2.0, denomdf=10.0)', 1.25, 2.60417)
   benchmark('Fisher_F(numdf=2.0, denomdf=10.0, size=sz)',1.25,2.60417,buf)

   benchmark('gamma(mu=0.5)', 0.5, 0.5)
   benchmark('gamma(mu=0.5, size=sz)', 0.5, 0.5, buf)
   benchmark('gamma(mu=10.0)', 10.0, 10.0)
   benchmark('gamma(mu=10.0, size=sz)', 10.0, 10.0, buf)

   benchmark('geom(pr=0.5)', 1.0, 2.0)
   benchmark('geom(pr=0.5, size=sz)', 1.0, 2.0, buf)
   benchmark('geom(pr=0.95)', 19.0, 380.0)
   benchmark('geom(pr=0.95, size=sz)', 19.0, 380.0, buf)

   benchmark('Gumbel(mode=1.0,scale=1.0)', 1.57721, 1.644934)
   benchmark('Gumbel(mode=1.0,scale=1.0,size=sz)', 1.57721,1.644934, buf)

   benchmark('hypergeom(bad=10,good=25,sample=10)', 2.85714, 1.500600)
   benchmark('hypergeom(10,25,10,size=sz)', 2.85714, 1.500600, buf)
   benchmark('hypergeom(bad=10,good=25,sample=15)', 4.285714, 1.800720)
   benchmark('hypergeom(10,25,15,size=sz)', 4.285714, 1.800720, buf)
   benchmark('hypergeom(bad=10, good=25, sample=20)', 5.71429, 1.800720)
   benchmark('hypergeom(bad=25, good=10, sample=10)', 7.142857, 1.500600)
   benchmark('hypergeom(bad=10, good=25, sample=34)', 9.714286, 0.204082)
   benchmark('hypergeom(bad=10, good=25, sample=30)', 8.571429, 0.900360)
   benchmark('hypergeom(bad=10, good=25, sample=25)', 7.142857, 1.500600)
   benchmark('hypergeom(bad=25, good=10, sample=30)', 21.42857, 0.900360)
   benchmark('hypergeom(bad=25, good=475,sample=50)', 2.50, 2.14178)

   benchmark('Laplace(mu=0.0, scale=1.0)', 0.0, 2.0)
   benchmark('Laplace(mu=0.0, scale=1.0, size=sz)', 0.0, 2.0, buf)

   benchmark('logarithmic(p=0.5)', 1.442695, 0.804021)
   benchmark('logarithmic(p=0.5, size=sz)', 1.442695, 0.804021, buf)
   benchmark('logarithmic(p=0.975)',10.57232, 311.1188)
   benchmark('logarithmic(p=0.975, size=sz)', 10.57232, 311.1188, buf)

   benchmark('logistic(mu=0.0, scale=1.0)', 0.0, 3.289868)
   benchmark('logistic(mu=0.0, scale=1.0, size=sz)', 0.0, 3.289868, buf)

   benchmark('lognormal(mean=0.0, sigma=1.0)', 1.64872, 4.670774)
   benchmark('lognormal(mean=0.0, sigma=1.0, size=sz)',1.6487, 4.6708, buf)

   benchmark('negative_binomial(r=0.5, pr=0.5)', 0.5, 1.0)
   benchmark('negative_binomial(r=0.5,pr=0.5,size=sz)',0.5,1.0,buf)
   benchmark('negative_binomial(r=2.0, pr=0.5)', 2.0, 4.0)
   benchmark('negative_binomial(r=2.0,pr=0.5,size=sz)',2.0,4.0,buf)
   benchmark('negative_binomial(r=0.5, pr=0.9)', 4.5, 45.0)
   benchmark('negative_binomial(r=2.0, pr=0.1)', 0.22222, 0.24691)

   benchmark('normal(mu=0.0, sigma=1.0)', 0.0, 1.0)
   benchmark('normal(mu=0.0, sigma=1.0, size=sz)', 0.0, 1.0, buf)

   benchmark('Pareto(mode=1.0, shape=4.0)', 1.33333, 0.222222)
   benchmark('Pareto(mode=1.0, shape=4.0, size=sz)',1.33333,0.222222,buf)

   benchmark('Poisson(rate=5.0)', 5.0, 5.0)
   benchmark('Poisson(rate=5.0, size=sz)', 5.0, 5.0,buf)
   benchmark('Poisson(rate=20.0)', 20.0, 20.0)
   benchmark('Poisson(rate=20.0, size=sz)', 20.0, 20.0,buf)
   benchmark('Poisson(rate=200.0)', 200.0, 200.0)

   benchmark('randint(lowint=-1, upint=+1)', 0.0, 0.66667)
   benchmark('randint(lowint=-1, upint=+1, size=sz)', 0.0, 0.66667 ,buf)
   benchmark('randint(lowint=0, upint=1)', 0.5, 0.25)
   benchmark('randint(lowint=0, upint=1, size=sz)', 0.5, 0.25,buf)

   benchmark('Rayleigh(mode=1.0)', 1.253314, 0.429204)
   benchmark('Rayleigh(mode=1.0, size=sz)', 1.253314, 0.429204,buf)

   benchmark('Student_t(df=1.0)', 0.0, 1e100)
   benchmark('Student_t(df=1.0, size=sz)', 0.0, 1e100, buf)
   benchmark('Student_t(df=100.0)', 0.0, 1.0204082)
   benchmark('Student_t(df=100.0, size=sz)', 0.0, 1.0204082, buf)
   benchmark('Student_t(df=3.5)', 0.0, 2.333333)

   benchmark('triangular(left=0.0, mode=0.5, right=1.0)', 0.5, 0.04166667)
   benchmark('triangular(0.0, 0.5, 1.0, buf)', 0.5, 0.04166667, buf)
   benchmark('triangular(left=0.0, mode=0.0, right=1.0)', 0.33333, 0.055556)
   benchmark('triangular(0.0, 0.0, 1.0, buf)', 0.33333, 0.055556, buf)
   benchmark('triangular(left=0.0, mode=1.0, right=1.0)', 0.66667, 0.055556)
   benchmark('triangular(0.0, 1.0, 1.0, buf)', 0.66667, 0.055556, buf)

   benchmark('uniform(lower=-0.5, upper=+0.5)', 0.0, 0.0833333)
   benchmark('uniform(lower=-0.5, upper=+0.5, size=sz)',0.0,0.0833333,buf)

   benchmark('von_Mises(mean=0.0, shape=1.0)', 0.0, 1.61)
   benchmark('von_Mises(mean=0.0, shape=1.0, size=sz)', 0.0, 1.61,buf)

   benchmark('Wald(mean=1.0, scale=1.0)', 1.0, 1.0)
   benchmark('Wald(mean=1.0, scale=1.0, size=sz)', 1.0, 1.0,buf)

   benchmark('Weibull(scale=1.0, shape=0.5)', 2.01, 20.16)
   benchmark('Weibull(scale=1.0, shape=0.5, size=sz)', 2.01, 20.16,buf)

   benchmark('Zipf(a=4.0)', 1.110630, 0.28632)
   benchmark('Zipf(a=4.0, size=sz)', 1.110630, 0.28632, buf)

   save_state()
   print
   print  "   End of Module rv Test"

#if __name__ == '__main__': _rvtest()
