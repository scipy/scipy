# Copyright (c) Gary Strangman.  All rights reserved
#
# Disclaimer
# 
# This software is provided "as-is".  There are no expressed or implied
# warranties of any kind, including, but not limited to, the warranties
# of merchantability and fittness for a given application.  In no event
# shall Gary Strangman be liable for any direct, indirect, incidental,
# special, exemplary or consequential damages (including, but not limited
# to, loss of use, data or profits, or business interruption) however
# caused and on any theory of liability, whether in contract, strict
# liability or tort (including negligence or otherwise) arising in any way
# out of the use of this software, even if advised of the possibility of
# such damage.
#

#
# Heavily adapted for use by SciPy 2002 by Travis Oliphant
"""
stats.py module

(Requires pstat.py modules.)

#################################################
#######  Written by:  Gary Strangman  ###########
#######  Last modified:  Apr 13, 2000 ###########
#################################################

A collection of basic statistical functions for python.  The function
names appear below.

 *** Some scalar functions defined here are also available in the scipy.special 
     package where they work on arbitrary sized arrays. ****

Disclaimers:  The function list is obviously incomplete and, worse, the
functions are not optimized.  All functions have been tested (some more
so than others), but they are far from bulletproof.  Thus, as with any
free software, no warranty or guarantee is expressed or implied. :-)  A
few extra functions that don't appear in the list below can be found by
interested treasure-hunters.  These functions don't necessarily have
both list and array versions but were deemed useful

CENTRAL TENDENCY:  gmean    (geometric mean)
                   hmean    (harmonic mean)
                   mean
                   median
                   medianscore
                   mode

MOMENTS:  moment
          variation
          skew
          kurtosis
          normaltest (for arrays only)

ALTERED VERSIONS:  tmean 
                   tvar  
                   tstd
                   tsem  
                   describe

FREQUENCY STATS:  freqtable
                  itemfreq
                  scoreatpercentile
                  percentileofscore
                  histogram
                  cumfreq
                  relfreq

VARIABILITY:  obrientransform
              samplevar
              samplestd
              signaltonoise (for arrays only)
              var
              std
              stderr
              sem
              z
              zs

TRIMMING FCNS:  threshold (for arrays only)
                trimboth
                trim1
                around (round all vals to 'n' decimals)

CORRELATION FCNS:  paired
                   pearsonr
                   spearmanr
                   pointbiserialr
                   kendalltau
                   linregress

INFERENTIAL STATS:  ttest_1samp
                    ttest_ind
                    ttest_rel
                    chisquare
                    ks_2samp
                    mannwhitneyu
                    ranksums
                    wilcoxon
                    kruskal
                    friedmanchisquare

PROBABILITY CALCS:  chisqprob
                    erfcc
                    zprob
                    fprob
                    betai

## Note that scipy.stats.distributions has many more statistical probability functions defined.


ANOVA FUNCTIONS:  anova (NumPy required)
                  f_oneway
                  f_value

SUPPORT FUNCTIONS:  writecc
                    incr
                    sum
                    cumsum
                    ss
                    summult
                    square_of_sums
                    sumdiffsquared
                    shellsort
                    rankdata
                    outputpairedstats
                    findwithin
"""
## CHANGE LOG:
## ===========
## 02-02-10 ... require Numeric, eliminate "list-only" functions 
##              (only 1 set of functions now and no Dispatch class),
##              removed all references to aXXXX functions.
## 00-04-13 ... pulled all "global" statements, except from aanova()
##              added/fixed lots of documentation, removed io.py dependency
##              changed to version 0.5
## 99-11-13 ... added asign() function
## 99-11-01 ... changed version to 0.4 ... enough incremental changes now
## 99-10-25 ... added acovariance and acorrelation functions
## 99-10-10 ... fixed askew/akurtosis to avoid divide-by-zero errors
##              added aglm function (crude, but will be improved)
## 99-10-04 ... upgraded acumsum, ass, asummult, asamplevar, var, etc. to
##                   all handle lists of 'dimension's and keepdims
##              REMOVED ar0, ar2, ar3, ar4 and replaced them with around
##              reinserted fixes for abetai to avoid math overflows
## 99-09-05 ... rewrote achisqprob/aerfcc/aksprob/afprob/abetacf/abetai to
##                   handle multi-dimensional arrays (whew!)
## 99-08-30 ... fixed l/amoment, l/askew, l/akurtosis per D'Agostino (1990)
##              added anormaltest per same reference
##              re-wrote azprob to calc arrays of probs all at once
## 99-08-22 ... edited attest_ind printing section so arrays could be rounded
## 99-08-19 ... fixed amean and aharmonicmean for non-error(!) overflow on
##                   short/byte arrays (mean of #s btw 100-300 = -150??)
## 99-08-09 ... fixed asum so that the None case works for Byte arrays
## 99-08-08 ... fixed 7/3 'improvement' to handle t-calcs on N-D arrays
## 99-07-03 ... improved attest_ind, attest_rel (zero-division errortrap)
## 99-06-24 ... fixed bug(?) in attest_ind (n1=a.shape[0])
## 04/11/99 ... added asignaltonoise, athreshold functions, changed all
##                   max/min in array section to maximum/minimum,
##                   fixed square_of_sums to prevent integer overflow
## 04/10/99 ... !!! Changed function name ... sumsquared ==> square_of_sums
## 03/18/99 ... Added ar0, ar2, ar3 and ar4 rounding functions
## 02/28/99 ... Fixed aobrientransform to return an array rather than a list
## 01/15/99 ... Essentially ceased updating list-versions of functions (!!!)
## 01/13/99 ... CHANGED TO VERSION 0.3
##              fixed bug in a/lmannwhitneyu p-value calculation
## 12/31/98 ... fixed variable-name bug in ldescribe
## 12/19/98 ... fixed bug in findwithin (fcns needed pstat. prefix)
## 12/16/98 ... changed amedianscore to return float (not array) for 1 score
## 12/14/98 ... added atmin and atmax functions
##              removed umath from import line (not needed)
##              l/ageometricmean modified to reduce chance of overflows (take
##                   nth root first, then multiply)
## 12/07/98 ... added __version__variable (now 0.2)
##              removed all 'stats.' from anova() fcn
## 12/06/98 ... changed those functions (except shellsort) that altered
##                   arguments in-place ... cumsum, ranksort, ...
##              updated (and fixed some) doc-strings
## 12/01/98 ... added anova() function (requires NumPy)
##              incorporated Dispatch class
## 11/12/98 ... added functionality to amean, aharmonicmean, ageometricmean
##              added 'asum' function (added functionality to add.reduce)
##              fixed both moment and amoment (two errors)
##              changed name of skewness and askewness to skew and askew
##              fixed (a)histogram (which sometimes counted points <lowerlimit)



import string, sys, _support
from types import *

__version__ = 0.8

from scipy_base import *
import scipy_base as sb
import scipy_base.fastumath as math
from scipy_base.fastumath import *
import Numeric
N = Numeric
import scipy.special as special
import scipy.linalg as linalg
import scipy

SequenceType = [ListType, TupleType, ArrayType]


# These two functions replace letting axis be a sequence and the
#  keepdims features used throughout.  These ideas
#  did not match the rest of Scipy.
from scipy_base import expand_dims, apply_over_axes

def _chk_asarray(a, axis):
    if axis is None:
        a = ravel(a)
        outaxis = 0
    else:
        a = asarray(a)
        outaxis = axis
    return a, outaxis

def _chk2_asarray(a, b,  axis):
    if axis is None:
        a = ravel(a)
        b = ravel(b)
        outaxis = 0
    else:
        a = asarray(a)
        b = asarray(b)
        outaxis = axis
    return a, b, outaxis

#######
### NAN friendly functions
########

def nanmean(x,axis=-1):
    """Compute the mean over the given axis ignoring nans.
    """
    x, axis = _chk_asarray(x,axis)
    x = x.copy()
    Norig = x.shape[axis]
    factor = 1.0-sum(isnan(x),axis)*1.0/Norig
    n = N-sum(isnan(x),axis)
    N.putmask(x,isnan(x),0)
    return stats.mean(x,axis)/factor
    
def nanstd(x,axis=-1,bias=0):
    """Compute the standard deviation over the given axis ignoring nans
    """
    x, axis = _chk_asarray(x,axis)
    x = x.copy()
    Norig = x.shape[axis]
    n = Norig - sum(isnan(x),axis)*1.0
    factor = n/Norig
    n = N-sum(isnan(x),axis)
    N.putmask(x,isnan(x),0)
    m1 = stats.mean(x,axis)
    m1c = m1/factor
    m2 = stats.mean((x-m1c)**2.0,axis)
    if bias:
        m2c = m2/factor
    else:
        m2c = m2*Norig/(n-1.0)
    return m2c

def _nanmedian(arr1d):  # This only works on 1d arrays
   cond = 1-isnan(arr1d)
   x = sb.sort(compress(cond,arr1d))
   return median(x)
   
def nanmedian(x, axis=-1):
    """ Compute the median along the given axis ignoring nan values
    """
    x, axis = _chk_asarray(x,axis)
    x = x.copy()
    return sb.apply_along_axis(_nanmedian,axis,x)


#####################################
########  ACENTRAL TENDENCY  ########
#####################################

def gmean (a,axis=-1):
    """Calculates the geometric mean of the values in the passed array.
    
    That is:  n-th root of (x1 * x2 * ... * xn).
    
    If a is 1D, a single value is returned.  If a is multi-dimensional, 
    the geometric mean along the dimension specified is calculated.  The
    returned array has one less dimension than a.  dimension defaults
    to the last dimension of the array.  This means that, for a two
    dimensional array, the default is to calculate the geometric mean
    of each row.
    """
    a, axis = _chk_asarray(a, axis)
    size = a.shape[axis]
    prod = product(a,axis)
    return power(prod,1./size)


def hmean(a,axis=-1):
    """Calculates the harmonic mean of the values in the passed array.

    That is:  n / (1/x1 + 1/x2 + ... + 1/xn).  Defaults to ALL values in 
    the passed array.  REMEMBER: if axis=0, it collapses over 
    axis 0 ('rows' in a 2D array) only, and if axis is a 
    sequence, it collapses over all specified axes.
    
    Returns: harmonic mean computed over dim(s) in axis     
    """
    a, axis = _chk_asarray(a, axis)
    size = a.shape[axis]
    return size / add.reduce(1.0/a, axis)

def mean(a,axis=-1):
    """Returns the mean of m along the given dimension.

       If m is of integer type, returns a floating point answer.
    """
    a, axis = _chk_asarray(a, axis)
    return add.reduce(a,axis)/float(a.shape[axis])

def cmedian (a,numbins=1000):
    """Calculates the COMPUTED median value of an array of numbers, given the
    number of bins to use for the histogram (more bins approaches finding the
    precise median value of the array; default number of bins = 1000).  From
    G.W. Heiman's Basic Stats, or CRC Probability & Statistics.
    NOTE:  THIS ROUTINE ALWAYS uses the entire passed array (flattens it first).
    
    Returns: median calculated over ALL values in the array
"""
    a = ravel(a)
    (hist, smallest, binsize, extras) = histogram(a,numbins)
    cumhist = cumsum(hist)            # make cumulative histogram
    otherbins = greater_equal(cumhist,len(a)/2.0)
    otherbins = list(otherbins)         # list of 0/1s, 1s start at median bin
    cfbin = otherbins.index(1)          # get 1st(!) index holding 50%ile score
    LRL = smallest + binsize*cfbin      # get lower read limit of that bin
    cfbelow = add.reduce(hist[0:cfbin])       # cum. freq. below bin
    freq = hist[cfbin]                  # frequency IN the 50%ile bin
    median = LRL + ((len(a)/2.0-cfbelow)/float(freq))*binsize # MEDIAN
    return median

def median (a,axis=-1):
    """ Returns the median of the passed array along the given axis.
    If there is an even number of entires, the mean of the
    2 middle values is returned.

    Returns: median of the array, or the mean of the 2 middle values.
    """
    a, axis = _chk_asarray(a, axis)
    a = sort(a,axis)
    indx = int(a.shape[axis]/2.0)
    if a.shape[axis] % 2 == 0:   # if even number of elements
        median = asarray(take(a,[indx],axis)+take(a,[indx-1],axis)) / 2.0
    else:
        median = take(a,[indx],axis)
    return sb.squeeze(median)

def mode(a, axis=-1):
    """Returns an array of the modal (most common) value in the passed array.
    If there is more than one such value, ONLY THE FIRST is returned.
    The bin-count for the modal bins is also returned.  Operates on whole
    array (axis=None), or on a given axis.

    Returns: array of bin-counts for mode(s), array of corresponding modal value
    """

    a, axis = _chk_asarray(a, axis)
    scores = sb.unique(ravel(a))       # get ALL unique values
    testshape = list(a.shape)
    testshape[axis] = 1
    oldmostfreq = zeros(testshape)
    oldcounts = zeros(testshape)
    for score in scores:
        template = equal(a,score)
        counts = expand_dims(sum(template,axis),axis)
        mostfrequent = where(greater(counts,oldcounts),score,oldmostfreq)
        oldcounts = where(greater(counts,oldcounts),counts,oldcounts)
        oldmostfreq = mostfrequent
    return mostfrequent, oldcounts

def tmean(a,limits=None,inclusive=(1,1)):
     """Returns the arithmetic mean of all values in an array, ignoring values
     strictly outside the sequence passed to 'limits'.   Note: either limit
     in the sequence, or the value of limits itself, can be set to None.  The
     inclusive list/tuple determines whether the lower and upper limiting bounds
     (respectively) are open/exclusive (0) or closed/inclusive (1).

     """
     a = asarray(a)
     if a.typecode() in ['l','s','b']:
         a = a.astype(Float)
     if limits is None:
         return mean(a,None)
     assert type(limits) in SequenceType, "Wrong type for limits in tmean"
     if inclusive[0]:    lowerfcn = greater_equal
     else:               lowerfcn = greater
     if inclusive[1]:    upperfcn = less_equal
     else:               upperfcn = less
     if limits[0] > maximum.reduce(ravel(a)) or limits[1] < minimum.reduce(ravel(a)):
         raise ValueError, "No array values within given limits (tmean)."
     elif limits[0] is None and limits[1] is not None:
         mask = upperfcn(a,limits[1])
     elif limits[0] is not None and limits[1] is None:
         mask = lowerfcn(a,limits[0])
     elif limits[0] is not None and limits[1] is not None:
         mask = lowerfcn(a,limits[0])*upperfcn(a,limits[1])
     s = float(add.reduce(ravel(a*mask)))
     n = float(add.reduce(ravel(mask)))
     return s/n


def tvar(a,limits=None,inclusive=(1,1)):
     """Returns the sample variance of values in an array, (i.e., using N-1),
     ignoring values strictly outside the sequence passed to 'limits'.  
     Note: either limit in the sequence, or the value of limits itself,
     can be set to None.  The inclusive list/tuple determines whether the lower
     and upper limiting bounds (respectively) are open/exclusive (0) or
     closed/inclusive (1).
     """
     a = asarray(a)
     a = a.astype(Float)
     if limits is None or limits == [None,None]:
         term1 = add.reduce(ravel(a*a))
         n = float(len(ravel(a))) - 1
         term2 = add.reduce(ravel(a))**2 / n
         return (term1 - term2) / n
     assert type(limits) in SequenceType, "Wrong type for limits in tvar"
     if inclusive[0]:    lowerfcn = greater_equal
     else:               lowerfcn = greater
     if inclusive[1]:    upperfcn = less_equal
     else:               upperfcn = less
     if limits[0] > maximum.reduce(ravel(a)) or limits[1] < minimum.reduce(ravel(a)):
         raise ValueError, "No array values within given limits (tvar)."
     elif limits[0] is None and limits[1] is not None:
         mask = upperfcn(a,limits[1])
     elif limits[0] is not None and limits[1] is None:
         mask = lowerfcn(a,limits[0])
     elif limits[0] is not None and limits[1] is not None:
         mask = lowerfcn(a,limits[0])*upperfcn(a,limits[1])
     term1 = add.reduce(ravel(a*a*mask))
     n = float(add.reduce(ravel(mask))) - 1
     term2 = add.reduce(ravel(a*mask))**2 / n
     return (term1 - term2) / n

def tmin(a,lowerlimit=None,axis=-1,inclusive=1):
     """Returns the minimum value of a, along axis, including only values less
     than (or equal to, if inclusive=1) lowerlimit.  If the limit is set to None,
     all values in the array are used.
     """
     a, axis = _chk_asarray(a, axis)
     if inclusive:       lowerfcn = greater
     else:               lowerfcn = greater_equal
     if lowerlimit is None:
         lowerlimit = minimum.reduce(ravel(a))-11
     biggest = maximum.reduce(ravel(a))
     ta = where(lowerfcn(a,lowerlimit),a,biggest)
     return minimum.reduce(ta,axis)

def tmax(a,upperlimit,axis=-1,inclusive=1):
     """Returns the maximum value of a, along axis, including only values greater
     than (or equal to, if inclusive=1) upperlimit.  If the limit is set to None,
     a limit larger than the max value in the array is used.
     """
     a, axis = asarray(a, axis)
     if inclusive:       upperfcn = less
     else:               upperfcn = less_equal
     if upperlimit is None:
         upperlimit = maximum.reduce(ravel(a))+1
     smallest = minimum.reduce(ravel(a))
     ta = where(upperfcn(a,upperlimit),a,smallest)
     return maximum.reduce(ta,axis)


def tstd(a,limits=None,inclusive=(1,1)):
     """Returns the standard deviation of all values in an array, ignoring values
     strictly outside the sequence passed to 'limits'.   Note: either limit
     in the sequence, or the value of limits itself, can be set to None.  The
     inclusive list/tuple determines whether the lower and upper limiting bounds
     (respectively) are open/exclusive (0) or closed/inclusive (1).     
     """
     return sqrt(tvar(a,limits,inclusive))


def tsem(a,limits=None,inclusive=(1,1)):
    """Returns the standard error of the mean for the values in an array,
    (i.e., using N for the denominator), ignoring values strictly outside
    the sequence passed to 'limits'.   Note: either limit in the sequence,
    or the value of limits itself, can be set to None.  The inclusive list/tuple
    determines whether the lower and upper limiting bounds (respectively) are
    open/exclusive (0) or closed/inclusive (1).
    """
    a = asarray(a)
    sd = tstd(a,limits,inclusive)
    if limits is None or limits == [None,None]:
        n = float(len(ravel(a)))
    assert type(limits) in SequenceType, "Wrong type for limits in tsem"
    if inclusive[0]:    lowerfcn = greater_equal
    else:               lowerfcn = greater
    if inclusive[1]:    upperfcn = less_equal
    else:               upperfcn = less
    if limits[0] > maximum.reduce(ravel(a)) or limits[1] < minimum.reduce(ravel(a)):
        raise ValueError, "No array values within given limits (tsem)."
    elif limits[0] is None and limits[1] is not None:
        mask = upperfcn(a,limits[1])
    elif limits[0] is not None and limits[1] is None:
        mask = lowerfcn(a,limits[0])
    elif limits[0] is not None and limits[1] is not None:
        mask = lowerfcn(a,limits[0])*upperfcn(a,limits[1])
    #term1 = add.reduce(ravel(a*a*mask))
    n = float(add.reduce(ravel(mask)))
    return sd/math.sqrt(n)


#####################################
############  AMOMENTS  #############
#####################################

def moment(a,moment=1,axis=-1):
    """Calculates the nth moment about the mean for a sample (defaults to the
    1st moment).  Generally used to calculate coefficients of skewness and
    kurtosis.  Axis can equal None (ravel array first), or an integer
    (the axis over which to operate).
        
    Returns: appropriate moment along given axis
    """
    a, axis = _chk_asarray(a, axis)
    if moment == 1:
        return 0.0
    else:
        mn = expand_dims(mean(a,axis),axis)
        s = power((a-mn),moment)
        return mean(s,axis)


def variation(a,axis=-1):
    """Returns the coefficient of variation, as defined in CRC Standard
    Probability and Statistics, p.6. Axis can equal None (ravel array
    first), or an integer (the axis over which to operate)
    """
    return samplestd(a,axis)/mean(a,axis)


def skew(a,axis=-1,bias=1):
    """Returns the skewness of a distribution (normal ==> 0.0; >0 means extra
    weight in left tail).  Use skewtest() to see if it's close enough.
    Axis can equal None (ravel array first), or an integer (the
    axis over which to operate).
    
    Returns: skew of vals in a along axis, returning ZERO where all vals equal
    """
    a, axis = _chk_asarray(a,axis)
    n = a.shape[axis]
    m2 = moment(a,2,axis)
    m3 = moment(a,3,axis)
    zero = (m2 == 0)
    vals = where(zero, 0, m3/ m2**1.5)
    if not bias:
        can_correct = (n > 2) & (m2 > 0)
        if sb.any(can_correct):
            m2 = sb.extract(can_correct,m2)
            m3 = sb.extract(can_correct,m3)
            nval = sqrt((n-1.0)*n)/(n-2.0)*m3/m2**1.5
            sb.insert(vals, can_correct, nval)
    return vals

def kurtosis(a,axis=-1,fisher=1,bias=1):
    """Returns the kurtosis (fisher or pearson) of a distribution

    Kurtosis is the fourth central moment divided by the square of the
      variance.  If Fisher's definition is used, then 3.0 is subtracted
      from the result to give 0.0 for a normal distribution. 

    Axis can equal None (ravel arrayfirst), or an integer
    (the axis over which to operate)

    If bias is 0 then the kurtosis is calculated using k statistics to
       eliminate bias comming from biased moment estimators

    Use kurtosistest() to see if result is close enough to normal.

    Returns: kurtosis of values in a along axis, and ZERO where all vals equal
    """
    a, axis = _chk_asarray(a, axis)
    n = a.shape[axis]
    m2 = moment(a,2,axis)
    m4 = moment(a,4,axis)
    zero = (m2 == 0)
    vals = where(zero, 0, m4/ m2**2.0)
    if not bias:
        can_correct = (n > 3) & (m2 > 0)
        if sb.any(can_correct):
            m2 = sb.extract(can_correct,m2)
            m4 = sb.extract(can_correct,m4)
            nval = 1.0/(n-2)/(n-3)*((n*n-1.0)*m4/m2**2.0-3*(n-1)**2.0)
            sb.insert(vals, can_correct, nval+3.0)
    if fisher:
        return vals - 3
    else:
        return vals

def describe(a,axis=-1):
     """Returns several descriptive statistics of the passed array.  Axis
     can equal None (ravel array first), or an integer (the axis over
     which to operate)
     
     Returns: n, (min,max), mean, standard deviation, skew, kurtosis
     """
     a, axis = _chk_asarray(a, axis)
     n = a.shape[axis]
     mm = (minimum.reduce(a),maximum.reduce(a))
     m = mean(a,axis)
     v = var(a,axis)
     sk = skew(a,axis)
     kurt = kurtosis(a,axis)
     return n, mm, m, v, sk, kurt

#####################################
########  NORMALITY TESTS  ##########
#####################################

def skewtest(a,axis=-1):
    """Tests whether the skew is significantly different from a normal
    distribution.  Axis can equal None (ravel array first), an integer
    (the axis over which to operate), or a sequence (operate over
    multiple axes).

    Returns: z-score and 2-tail z-probability
    """
    a, axis = _chk_asarray(a, axis)
    if axis is None:
        a = ravel(a)
        axis = 0
    b2 = skew(a,axis)
    n = float(a.shape[axis])
    if n<8:
        print "kurtosistest only valid for n>=8 ... continuing anyway, n=",n
    y = b2 * sqrt(((n+1)*(n+3)) / (6.0*(n-2)) )
    beta2 = ( 3.0*(n*n+27*n-70)*(n+1)*(n+3) ) / ( (n-2.0)*(n+5)*(n+7)*(n+9) )
    W2 = -1 + sqrt(2*(beta2-1))
    delta = 1/sqrt(log(sqrt(W2)))
    alpha = sqrt(2.0/(W2-1))
    y = where(equal(y,0),1,y)
    Z = delta*log(y/alpha + sqrt((y/alpha)**2+1))
    return Z, (1.0-zprob(Z))*2


def kurtosistest(a,axis=-1):
    """
Tests whether a dataset has normal kurtosis (i.e.,
kurtosis=3(n-1)/(n+1)) Valid only for n>20.  Axis can equal None
(ravel array first), an integer (the axis over which to operate),
or a sequence (operate over multiple axes).

Returns: z-score and 2-tail z-probability, returns 0 for bad pixels
"""
    a, axis = _chk_asarray(a, axis)
    n = float(a.shape[axis])
    if n<20:
        print "kurtosistest only valid for n>=20 ... continuing anyway, n=",n
    b2 = kurtosis(a,axis)
    E = 3.0*(n-1) /(n+1)
    varb2 = 24.0*n*(n-2)*(n-3) / ((n+1)*(n+1)*(n+3)*(n+5))
    x = (b2-E)/sqrt(varb2)
    sqrtbeta1 = 6.0*(n*n-5*n+2)/((n+7)*(n+9)) * sqrt((6.0*(n+3)*(n+5))/
                                                       (n*(n-2)*(n-3)))
    A = 6.0 + 8.0/sqrtbeta1 *(2.0/sqrtbeta1 + sqrt(1+4.0/(sqrtbeta1**2)))
    term1 = 1 -2/(9.0*A)
    denom = 1 +x*sqrt(2/(A-4.0))
    denom = where(less(denom,0), 99, denom)
    term2 = where(equal(denom,0), term1, power((1-2.0/A)/denom,1/3.0))
    Z = ( term1 - term2 ) / sqrt(2/(9.0*A))
    Z = where(equal(denom,99), 0, Z)
    return Z, (1.0-zprob(Z))*2


def normaltest(a,axis=-1):
    """
Tests whether skew and/OR kurtosis of dataset differs from normal
curve.  Can operate over multiple axes.  Axis can equal
None (ravel array first), an integer (the axis over which to
operate), or a sequence (operate over multiple axes).

Omnibus test of D'Agostino and Pearson, 1973

Returns: Score and 2-tail probability
"""
    a, axis = _chk_asarray(a, axis)
    s,p = skewtest(a,axis)
    k,p = kurtosistest(a,axis)
    k2 = s*s + k*k
    return k2, chisqprob(k2,2)

# Martinez-Iglewicz test
# K-S test

#####################################
######  AFREQUENCY FUNCTIONS  #######
#####################################

def itemfreq(a):
    """
Returns a 2D array of item frequencies.  Column 1 contains item values,
column 2 contains their respective counts.  Assumes a 1D array is passed.

Returns: a 2D frequency table (col [0:n-1]=scores, col n=frequencies)
"""
    scores = _support.unique(a)
    scores = sort(scores)
    freq = zeros(len(scores))
    for i in range(len(scores)):
        freq[i] = add.reduce(equal(a,scores[i]))
    return array(_support.abut(scores, freq))


def scoreatpercentile (a, percent):
    """
Returns: score at given percentile, relative to a distribution
"""
    percent = percent / 100.0
    targetcf = percent*len(a)
    h, lrl, binsize, extras = histogram(a)
    cumhist = cumsum(h*1)
    for i in range(len(cumhist)):
        if cumhist[i] >= targetcf:
            break
    score = binsize * ((targetcf - cumhist[i-1]) / float(h[i])) + (lrl+binsize*i)
    return score


def percentileofscore (a,score,histbins=10,defaultlimits=None):
    """
Note: result of this function depends on the values used to histogram
the data(!).

Returns: percentile-position of score (0-100) relative to a
"""
    h, lrl, binsize, extras = histogram(a,histbins,defaultlimits)
    cumhist = cumsum(h*1)
    i = int((score - lrl)/float(binsize))
    pct = (cumhist[i-1]+((score-(lrl+binsize*i))/float(binsize))*h[i])/float(len(a)) * 100
    return pct


def histogram2(a, bins):
    """ histogram2(a,bins) -- Compute histogram of a using divisions in bins

         Description:
            Count the number of times values from array a fall into
            numerical ranges defined by bins.  Range x is given by
            bins[x] <= range_x < bins[x+1] where x =0,N and N is the
            length of the bins array.  The last range is given by
            bins[N] <= range_N < infinity.  Values less than bins[0] are
            not included in the histogram.
         Arguments:
            a -- 1D array.  The array of values to be divied into bins
            bins -- 1D array.  Defines the ranges of values to use during
                    histogramming.
         Returns:
            1D array.  Each value represents the occurences for a given
            bin (range) of values.

         Caveat:
            This should probably have an axis argument that would histogram
            along a specific axis (kinda like matlab)

    """
    n = Numeric.searchsorted(Numeric.sort(a), bins)
    n = Numeric.concatenate([ n, [len(a)]])
    return n[ 1:]-n[:-1]




def histogram (a,numbins=10,defaultlimits=None,printextras=1):
    """
Returns (i) an array of histogram bin counts, (ii) the smallest value
of the histogram binning, and (iii) the bin width (the last 2 are not
necessarily integers).  Default number of bins is 10.  Defaultlimits
can be None (the routine picks bins spanning all the numbers in the
a) or a 2-sequence (lowerlimit, upperlimit).  Returns all of the
following: array of bin values, lowerreallimit, binsize, extrapoints.

Returns: (array of bin counts, bin-minimum, min-width, #-points-outside-range)
"""
    a = ravel(a)               # flatten any >1D arrays
    if (defaultlimits <> None):
        lowerreallimit = defaultlimits[0]
        upperreallimit = defaultlimits[1]
        binsize = (upperreallimit-lowerreallimit) / float(numbins)
    else:
        Min = minimum.reduce(a)
        Max = maximum.reduce(a)
        estbinwidth = float(Max - Min)/float(numbins) + 1
        binsize = (Max-Min+estbinwidth)/float(numbins)
        lowerreallimit = Min - binsize/2.0  #lower real limit,1st bin
    bins = zeros(numbins)
    extrapoints = 0
    for num in a:
        try:
            if (num-lowerreallimit) < 0:
                extrapoints = extrapoints + 1
            else:
                bintoincrement = int((num-lowerreallimit) / float(binsize))
                bins[bintoincrement] = bins[bintoincrement] + 1
        except:                           # point outside lower/upper limits
            extrapoints = extrapoints + 1
    if (extrapoints > 0 and printextras == 1):
        print '\nPoints outside given histogram range =',extrapoints
    return (bins, lowerreallimit, binsize, extrapoints)


def cumfreq(a,numbins=10,defaultreallimits=None):
    """
Returns a cumulative frequency histogram, using the histogram function.
Defaultreallimits can be None (use all data), or a 2-sequence containing
lower and upper limits on values to include.

Returns: array of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(a,numbins,defaultreallimits)
    cumhist = cumsum(h*1)
    return cumhist,l,b,e


def relfreq(a,numbins=10,defaultreallimits=None):
    """
Returns a relative frequency histogram, using the histogram function.
Defaultreallimits can be None (use all data), or a 2-sequence containing
lower and upper limits on values to include.

Returns: array of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(a,numbins,defaultreallimits)
    h = array(h/float(a.shape[0]))
    return h,l,b,e


#####################################
######  AVARIABILITY FUNCTIONS  #####
#####################################

def obrientransform(*args):
    """
Computes a transform on input data (any number of columns).  Used to
test for homogeneity of variance prior to running one-way stats.  Each
array in *args is one level of a factor.  If an F_oneway() run on the
transformed data and found significant, variances are unequal.   From
Maxwell and Delaney, p.112.

Returns: transformed data for use in an ANOVA
"""
    TINY = 1e-10
    k = len(args)
    n = zeros(k,Float)
    v = zeros(k,Float)
    m = zeros(k,Float)
    nargs = []
    for i in range(k):
        nargs.append(args[i].astype(Float))
        n[i] = float(len(nargs[i]))
        v[i] = var(nargs[i])
        m[i] = mean(nargs[i],None)
    for j in range(k):
        for i in range(n[j]):
            t1 = (n[j]-1.5)*n[j]*(nargs[j][i]-m[j])**2
            t2 = 0.5*v[j]*(n[j]-1.0)
            t3 = (n[j]-1.0)*(n[j]-2.0)
            nargs[j][i] = (t1-t2) / float(t3)
    check = 1
    for j in range(k):
        if v[j] - mean(nargs[j],None) > TINY:
            check = 0
    if check <> 1:
        raise ValueError, 'Lack of convergence in obrientransform.'
    else:
        return array(nargs)


def samplevar (a,axis=-1):
    """
Returns the sample standard deviation of the values in the passed
array (i.e., using N).  Axis can equal None (ravel array first),
an integer (the axis over which to operate)
"""
    a, axis = _chk_asarray(a, axis)
    mn = expand_dims(mean(a, axis), axis)
    deviations = a - mn 
    n = a.shape[axis]
    svar = ss(deviations,axis) / float(n)
    return svar


def samplestd (a, axis=-1):
    """Returns the sample standard deviation of the values in the passed
array (i.e., using N).  Axis can equal None (ravel array first),
an integer (the axis over which to operate).
"""
    return sqrt(samplevar(a,axis))


def signaltonoise(instack,axis=-1):
    """
Calculates signal-to-noise.  Axis can equal None (ravel array
first), an integer (the axis over which to operate).

Returns: array containing the value of (mean/stdev) along axis,
         or 0 when stdev=0
"""
    m = mean(instack,axis)
    sd = samplestd(instack,axis)
    return where(equal(sd,0),0,m/sd)


def var(a, axis=-1, bias=0):
    """
Returns the estimated population variance of the values in the passed
array (i.e., N-1).  Axis can equal None (ravel array first), or an
integer (the axis over which to operate).
"""
    a, axis = _chk_asarray(a, axis)
    mn = expand_dims(mean(a,axis),axis)
    deviations = a - mn
    n = a.shape[axis]
    vals = ss(deviations,axis)/(n-1.0)
    if bias:
        return vals * (n-1.0)/n
    else:
        return vals

def std (a, axis=-1, bias=0):
    """
Returns the estimated population standard deviation of the values in
the passed array (i.e., N-1).  Axis can equal None (ravel array
first), or an integer (the axis over which to operate).
"""
    return sqrt(var(a,axis,bias))


def stderr (a, axis=-1):
    """
Returns the estimated population standard error of the values in the
passed array (i.e., N-1).  Axis can equal None (ravel array
first), or an integer (the axis over which to operate).
"""
    a, axis = _chk_asarray(a, axis)
    return std(a,axis) / float(sqrt(a.shape[axis]))


def sem (a, axis=-1):
    """
Returns the standard error of the mean (i.e., using N) of the values
in the passed array.  Axis can equal None (ravel array first), or an
integer (the axis over which to operate)
"""
    a, axis = _chk_asarray(a, axis)
    n = a.shape[axis]
    s = samplestd(a,axis) / sqrt(n-1)
    return s


def z (a, score):
    """
Returns the z-score of a given input score, given thearray from which
that score came.  Not appropriate for population calculations, nor for
arrays > 1D.

"""
    z = (score-mean(a,None)) / samplestd(a)
    return z


def zs (a):
    """
Returns a 1D array of z-scores, one for each score in the passed array,
computed relative to the passed array.

"""
    zscores = []
    for item in a:
        zscores.append(z(a,item))
    return array(zscores)


def zmap (scores, compare, axis=-1):
    """
Returns an array of z-scores the shape of scores (e.g., [x,y]), compared to
array passed to compare (e.g., [time,x,y]).  Assumes collapsing over dim 0
of the compare array.

"""
    mns = mean(compare,axis)
    sstd = samplestd(compare,0)
    return (scores - mns) / sstd


#####################################
#######  ATRIMMING FUNCTIONS  #######
#####################################

# Removed round --- same as Numeric.around 

def threshold(a,threshmin=None,threshmax=None,newval=0):
    """
Like Numeric.clip() except that values <threshmid or >threshmax are replaced
by newval instead of by threshmin/threshmax (respectively).

Returns: a, with values <threshmin or >threshmax replaced with newval
"""
    a = asarray(a)
    mask = zeros(a.shape)
    if threshmin <> None:
        mask = mask + where(less(a,threshmin),1,0)
    if threshmax <> None:
        mask = mask + where(greater(a,threshmax),1,0)
    mask = clip(mask,0,1)
    return where(mask,newval,a)


def trimboth (a,proportiontocut):
    """
Slices off the passed proportion of items from BOTH ends of the passed
array (i.e., with proportiontocut=0.1, slices 'leftmost' 10% AND
'rightmost' 10% of scores.  You must pre-sort the array if you want
"proper" trimming.  Slices off LESS if proportion results in a
non-integer slice index (i.e., conservatively slices off
proportiontocut).

Returns: trimmed version of array a
"""
    a = asarray(a)
    lowercut = int(proportiontocut*len(a))
    uppercut = len(a) - lowercut
    if (lowercut >= uppercut):
        raise ValueError, "Proportion too big."
    return a[lowercut:uppercut]


def trim1 (a,proportiontocut,tail='right'):
    """
    Slices off the passed proportion of items from ONE end of the passed
    array (i.e., if proportiontocut=0.1, slices off 'leftmost' or 'rightmost'
    10% of scores).  Slices off LESS if proportion results in a non-integer
    slice index (i.e., conservatively slices off proportiontocut).
    
    Returns: trimmed version of array a
    """
    a = asarray(a)
    if string.lower(tail) == 'right':
        lowercut = 0
        uppercut = len(a) - int(proportiontocut*len(a))
    elif string.lower(tail) == 'left':
        lowercut = int(proportiontocut*len(a))
        uppercut = len(a)
    return a[lowercut:uppercut]

def trim_mean(a,proportiontocut):
    """Return mean with proportiontocut chopped from each of the lower and
    upper tails.
    """
    newa = trimboth(N.sort(a),proportiontocut)
    return mean(newa)
        


#####################################
#####  ACORRELATION FUNCTIONS  ######
#####################################

#  Cov is more flexible than the original
#    covariance and computes an unbiased covariance matrix
#    by default. 
def cov(m,y=None, rowvar=0, bias=0):
    """Estimate the covariance matrix.

    If m is a vector, return the variance.  For matrices where each row
    is an observation, and each column a variable, return the covariance
    matrix.  Note that in this case diag(cov(m)) is a vector of
    variances for each column.

    cov(m) is the same as cov(m, m)

    Normalization is by (N-1) where N is the number of observations
    (unbiased estimate).  If bias is 1 then normalization is by N.

    If rowvar is zero, then each row is a variables with
    observations in the columns.
    """
    m = asarray(m)
    if y is None:
        y = m
    else:
        y = asarray(y)
    if rowvar:
        m = transpose(m)
        y = transpose(y)
    N = m.shape[0]
    if (y.shape[0] != N):
        raise ValueError, "x and y must have the same number of observations."    
    m = m - mean(m,axis=0)
    y = y - mean(y,axis=0)
    if bias:
        fact = N*1.0
    else:
        fact = N-1.0
    val = sb.squeeze(dot(transpose(m),conjugate(y))) / fact
    return val

def corrcoef(x, y=None, rowvar=0, bias=0):
    """The correlation coefficients formed from 2-d array x, where the
    rows are the observations, and the columns are variables.

    corrcoef(x,y) where x and y are 1d arrays is the same as
    corrcoef(transpose([x,y]))
    
    If rowvar is zero, then each row is a variables with
    observations in the columns.
    """
    if y is not None:
        x = transpose([x,y])
        y = None
    c = cov(x, y, rowvar=rowvar, bias=bias)
    d = scipy.diag(c)
    return c/sqrt(multiply.outer(d,d))



def f_oneway(*args):
    """
Performs a 1-way ANOVA, returning an F-value and probability given
any number of groups.  From Heiman, pp.394-7.

Usage:   f_oneway (*args)    where *args is 2 or more arrays, one per
                                  treatment group
Returns: f-value, probability
"""
    na = len(args)            # ANOVA on 'na' groups, each in it's own array
    tmp = map(N.array,args)
    #means = map(mean,tmp)
    #vars = map(var,tmp)
    #ns = map(len,args)
    alldata = N.concatenate(args)
    bign = len(alldata)
    sstot = ss(alldata)-(square_of_sums(alldata)/float(bign))
    ssbn = 0
    for a in args:
        ssbn = ssbn + square_of_sums(N.array(a))/float(len(a))
    ssbn = ssbn - (square_of_sums(alldata)/float(bign))
    sswn = sstot-ssbn
    dfbn = na-1
    dfwn = bign - na
    msb = ssbn/float(dfbn)
    msw = sswn/float(dfwn)
    f = msb/msw
    prob = fprob(dfbn,dfwn,f)
    return f, prob


def paired(x,y):
    """
Interactively determines the type of data in x and y, and then runs the
appropriated statistic for paired group data.

Returns: appropriate statistic name, value, and probability
"""
    x,y = map(asarray, (x,y))
    samples = ''
    while samples not in ['i','r','I','R','c','C']:
        print '\nIndependent or related samples, or correlation (i,r,c): ',
        samples = raw_input()

    if samples in ['i','I','r','R']:
        print '\nComparing variances ...',
# USE O'BRIEN'S TEST FOR HOMOGENEITY OF VARIANCE, Maxwell & delaney, p.112
        r = obrientransform(x,y)
        f,p = f_oneway(_support.colex(r,0),_support.colex(r,1))
        if p<0.05:
            vartype='unequal, p='+str(around(p,4))
        else:
            vartype='equal'
        print vartype
        if samples in ['i','I']:
            if vartype[0]=='e':
                t,p = ttest_ind(x,y,None,0)
                print '\nIndependent samples t-test:  ', around(t,4),around(p,4)
            else:
                if len(x)>20 or len(y)>20:
                    z,p = ranksums(x,y)
                    print '\nRank Sums test (NONparametric, n>20):  ', around(z,4),around(p,4)
                else:
                    u,p = mannwhitneyu(x,y)
                    print '\nMann-Whitney U-test (NONparametric, ns<20):  ', around(u,4),around(p,4)

        else:  # RELATED SAMPLES
            if vartype[0]=='e':
                t,p = ttest_rel(x,y,0)
                print '\nRelated samples t-test:  ', around(t,4),around(p,4)
            else:
                t,p = ranksums(x,y)
                print '\nWilcoxon T-test (NONparametric):  ', around(t,4),around(p,4)
    else:  # CORRELATION ANALYSIS
        corrtype = ''
        while corrtype not in ['c','C','r','R','d','D']:
            print '\nIs the data Continuous, Ranked, or Dichotomous (c,r,d): ',
            corrtype = raw_input()
        if corrtype in ['c','C']:
            m,b,r,p,see = linregress(x,y)
            print '\nLinear regression for continuous variables ...'
            lol = [['Slope','Intercept','r','Prob','SEestimate'],[around(m,4),around(b,4),around(r,4),around(p,4),around(see,4)]]
            _support.printcc(lol)
        elif corrtype in ['r','R']:
            r,p = spearmanr(x,y)
            print '\nCorrelation for ranked variables ...'
            print "Spearman's r: ",around(r,4),around(p,4)
        else: # DICHOTOMOUS
            r,p = pointbiserialr(x,y)
            print '\nAssuming x contains a dichotomous variable ...'
            print 'Point Biserial r: ',around(r,4),around(p,4)
    print '\n\n'
    return None


def pearsonr(x,y):
    """
Calculates a Pearson correlation coefficient and returns p.  Taken
from Heiman's Basic Statistics for the Behav. Sci (2nd), p.195.

Returns: Pearson's r, two-tailed p-value
"""
    # x and y should have same length.
    x,y = map(asarray, (x,y))
    TINY = 1.0e-20
    n = len(x)
    mx,my = mean(x), mean(y)
    xm,ym = x-mx, y-my
    r_num = n*(add.reduce(xm*ym))
    r_den = n*math.sqrt(ss(xm)*ss(ym))
    r = (r_num / r_den)
    if (r > 1.0): r = 1.0  # numerical error caused this
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    return r,prob


def spearmanr(x,y):
    """
Calculates a Spearman rank-order correlation coefficient.  Taken
from Heiman's Basic Statistics for the Behav. Sci (1st), p.192.

Returns: Spearman's r, two-tailed p-value
"""
    n = len(x)
    assert (n>2), "Length of array must be > 2."
    rankx = rankdata(x)
    ranky = rankdata(y)
    dsq = add.reduce((rankx-ranky)**2)
    rs = 1 - 6*dsq / float(n*(n**2-1))
    df = n-2
    try:
        t = rs * math.sqrt((n-2) / ((rs+1.0)*(1.0-rs)))
        probrs = betai(0.5*df,0.5,df/(df+t*t))
    except ZeroDivisionError:
        probrs = 0.0

# probability values for rs are from part 2 of the spearman function in
# Numerical Recipies, p.510.  They close to tables, but not exact.(?)
    return rs, probrs


def pointbiserialr(x,y):
    """
Calculates a point-biserial correlation coefficient and the associated
probability value.  Taken from Heiman's Basic Statistics for the Behav.
Sci (1st), p.194.

Returns: Point-biserial r, two-tailed p-value
"""
    x,y = map(asarray, (x,y))
    TINY = 1e-30
    categories = _support.unique(x)
    data = _support.abut(x,y)
    if len(categories) <> 2:
        raise ValueError, "Exactly 2 categories required (in x) for pointbiserialr()."
    else:   # there are 2 categories, continue
        #codemap = _support.abut(categories,arange(2))
        #recoded = _support.recode(data,codemap,0)
        x = _support.linexand(data,0,categories[0])
        y = _support.linexand(data,0,categories[1])
        xmean = mean(_support.colex(x,1),None)
        ymean = mean(_support.colex(y,1),None)
        n = len(data)
        adjust = math.sqrt((len(x)/float(n))*(len(y)/float(n)))
        rpb = (ymean - xmean)/samplestd(_support.colex(data,1))*adjust
        df = n-2
        t = rpb*math.sqrt(df/((1.0-rpb+TINY)*(1.0+rpb+TINY)))
        prob = betai(0.5*df,0.5,df/(df+t*t))
        return rpb, prob


def kendalltau(x,y):
    """
Calculates Kendall's tau ... correlation of ordinal data.  Adapted
from function kendl1 in Numerical Recipies.  Needs good test-cases.@@@

Returns: Kendall's tau, two-tailed p-value
"""
    n1 = 0
    n2 = 0
    iss = 0
    for j in range(len(x)-1):
        for k in range(j+1,len(y)):
            a1 = x[j] - x[k]
            a2 = y[j] - y[k]
            aa = a1 * a2
            if (aa):             # neither array has a tie
                n1 = n1 + 1
                n2 = n2 + 1
                if aa > 0:
                    iss = iss + 1
                else:
                    iss = iss -1
            else:
                if a1:
                    n1 = n1 + 1
                if a2:
                    n2 = n2 + 1
    tau = iss / math.sqrt(float(n1*n2))
    svar = (4.0*len(x)+10.0) / (9.0*len(x)*(len(x)-1))
    z = tau / math.sqrt(svar)
    prob = erfc(abs(z)/1.4142136)
    return tau, prob


def linregress(*args):
    """
Calculates a regression line on two arrays, x and y, corresponding to x,y
pairs.  If a single 2D array is passed, alinregress finds dim with 2 levels
and splits data into x,y pairs along that dim.

Returns: slope, intercept, r, two-tailed prob, stderr-of-the-estimate
"""
    TINY = 1.0e-20
    if len(args) == 1:  # more than 1D array?
        args = asarray(args[0])
        if len(args) == 2:
            x = args[0]
            y = args[1]
        else:
            x = args[:,0]
            y = args[:,1]
    else:
        x = asarray(args[0])
        y = asarray(args[1])
    n = len(x)
    xmean = mean(x,None)
    ymean = mean(y,None)
    xm,ym = x-xmean, y-ymean
    r_num = add.reduce(xm*ym)
    r_den = math.sqrt(ss(xm)*ss(ym))
    if r_den == 0.0:
        r = 0.0
    else:
        r = r_num / r_den
        if (r > 1.0): r = 1.0 # from numerical error
    #z = 0.5*math.log((1.0+r+TINY)/(1.0-r+TINY))
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    slope = r_num / ss(xm)
    intercept = ymean - slope*xmean
    sterrest = math.sqrt(1-r*r)*samplestd(y)
    return slope, intercept, r, prob, sterrest


#####################################
#####  AINFERENTIAL STATISTICS  #####
#####################################

def ttest_1samp(a,popmean,printit=0,name='Sample',writemode='a'):
    """
Calculates the t-obtained for the independent samples T-test on ONE group
of scores a, given a population mean.  If printit=1, results are printed
to the screen.  If printit='filename', the results are output to 'filename'
using the given writemode (default=append).  Returns t-value, and prob.

Returns: t-value, two-tailed prob
"""
    a = asarray(a)
    x = mean(a,None)
    v = var(a)
    n = len(a)
    df = n-1
    svar = ((n-1)*v) / float(df)
    t = (x-popmean)/math.sqrt(svar*(1.0/n))
    prob = betai(0.5*df,0.5,df/(df+t*t))

    if printit <> 0:
        statname = 'Single-sample T-test.'
        outputpairedstats(printit,writemode,
                          'Population','--',popmean,0,0,0,
                          name,n,x,v,minimum.reduce(ravel(a)),
                          maximum.reduce(ravel(a)),
                          statname,t,prob)
    return t,prob


def ttest_ind (a, b, axis=-1, printit=0, name1='Samp1', name2='Samp2',writemode='a'):
    """
Calculates the t-obtained T-test on TWO INDEPENDENT samples of scores
a, and b.  From Numerical Recipies, p.483.  If printit=1, results are
printed to the screen.  If printit='filename', the results are output
to 'filename' using the given writemode (default=append).  Axis
can equal None (ravel array first), or an integer (the axis over
which to operate on a and b).

Returns: t-value, two-tailed p-value
"""
    a, b, axis = _chk2_asarray(a, b, axis)
    x1 = mean(a,axis)
    x2 = mean(b,axis)
    v1 = var(a,axis)
    v2 = var(b,axis)
    n1 = a.shape[axis]
    n2 = b.shape[axis]
    df = n1+n2-2
    svar = ((n1-1)*v1+(n2-1)*v2) / float(df)
    zerodivproblem = equal(svar,0)
    t = (x1-x2)/sqrt(svar*(1.0/n1 + 1.0/n2))  # N-D COMPUTATION HERE!!!!!!
    t = where(zerodivproblem,1.0,t)           # replace NaN t-values with 1.0
    probs = betai(0.5*df,0.5,float(df)/(df+t*t))

    if type(t) == ArrayType:
        probs = reshape(probs,t.shape)
    if len(probs) == 1:
        probs = probs[0]
        
    if printit <> 0:
        if type(t) == ArrayType:
            t = t[0]
        if type(probs) == ArrayType:
            probs = probs[0]
        statname = 'Independent samples T-test.'
        outputpairedstats(printit,writemode,
                          name1,n1,x1,v1,minimum.reduce(ravel(a)),
                          maximum.reduce(ravel(a)),
                          name2,n2,x2,v2,minimum.reduce(ravel(b)),
                          maximum.reduce(ravel(b)),
                          statname,t,probs)
        return
    return t, probs


def ttest_rel (a,b,axis=None,printit=0,name1='Samp1',name2='Samp2',writemode='a'):
    """
Calculates the t-obtained T-test on TWO RELATED samples of scores, a
and b.  From Numerical Recipies, p.483.  If printit=1, results are
printed to the screen.  If printit='filename', the results are output
to 'filename' using the given writemode (default=append).  Axis
can equal None (ravel array first), or an integer (the axis over
which to operate on a and b).

Returns: t-value, two-tailed p-value
"""
    a, b, axis = _chk2_asarray(a, b, axis)
    if len(a)<>len(b):
        raise ValueError, 'Unequal length arrays.'
    x1 = mean(a,axis)
    x2 = mean(b,axis)
    v1 = var(a,axis)
    v2 = var(b,axis)
    n = a.shape[axis]
    df = float(n-1)
    d = (a-b).astype('d')

    denom = sqrt((n*add.reduce(d*d,axis) - add.reduce(d,axis)**2) /df)
    zerodivproblem = equal(denom,0)
    t = add.reduce(d,axis) / denom      # N-D COMPUTATION HERE!!!!!!
    t = where(zerodivproblem,1.0,t)          # replace NaN t-values with 1.0
    t = where(zerodivproblem,1.0,t)           # replace NaN t-values with 1.0
    probs = betai(0.5*df,0.5,float(df)/(df+t*t))
    if type(t) == ArrayType:
        probs = reshape(probs,t.shape)
    if len(probs) == 1:
        probs = probs[0]

    if printit <> 0:
        statname = 'Related samples T-test.'
        outputpairedstats(printit,writemode,
                          name1,n,x1,v1,minimum.reduce(ravel(a)),
                          maximum.reduce(ravel(a)),
                          name2,n,x2,v2,minimum.reduce(ravel(b)),
                          maximum.reduce(ravel(b)),
                          statname,t,probs)
        return
    return t, probs


import scipy.stats
import distributions
def kstest(rvs,cdf,args=(),N=20):
    """Return the D-value and the p-value for a
    Kolmogorov-Smirnov test of the null that N RV's generated by
    the rvs fits the cdf given the extra arguments.  rvs
    needs to accept the size= keyword if a function.  rvs can also
    be a vector of RVs.

    cdf can be a function or a string indicating the distriubtion type.

    if the p-value is greater than the significance level (say 5%), then we
    cannot reject the hypothesis that the data come from the given
    distribution.
    """
    if type(rvs) is StringType:
        cdf = eval("scipy.stats."+rvs+".cdf")    
        rvs = eval("scipy.stats."+rvs+".rvs")
    if type(cdf) is StringType:
        cdf = eval("scipy.stats."+cdf+".cdf")        
    if callable(rvs):
        kwds = {'size':N}
        vals = sb.sort(rvs(*args,**kwds))
    else:
        vals = sb.sort(rvs)
        N = len(vals)
    cdfvals = cdf(vals, *args)
    D1 = sb.amax(abs(cdfvals - sb.arange(1.0,N+1)/N))
#    D2 = sb.amax(abs(cdfvals - sb.arange(0.0,N)/N))
#    D = max(D1,D2)
    D = D1
    return D, distributions.ksone.sf(D,N)

def chisquare(f_obs,f_exp=None):
    """
Calculates a one-way chi square for array of observed frequencies and returns
the result.  If no expected frequencies are given, the total N is assumed to
be equally distributed across all groups.

Returns: chisquare-statistic, associated p-value
"""

    f_obs = asarray(f_obs)
    k = len(f_obs)
    if f_exp is None:
        f_exp = array([sum(f_obs)/float(k)] * len(f_obs),Float)
    f_exp = f_exp.astype(Float)
    chisq = add.reduce((f_obs-f_exp)**2 / f_exp)
    return chisq, chisqprob(chisq, k-1)


def ks_2samp (data1,data2):
    """
Computes the Kolmogorov-Smirnof statistic on 2 samples.  Modified from
Numerical Recipies in C, page 493.  Returns KS D-value, prob.  Not ufunc-
like.

Returns: KS D-value, p-value
"""
    data1, data2 = map(asarray, (data1, data2))
    j1 = 0    # zeros(data1.shape[1:]) TRIED TO MAKE THIS UFUNC-LIKE
    j2 = 0    # zeros(data2.shape[1:])
    fn1 = 0.0 # zeros(data1.shape[1:],Float)
    fn2 = 0.0 # zeros(data2.shape[1:],Float)
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    en1 = n1*1
    en2 = n2*1
    d = zeros(data1.shape[1:],Float)
    data1 = sort(data1,0)
    data2 = sort(data2,0)
    while j1 < n1 and j2 < n2:
        d1=data1[j1]
        d2=data2[j2]
        if d1 <= d2:
            fn1 = (j1)/float(en1)
            j1 = j1 + 1
        if d2 <= d1:
            fn2 = (j2)/float(en2)
            j2 = j2 + 1
        dt = (fn2-fn1)
        if abs(dt) > abs(d):
            d = dt
    try:
        en = math.sqrt(en1*en2/float(en1+en2))
        prob = ksprob((en+0.12+0.11/en)*fabs(d))
    except:
        prob = 1.0
    return d, prob


def mannwhitneyu(x,y):
    """
Calculates a Mann-Whitney U statistic on the provided scores and
returns the result.  Use only when the n in each condition is < 20 and
you have 2 independent samples of ranks.  REMEMBER: Mann-Whitney U is
significant if the u-obtained is LESS THAN or equal to the critical
value of U.

Returns: u-statistic, one-tailed p-value (i.e., p(z(U)))
"""
    x,y = asarray(x, y)
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(concatenate((x,y)))
    rankx = ranked[0:n1]       # get the x-ranks
    #ranky = ranked[n1:]        # the rest are y-ranks
    u1 = n1*n2 + (n1*(n1+1))/2.0 - sum(rankx)  # calc U for x
    u2 = n1*n2 - u1                            # remainder is U for y
    bigu = max(u1,u2)
    smallu = min(u1,u2)
    T = math.sqrt(tiecorrect(ranked))  # correction factor for tied scores
    if T == 0:
        raise ValueError, 'All numbers are identical in amannwhitneyu'
    sd = math.sqrt(T*n1*n2*(n1+n2+1)/12.0)
    z = abs((bigu-n1*n2/2.0) / sd)  # normal approximation for prob calc
    return smallu, 1.0 - zprob(z)


def tiecorrect(rankvals):
    """
Tie-corrector for ties in Mann Whitney U and Kruskal Wallis H tests.
See Siegel, S. (1956) Nonparametric Statistics for the Behavioral
Sciences.  New York: McGraw-Hill.  Code adapted from |Stat rankind.c
code.

Returns: T correction factor for U or H
"""
    sorted,posn = fastsort(asarray(rankvals))
    n = len(sorted)
    T = 0.0
    i = 0
    while (i<n-1):
        if sorted[i] == sorted[i+1]:
            nties = 1
            while (i<n-1) and (sorted[i] == sorted[i+1]):
                nties = nties +1
                i = i +1
            T = T + nties**3 - nties
        i = i+1
    T = T / float(n**3-n)
    return 1.0 - T


def ranksums(x,y):
    """
Calculates the rank sums statistic on the provided scores and returns
the result.

Returns: z-statistic, two-tailed p-value
"""
    x,y = map(asarray, (x, y))
    n1 = len(x)
    n2 = len(y)
    alldata = concatenate((x,y))
    ranked = rankdata(alldata)
    x = ranked[:n1]
    y = ranked[n1:]
    s = sum(x)
    expected = n1*(n1+n2+1) / 2.0
    z = (s - expected) / math.sqrt(n1*n2*(n1+n2+1)/12.0)
    prob = 2*(1.0 -zprob(abs(z)))
    return z, prob


    
def kruskal(*args):
    """
The Kruskal-Wallis H-test is a non-parametric ANOVA for 2 or more
groups, requiring at least 5 subjects in each group.  This function
calculates the Kruskal-Wallis H and associated p-value for 2 or more
independent samples.

Returns: H-statistic (corrected for ties), associated p-value
"""
    assert len(args) >= 2, "Need at least 2 groups in stats.kruskal()"
    n = map(len,args)
    all = []
    for i in range(len(args)):
        all.extend(args[i].tolist())
    ranked = list(rankdata(all))
    T = tiecorrect(ranked)
    for i in range(len(args)):
        args[i] = ranked[0:n[i]]
        del ranked[0:n[i]]
    rsums = []
    for i in range(len(args)):
        rsums.append(sum(args[i])**2)
        rsums[i] = rsums[i] / float(n[i])
    ssbn = sum(rsums)
    totaln = sum(n)
    h = 12.0 / (totaln*(totaln+1)) * ssbn - 3*(totaln+1)
    df = len(args) - 1
    if T == 0:
        raise ValueError, 'All numbers are identical in kruskal'
    h = h / float(T)
    return h, chisqprob(h,df)


def friedmanchisquare(*args):
    """
Friedman Chi-Square is a non-parametric, one-way within-subjects
ANOVA.  This function calculates the Friedman Chi-square test for
repeated measures and returns the result, along with the associated
probability value.  It assumes 3 or more repeated measures.  Only 3
levels requires a minimum of 10 subjects in the study.  Four levels
requires 5 subjects per level(??).

Returns: chi-square statistic, associated p-value
"""
    k = len(args)
    if k < 3:
        raise ValueError, '\nLess than 3 levels.  Friedman test not appropriate.\n'
    n = len(args[0])
    data = apply(_support.abut,args)
    data = data.astype(Float)
    for i in range(len(data)):
        data[i] = rankdata(data[i])
    ssbn = sum(sum(args,1)**2)
    chisq = 12.0 / (k*n*(k+1)) * ssbn - 3*n*(k+1)
    return chisq, chisqprob(chisq,k-1)


#####################################
####  PROBABILITY CALCULATIONS  ####
#####################################

zprob = special.ndtr
erfc = special.erfc

def chisqprob(chisq,df):
    """Returns the (1-tail) probability value associated with the provided chi-square
    value and df.
    """
    return special.chdtrc(df,chisq)

ksprob = special.kolmogorov
fprob = special.fdtrc

def betai(a,b,x):
    """
    Returns the incomplete beta function:
    
    I-sub-x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)

    where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
    function of a.
    """
    x = asarray(x)
    x = where(x < 1.0, x, 1.0)  # if x > 1 then return 1.0
    return special.betainc(a,b,x)

#####################################
#######  AANOVA CALCULATIONS  #######
#####################################

def glm(data,para):
    """
Calculates a linear model fit ... anova/ancova/lin-regress/t-test/etc. Taken
from:
    Peterson et al. Statistical limitations in functional neuroimaging
    I. Non-inferential methods and statistical models.  Phil Trans Royal Soc
    Lond B 354: 1239-1260.

Returns: statistic, p-value ???
"""
    if len(para) <> len(data):
        print "data and para must be same length in aglm"
        return
    n = len(para)
    p = _support.unique(para)
    x = zeros((n,len(p)))  # design matrix
    for l in range(len(p)):
        x[:,l] = equal(para,p[l])
    b = dot(dot(linalg.inv(dot(transpose(x),x)),  # i.e., b=inv(X'X)X'Y
                    transpose(x)),data)
    diffs = (data - dot(x,b))
    s_sq = 1./(n-len(p)) * dot(transpose(diffs), diffs)

    if len(p) == 2:  # ttest_ind
        c = array([1,-1])
        df = n-2
        fact = sum(1.0/sum(x,0))  # i.e., 1/n1 + 1/n2 + 1/n3 ...
        t = dot(c,b) / sqrt(s_sq*fact)
        probs = betai(0.5*df,0.5,float(df)/(df+t*t))
        return t, probs

def anova(data,effects=['A','B','C','D','E','F','G','H','I','J','K']):
    """
Prints the results of single-variable between- and within-subject ANOVA
designs.  The function can only handle univariate ANOVAs with a single
random factor.  The random factor is coded in column 0 of the input
list/array (see below) and the measured variable is coded in the last
column of the input list/array. The following were used as references
when writing the code:

Maxwell, SE, Delaney HD (1990)  Designing Experiments and Analyzing
    Data, Wadsworth: Belmont, CA.
Lindman, HR (1992) Analysis of Variance in Experimental Design,
    Springer-Verlag: New York.

TO DO:  Increase Current Max Of 10 Levels Per W/I-Subject Factor
        Consolidate Between-Subj Analyses For Between And Within/Between
        Front-end for different input data-array shapes/organization
        Axe mess of 'global' statements (particularly for d_restrict fcns)

Usage:   anova(data,                         data = |Stat format
               effects=['A','B','C','D','E','F','G','H','I','J','K'])

Note: |Stat format is as follows ... one datum per row, first element of
row is the subject identifier, followed by all within/between subject
variable designators, and the measured data point as the last element in the
row.  Thus, [1, 'short', 'drugY', 2, 14.7] represents subject 1 when measured
in the short / drugY / 2 condition, and subject 1 gave a measured value of
14.7 in this combination of conditions.  Thus, all input lists are '2D'
lists-of-lists.
"""
    global alluniqueslist, Nlevels, Nfactors, Nsubjects, Nblevels, Nallsources
    global Bscols, Bbetweens, SSlist, SSsources, Bwonly_sources, D
    global alleffects, alleffsources
    outputlist = []
    #SSbtw = []
    #SSbtwsources = []
    #SSwb = []
    #SSwbsources = []
    alleffects = []
    alleffsources = []
    SSlist = []
    SSsources = []

    print
    variables = 1       # this function only handles one measured variable
    if type(data)==ArrayType:
        data = data.tolist()

## Create a list of all unique values in each column, and a list of these Ns
    alluniqueslist = [0]*(len(data[0])-variables) # all cols but data cols
    Nlevels = [0]*(len(data[0])-variables)        # (as above)
    for column in range(len(Nlevels)): 
        alluniqueslist[column] = _support.unique(_support.colex(data,column))
        Nlevels[column] = len(alluniqueslist[column])

    #Ncells = multiply.reduce(Nlevels[1:]) # total num cells (w/i AND btw)
    Nfactors = len(Nlevels[1:])             # total num factors
    Nallsources = 2**(Nfactors+1)  # total no. possible sources (factor-combos)
    Nsubjects = len(alluniqueslist[0])  # total # subj in study (# of diff. subj numbers in column 0)

## Within-subj factors defined as those where there are fewer subj than
## scores in the first level of a factor (quick and dirty; findwithin() below)
    Bwithins = findwithin(data)         # binary w/i subj factors (excl. col 0)
    Bbetweens = ~Bwithins & (Nallsources-1) - 1

    Wcolumns = makelist(Bwithins,Nfactors+1)  # get list of cols of w/i factors
    Wscols = [0] + Wcolumns                   # w/i subj columns INCL col 0
    Bscols = makelist(Bbetweens+1,Nfactors+1) #list of btw-subj cols,INCL col 0
    Nwifactors = len(Wscols) - 1 # WAS len(Wcolumns)
    #Nwlevels = take(array(Nlevels),Wscols) # no.lvls for each w/i subj fact
    #Nbtwfactors = len(Bscols) - 1 # WASNfactors - Nwifactors + 1
    Nblevels = take(array(Nlevels),Bscols)

    Nwsources = 2**Nwifactors - 1 # num within-subject factor-combos
    #Nbsources = Nallsources - Nwsources

    #
    # CALC M-VARIABLE (LIST) and Marray/Narray VARIABLES (ARRAY OF CELL MNS/NS)
    #
    # Eliminate replications for the same subject in same condition as well as
    # within-subject repetitions, keep as list
    M = _support.collapse(data,Bscols,-1,0,0)
    # Create an arrays of Nblevels shape (excl. subj dim)
    Marray = zeros(Nblevels[1:],'f')
    Narray = zeros(Nblevels[1:],'f')
    # Fill arrays by looping through all scores in the (collapsed) M
    for row in M:
        idx = []
        for i in range(len(row[:-1])):
            idx.append(alluniqueslist[Bscols[i]].index(row[i]))
        idx = idx[1:]
        Marray[idx] = Marray[idx] + row[-1]
        Narray[idx] = Narray[idx] + 1
    Marray = Marray / Narray

    #
    # CREATE DATA ARRAY, DA, FROM ORIGINAL INPUT DATA
    # (this is an unbelievably bad, wasteful data structure, but it makes lots
    # of tasks much easier; should nevertheless be fixed someday)

    # This limits the within-subject level count to 10!
    coefflist =[[[1]],
                [[-1,1]],
                [[-1,0,1],[1,-2,1]],
                [[-3,-1,1,3],[1,-1,-1,1],[-1,3,-3,1]],
                [[-2,-1,0,1,2],[2,-1,-2,-1,2],[-1,2,0,-2,1],[1,-4,6,-4,1]],
                [[-5,-3,-1,1,3,5],[5,-1,-4,-4,-1,5],[-5,7,4,-4,-7,5],
                 [1,-3,2,2,-3,1],[-1,5,-10,10,-5,1]],
                [[-3,-2,-1,0,1,2,3],[5,0,-3,-4,-3,0,5],[-1,1,1,0,-1,-1,1],
                 [3,-7,1,6,1,-7,3],[-1,4,-5,0,5,-4,1],[1,-6,15,-20,15,-6,1]],
                [[-7,-5,-3,-1,1,3,5,7],[7,1,-3,-5,-5,-3,1,7],
                 [-7,5,7,3,-3,-7,-5,7],[7,-13,-3,9,9,-3,-13,7],
                 [-7,23,-17,-15,15,17,-23,7],[1,-5,9,-5,-5,9,-5,1],
                 [-1,7,-21,35,-35,21,-7,1]],
                [[-4,-3,-2,-1,0,1,2,3,4],[28,7,-8,-17,-20,-17,-8,7,28],
                 [-14,7,13,9,0,-9,-13,-7,14],[14,-21,-11,9,18,9,-11,-21,14],
                 [-4,11,-4,-9,0,9,4,-11,4],[4,-17,22,1,-20,1,22,-17,4],
                 [-1,6,-14,14,0,-14,14,-6,1],[1,-8,28,-56,70,-56,28,-8,1]],
                [[-9,-7,-5,-3,-1,1,3,5,7,9],[6,2,-1,-3,-4,-4,-3,-1,2,6],
                 [-42,14,35,31,12,-12,-31,-35,-14,42],
                 [18,-22,-17,3,18,18,3,-17,-22,18],
                 [-6,14,-1,-11,-6,6,11,1,-14,6],[3,-11,10,6,-8,-8,6,10,-11,3],
                 [9,-47,86,-42,-56,56,42,-86,47,-9],
                 [1,-7,20,-28,14,14,-28,20,-7,1],
                 [-1,9,-36,84,-126,126,-84,36,-9,1]]]

    dindex = 0
    # Prepare a list to be filled with arrays of D-variables, array per within-
    # subject combo (i.e., for 2 w/i subj factors E and F ... E, F, ExF)
    NDs = [0]* Nwsources
    for source in range(Nwsources):
        if subset(source,Bwithins):
            NDs[dindex] = numlevels(source,Nlevels)
            dindex = dindex + 1

    # Collapse multiple repetitions on the same subject and same condition
    cdata = _support.collapse(data,range(Nfactors+1),-1,0,0)

    # Find a value that's not a data score with which to fill the array DA
    dummyval = -1
    datavals = _support.colex(data,-1)
    while dummyval in datavals:  # find a value that's not a data score
        dummyval = dummyval - 1
    DA = ones(Nlevels,'f')*dummyval # create plenty of data-slots to fill

    if len(Bscols) == 1: # ie., if no btw-subj factors
        # 1 (below) needed because we need 2D array even w/ only 1 group of subjects
        subjslots = ones((Nsubjects,1)) 
    else: # create array to hold 1s (subj present) and 0s (subj absent)
        subjslots = zeros(Nblevels)
    for i in range(len(data)): # for every datapoint given as input
        idx = []
        for j in range(Nfactors+1): # get n-D bin idx for this datapoint
            new = alluniqueslist[j].index(data[i][j])
            idx.append(new)
        DA[idx] = data[i][-1] # put this data point in proper place in DA
        btwidx = take(idx,array(Bscols))
        subjslots[btwidx] = 1
    # DONE CREATING DATA ARRAY, DA ... #dims = numfactors+1, dim 0=subjects
    # dim -1=measured values, dummyval = values used to fill empty slots in DA

    # PREPARE FOR MAIN LOOP
    dcount = -1     # prepare for pre-increment of D-variable pointer
    Bwsources = []  # binary #s, each=source containing w/i subj factors
    Bwonly_sources = [] # binary #s, each=source of w/i-subj-ONLY factors
    D = zeros(Nwsources,PyObject) # one slot for each Dx,2**Nwifactors
    DM = [0] *Nwsources # Holds arrays of cell-means
    DN = [0] *Nwsources # Holds arrays of cell-ns

    # BEGIN MAIN LOOP!!!!!
    # BEGIN MAIN LOOP!!!!!
    # BEGIN MAIN LOOP!!!!!
    for source in range(3,Nallsources,2): # all sources that incl. subjects
        if ((source-1) & Bwithins) != 0: # 1 or more w/i subj sources?
            Bwsources.append(source-1)   # add it to a list
        #
        # WITHIN-SUBJECT-ONLY TERM?  IF SO ... NEED TO CALCULATE NEW D-VARIABLE
        # (per Maxwell & Delaney pp.622-4)
        if subset((source-1),Bwithins):
            # Keep track of which D-var set we're working with (De, Df, Def, etc.)
            dcount = dcount + 1
            Bwonly_sources.append(source-1) #add source, minus subj,to list
            dwsc = 1.0 * DA       # get COPY of w/i-subj data array
            # Find all non-source columns, note ~source alone (below) -> negative number
            Bnonsource = (Nallsources-1) & ~source
            Bwscols = makebin(Wscols) # make a binary version of Wscols
            # Figure out which cols from the ORIGINAL (input) data matrix are both non-
            # source and also within-subj vars (excluding subjects col)
            Bwithinnonsource = Bnonsource & Bwscols 

            # Next, make a list of the above.  The list is a list of axes in DA
            # because DA has the same number of axes as there are factors
            # (including subjects), but with extra dummyval='-1' values the original
            # data array (assuming between-subj vars exist)
            Lwithinnonsource = makelist(Bwithinnonsource,Nfactors+1)

            # Collapse all non-source, w/i subj dims, FROM THE END (otherwise the
            # dim-numbers change as you collapse).  THIS WORKS BECAUSE WE'RE
            # COLLAPSING ACROSS W/I SUBJECT AXES, WHICH WILL ALL HAVE THE
            # SAME SUBJ IN THE SAME ARRAY LOCATIONS (i.e., dummyvals will still exist
            # but should remain the same value through the mean() function
            for i in range(len(Lwithinnonsource)-1,-1,-1):
                dwsc = mean(dwsc,Lwithinnonsource[i])
            mns = dwsc

            # NOW, ACTUALLY COMPUTE THE D-VARIABLE ENTRIES FROM DA
            # CREATE LIST OF COEFF-COMBINATIONS TO DO (len=e-1, f-1, (e-1)*(f-1), etc...)
            #
            # Figure out which cols are both source and within-subjects, including col 0
            Bwithinsource = source & Bwscols 
            # Make a list of within-subj cols, incl subjects col (0)
            Lwithinsourcecol = makelist(Bwithinsource, Nfactors+1)
            # Make a list of cols that are source within-subj OR btw-subj 
            Lsourceandbtws = makelist(source | Bbetweens, Nfactors+1)
            if Lwithinnonsource <> []:
                Lwithinsourcecol = map(Lsourceandbtws.index,Lwithinsourcecol)
                # Now indxlist should hold a list of indices into the list of possible
                # coefficients, one row per combo of coefficient. Next line PRESERVES dummyval
            dvarshape = array(take(mns.shape,Lwithinsourcecol[1:])) -1
            idxarray = indices(dvarshape)
            newshape = array([idxarray.shape[0],
                                multiply.reduce(idxarray.shape[1:])])
            indxlist = swapaxes(reshape(idxarray,newshape),0,1)

            # The following is what makes the D-vars 2D.  It takes an n-dim array
            # and retains the first (num of factors) dim while making the 2nd dim
            # equal to the total number of source within-subject cells.

            #
            # CREATE ALL D-VARIABLES FOR THIS COMBINATION OF FACTORS
            #
            for i in range(len(indxlist)):
                #
                # FILL UP COEFFMATRIX (OF SHAPE = MNS) WITH CORRECT COEFFS FOR 1 D-VAR
                #
                coeffmatrix = ones(mns.shape,Float) # fewer dims than DA (!!)
                # Make a list of dim #s that are both in source AND w/i subj fact, incl subj
                #Wsourcecol = makelist(Bwscols&source,Nfactors+1)
                # Fill coeffmatrix with a complete set of coeffs (1 per w/i-source factor)
                for wfactor in range(len(Lwithinsourcecol[1:])):
                    #put correct coeff. axis as first axis, or "swap it up"
                    coeffmatrix = swapaxes(coeffmatrix,0,
                                             Lwithinsourcecol[wfactor+1])
                    # Find appropriate ROW of (static) coefflist we need
                    nlevels = coeffmatrix.shape[0]
                    # Get the next coeff in that row
                    try:
                        nextcoeff = coefflist[nlevels-1][indxlist[i,wfactor]]
                    except IndexError:
                        raise IndexError, "anova() can only handle up to 10 levels on a within-subject factors"
                    for j in range(nlevels):
                        coeffmatrix[j] = coeffmatrix[j] * nextcoeff[j]
                    # Swap it back to where it came from
                    coeffmatrix = swapaxes(coeffmatrix,0,
                                             Lwithinsourcecol[wfactor+1])

                # CALCULATE D VARIABLE
                scratch = coeffmatrix * mns
                # Collapse all axes EXCEPT subjects dim (dim 0)
                for j in range(len(coeffmatrix.shape[1:])):
                    scratch = add.reduce(scratch,1)
                if len(scratch.shape) == 1:
                    scratch.shape = list(scratch.shape)+[1]
                try:
                    # Tack this column onto existing ones
                    #tmp = D[dcount].shape
                    D[dcount] = _support.abut(D[dcount],scratch)
                except AttributeError: # i.e., D[dcount]=integer/float
                    # If this is the first, plug it in
                    D[dcount] = scratch


            # Big long thing to create DMarray (list of DM variables) for this source
            variables = D[dcount].shape[1] # Num variables for this source
            tidx = range(1,len(subjslots.shape)) + [0] # [0] = Ss dim
            tsubjslots = transpose(subjslots,tidx) # put Ss in last dim
            DMarray = zeros(list(tsubjslots.shape[0:-1]) +
                              [variables],'f') # btw-subj dims, then vars
            DNarray = zeros(list(tsubjslots.shape[0:-1]) +
                              [variables],'f') # btw-subj dims, then vars
            idx = [0] *len(tsubjslots.shape[0:-1])
            idx[0] = -1
            loopcap = array(tsubjslots.shape[0:-1]) -1
            while incr(idx,loopcap) <> -1:
                DNarray[idx] = float(sum(tsubjslots[idx]))
                thismean =  (add.reduce(tsubjslots[idx] * # 1=subj dim
                                          transpose(D[dcount]),1) /
                             DNarray[idx])
                thismean = array(thismean,PyObject)
                DMarray[idx] = thismean
            DM[dcount] = DMarray
            DN[dcount] = DNarray

        #
        # DONE CREATING M AND D VARIABLES ... TIME FOR SOME SS WORK
        # DONE CREATING M AND D VARIABLES ... TIME FOR SOME SS WORK
        #
        if Bscols[1:] <> []:
            BNs = _support.colex([Nlevels],Bscols[1:])
        else:
            BNs = [1]
            #
            # FIGURE OUT WHICH VARS TO RESTRICT, see p.680 (Maxwell&Delaney)
            #
            # BETWEEN-SUBJECTS VARIABLES ONLY, use M variable for analysis
            #
        if ((source-1) & Bwithins) == 0:  # btw-subjects vars only?
            #sourcecols = makelist(source-1,Nfactors+1)

            # Determine cols (from input list) required for n-way interaction
            Lsource = makelist((Nallsources-1)&Bbetweens,Nfactors+1)
            # NOW convert this list of between-subject column numbers to a list of
            # AXES in M, since M has fewer dims than the original data array
            # (assuming within-subj vars exist); Bscols has list of between-subj cols
            # from input list, the indices of which correspond to that var's loc'n in M
            btwcols = map(Bscols.index,Lsource)
            # Obviously-needed loop to get cell means is embedded in the collapse fcn, -1
            # represents last (measured-variable) column, None=std, 1=retain Ns

            #hn = ahmean(Narray,-1) # -1=unravel first

            # CALCULATE SSw ... SUBTRACT APPROPRIATE CELL MEAN FROM EACH SUBJ SCORE
            SSw = 0.0
            #idxlist = _support.unique(_support.colex(M,btwcols))
            for row in M:
                idx = []
                for i in range(len(row[:-1])):
                    idx.append(alluniqueslist[Bscols[i]].index(row[i]))
                idx = idx[1:]   # Strop off Ss col/dim
                newval = row[-1] - Marray[idx]
                SSw = SSw + (newval)**2

            # Determine which cols from input are required for this source
            Lsource = makelist(source-1,Nfactors+1)
            # NOW convert this list of between-subject column numbers to a list of
            # AXES in M, since M has fewer dims than the original data array
            # (assuming within-subj vars exist); Bscols has list of between-subj cols
            # from input list, the indices of which correspond to that var's loc'n in M
            #btwsourcecols = (array(map(Bscols.index,Lsource))-1).tolist()

            # Average Marray and get harmonic means of Narray OVER NON-SOURCE DIMS
            Bbtwnonsourcedims = ~source & Bbetweens
            Lbtwnonsourcedims = makelist(Bbtwnonsourcedims,Nfactors+1)
            btwnonsourcedims = (array(map(Bscols.index,Lbtwnonsourcedims))-1).tolist()

    ## Average Marray over non-source axes (1=keep squashed dims)
            sourceMarray = apply_over_axes(mean, Marray,btwnonsourcedims)

    ## Calculate harmonic means for each level in source
            sourceNarray = apply_over_axes(hmean, Narray,btwnonsourcedims)

    ## Calc grand average (ga), used for ALL effects
            ga = sum((sourceMarray*sourceNarray)/
                            sum(sourceNarray))
            ga = reshape(ga,ones(len(Marray.shape)))

    ## If GRAND interaction, use harmonic mean of ALL cell Ns
            if source == Nallsources-1:
                sourceNarray = hmean(Narray, None)

    ## Calc all SUBSOURCES to be subtracted from sourceMarray (M&D p.320)
            sub_effects = 1.0 * ga # start with grand mean
            for subsource in range(3,source,2):
        ## Make a list of the non-subsource axes
                if subset(subsource-1,source-1):
                    sub_effects = (sub_effects +
                                   alleffects[alleffsources.index(subsource)])
        ## Calc this effect (a(j)'s, b(k)'s, ab(j,k)'s, whatever)
            effect = sourceMarray - sub_effects

        ## Save it so you don't have to calculate it again next time
            alleffects.append(effect)
            alleffsources.append(source)

    ## Calc and save sums of squares for this source
            SS = sum((effect**2 *sourceNarray) *
                      multiply.reduce(take(Marray.shape,btwnonsourcedims)))
        ## Save it so you don't have to calculate it again next time
            SSlist.append(SS)
            SSsources.append(source)

            collapsed = _support.collapse(M,btwcols,-1,0,1)
            # Obviously needed for-loop to get source cell-means embedded in collapse fcns
            #contrastmns = _support.collapse(collapsed,btwsourcecols,-2,1,1)
            # Collapse again, this time SUMMING instead of averaging (to get cell Ns)
            #contrastns = _support.collapse(collapsed,btwsourcecols,-1,0,0,
            #                            sum)
                
            # Collapse again, this time calculating hmeans (for hns)
            #contrasthns = _support.collapse(collapsed,btwsourcecols,-1,0,0,
            #                             hmean)
            # CALCULATE *BTW-SUBJ* dfnum, dfden
            sourceNs = _support.colex([Nlevels],makelist(source-1,Nfactors+1))
            dfnum = multiply.reduce(ravel(array(sourceNs)-1))
            dfden = Nsubjects - multiply.reduce(ravel(BNs))

            # CALCULATE MS, MSw, F AND PROB FOR ALL-BETWEEN-SUBJ SOURCES ONLY
            MS = SS / dfnum
            MSw = SSw / dfden
            if MSw <> 0:
                f = MS / MSw
            else:
                f = 0  # i.e., absolutely NO error in the full model

            if f >= 0:
                prob = fprob(dfnum, dfden, f)
            else:
                prob = 1.0
    # Now this falls thru to output stage

    #
    # SOME WITHIN-SUBJECTS FACTORS TO DEAL WITH ... use appropriate D variable
    #
        else:  # Source has some w/i subj factors
            # FIGURE OUT WHICH D-VAR TO USE BASED ON WHICH W/I-SUBJ FACTORS ARE IN SOURCE
            # Determine which w/i-subj factors are in this source
            sourcewithins = (source-1) & Bwithins
            # Use D-var that was created for that w/i subj combo (the position of that
            # source within Bwsources determines the index of that D-var in D)
            workd = asarray(D[Bwonly_sources.index(sourcewithins)])

            # CALCULATE Er, Ef
    ## Set up workd and subjslots for upcoming calcs
            if len(workd.shape)==1:
                workd = workd[:,NewAxis]
            if len(subjslots.shape)==1:
                subjslots = subjslots[:,NewAxis]

    ## Calculate full-model sums of squares
            ef = d_full_model(workd,subjslots) # Uses cell-means model

            #
            # **ONLY** WITHIN-SUBJECT VARIABLES TO CONSIDER
            #
            if subset((source-1),Bwithins):
                # restrict grand mean, as per M&D p.680
                er = d_restrict_mean(workd,subjslots) 
        #
        # **BOTH** WITHIN- AND BETWEEN-SUBJECTS VARIABLES TO CONSIDER
        #
            else:
                er = d_restrict_source(workd,subjslots,source) + ef
            SSw = linalg.det(ef)
            SS = linalg.det(er) - SSw

        # CALCULATE *W/I-SUBJ* dfnum, dfden
            sourceNs = _support.colex([Nlevels],makelist(source,Nfactors+1))
            # Calculation of dfnum is straightforward regardless
            dfnum = multiply.reduce(ravel(array(sourceNs)-1)[1:])
            # If only within-subject factors are involved, dfden is straightforward
            if subset(source-1,Bwithins):
                dfden = Nsubjects -multiply.reduce(ravel(BNs))-dfnum +1
                MS = SS / dfnum
                MSw = SSw / dfden
                if MSw <> 0:
                    f = MS / MSw
                else:
                    f = 0  # i.e., absolutely NO error in full model

                if f >= 0:
                    prob = fprob(dfnum, dfden, f)
                else:
                    prob = 1.0

            # If combined within-between source, must use Rao's approximation for dfden
            # from Tatsuoka, MM (1988) Multivariate Analysis (2nd Ed), MacMillan: NY p93
            else: # it's a within-between combo source
                try:
                    p = workd.shape[1]
                except IndexError:
                    p = 1
                k = multiply.reduce(ravel(BNs))
                m = Nsubjects -1 -(p+k)/2.0
                d_en = float(p**2 + (k-1)**2 - 5)
                if d_en == 0.0:
                    s = 1.0
                else:
                    s = math.sqrt(((p*(k-1))**2-4) / d_en)
                dfden = m*s - dfnum/2.0 + 1

                # Given a within-between combined source, Wilk's Lambda is appropriate
                if linalg.det(er) <> 0:
                    lmbda = linalg.det(ef) / linalg.det(er)
                    W = math.pow(lmbda,(1.0/s))
                    f = ((1.0-W)/W) * (dfden/dfnum)
                else:
                    f = 0  # i.e., absolutely NO error in restricted model

                if f >= 0:
                    prob = fprob(dfnum,dfden,f)
                else:
                    prob = 1.0

        #
        # CREATE STRING-LIST FOR RESULTS FROM THIS PARTICULAR SOURCE
        #
        suffix = ''                       # for *s after the p-value
        if  prob < 0.001:  suffix = '***'
        elif prob < 0.01:  suffix = '**'
        elif prob < 0.05:  suffix = '*'
        adjsourcecols = array(makelist(source-1,Nfactors+1)) -1
        thiseffect = ''
        for col in adjsourcecols:
            if len(adjsourcecols) > 1:
                thiseffect = thiseffect + effects[col][0]
            else:
                thiseffect = thiseffect + (effects[col])
        outputlist = (outputlist
        # These terms are for the numerator of the current effect/source
                      + [[thiseffect, round4(SS),dfnum,
                          round4(SS/float(dfnum)),round4(f),
                          round4(prob),suffix]] 
        # These terms are for the denominator for the current effect/source
                      + [[thiseffect+'/w', round4(SSw),dfden,
                          round4(SSw/float(dfden)),'','','']]
                      + [['\n']])

        #
        # PRINT OUT ALL MEANS AND Ns FOR THIS SOURCE (i.e., this combo of factors)
        #
        Lsource = makelist(source-1,Nfactors+1)
        collapsed = _support.collapse(cdata,Lsource,-1,1,1)

        # First, get the list of level-combos for source cells
        prefixcols = range(len(collapsed[0][:-3]))
        outlist = _support.colex(collapsed,prefixcols)
        # Start w/ factor names (A,B,C, or ones input to anova())
        eff = []
        for col in Lsource:
            eff.append(effects[col-1])
        # Add in the mean and N labels for printout
        for item in ['MEAN','STDERR','N']:
            eff.append(item)
        # To the list of level-combos, abut the corresp. means and Ns
        outlist = _support.abut(outlist,
                             map(round4,_support.colex(collapsed,-3)),
                             map(round4,_support.colex(collapsed,-2)),
                             map(round4,_support.colex(collapsed,-1)))
        outlist = [eff] + outlist # add titles to the top of the list
        _support.printcc(outlist)    # print it in customized columns
        print


###
### OUTPUT FINAL RESULTS (ALL SOURCES TOGETHER)
### Note: All 3 types of source-calcs fall through to here
###
    print
    title = [['FACTORS: ','RANDOM'] + effects[:Nfactors]]
    title = title + [['LEVELS:  ']+Nlevels]
    facttypes = ['BETWEEN']*Nfactors
    for i in range(len(Wscols[1:])):
        facttypes[Wscols[i+1]-1] = 'WITHIN'
    title = title + [['TYPE:    ','RANDOM']+facttypes]
    _support.printcc(title)
    print

    title = [['Effect','SS','DF','MS','F','p','sig']] + ['dashes']
    outputlist = title + outputlist
    _support.printcc(outputlist)
    return


def d_full_model(workd,subjslots):
     """
     RESTRICTS NOTHING (i.e., FULL MODEL CALCULATION).  Subtracts D-variable
cell-mean for each between-subj group and then calculates the SS array.
     """
     workd = subtr_cellmeans(workd,subjslots)
     sserr = multivar_SScalc(workd)
     return sserr


def d_restrict_mean(workd,subjslots):
     """
     RESTRICTS GRAND MEA  Subtracts D-variable cell-mean for each between-
subj group, and then adds back each D-variable's grand mean.
     """
     # subtract D-variable cell-mean for each (btw-subj) group
     errors = subtr_cellmeans(workd,subjslots)

     # add back in appropriate grand mean from individual scores
     grandDmeans = expand_dims(mean(workd,0),0)
     errors = errors + transpose(grandDmeans) # errors has reversed dims!!
     # SS for mean-restricted model is calculated below.  Note: already put
     # subj as last dim because later code expects this code here to leave
     # workd that way
     sserr = multivar_SScalc(errors)
     return sserr


def d_restrict_source(workd,subjslots,source):
     """
Calculates error for a given model on array workd.  Subjslots is an
array of 1s and 0s corresponding to whether or not the subject is a
member of that between-subjects variable combo.  source is the code
for the type of model to calculate.  source=-1 means no restriction;
source=0 means to restrict workd's grand mean; source>0 means to
restrict the columns of the main data array, DA, specified (in binary)
by the source-value.

Returns: SS array for multivariate F calculation
"""
###
### RESTRICT COLUMNS/AXES SPECIFIED IN source (BINARY)
### (i.e., is the value of source not equal to 0 or -1?)
###
     global D
     if source > 0:
         sourcewithins = (source-1) & Bwithins
         sourcebetweens = (source-1) & Bbetweens
         dindex = Bwonly_sources.index(sourcewithins)
         all_cellmeans = transpose(DM[dindex],[-1]+range(0,len(DM[dindex].shape)-1))
         all_cellns = transpose(DN[dindex],[-1]+range(0,len(DN[dindex].shape)-1))
         hn = hmean(all_cellns, None)

         levels = D[dindex].shape[1]  # GENERAL, 'cause each workd is always 2D
         SSm = zeros((levels,levels),'f') #called RCm=SCm in Lindman,p.317-8
         tworkd = transpose(D[dindex])

     ## Calculate SSw, within-subj variance (Lindman approach)
         RSw = zeros((levels,levels),'f')
         RSinter = zeros((levels,levels),PyObject)  
         for i in range(levels):
             for j in range(i,levels):
                 RSw[i,j] = RSw[j,i] = sum(tworkd[i]*tworkd[j])
                 cross = all_cellmeans[i] * all_cellmeans[j]
                 multfirst = sum(cross*all_cellns[i])
                 RSinter[i,j] = RSinter[j,i] = asarray(multfirst)
                 SSm[i,j] = SSm[j,i] = (mean(all_cellmeans[i],None) *
                                        mean(all_cellmeans[j],None) *
                                        len(all_cellmeans[i]) *hn)
         #SSw = RSw - RSinter

### HERE BEGINS THE MAXWELL & DELANEY APPROACH TO CALCULATING SS
         Lsource = makelist(sourcebetweens,Nfactors+1)
         #btwsourcecols = (array(map(Bscols.index,Lsource))-1).tolist()
         Bbtwnonsourcedims = ~source & Bbetweens
         Lbtwnonsourcedims = makelist(Bbtwnonsourcedims,Nfactors+1)
         btwnonsourcedims = (array(map(Bscols.index,Lbtwnonsourcedims))-1).tolist()

       ## Average Marray over non-source axes
         sourceDMarray = DM[dindex] *1.0
         for dim in btwnonsourcedims: # collapse all non-source dims
             if dim == len(DM[dindex].shape)-1:
                 raise ValueError, "Crashing ... shouldn't ever collapse ACROSS variables"
             sourceDMarray = expand_dims(mean(sourceDMarray,dim),dim)

       ## Calculate harmonic means for each level in source
         sourceDNarray = apply_over_axes(hmean, DN[dindex],btwnonsourcedims)

       ## Calc grand average (ga), used for ALL effects
         variableNs = apply_over_axes(sum, sourceDNarray,
                                      range(len(sourceDMarray.shape)-1))
         ga = apply_over_axes(sum, (sourceDMarray*sourceDNarray) / \
                              variableNs,
                              range(len(sourceDMarray.shape)-1))

       ## If GRAND interaction, use harmonic mean of ALL cell Ns
         if source == Nallsources-1:
             sourceDNarray = hmean(DN[dindex],
                                          range(len(sourceDMarray.shape)-1))
                
       ## Calc all SUBSOURCES to be subtracted from sourceMarray (M&D p.320)
         sub_effects = ga *1.0   # start with grand mean
         for subsource in range(3,source-2,2):
       ## Make a list of the non-subsource axes
             #subsourcebtw = (subsource-1) & Bbetweens
             if (propersubset(subsource-1,source-1) and
                 (subsource-1)&Bwithins == (source-1)&Bwithins and
                 (subsource-1) <> (source-1)&Bwithins):
                 sub_effects = (sub_effects +
                                alleffects[alleffsources.index(subsource)])

       ## Calc this effect (a(j)'s, b(k)'s, ab(j,k)'s, whatever)
         effect = sourceDMarray - sub_effects

       ## Save it so you don't have to calculate it again next time
         alleffects.append(effect)
         alleffsources.append(source)

       ## Calc and save sums of squares for this source
         SS = zeros((levels,levels),'f')
         SS = sum((effect**2 *sourceDNarray) *
                   multiply.reduce(take(DM[dindex].shape,btwnonsourcedims)),
                         range(len(sourceDMarray.shape)-1))
       ## Save it so you don't have to calculate it again next time
         SSlist.append(SS)
         SSsources.append(source)

         return SS


def multivar_sscalc(workd):
###
### DO SS CALCS ON THE OUTPUT FROM THE SOURCE=0 AND SOURCE=-1 CASES
###
     # this section expects workd to have subj. in LAST axis!!!!!!
     if len(workd.shape) == 1:
         levels = 1
     else:
         levels = workd.shape[0] # works because workd is always 2D

     sserr = zeros((levels,levels),'f')
     for i in range(levels):
         for j in range(i,levels):
             ssval = add.reduce(workd[i]*workd[j])
             sserr[i,j] = ssval
             sserr[j,i] = ssval
     return sserr


def subtr_cellmeans(workd,subjslots):
     """
Subtract all cell means when within-subjects factors are present ...
i.e., calculate full-model using a D-variable.
"""
     # Get a list of all dims that are source and between-subj
     sourcedims = makelist(Bbetweens,Nfactors+1)

     # Now, fix this list by mapping the dims from the original source
     # to dims for a between-subjects variable (namely, subjslots)
     transidx = range(len(subjslots.shape))[1:] + [0] # put subj dim at end
     tsubjslots = transpose(subjslots,transidx) # get all Ss for this idx
     tworkd = transpose(workd) # swap subj. and variable dims
     errors = 1.0 * tworkd

     if len(sourcedims) == 0:
         idx = [-1]
         loopcap = [0]
     if len(sourcedims) <> 0:
         btwsourcedims = map(Bscols.index,sourcedims)
         idx = [0] * len(btwsourcedims)
         idx[0] = -1 # compensate for pre-increment of 1st slot in incr()
             
         # Get a list of the maximum values each factor can handle
         loopcap = take(array(Nlevels),sourcedims)-1

### WHILE STILL MORE GROUPS, CALCULATE GROUP MEAN FOR EACH D-VAR
     while incr(idx,loopcap) <> -1:  # loop through source btw level-combos
         mask = tsubjslots[idx]
         thisgroup = tworkd*mask[NewAxis,:]
         groupmns = mean(compress(mask,thisgroup),1)

### THEN SUBTRACT THEM FROM APPROPRIATE SUBJECTS
         errors = errors - multiply.outer(groupmns,mask)
     return errors


def f_value_wilks_lambda(ER, EF, dfnum, dfden, a, b):
     """
Calculation of Wilks lambda F-statistic for multivarite data, per
Maxwell & Delaney p.657.

"""
     if type(ER) in [IntType, FloatType]:
         ER = array([[ER]])
     if type(EF) in [IntType, FloatType]:
         EF = array([[EF]])
     lmbda = linalg.det(EF) / linalg.det(ER)
     if (a-1)**2 + (b-1)**2 == 5:
         q = 1
     else:
         q = math.sqrt( ((a-1)**2*(b-1)**2 - 2) / ((a-1)**2 + (b-1)**2 -5) )
     n_um = (1 - lmbda**(1.0/q))*(a-1)*(b-1)
     d_en = lmbda**(1.0/q) / (n_um*q - 0.5*(a-1)*(b-1) + 1)
     return n_um / d_en

def member(factor,source):
     return (1 << factor) & source != 0

def setsize(source):
     size = 0
     for bit in source:
         if bit == 1:
             size = size + 1
     return size

def subset (a,b):
     return (a&b)==a

def propersubset (a,b):
     sub = ((a&b)==a)
     if a==b:
         sub = 0
     return sub

def numlevels(source,Nlevels):
     for i in range(30): # find the biggest i such that 2**i >= source
         if 1<<i >= source:
             break
     levelcount = 1
     for j in range(i): # loop up through each bit
         if subset(1<<j,source):
             levelcount = levelcount * Nlevels[j] - 1
     return levelcount

def numbitson(a):
     numon = 0
     while a>0:
         numon = numon + a%2
         a = a>>1
     return numon

def makebin(sourcelist):
     outbin = 0
     for item in sourcelist:
         outbin = outbin + 2**item
     return outbin

def makelist(source,ncols):
     levellist = []
     for j in range(ncols):
         if subset(1<<j,source):
             levellist.append(j)
     return levellist

def round4(num):
     try:
         return around(num,4)
     except:
         return 'N/A'


def f_value (ER,EF,dfR,dfF):
    """
Returns an F-statistic given the following:
        ER  = error associated with the null hypothesis (the Restricted model)
        EF  = error associated with the alternate hypothesis (the Full model)
        dfR = degrees of freedom the Restricted model
        dfF = degrees of freedom associated with the Restricted model
"""
    return ((ER-EF)/float(dfR-dfF) / (EF/float(dfF)))


def outputfstats(Enum, Eden, dfnum, dfden, f, prob):
     Enum = around(Enum,3)
     Eden = around(Eden,3)
     dfnum = around(Enum,3)
     dfden = around(dfden,3)
     f = around(f,3)
     prob = around(prob,3)
     suffix = ''                       # for *s after the p-value
     if  prob < 0.001:  suffix = '  ***'
     elif prob < 0.01:  suffix = '  **'
     elif prob < 0.05:  suffix = '  *'
     title = [['EF/ER','DF','Mean Square','F-value','prob','']]
     lofl = title+[[Enum, dfnum, around(Enum/float(dfnum),3), f, prob, suffix],
                   [Eden, dfden, around(Eden/float(dfden),3),'','','']]
     _support.printcc(lofl)
     return


def f_value_multivariate(ER, EF, dfnum, dfden):
     """
Returns an F-statistic given the following:
        ER  = error associated with the null hypothesis (the Restricted model)
        EF  = error associated with the alternate hypothesis (the Full model)
        dfR = degrees of freedom the Restricted model
        dfF = degrees of freedom associated with the Restricted model
where ER and EF are matrices from a multivariate F calculation.
"""
     if type(ER) in [IntType, FloatType]:
         ER = array([[ER]])
     if type(EF) in [IntType, FloatType]:
         EF = array([[EF]])
     n_um = (linalg.det(ER) - linalg.det(EF)) / float(dfnum)
     d_en = linalg.det(EF) / float(dfden)
     return n_um / d_en


#####################################
#######  SUPPORT FUNCTIONS  ########
#####################################

# sign is in Numeric
sum = scipy.sum
cumsum = scipy.cumsum

def ss(a, axis=-1):
    """
Squares each value in the passed array, adds these squares & returns
the result.  Axis can equal None (ravel array first), an integer
(the axis over which to operate), or a sequence (operate over
multiple axes).

Returns: sum-along-'axis' for (a*a)
"""
    a, axis = _chk_asarray(a, axis)
    return sum(a*a,axis)

def summult (array1,array2,axis=-1):
    """
Multiplies elements in array1 and array2, element by element, and
returns the sum (along 'axis') of all resulting multiplications.
Axis can equal None (ravel array first), or an integer (the
axis over which to operate),
"""
    array1, array2 = map(asarray, (array1, array2))
    if axis is None:
        array1 = ravel(array1)
        array2 = ravel(array2)
        axis = 0
    return sum(array1*array2,axis)


def square_of_sums(a, axis=-1):    
    """Adds the values in the passed array, squares that sum, and returns the
result.

Returns: the square of the sum over axis.
"""
    a, axis = _chk_asarray(a, axis)
    s = sum(a,axis)
    if type(s) == ArrayType:
        return s.astype(Float)*s
    else:
        return float(s)*s


def sumdiffsquared(a, b, axis=-1):
    """
Takes pairwise differences of the values in arrays a and b, squares
these differences, and returns the sum of these squares.  Axis
can equal None (ravel array first), an integer (the axis over
which to operate).

Returns: sum[(a-b)**2]
"""

    a, b = _chk2_asarray(a, b, axis)
    return sum((a-b)**2,axis)

def fastsort(a):
    it = argsort(a)
    as = take(a, it)
    return as, it

def shellsort(a):
    """
Shellsort algorithm.  Sorts a 1D-array.

Returns: sorted-a, sorting-index-vector (for original array)

Use fastsort for speed.
"""
    a = asarray(a)
    n = len(a)
    svec = a*1.0
    ivec = range(n)
    gap = n/2   # integer division needed
    while gap >0:
        for i in range(gap,n):
            for j in range(i-gap,-1,-gap):
                while j>=0 and svec[j]>svec[j+gap]:
                    temp        = svec[j]
                    svec[j]     = svec[j+gap]
                    svec[j+gap] = temp
                    itemp       = ivec[j]
                    ivec[j]     = ivec[j+gap]
                    ivec[j+gap] = itemp
        gap = gap / 2  # integer division needed
#    svec is now sorted input vector, ivec has the order svec[i] = vec[ivec[i]]
    return array(svec), array(ivec)

def rankdata(a):
    """
Ranks the data in a, dealing with ties appropritely.  First ravels 
a.  Adapted from Gary Perlman's |Stat ranksort.

Returns: array of length equal to a, containing rank scores
"""
    a = ravel(a)
    n = len(a)
    svec, ivec = fastsort(a)
    sumranks = 0
    dupcount = 0
    newarray = zeros(n,Float)
    for i in range(n):
        sumranks = sumranks + i
        dupcount = dupcount + 1
        if i==n-1 or svec[i] <> svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newarray


def writecc (listoflists,file,writetype='w',extra=2):
    """
Writes a list of lists to a file in columns, customized by the max
size of items within the columns (max size of items in col, +2 characters)
to specified file.  File-overwrite is the default.

Usage:   writecc (listoflists,file,writetype='w',extra=2)
Returns: None
"""
    if type(listoflists[0]) not in [ListType,TupleType]:
        listoflists = [listoflists]
    outfile = open(file,writetype)
    rowstokill = []
    list2print = copy.deepcopy(listoflists)
    for i in range(len(listoflists)):
        if listoflists[i] == ['\n'] or listoflists[i]=='\n' or listoflists[i]=='dashes':
            rowstokill = rowstokill + [i]
    rowstokill.reverse()
    for row in rowstokill:
        del list2print[row]
    maxsize = [0]*len(list2print[0])
    for col in range(len(list2print[0])):
        items = _support.colex(list2print,col)
        items = map(_support.makestr,items)
        maxsize[col] = max(map(len,items)) + extra
    for row in listoflists:
        if row == ['\n'] or row == '\n':
            outfile.write('\n')
        elif row == ['dashes'] or row == 'dashes':
            dashes = [0]*len(maxsize)
            for j in range(len(maxsize)):
                dashes[j] = '-'*(maxsize[j]-2)
            outfile.write(_support.lineincustcols(dashes,maxsize))
        else:
            outfile.write(_support.lineincustcols(row,maxsize))
        outfile.write('\n')
    outfile.close()
    return None


def outputpairedstats(fname,writemode,name1,n1,m1,se1,min1,max1,name2,n2,m2,se2,min2,max2,statname,stat,prob):
    """
Prints or write to a file stats for two groups, using the name, n,
mean, sterr, min and max for each group, as well as the statistic name,
its value, and the associated p-value.

Usage:   outputpairedstats(fname,writemode,
                           name1,n1,mean1,stderr1,min1,max1,
                           name2,n2,mean2,stderr2,min2,max2,
                           statname,stat,prob)
Returns: None
"""
    suffix = ''                       # for *s after the p-value
    try:
        x = prob.shape
        prob = prob[0]
    except:
        pass
    if  prob < 0.001:  suffix = '  ***'
    elif prob < 0.01:  suffix = '  **'
    elif prob < 0.05:  suffix = '  *'
    title = [['Name','N','Mean','SD','Min','Max']]
    lofl = title+[[name1,n1,round(m1,3),round(math.sqrt(se1),3),min1,max1],
                  [name2,n2,round(m2,3),round(math.sqrt(se2),3),min2,max2]]
    if type(fname)<>StringType or len(fname)==0:
        print
        print statname
        print
        _support.printcc(lofl)
        print
        try:
            if stat.shape == ():
                stat = stat[0]
            if prob.shape == ():
                prob = prob[0]
        except:
            pass
        print 'Test statistic = ',round(stat,3),'   p = ',round(prob,3),suffix
        print
    else:
        file = open(fname,writemode)
        file.write('\n'+statname+'\n\n')
        file.close()
        writecc(lofl,fname,'a')
        file = open(fname,'a')
        try:
            if stat.shape == ():
                stat = stat[0]
            if prob.shape == ():
                prob = prob[0]
        except:
            pass
        file.write(_support.list2string(['\nTest statistic = ',round(stat,4),'   p = ',round(prob,4),suffix,'\n\n']))
        file.close()
    return None


def findwithin(data):
    """
Returns a binary vector, 1=within-subject factor, 0=between.  Input
equals the entire data array (i.e., column 1=random factor, last
column = measured values.

"""
    numfact = len(data[0])-2
    withinvec = [0]*numfact
    for col in range(1,numfact+1):
        rows = _support.linexand(data,col,_support.unique(_support.colex(data,1))[0])  # get 1 level of this factor
        if len(_support.unique(_support.colex(rows,0))) < len(rows):   # if fewer subjects than scores on this factor
            withinvec[col-1] = 1
    return withinvec
