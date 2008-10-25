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

#################################################
#######  Written by:  Gary Strangman  ###########
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

MOMENTS HANDLING NAN: nanmean
                      nanmedian
                      nanstd

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

## Note that scipy.stats.distributions has many more statistical probability
## functions defined.


ANOVA FUNCTIONS:  f_oneway
                  f_value

SUPPORT FUNCTIONS:  ss
                    square_of_sums
                    shellsort
                    rankdata

References
----------
[CRCProbStat2000] Zwillinger, D. and Kokoska, S. _CRC Standard Probablity and
Statistics Tables and Formulae_. Chapman & Hall: New York. 2000.
"""
## CHANGE LOG:
## ===========
## since 2001-06-25 ... see scipy SVN changelog
## 05-11-29 ... fixed default axis to be 0 for consistency with scipy;
##              cleanup of redundant imports, dead code, {0,1} -> booleans
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

# Standard library imports.
import warnings
import math

# Scipy imports.
from numpy import array, asarray, dot, ma, zeros, sum
import scipy.special as special
import scipy.linalg as linalg
import numpy as np

# Local imports.
import _support

__all__ = ['gmean', 'hmean', 'mean', 'cmedian', 'median', 'mode',
           'tmean', 'tvar', 'tmin', 'tmax', 'tstd', 'tsem',
           'moment', 'variation', 'skew', 'kurtosis', 'describe',
           'skewtest', 'kurtosistest', 'normaltest',
           'itemfreq', 'scoreatpercentile', 'percentileofscore',
           'histogram', 'histogram2', 'cumfreq', 'relfreq',
           'obrientransform', 'samplevar', 'samplestd', 'signaltonoise',
           'var', 'std', 'stderr', 'sem', 'z', 'zs', 'zmap',
           'threshold', 'trimboth', 'trim1', 'trim_mean',
           'cov', 'corrcoef', 'f_oneway', 'pearsonr', 'spearmanr',
           'pointbiserialr', 'kendalltau', 'linregress',
           'ttest_1samp', 'ttest_ind', 'ttest_rel',
           'kstest', 'chisquare', 'ks_2samp', 'mannwhitneyu',
           'tiecorrect', 'ranksums', 'kruskal', 'friedmanchisquare',
           'zprob', 'erfc', 'chisqprob', 'ksprob', 'fprob', 'betai',
           'glm', 'f_value_wilks_lambda',
           'f_value', 'f_value_multivariate',
           'ss', 'square_of_sums',
           'fastsort', 'rankdata',
           'nanmean', 'nanstd', 'nanmedian',
          ]


def _chk_asarray(a, axis):
    if axis is None:
        a = np.ravel(a)
        outaxis = 0
    else:
        a = np.asarray(a)
        outaxis = axis
    return a, outaxis

def _chk2_asarray(a, b, axis):
    if axis is None:
        a = np.ravel(a)
        b = np.ravel(b)
        outaxis = 0
    else:
        a = np.asarray(a)
        b = np.asarray(b)
        outaxis = axis
    return a, b, outaxis

#######
### NAN friendly functions
########

def nanmean(x, axis=0):
    """Compute the mean over the given axis ignoring nans.

    :Parameters:
        x : ndarray
            input array
        axis : int
            axis along which the mean is computed.

    :Results:
        m : float
            the mean."""
    x, axis = _chk_asarray(x,axis)
    x = x.copy()
    Norig = x.shape[axis]
    factor = 1.0-np.sum(np.isnan(x),axis)*1.0/Norig

    x[np.isnan(x)] = 0
    return np.mean(x,axis)/factor

def nanstd(x, axis=0, bias=False):
    """Compute the standard deviation over the given axis ignoring nans

    :Parameters:
        x : ndarray
            input array
        axis : int
            axis along which the standard deviation is computed.
        bias : boolean
            If true, the biased (normalized by N) definition is used. If false,
            the unbiased is used (the default).

    :Results:
        s : float
            the standard deviation."""
    x, axis = _chk_asarray(x,axis)
    x = x.copy()
    Norig = x.shape[axis]

    Nnan = np.sum(np.isnan(x),axis)*1.0
    n = Norig - Nnan

    x[np.isnan(x)] = 0.
    m1 = np.sum(x,axis)/n

    # Kludge to subtract m1 from the correct axis
    if axis!=0:
        shape = np.arange(x.ndim).tolist()
        shape.remove(axis)
        shape.insert(0,axis)
        x = x.transpose(tuple(shape))
        d = (x-m1)**2.0
        shape = tuple(array(shape).argsort())
        d = d.transpose(shape)
    else:
        d = (x-m1)**2.0
    m2 = np.sum(d,axis)-(m1*m1)*Nnan
    if bias:
        m2c = m2 / n
    else:
        m2c = m2 / (n - 1.)
    return np.sqrt(m2c)

def _nanmedian(arr1d):  # This only works on 1d arrays
    """Private function for rank a arrays. Compute the median ignoring Nan.

    :Parameters:
        arr1d : rank 1 ndarray
            input array

    :Results:
        m : float
            the median."""
    cond = 1-np.isnan(arr1d)
    x = np.sort(np.compress(cond,arr1d,axis=-1))
    if x.size == 0:
        return np.nan
    return median(x)

def nanmedian(x, axis=0):
    """ Compute the median along the given axis ignoring nan values

    :Parameters:
        x : ndarray
            input array
        axis : int
            axis along which the median is computed.

    :Results:
        m : float
            the median."""
    x, axis = _chk_asarray(x,axis)
    x = x.copy()
    return np.apply_along_axis(_nanmedian,axis,x)


#####################################
########  CENTRAL TENDENCY  ########
#####################################

def gmean(a, axis=0):
    """Calculates the geometric mean of the values in the passed array.

    That is:  n-th root of (x1 * x2 * ... * xn)

    Parameters
    ----------
    a : array of positive values
    axis : int or None
    zero_sub : value to substitute for zero values. Default is 0.

    Returns
    -------
    The geometric mean computed over a single dimension of the input array or
    all values in the array if axis==None.
    """
    a, axis = _chk_asarray(a, axis)
    log_a = np.log(a)
    return np.exp(log_a.mean(axis=axis))


def hmean(a, axis=0, zero_sub=0):
    """Calculates the harmonic mean of the values in the passed array.

    That is:  n / (1/x1 + 1/x2 + ... + 1/xn)

    Parameters
    ----------
    a : array
    axis : int or None

    Returns
    -------
    The harmonic mean computed over a single dimension of the input array or all
    values in the array if axis=None.
    """
    a, axis = _chk_asarray(a, axis)
    size = a.shape[axis]
    return size / np.sum(1.0/a, axis)

def mean(a, axis=0):
    # fixme: This seems to be redundant with numpy.mean(,axis=0) or even
    # the ndarray.mean() method.
    """Returns the arithmetic mean of m along the given dimension.

    That is: (x1 + x2 + .. + xn) / n

    Parameters
    ----------
    a : array
    axis : int or None

    Returns
    -------
    The arithmetic mean computed over a single dimension of the input array or
    all values in the array if axis=None. The return value will have a floating
    point dtype even if the input data are integers.
    """
    warnings.warn("""\
scipy.stats.mean is deprecated; please update your code to use numpy.mean.
Please note that:
    - numpy.mean axis argument defaults to None, not 0
    - numpy.mean has a ddof argument to replace bias in a more general manner.
      scipy.stats.mean(a, bias=True) can be replaced by numpy.mean(x,
axis=0, ddof=1).""", DeprecationWarning)
    a, axis = _chk_asarray(a, axis)
    return a.mean(axis)

def cmedian(a, numbins=1000):
    # fixme: numpy.median() always seems to be a better choice.
    # A better version of this function would take already-histogrammed data
    # and compute the median from that.
    # fixme: the wording of the docstring is a bit wonky.
    """Returns the computed median value of an array.

    All of the values in the input array are used. The input array is first
    histogrammed using numbins bins. The bin containing the median is
    selected by searching for the halfway point in the cumulative histogram.
    The median value is then computed by linearly interpolating across that bin.

    Parameters
    ----------
    a : array
    numbins : int
        The number of bins used to histogram the data. More bins give greater
        accuracy to the approximation of the median.

    Returns
    -------
    A floating point value approximating the median.

    References
    ----------
    [CRCProbStat2000] Section 2.2.6
    """
    a = np.ravel(a)
    n = float(len(a))

    # We will emulate the (fixed!) bounds selection scheme used by
    # scipy.stats.histogram(), but use numpy.histogram() since it is faster.
    amin = a.min()
    amax = a.max()
    estbinwidth = (amax - amin)/float(numbins - 1)
    binsize = (amax - amin + estbinwidth) / float(numbins)
    (hist, bins) = np.histogram(a, numbins,
        range=(amin-binsize*0.5, amax+binsize*0.5))
    binsize = bins[1] - bins[0]
    cumhist = np.cumsum(hist)           # make cumulative histogram
    cfbin = np.searchsorted(cumhist, n/2.0)
    LRL = bins[cfbin]      # get lower read limit of that bin
    if cfbin == 0:
        cfbelow = 0.0
    else:
        cfbelow = cumhist[cfbin-1]       # cum. freq. below bin
    freq = hist[cfbin]                  # frequency IN the 50%ile bin
    median = LRL + ((n/2.0-cfbelow)/float(freq))*binsize # MEDIAN
    return median

def median(a, axis=0):
    # fixme: This would be redundant with numpy.median() except that the latter
    # does not deal with arbitrary axes.
    """Returns the median of the passed array along the given axis.

    If there is an even number of entries, the mean of the
    2 middle values is returned.

    Parameters
    ----------
    a : array
    axis=0 : int

    Returns
    -------
    The median of each remaining axis, or of all of the values in the array
    if axis is None.
    """
    warnings.warn("""\
scipy.stats.median is deprecated; please update your code to use numpy.median.
Please note that:
    - numpy.median axis argument defaults to None, not 0
    - numpy.median has a ddof argument to replace bias in a more general manner.
      scipy.stats.median(a, bias=True) can be replaced by numpy.median(x,
axis=0, ddof=1).""", DeprecationWarning)
    return np.median(a, axis)

def mode(a, axis=0):
    """Returns an array of the modal (most common) value in the passed array.

    If there is more than one such value, only the first is returned.
    The bin-count for the modal bins is also returned.

    Parameters
    ----------
    a : array
    axis=0 : int

    Returns
    -------
    (array of modal values, array of counts for each mode)
    """
    a, axis = _chk_asarray(a, axis)
    scores = np.unique(np.ravel(a))       # get ALL unique values
    testshape = list(a.shape)
    testshape[axis] = 1
    oldmostfreq = np.zeros(testshape)
    oldcounts = np.zeros(testshape)
    for score in scores:
        template = (a == score)
        counts = np.expand_dims(np.sum(template, axis),axis)
        mostfrequent = np.where(counts > oldcounts, score, oldmostfreq)
        oldcounts = np.maximum(counts, oldcounts)
        oldmostfreq = mostfrequent
    return mostfrequent, oldcounts

def mask_to_limits(a, limits, inclusive):
    """Mask an array for values outside of given limits.

    This is primarily a utility function.

    Parameters
    ----------
    a : array
    limits : (float or None, float or None)
        A tuple consisting of the (lower limit, upper limit).  Values in the
        input array less than the lower limit or greater than the upper limit
        will be masked out. None implies no limit.
    inclusive : (bool, bool)
        A tuple consisting of the (lower flag, upper flag).  These flags
        determine whether values exactly equal to lower or upper are allowed.

    Returns
    -------
    A MaskedArray.

    Raises
    ------
    A ValueError if there are no values within the given limits.
    """
    lower_limit, upper_limit = limits
    lower_include, upper_include = inclusive
    am = ma.MaskedArray(a)
    if lower_limit is not None:
        if lower_include:
            am = ma.masked_less(am, lower_limit)
        else:
            am = ma.masked_less_equal(am, lower_limit)
    if upper_limit is not None:
        if upper_include:
            am = ma.masked_greater(am, upper_limit)
        else:
            am = ma.masked_greater_equal(am, upper_limit)
    if am.count() == 0:
        raise ValueError("No array values within given limits")
    return am

def tmean(a, limits=None, inclusive=(True, True)):
    """Returns the arithmetic mean of all values in an array, ignoring values
    strictly outside given limits.

    Parameters
    ----------
    a : array
    limits : None or (lower limit, upper limit)
        Values in the input array less than the lower limit or greater than the
        upper limit will be masked out. When limits is None, then all values are
        used. Either of the limit values in the tuple can also be None
        representing a half-open interval.
    inclusive : (bool, bool)
        A tuple consisting of the (lower flag, upper flag).  These flags
        determine whether values exactly equal to lower or upper are allowed.

    Returns
    -------
    A float.
    """
    a = asarray(a)

    # Cast to a float if this is an integer array. If it is already a float
    # array, leave it as is to preserve its precision.
    if issubclass(a.dtype.type, np.integer):
        a = a.astype(float)

    # No trimming.
    if limits is None:
        return mean(a,None)

    am = mask_to_limits(a.ravel(), limits, inclusive)
    return am.mean()

def masked_var(am):
    m = am.mean()
    s = ma.add.reduce((am - m)**2)
    n = am.count() - 1.0
    return s / n

def tvar(a, limits=None, inclusive=(1,1)):
    """Returns the sample variance of values in an array, (i.e., using
    N-1), ignoring values strictly outside the sequence passed to
    'limits'.  Note: either limit in the sequence, or the value of
    limits itself, can be set to None.  The inclusive list/tuple
    determines whether the lower and upper limiting bounds
    (respectively) are open/exclusive (0) or closed/inclusive (1).
    """
    a = asarray(a)
    a = a.astype(float).ravel()
    if limits is None:
        n = len(a)
        return a.var()*(n/(n-1.))
    am = mask_to_limits(a, limits, inclusive)
    return masked_var(am)

def tmin(a, lowerlimit=None, axis=0, inclusive=True):
    """Returns the minimum value of a, along axis, including only values
    less than (or equal to, if inclusive is True) lowerlimit.  If the
    limit is set to None, all values in the array are used.
    """
    a, axis = _chk_asarray(a, axis)
    am = mask_to_limits(a, (lowerlimit, None), (inclusive, False))
    return ma.minimum.reduce(am, axis)

def tmax(a, upperlimit, axis=0, inclusive=True):
    """Returns the maximum value of a, along axis, including only values
    greater than (or equal to, if inclusive is True) upperlimit.  If the limit
    is set to None, a limit larger than the max value in the array is
    used.
    """
    a, axis = _chk_asarray(a, axis)
    am = mask_to_limits(a, (None, upperlimit), (False, inclusive))
    return ma.maximum.reduce(am, axis)

def tstd(a, limits=None, inclusive=(1,1)):
    """Returns the standard deviation of all values in an array,
    ignoring values strictly outside the sequence passed to 'limits'.
    Note: either limit in the sequence, or the value of limits itself,
    can be set to None.  The inclusive list/tuple determines whether the
    lower and upper limiting bounds (respectively) are open/exclusive
    (0) or closed/inclusive (1).
    """
    return np.sqrt(tvar(a,limits,inclusive))


def tsem(a, limits=None, inclusive=(True,True)):
    """Returns the standard error of the mean for the values in an array,
    (i.e., using N for the denominator), ignoring values strictly outside
    the sequence passed to 'limits'.   Note: either limit in the
    sequence, or the value of limits itself, can be set to None.  The
    inclusive list/tuple determines whether the lower and upper limiting
    bounds (respectively) are open/exclusive (0) or closed/inclusive (1).
    """
    a = np.asarray(a).ravel()
    if limits is None:
        n = float(len(a))
        return a.std()/np.sqrt(n)
    am = mask_to_limits(a.ravel(), limits, inclusive)
    sd = np.sqrt(masked_var(am))
    return sd / am.count()


#####################################
############  MOMENTS  #############
#####################################

def moment(a, moment=1, axis=0):
    """Calculates the nth moment about the mean for a sample.

    Generally used to calculate coefficients of skewness and
    kurtosis.

    Parameters
    ----------
    a : array
    moment : int
    axis : int or None

    Returns
    -------
    The appropriate moment along the given axis or over all values if axis is
    None.
    """
    a, axis = _chk_asarray(a, axis)
    if moment == 1:
        # By definition the first moment about the mean is 0.
        shape = list(a.shape)
        del shape[axis]
        if shape:
            # return an actual array of the appropriate shape
            return np.zeros(shape, dtype=float)
        else:
            # the input was 1D, so return a scalar instead of a rank-0 array
            return np.float64(0.0)
    else:
        mn = np.expand_dims(np.mean(a,axis), axis)
        s = np.power((a-mn), moment)
        return np.mean(s, axis)


def variation(a, axis=0):
    """Computes the coefficient of variation, the ratio of the biased standard
    deviation to the mean.

    Parameters
    ----------
    a : array
    axis : int or None

    References
    ----------
    [CRCProbStat2000] section 2.2.20
    """
    a, axis = _chk_asarray(a, axis)
    n = a.shape[axis]
    return a.std(axis)/a.mean(axis)


def skew(a, axis=0, bias=True):
    """Computes the skewness of a data set.

    For normally distributed data, the skewness should be about 0. A skewness
    value > 0 means that there is more weight in the left tail of the
    distribution. The function skewtest() can be used to determine if the
    skewness value is close enough to 0, statistically speaking.

    Parameters
    ----------
    a : array
    axis : int or None
    bias : bool
        If False, then the calculations are corrected for statistical bias.

    Returns
    -------
    The skewness of values along an axis, returning 0 where all values are
    equal.

    References
    ----------
    [CRCProbStat2000] section 2.2.24.1
    """
    a, axis = _chk_asarray(a,axis)
    n = a.shape[axis]
    m2 = moment(a, 2, axis)
    m3 = moment(a, 3, axis)
    zero = (m2 == 0)
    vals = np.where(zero, 0, m3 / m2**1.5)
    if not bias:
        can_correct = (n > 2) & (m2 > 0)
        if can_correct.any():
            m2 = np.extract(can_correct, m2)
            m3 = np.extract(can_correct, m3)
            nval = np.sqrt((n-1.0)*n)/(n-2.0)*m3/m2**1.5
            np.place(vals, can_correct, nval)
    return vals

def kurtosis(a, axis=0, fisher=True, bias=True):
    """Computes the kurtosis (Fisher or Pearson) of a dataset.

    Kurtosis is the fourth central moment divided by the square of the variance.
    If Fisher's definition is used, then 3.0 is subtracted from the result to
    give 0.0 for a normal distribution.

    If bias is False then the kurtosis is calculated using k statistics to
    eliminate bias comming from biased moment estimators

    Use kurtosistest() to see if result is close enough to normal.

    Parameters
    ----------
    a : array
    axis : int or None
    fisher : bool
        If True, Fisher's definition is used (normal ==> 0.0). If False,
        Pearson's definition is used (normal ==> 3.0).
    bias : bool
        If False, then the calculations are corrected for statistical bias.

    Returns
    -------
    The kurtosis of values along an axis. If all values are equal, return -3 for Fisher's
    definition and 0 for Pearson's definition.


    References
    ----------
    [CRCProbStat2000] section 2.2.25
    """
    a, axis = _chk_asarray(a, axis)
    n = a.shape[axis]
    m2 = moment(a,2,axis)
    m4 = moment(a,4,axis)
    zero = (m2 == 0)
    vals = np.where(zero, 0, m4/ m2**2.0)
    if not bias:
        can_correct = (n > 3) & (m2 > 0)
        if can_correct.any():
            m2 = np.extract(can_correct, m2)
            m4 = np.extract(can_correct, m4)
            nval = 1.0/(n-2)/(n-3)*((n*n-1.0)*m4/m2**2.0-3*(n-1)**2.0)
            np.place(vals, can_correct, nval+3.0)
    if fisher:
        return vals - 3
    else:
        return vals

def describe(a, axis=0):
    """Computes several descriptive statistics of the passed array.

    Parameters
    ----------
    a : array
    axis : int or None

    Returns
    -------
    (size of the data,
     (min, max),
     arithmetic mean,
     unbiased variance,
     biased skewness,
     biased kurtosis)
    """
    a, axis = _chk_asarray(a, axis)
    n = a.shape[axis]
    mm = (np.minimum.reduce(a), np.maximum.reduce(a))
    m = mean(a, axis)
    v = var(a, axis)
    sk = skew(a, axis)
    kurt = kurtosis(a, axis)
    return n, mm, m, v, sk, kurt

#####################################
########  NORMALITY TESTS  ##########
#####################################

def skewtest(a, axis=0):
    """Tests whether the skew is significantly different from a normal
    distribution.

    The size of the dataset should be >= 8.

    Parameters
    ----------
    a : array
    axis : int or None

    Returns
    -------
    (Z-score,
     2-tail Z-probability,
    )
    """
    a, axis = _chk_asarray(a, axis)
    if axis is None:
        a = np.ravel(a)
        axis = 0
    b2 = skew(a,axis)
    n = float(a.shape[axis])
    if n < 8:
        warnings.warn(
            "skewtest only valid for n>=8 ... continuing anyway, n=%i" %
            int(n))
    y = b2 * math.sqrt(((n+1)*(n+3)) / (6.0*(n-2)) )
    beta2 = ( 3.0*(n*n+27*n-70)*(n+1)*(n+3) ) / ( (n-2.0)*(n+5)*(n+7)*(n+9) )
    W2 = -1 + math.sqrt(2*(beta2-1))
    delta = 1/math.sqrt(0.5*math.log(W2))
    alpha = math.sqrt(2.0/(W2-1))
    y = np.where(y==0, 1, y)
    Z = delta*np.log(y/alpha + np.sqrt((y/alpha)**2+1))
    return Z, (1.0 - zprob(Z))*2


def kurtosistest(a, axis=0):
    """Tests whether a dataset has normal kurtosis (i.e.,
    kurtosis=3(n-1)/(n+1)).

    Valid only for n>20.

    Parameters
    ----------
    a : array
    axis : int or None

    Returns
    -------
    (Z-score,
     2-tail Z-probability)
    The Z-score is set to 0 for bad entries.
    """
    a, axis = _chk_asarray(a, axis)
    n = float(a.shape[axis])
    if n < 20:
        warnings.warn(
            "kurtosistest only valid for n>=20 ... continuing anyway, n=%i" %
            int(n))
    b2 = kurtosis(a, axis, fisher=False)
    E = 3.0*(n-1) /(n+1)
    varb2 = 24.0*n*(n-2)*(n-3) / ((n+1)*(n+1)*(n+3)*(n+5))
    x = (b2-E)/np.sqrt(varb2)
    sqrtbeta1 = 6.0*(n*n-5*n+2)/((n+7)*(n+9)) * np.sqrt((6.0*(n+3)*(n+5))/
                                                       (n*(n-2)*(n-3)))
    A = 6.0 + 8.0/sqrtbeta1 *(2.0/sqrtbeta1 + np.sqrt(1+4.0/(sqrtbeta1**2)))
    term1 = 1 -2/(9.0*A)
    denom = 1 +x*np.sqrt(2/(A-4.0))
    denom = np.where(denom < 0, 99, denom)
    term2 = np.where(denom < 0, term1, np.power((1-2.0/A)/denom,1/3.0))
    Z = ( term1 - term2 ) / np.sqrt(2/(9.0*A))
    Z = np.where(denom == 99, 0, Z)
    return Z, (1.0-zprob(Z))*2


def normaltest(a, axis=0):
    """Tests whether skew and/or kurtosis of dataset differs from normal curve.

    Parameters
    ----------
    a : array
    axis : int or None

    Returns
    -------
    (Chi^2 score,
     2-tail probability)

    Based on the D'Agostino and Pearson's test that combines skew and
    kurtosis to produce an omnibus test of normality.

    D'Agostino, R. B. and Pearson, E. S. (1971), "An Omnibus Test of
    Normality for Moderate and Large Sample Size," Biometrika, 58, 341-348

    D'Agostino, R. B. and Pearson, E. S. (1973), "Testing for departures from
    Normality," Biometrika, 60, 613-622

    """
    a, axis = _chk_asarray(a, axis)
    s,p = skewtest(a,axis)
    k,p = kurtosistest(a,axis)
    k2 = s*s + k*k
    return k2, chisqprob(k2,2)

# Martinez-Iglewicz test
# K-S test

#####################################
######  FREQUENCY FUNCTIONS  #######
#####################################

def itemfreq(a):
    # fixme: I'm not sure I understand what this does. The docstring is
    # internally inconsistent.
    # comment: fortunately, this function doesn't appear to be used elsewhere
    """Returns a 2D array of item frequencies.

    Column 1 contains item values, column 2 contains their respective counts.
    Assumes a 1D array is passed.

    Parameters
    ----------
    a : array

    Returns
    -------
    A 2D frequency table (col [0:n-1]=scores, col n=frequencies)
    """
    scores = _support.unique(a)
    scores = np.sort(scores)
    freq = zeros(len(scores))
    for i in range(len(scores)):
        freq[i] = np.add.reduce(np.equal(a,scores[i]))
    return array(_support.abut(scores, freq))


def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a)*fraction;

def scoreatpercentile(a, per, limit=()):
    """Calculate the score at the given 'per' percentile of the
    sequence a.  For example, the score at per=50 is the median.

    If the desired quantile lies between two data points, we
    interpolate between them.

    If the parameter 'limit' is provided, it should be a tuple (lower,
    upper) of two values.  Values of 'a' outside this (closed)
    interval will be ignored.

    """
    # TODO: this should be a simple wrapper around a well-written quantile
    # function.  GNU R provides 9 quantile algorithms (!), with differing
    # behaviour at, for example, discontinuities.
    values = np.sort(a,axis=0)
    if limit:
        values = values[(limit[0] < a) & (a < limit[1])]

    idx = per /100. * (values.shape[0] - 1)
    if (idx % 1 == 0):
        return values[idx]
    else:
        return _interpolate(values[int(idx)], values[int(idx) + 1], idx % 1)


def percentileofscore(a, score, histbins=10, defaultlimits=None):
    # fixme: Again with the histogramming. This probably should be replaced by
    # an empirical CDF approach.
    """
Note: result of this function depends on the values used to histogram
the data(!).

Returns: percentile-position of score (0-100) relative to a
"""
    h, lrl, binsize, extras = histogram(a,histbins,defaultlimits)
    cumhist = np.cumsum(h*1, axis=0)
    i = int((score - lrl)/float(binsize))
    pct = (cumhist[i-1]+((score-(lrl+binsize*i))/float(binsize))*h[i])/float(len(a)) * 100
    return pct


def histogram2(a, bins):
    # comment: probably obsoleted by numpy.histogram()
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
    n = np.searchsorted(np.sort(a), bins)
    n = np.concatenate([ n, [len(a)]])
    return n[ 1:]-n[:-1]




def histogram(a, numbins=10, defaultlimits=None, printextras=True):
    # fixme: use numpy.histogram() to implement
    """
Returns (i) an array of histogram bin counts, (ii) the smallest value
of the histogram binning, and (iii) the bin width (the last 2 are not
necessarily integers).  Default number of bins is 10.  Defaultlimits
can be None (the routine picks bins spanning all the numbers in the
a) or a 2-sequence (lowerlimit, upperlimit).  Returns all of the
following: array of bin values, lowerreallimit, binsize, extrapoints.

Returns: (array of bin counts, bin-minimum, min-width, #-points-outside-range)
"""
    a = np.ravel(a)               # flatten any >1D arrays
    if (defaultlimits is not None):
        lowerreallimit = defaultlimits[0]
        upperreallimit = defaultlimits[1]
        binsize = (upperreallimit-lowerreallimit) / float(numbins)
    else:
        Min = a.min()
        Max = a.max()
        estbinwidth = float(Max - Min)/float(numbins - 1)
        binsize = (Max-Min+estbinwidth)/float(numbins)
        lowerreallimit = Min - binsize/2.0  #lower real limit,1st bin
    bins = zeros(numbins)
    extrapoints = 0
    for num in a:
        try:
            if (num-lowerreallimit) < 0:
                extrapoints += 1
            else:
                bintoincrement = int((num-lowerreallimit) / float(binsize))
                bins[bintoincrement] = bins[bintoincrement] + 1
        except:                           # point outside lower/upper limits
            extrapoints += 1
    if extrapoints > 0 and printextras:
        # fixme: warnings.warn()
        print '\nPoints outside given histogram range =',extrapoints
    return (bins, lowerreallimit, binsize, extrapoints)


def cumfreq(a, numbins=10, defaultreallimits=None):
    """
Returns a cumulative frequency histogram, using the histogram function.
Defaultreallimits can be None (use all data), or a 2-sequence containing
lower and upper limits on values to include.

Returns: array of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(a,numbins,defaultreallimits)
    cumhist = np.cumsum(h*1, axis=0)
    return cumhist,l,b,e


def relfreq(a, numbins=10, defaultreallimits=None):
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
######  VARIABILITY FUNCTIONS  #####
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
    n = zeros(k)
    v = zeros(k)
    m = zeros(k)
    nargs = []
    for i in range(k):
        nargs.append(args[i].astype(float))
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
    if check != 1:
        raise ValueError, 'Lack of convergence in obrientransform.'
    else:
        return array(nargs)


def samplevar(a, axis=0):
    """
Returns the sample standard deviation of the values in the passed
array (i.e., using N).  Axis can equal None (ravel array first),
an integer (the axis over which to operate)
"""
    a, axis = _chk_asarray(a, axis)
    mn = np.expand_dims(mean(a, axis), axis)
    deviations = a - mn
    n = a.shape[axis]
    svar = ss(deviations,axis) / float(n)
    return svar


def samplestd(a, axis=0):
    """Returns the sample standard deviation of the values in the passed
array (i.e., using N).  Axis can equal None (ravel array first),
an integer (the axis over which to operate).
"""
    return np.sqrt(samplevar(a,axis))


def signaltonoise(instack, axis=0):
    """
Calculates signal-to-noise.  Axis can equal None (ravel array
first), an integer (the axis over which to operate).

Returns: array containing the value of (mean/stdev) along axis,
         or 0 when stdev=0
"""
    m = mean(instack,axis)
    sd = samplestd(instack,axis)
    return np.where(sd == 0, 0, m/sd)

def var(a, axis=0, bias=False):
    """
Returns the estimated population variance of the values in the passed
array (i.e., N-1).  Axis can equal None (ravel array first), or an
integer (the axis over which to operate).
"""
    warnings.warn("""\
scipy.stats.var is deprecated; please update your code to use numpy.var.
Please note that:
    - numpy.var axis argument defaults to None, not 0
    - numpy.var has a ddof argument to replace bias in a more general manner.
      scipy.stats.var(a, bias=True) can be replaced by numpy.var(x,
axis=0, ddof=1).""", DeprecationWarning)
    a, axis = _chk_asarray(a, axis)
    mn = np.expand_dims(mean(a,axis),axis)
    deviations = a - mn
    n = a.shape[axis]
    vals = sum(abs(deviations)**2,axis)/(n-1.0)
    if bias:
        return vals * (n-1.0)/n
    else:
        return vals

def std(a, axis=0, bias=False):
    """
Returns the estimated population standard deviation of the values in
the passed array (i.e., N-1).  Axis can equal None (ravel array
first), or an integer (the axis over which to operate).
"""
    warnings.warn("""\
scipy.stats.std is deprecated; please update your code to use numpy.std.
Please note that:
    - numpy.std axis argument defaults to None, not 0
    - numpy.std has a ddof argument to replace bias in a more general manner.
      scipy.stats.std(a, bias=True) can be replaced by numpy.std(x,
axis=0, ddof=1).""", DeprecationWarning)
    return np.sqrt(var(a,axis,bias))


def stderr(a, axis=0):
    """
Returns the estimated population standard error of the values in the
passed array (i.e., N-1).  Axis can equal None (ravel array
first), or an integer (the axis over which to operate).
"""
    a, axis = _chk_asarray(a, axis)
    return std(a,axis) / float(np.sqrt(a.shape[axis]))


def sem(a, axis=0):
    """
Returns the standard error of the mean (i.e., using N) of the values
in the passed array.  Axis can equal None (ravel array first), or an
integer (the axis over which to operate)
"""
    a, axis = _chk_asarray(a, axis)
    n = a.shape[axis]
    s = samplestd(a,axis) / np.sqrt(n-1)
    return s


def z(a, score):
    """
Returns the z-score of a given input score, given thearray from which
that score came.  Not appropriate for population calculations, nor for
arrays > 1D.

"""
    z = (score-mean(a,None)) / samplestd(a)
    return z


def zs(a):
    """
Returns a 1D array of z-scores, one for each score in the passed array,
computed relative to the passed array.

"""
    mu = mean(a,None)
    sigma = samplestd(a)
    return (array(a)-mu)/sigma

def zmap(scores, compare, axis=0):
    """
Returns an array of z-scores the shape of scores (e.g., [x,y]), compared to
array passed to compare (e.g., [time,x,y]).  Assumes collapsing over dim 0
of the compare array.

"""
    mns = mean(compare,axis)
    sstd = samplestd(compare,0)
    return (scores - mns) / sstd


#####################################
#######  TRIMMING FUNCTIONS  #######
#####################################

def threshold(a, threshmin=None, threshmax=None, newval=0):
    """Clip array to a given value.

Similar to numpy.clip(), except that values less than threshmin or
greater than threshmax are replaced by newval, instead of by
threshmin and threshmax respectively.

Returns: a, with values less than threshmin or greater than threshmax
         replaced with newval

"""
    a = asarray(a).copy()
    mask = zeros(a.shape, dtype=bool)
    if threshmin is not None:
        mask |= (a < threshmin)
    if threshmax is not None:
        mask |= (a > threshmax)
    a[mask] = newval
    return a


def trimboth(a, proportiontocut):
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


def trim1(a, proportiontocut, tail='right'):
    """
    Slices off the passed proportion of items from ONE end of the passed
    array (i.e., if proportiontocut=0.1, slices off 'leftmost' or 'rightmost'
    10% of scores).  Slices off LESS if proportion results in a non-integer
    slice index (i.e., conservatively slices off proportiontocut).

    Returns: trimmed version of array a
    """
    a = asarray(a)
    if tail.lower() == 'right':
        lowercut = 0
        uppercut = len(a) - int(proportiontocut*len(a))
    elif tail.lower() == 'left':
        lowercut = int(proportiontocut*len(a))
        uppercut = len(a)
    return a[lowercut:uppercut]

def trim_mean(a, proportiontocut):
    """Return mean with proportiontocut chopped from each of the lower and
    upper tails.
    """
    newa = trimboth(np.sort(a),proportiontocut)
    return mean(newa,axis=0)



#####################################
#####  CORRELATION FUNCTIONS  ######
#####################################

#  Cov is more flexible than the original
#    covariance and computes an unbiased covariance matrix
#    by default.
def cov(m, y=None, rowvar=False, bias=False):
    """Estimate the covariance matrix.

    If m is a vector, return the variance.  For matrices where each row
    is an observation, and each column a variable, return the covariance
    matrix.  Note that in this case diag(cov(m)) is a vector of
    variances for each column.

    cov(m) is the same as cov(m, m)

    Normalization is by (N-1) where N is the number of observations
    (unbiased estimate).  If bias is True then normalization is by N.

    If rowvar is False, then each row is a variable with
    observations in the columns.
    """
    warnings.warn("""\
scipy.stats.cov is deprecated; please update your code to use numpy.cov.
Please note that:
    - numpy.cov rowvar argument defaults to true, not false
    - numpy.cov bias argument defaults to false, not true
""", DeprecationWarning)
    m = asarray(m)
    if y is None:
        y = m
    else:
        y = asarray(y)
    if rowvar:
        m = np.transpose(m)
        y = np.transpose(y)
    N = m.shape[0]
    if (y.shape[0] != N):
        raise ValueError, "x and y must have the same number of observations."
    m = m - mean(m,axis=0)
    y = y - mean(y,axis=0)
    if bias:
        fact = N*1.0
    else:
        fact = N-1.0
    val = np.squeeze(np.dot(np.transpose(m),np.conjugate(y))) / fact
    return val

def corrcoef(x, y=None, rowvar=False, bias=True):
    """The correlation coefficients formed from 2-d array x, where the
    rows are the observations, and the columns are variables.

    corrcoef(x,y) where x and y are 1d arrays is the same as
    corrcoef(transpose([x,y]))

    If rowvar is True, then each row is a variables with
    observations in the columns.
    """
    warnings.warn("""\
scipy.stats.corrcoef is deprecated; please update your code to use numpy.corrcoef.
Please note that:
    - numpy.corrcoef rowvar argument defaults to true, not false
    - numpy.corrcoef bias argument defaults to false, not true
""", DeprecationWarning)
    if y is not None:
        x = np.transpose([x,y])
        y = None
    c = cov(x, y, rowvar=rowvar, bias=bias)
    d = np.diag(c)
    return c/np.sqrt(np.multiply.outer(d,d))



def f_oneway(*args):
    """
Performs a 1-way ANOVA, returning an F-value and probability given
any number of groups.  From Heiman, pp.394-7.

Usage:   f_oneway (*args)    where *args is 2 or more arrays, one per
                                  treatment group
Returns: f-value, probability
"""
    na = len(args)            # ANOVA on 'na' groups, each in it's own array
    tmp = map(np.array,args)
    alldata = np.concatenate(args)
    bign = len(alldata)
    sstot = ss(alldata)-(square_of_sums(alldata)/float(bign))
    ssbn = 0
    for a in args:
        ssbn = ssbn + square_of_sums(array(a))/float(len(a))
    ssbn = ssbn - (square_of_sums(alldata)/float(bign))
    sswn = sstot-ssbn
    dfbn = na-1
    dfwn = bign - na
    msb = ssbn/float(dfbn)
    msw = sswn/float(dfwn)
    f = msb/msw
    prob = fprob(dfbn,dfwn,f)
    return f, prob



def pearsonr(x, y):
    """Calculates a Pearson correlation coefficient and the p-value for testing
    non-correlation.

    The Pearson correlation coefficient measures the linear relationship
    between two datasets. Strictly speaking, Pearson's correlation requires
    that each dataset be normally distributed. Like other correlation
    coefficients, this one varies between -1 and +1 with 0 implying no
    correlation. Correlations of -1 or +1 imply an exact linear
    relationship. Positive correlations imply that as x increases, so does
    y. Negative correlations imply that as x increases, y decreases.

    The p-value roughly indicates the probability of an uncorrelated system
    producing datasets that have a Pearson correlation at least as extreme
    as the one computed from these datasets. The p-values are not entirely
    reliable but are probably reasonable for datasets larger than 500 or so.

    Parameters
    ----------
    x : 1D array
    y : 1D array the same length as x

    Returns
    -------
    (Pearson's correlation coefficient,
     2-tailed p-value)

    References
    ----------
    http://www.statsoft.com/textbook/glosp.html#Pearson%20Correlation
    """
    # x and y should have same length.
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    mx = x.mean()
    my = y.mean()
    xm, ym = x-mx, y-my
    r_num = n*(np.add.reduce(xm*ym))
    r_den = n*np.sqrt(ss(xm)*ss(ym))
    r = (r_num / r_den)

    # Presumably, if r > 1, then it is only some small artifact of floating
    # point arithmetic.
    r = min(r, 1.0)
    df = n-2

    # Use a small floating point value to prevent divide-by-zero nonsense
    # fixme: TINY is probably not the right value and this is probably not
    # the way to be robust. The scheme used in spearmanr is probably better.
    TINY = 1.0e-20
    t = r*np.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    return r,prob


def spearmanr(x, y):
    """Calculates a Spearman rank-order correlation coefficient and the p-value
    to test for non-correlation.

    The Spearman correlation is a nonparametric measure of the linear
    relationship between two datasets. Unlike the Pearson correlation, the
    Spearman correlation does not assume that both datasets are normally
    distributed. Like other correlation coefficients, this one varies
    between -1 and +1 with 0 implying no correlation. Correlations of -1 or
    +1 imply an exact linear relationship. Positive correlations imply that
    as x increases, so does y. Negative correlations imply that as x
    increases, y decreases.

    The p-value roughly indicates the probability of an uncorrelated system
    producing datasets that have a Spearman correlation at least as extreme
    as the one computed from these datasets. The p-values are not entirely
    reliable but are probably reasonable for datasets larger than 500 or so.

    Parameters
    ----------
    x : 1D array
    y : 1D array the same length as x
        The lengths of both arrays must be > 2.

    Returns
    -------
    (Spearman correlation coefficient,
     2-tailed p-value)

    References
    ----------
    [CRCProbStat2000] section 14.7
    """
    x = np.asanyarray(x)
    y = np.asanyarray(y)
    n = len(x)
    m = len(y)
    if n != m:
        raise ValueError("lengths of x and y must match: %s != %s" % (n, m))
    if n <= 2:
        raise ValueError("length must be > 2")
    rankx = rankdata(x)
    ranky = rankdata(y)
    dsq = np.add.reduce((rankx-ranky)**2)
    rs = 1 - 6*dsq / float(n*(n**2-1))
    df = n-2

    try:
        t = rs * np.sqrt((n-2) / ((rs+1.0)*(1.0-rs)))
        probrs = betai(0.5*df, 0.5, df/(df+t*t))
    except ZeroDivisionError:
        probrs = 0.0

    return rs, probrs


def pointbiserialr(x, y):
    # comment: I am changing the semantics somewhat. The original function is
    # fairly general and accepts an x sequence that has any type of thing in it as
    # along as there are only two unique items. I am going to restrict this to
    # a boolean array for my sanity.
    """Calculates a point biserial correlation coefficient and the associated
    p-value.

    The point biserial correlation is used to measure the relationship
    between a binary variable, x, and a continuous variable, y. Like other
    correlation coefficients, this one varies between -1 and +1 with 0
    implying no correlation. Correlations of -1 or +1 imply a determinative
    relationship.

    Parameters
    ----------
    x : array of bools
    y : array of floats

    Returns
    -------
    (point-biserial r,
     2-tailed p-value)

    References
    ----------
    http://www.childrens-mercy.org/stats/definitions/biserial.htm
    """

    ## Test data: http://support.sas.com/ctx/samples/index.jsp?sid=490&tab=output
    # x = [1,0,1,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1]
    # y = [14.8,13.8,12.4,10.1,7.1,6.1,5.8,4.6,4.3,3.5,3.3,3.2,3.0,2.8,2.8,2.5,
    #      2.4,2.3,2.1,1.7,1.7,1.5,1.3,1.3,1.2,1.2,1.1,0.8,0.7,0.6,0.5,0.2,0.2,
    #      0.1]
    # rpb = 0.36149

    x = np.asarray(x, dtype=bool)
    y = np.asarray(y, dtype=float)
    n = len(x)

    # phat is the fraction of x values that are True
    phat = x.sum() / float(len(x))
    y0 = y[~x]  # y-values where x is False
    y1 = y[x]  # y-values where x is True
    y0m = y0.mean()
    y1m = y1.mean()

    rpb = (y1m - y0m)*np.sqrt(phat * (1-phat)) / y.std()

    df = n-2
    # fixme: see comment about TINY in pearsonr()
    TINY = 1e-20
    t = rpb*np.sqrt(df/((1.0-rpb+TINY)*(1.0+rpb+TINY)))
    prob = betai(0.5*df, 0.5, df/(df+t*t))
    return rpb, prob


def kendalltau(x, y):
    """Calculates Kendall's tau, a correlation measure for ordinal data, and an
    associated p-value.

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
    tau = iss / np.sqrt(float(n1*n2))
    svar = (4.0*len(x)+10.0) / (9.0*len(x)*(len(x)-1))
    z = tau / np.sqrt(svar)
    prob = erfc(abs(z)/1.4142136)
    return tau, prob


def linregress(*args):
    """Calculates a regression line on two arrays, x and y, corresponding to
    x,y pairs.  If a single 2D array is passed, linregress finds dim with 2
    levels and splits data into x,y pairs along that dim.

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
    r_num = np.add.reduce(xm*ym)
    r_den = np.sqrt(ss(xm)*ss(ym))
    if r_den == 0.0:
        r = 0.0
    else:
        r = r_num / r_den
        if (r > 1.0): r = 1.0 # from numerical error
    #z = 0.5*log((1.0+r+TINY)/(1.0-r+TINY))
    df = n-2
    t = r*np.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    slope = r_num / ss(xm)
    intercept = ymean - slope*xmean
    sterrest = np.sqrt((1-r*r)*ss(y) / ss(x) / df)
    return slope, intercept, r, prob, sterrest


#####################################
#####  INFERENTIAL STATISTICS  #####
#####################################

def ttest_1samp(a, popmean):
    """
Calculates the t-obtained for the independent samples T-test on ONE group
of scores a, given a population mean.

Returns: t-value, two-tailed prob
"""
    a = asarray(a)
    x = mean(a,None)
    v = var(a,None)
    n = len(a)
    df = n-1
    svar = ((n-1)*v) / float(df)
    t = (x-popmean)/np.sqrt(svar*(1.0/n))
    prob = betai(0.5*df,0.5,df/(df+t*t))

    return t,prob


def ttest_ind(a, b, axis=0):
    """Calculates the t-obtained T-test on TWO INDEPENDENT samples of scores
    a, and b.  From Numerical Recipies, p.483. Axis can equal None (ravel
    array first), or an integer (the axis over which to operate on a and b).

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
    zerodivproblem = svar == 0
    t = (x1-x2)/np.sqrt(svar*(1.0/n1 + 1.0/n2))  # N-D COMPUTATION HERE!!!!!!
    t = np.where(zerodivproblem, 1.0, t)           # replace NaN t-values with 1.0
    probs = betai(0.5*df,0.5,float(df)/(df+t*t))

    if not np.isscalar(t):
        probs = probs.reshape(t.shape)
    if not np.isscalar(probs) and len(probs) == 1:
        probs = probs[0]
    return t, probs


def ttest_rel(a,b,axis=None):
    """Calculates the t-obtained T-test on TWO RELATED samples of scores, a
    and b.  From Numerical Recipies, p.483. Axis can equal None (ravel array
    first), or an integer (the axis over which to operate on a and b).

    Returns: t-value, two-tailed p-value
    """
    a, b, axis = _chk2_asarray(a, b, axis)
    if len(a)!=len(b):
        raise ValueError, 'unequal length arrays'
    x1 = mean(a,axis)
    x2 = mean(b,axis)
    v1 = var(a,axis)
    v2 = var(b,axis)
    n = a.shape[axis]
    df = float(n-1)
    d = (a-b).astype('d')

    denom = np.sqrt((n*np.add.reduce(d*d,axis) - np.add.reduce(d,axis)**2) /df)
    zerodivproblem = denom == 0
    t = np.add.reduce(d, axis) / denom      # N-D COMPUTATION HERE!!!!!!
    t = np.where(zerodivproblem, 1.0, t)    # replace NaN t-values with 1.0
    t = np.where(zerodivproblem, 1.0, t)    # replace NaN t-values with 1.0
    probs = betai(0.5*df,0.5,float(df)/(df+t*t))
    if not np.isscalar(t):
        probs = np.reshape(probs, t.shape)
    if not np.isscalar(probs) and len(probs) == 1:
        probs = probs[0]
    return t, probs


import scipy.stats
import distributions
def kstest(rvs, cdf, args=(), N=20):
    """Return the D-value and the p-value for a Kolmogorov-Smirnov test of
    the null that N RV's generated by the rvs fits the cdf given the extra
    arguments.  rvs needs to accept the size= keyword if a function.  rvs
    can also be a vector of RVs.

    cdf can be a function or a string indicating the distribution type.

    if the p-value is greater than the significance level (say 5%), then we
    cannot reject the hypothesis that the data come from the given
    distribution.
    """
    if isinstance(rvs, basestring):
        cdf = getattr(scipy.stats, rvs).cdf
        rvs = getattr(scipy.stats, rvs).rvs
    if isinstance(cdf, basestring):
        cdf = getattr(scipy.stats, cdf).cdf
    if callable(rvs):
        kwds = {'size':N}
        vals = np.sort(rvs(*args,**kwds))
    else:
        vals = np.sort(rvs)
        N = len(vals)
    cdfvals = cdf(vals, *args)
    D1 = np.absolute(cdfvals - np.arange(1.0, N+1)/N).max()
#    D2 = amax(abs(cdfvals - arange(0.0,N)/N))
#    D = max(D1,D2)
    D = D1
    return D, distributions.ksone.sf(D,N)

def chisquare(f_obs, f_exp=None):
    """ Calculates a one-way chi square for array of observed frequencies
    and returns the result.  If no expected frequencies are given, the total
    N is assumed to be equally distributed across all groups.

    Returns: chisquare-statistic, associated p-value
    """

    f_obs = asarray(f_obs)
    k = len(f_obs)
    if f_exp is None:
        f_exp = array([np.sum(f_obs,axis=0)/float(k)] * len(f_obs),float)
    f_exp = f_exp.astype(float)
    chisq = np.add.reduce((f_obs-f_exp)**2 / f_exp)
    return chisq, chisqprob(chisq, k-1)


def ks_2samp(data1, data2):
    """ Computes the Kolmogorov-Smirnof statistic on 2 samples.  Modified
    from Numerical Recipies in C, page 493.  Returns KS D-value, prob.  Not
    ufunc- like.

    Returns: KS D-value, p-value
    """
    data1, data2 = map(asarray, (data1, data2))
    j1 = 0    # zeros(data1.shape[1:]) TRIED TO MAKE THIS UFUNC-LIKE
    j2 = 0    # zeros(data2.shape[1:])
    fn1 = 0.0 # zeros(data1.shape[1:],float)
    fn2 = 0.0 # zeros(data2.shape[1:],float)
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    en1 = n1*1
    en2 = n2*1
    d = zeros(data1.shape[1:])
    data1 = np.sort(data1,0)
    data2 = np.sort(data2,0)
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
        en = np.sqrt(en1*en2/float(en1+en2))
        prob = ksprob((en+0.12+0.11/en)*np.fabs(d))
    except:
        prob = 1.0
    return d, prob


def mannwhitneyu(x, y):
    """Calculates a Mann-Whitney U statistic on the provided scores and
    returns the result.  Use only when the n in each condition is < 20 and
    you have 2 independent samples of ranks.  REMEMBER: Mann-Whitney U is
    significant if the u-obtained is LESS THAN or equal to the critical
    value of U.

    Returns: u-statistic, one-tailed p-value (i.e., p(z(U)))
    """
    x = asarray(x)
    y = asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(np.concatenate((x,y)))
    rankx = ranked[0:n1]       # get the x-ranks
    #ranky = ranked[n1:]        # the rest are y-ranks
    u1 = n1*n2 + (n1*(n1+1))/2.0 - np.sum(rankx,axis=0)  # calc U for x
    u2 = n1*n2 - u1                            # remainder is U for y
    bigu = max(u1,u2)
    smallu = min(u1,u2)
    T = np.sqrt(tiecorrect(ranked))  # correction factor for tied scores
    if T == 0:
        raise ValueError, 'All numbers are identical in amannwhitneyu'
    sd = np.sqrt(T*n1*n2*(n1+n2+1)/12.0)
    z = abs((bigu-n1*n2/2.0) / sd)  # normal approximation for prob calc
    return smallu, 1.0 - zprob(z)


def tiecorrect(rankvals):
    """Tie-corrector for ties in Mann Whitney U and Kruskal Wallis H tests.
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


def ranksums(x, y):
    """Calculates the rank sums statistic on the provided scores and
    returns the result.

    Returns: z-statistic, two-tailed p-value
    """
    x,y = map(np.asarray, (x, y))
    n1 = len(x)
    n2 = len(y)
    alldata = np.concatenate((x,y))
    ranked = rankdata(alldata)
    x = ranked[:n1]
    y = ranked[n1:]
    s = np.sum(x,axis=0)
    expected = n1*(n1+n2+1) / 2.0
    z = (s - expected) / np.sqrt(n1*n2*(n1+n2+1)/12.0)
    prob = 2*(1.0 -zprob(abs(z)))
    return z, prob



def kruskal(*args):
    """The Kruskal-Wallis H-test is a non-parametric ANOVA for 2 or more
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
    args = list(args)
    for i in range(len(args)):
        args[i] = ranked[0:n[i]]
        del ranked[0:n[i]]
    rsums = []
    for i in range(len(args)):
        rsums.append(np.sum(args[i],axis=0)**2)
        rsums[i] = rsums[i] / float(n[i])
    ssbn = np.sum(rsums,axis=0)
    totaln = np.sum(n,axis=0)
    h = 12.0 / (totaln*(totaln+1)) * ssbn - 3*(totaln+1)
    df = len(args) - 1
    if T == 0:
        raise ValueError, 'All numbers are identical in kruskal'
    h = h / float(T)
    return h, chisqprob(h,df)


def friedmanchisquare(*args):
    """Friedman Chi-Square is a non-parametric, one-way within-subjects
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
    data = data.astype(float)
    for i in range(len(data)):
        data[i] = rankdata(data[i])
    ssbn = np.sum(np.sum(args,1)**2,axis=0)
    chisq = 12.0 / (k*n*(k+1)) * ssbn - 3*n*(k+1)
    return chisq, chisqprob(chisq,k-1)


#####################################
####  PROBABILITY CALCULATIONS  ####
#####################################

zprob = special.ndtr
erfc = special.erfc

def chisqprob(chisq, df):
    """Returns the (1-tail) probability value associated with the provided
    chi-square value and degrees of freedom.

    Broadcasting rules apply.

    Parameters
    ----------
    chisq : array or float > 0
    df : array or float, probably int >= 1

    Returns
    -------
    The area from chisq to infinity under the Chi^2 probability distribution
    with degrees of freedom df.
    """
    return special.chdtrc(df,chisq)

ksprob = special.kolmogorov
fprob = special.fdtrc

def betai(a, b, x):
    """Returns the incomplete beta function.

    I_x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)

    where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
    function of a.

    The standard broadcasting rules apply to a, b, and x.

    Parameters
    ----------
    a : array or float > 0
    b : array or float > 0
    x : array or float
        x will be clipped to be no greater than 1.0 .

    Returns
    -------

    """
    x = np.asarray(x)
    x = np.where(x < 1.0, x, 1.0)  # if x > 1 then return 1.0
    return special.betainc(a, b, x)

#####################################
#######  ANOVA CALCULATIONS  #######
#####################################

def glm(data, para):
    """Calculates a linear model fit ...
    anova/ancova/lin-regress/t-test/etc. Taken from:

    Peterson et al. Statistical limitations in functional neuroimaging
    I. Non-inferential methods and statistical models.  Phil Trans Royal Soc
    Lond B 354: 1239-1260.

    Returns: statistic, p-value ???
    """
    if len(para) != len(data):
        raise ValueError("data and para must be same length in aglm")
    n = len(para)
    p = _support.unique(para)
    x = zeros((n,len(p)))  # design matrix
    for l in range(len(p)):
        x[:,l] = para == p[l]
    # fixme: normal equations are bad. Use linalg.lstsq instead.
    b = dot(dot(linalg.inv(dot(np.transpose(x),x)),  # i.e., b=inv(X'X)X'Y
                    np.transpose(x)),data)
    diffs = (data - dot(x,b))
    s_sq = 1./(n-len(p)) * dot(np.transpose(diffs), diffs)

    if len(p) == 2:  # ttest_ind
        c = array([1,-1])
        df = n-2
        fact = np.sum(1.0/np.sum(x,0),axis=0)  # i.e., 1/n1 + 1/n2 + 1/n3 ...
        t = dot(c,b) / np.sqrt(s_sq*fact)
        probs = betai(0.5*df,0.5,float(df)/(df+t*t))
        return t, probs
    else:
        raise ValueError("only ttest_ind implemented")


def f_value_wilks_lambda(ER, EF, dfnum, dfden, a, b):
    """Calculation of Wilks lambda F-statistic for multivarite data, per
    Maxwell & Delaney p.657.
    """
    if isinstance(ER, (int, float)):
        ER = array([[ER]])
    if isinstance(EF, (int, float)):
        EF = array([[EF]])
    lmbda = linalg.det(EF) / linalg.det(ER)
    if (a-1)**2 + (b-1)**2 == 5:
        q = 1
    else:
        q = np.sqrt( ((a-1)**2*(b-1)**2 - 2) / ((a-1)**2 + (b-1)**2 -5) )
    n_um = (1 - lmbda**(1.0/q))*(a-1)*(b-1)
    d_en = lmbda**(1.0/q) / (n_um*q - 0.5*(a-1)*(b-1) + 1)
    return n_um / d_en

def f_value(ER, EF, dfR, dfF):
    """Returns an F-statistic given the following:
        ER  = error associated with the null hypothesis (the Restricted model)
        EF  = error associated with the alternate hypothesis (the Full model)
        dfR = degrees of freedom the Restricted model
        dfF = degrees of freedom associated with the Restricted model
    """
    return ((ER-EF)/float(dfR-dfF) / (EF/float(dfF)))



def f_value_multivariate(ER, EF, dfnum, dfden):
    """Returns an F-statistic given the following:
        ER  = error associated with the null hypothesis (the Restricted model)
        EF  = error associated with the alternate hypothesis (the Full model)
        dfR = degrees of freedom the Restricted model
        dfF = degrees of freedom associated with the Restricted model
    where ER and EF are matrices from a multivariate F calculation.
    """
    if isinstance(ER, (int, float)):
        ER = array([[ER]])
    if isinstance(EF, (int, float)):
        EF = array([[EF]])
    n_um = (linalg.det(ER) - linalg.det(EF)) / float(dfnum)
    d_en = linalg.det(EF) / float(dfden)
    return n_um / d_en


#####################################
#######  SUPPORT FUNCTIONS  ########
#####################################

def ss(a, axis=0):
    """Squares each value in the passed array, adds these squares, and
    returns the result.

    Parameters
    ----------
    a : array
    axis : int or None

    Returns
    -------
    The sum along the given axis for (a*a).
    """
    a, axis = _chk_asarray(a, axis)
    return np.sum(a*a, axis)


def square_of_sums(a, axis=0):
    """Adds the values in the passed array, squares that sum, and returns the
result.

Returns: the square of the sum over axis.
"""
    a, axis = _chk_asarray(a, axis)
    s = np.sum(a,axis)
    if not np.isscalar(s):
        return s.astype(float)*s
    else:
        return float(s)*s


def fastsort(a):
    # fixme: the wording in the docstring is nonsense.
    """Sort an array and provide the argsort.

    Parameters
    ----------
    a : array

    Returns
    -------
    (sorted array,
     indices into the original array,
    )
    """
    it = np.argsort(a)
    as_ = a[it]
    return as_, it

def rankdata(a):
    """Ranks the data in a, dealing with ties appropriately.

    Equal values are assigned a rank that is the average of the ranks that
    would have been otherwise assigned to all of the values within that set.
    Ranks begin at 1, not 0.

    Example
    -------
    In [15]: stats.rankdata([0, 2, 2, 3])
    Out[15]: array([ 1. ,  2.5,  2.5,  4. ])

    Parameters
    ----------
    a : array
        This array is first flattened.

    Returns
    -------
    An array of length equal to the size of a, containing rank scores.
    """
    a = np.ravel(a)
    n = len(a)
    svec, ivec = fastsort(a)
    sumranks = 0
    dupcount = 0
    newarray = np.zeros(n, float)
    for i in xrange(n):
        sumranks += i
        dupcount += 1
        if i==n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in xrange(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newarray
