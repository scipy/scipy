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
# Comments and/or additions are welcome (send e-mail to:
# strang@nmr.mgh.harvard.edu).
# 
"""
stats.py module

(Requires pstat.py modules.)

#################################################
#######  Written by:  Gary Strangman  ###########
#######  Last modified:  Apr 13, 2000 ###########
#################################################

A collection of basic statistical functions for python.  The function
names appear below.


IMPORTANT:  There are really *3* sets of functions.  The first set has an 'l'
prefix, which can be used with list or tuple arguments.  The second set has
an 'a' prefix, which can accept NumPy array arguments.  These latter
functions are defined only when NumPy is available on the system.  The third
type has NO prefix (i.e., has the name that appears below).  Functions of
this set are members of a "Dispatch" class, c/o David Ascher.  This class
allows different functions to be called depending on the type of the passed
arguments.  Thus, stats.mean is a member of the Dispatch class and
stats.mean(range(20)) will call stats.lmean(range(20)) while
stats.mean(Numeric.arange(20)) will call stats.amean(Numeric.arange(20)).
This is a handy way to keep consistent function names when different
argument types require different functions to be called.  Having
implementated the Dispatch class, however, means that to get info on
a given function, you must use the REAL function name ... that is
"print stats.lmean.__doc__" or "print stats.amean.__doc__" work fine,
while "print stats.mean.__doc__" will print the doc for the Dispatch
class.  NUMPY FUNCTIONS ('a' prefix) generally have more argument options
but should otherwise be consistent with the corresponding list functions.

Disclaimers:  The function list is obviously incomplete and, worse, the
functions are not optimized.  All functions have been tested (some more
so than others), but they are far from bulletproof.  Thus, as with any
free software, no warranty or guarantee is expressed or implied. :-)  A
few extra functions that don't appear in the list below can be found by
interested treasure-hunters.  These functions don't necessarily have
both list and array versions but were deemed useful

CENTRAL TENDENCY:  geometricmean
                   harmonicmean
                   mean
                   median
                   medianscore
                   mode

MOMENTS:  moment
          variation
          skew
          kurtosis
          normaltest (for arrays only)

ALTERED VERSIONS:  tmean (for arrays only)
                   tvar  (for arrays only)
                   tstdev (for arrays only)
                   tsem  (for arrays only)
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
              samplestdev
              signaltonoise (for arrays only)
              var
              stdev
              sterr
              sem
              z
              zs

TRIMMING FCNS:  threshold (for arrays only)
                trimboth
                trim1
                around (round all vals to 'n' decimals; Numpy only)

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
                    wilcoxont
                    kruskalwallish
                    friedmanchisquare

PROBABILITY CALCS:  chisqprob
                    erfcc
                    zprob
                    fprob
                    betacf
                    gammln 
                    betai

  *** SEE ALSO the scipy.special package for probability calculation functions. ***

ANOVA FUNCTIONS:  anova (NumPy required)
                  F_oneway
                  F_value

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
## 00-04-13 ... pulled all "global" statements, except from aanova()
##              added/fixed lots of documentation, removed io.py dependency
##              changed to version 0.5
## 99-11-13 ... added asign() function
## 99-11-01 ... changed version to 0.4 ... enough incremental changes now
## 99-10-25 ... added acovariance and acorrelation functions
## 99-10-10 ... fixed askew/akurtosis to avoid divide-by-zero errors
##              added aglm function (crude, but will be improved)
## 99-10-04 ... upgraded acumsum, ass, asummult, asamplevar, avar, etc. to
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
##                   max/min in array section to N.maximum/N.minimum,
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
##              added 'asum' function (added functionality to N.add.reduce)
##              fixed both moment and amoment (two errors)
##              changed name of skewness and askewness to skew and askew
##              fixed (a)histogram (which sometimes counted points <lowerlimit)


import math, string, sys, pstat, copy
from types import *

__version__ = 0.5

############# DISPATCH CODE ##############


class Dispatch:
    """
The Dispatch class, care of David Ascher, allows different functions to
be called depending on the argument types.  This way, there can be one
function name regardless of the argument type.  To access function doc
in stats.py module, prefix the function with an 'l' or 'a' for list or
array arguments, respectively.  That is, print stats.lmean.__doc__ or
print stats.amean.__doc__ or whatever.
"""

    def __init__(self, *tuples):
        self._dispatch = {}
        for func, types in tuples:
            for t in types:
                if t in self._dispatch.keys():
                    raise ValueError, "can't have two dispatches on "+str(t)
                self._dispatch[t] = func
        self._types = self._dispatch.keys()

    def __call__(self, arg1, *args, **kw):
        if type(arg1) not in self._types:
            raise TypeError, "don't know how to dispatch %s arguments" %  type(arg1)
        return apply(self._dispatch[type(arg1)], (arg1,) + args, kw)


##########################################################################
########################   LIST-BASED FUNCTIONS   ########################
##########################################################################

### Define these regardless

####################################
#######  CENTRAL TENDENCY  #########
####################################

def lgeometricmean (inlist):
    """
Calculates the geometric mean of the values in the passed list.
That is:  n-th root of (x1 * x2 * ... * xn).  Assumes a '1D' list.

Usage:   lgeometricmean(inlist)
"""
    mult = 1.0
    one_over_n = 1.0/len(inlist)
    for item in inlist:
	mult = mult * pow(item,one_over_n)
    return mult


def lharmonicmean (inlist):
    """
Calculates the harmonic mean of the values in the passed list.
That is:  n / (1/x1 + 1/x2 + ... + 1/xn).  Assumes a '1D' list.

Usage:   lharmonicmean(inlist)
"""
    sum = 0
    for item in inlist:
	sum = sum + 1.0/item
    return len(inlist) / sum


def lmean (inlist):
    """
Returns the arithematic mean of the values in the passed list.
Assumes a '1D' list, but will function on the 1st dim of an array(!).

Usage:   lmean(inlist)
"""
    sum = 0
    for item in inlist:
	sum = sum + item
    return sum/float(len(inlist))


def lmedian (inlist,numbins=1000):
    """
Returns the computed median value of a list of numbers, given the
number of bins to use for the histogram (more bins brings the computed value
closer to the median score, default number of bins = 1000).  See G.W.
Heiman's Basic Stats (1st Edition), or CRC Probability & Statistics.

Usage:   lmedian (inlist, numbins=1000)
"""
    (hist, smallest, binsize, extras) = histogram(inlist,numbins) # make histog
    cumhist = cumsum(hist)              # make cumulative histogram
    for i in range(len(cumhist)):	# get 1st(!) index holding 50%ile score
	if cumhist[i]>=len(inlist)/2.0:
	    cfbin = i
	    break
    LRL = smallest + binsize*cfbin	# get lower read limit of that bin
    cfbelow = cumhist[cfbin-1]
    freq = float(hist[cfbin])		# frequency IN the 50%ile bin
    median = LRL + ((len(inlist)/2.0 - cfbelow)/float(freq))*binsize  # median formula
    return median


def lmedianscore (inlist):
    """
Returns the 'middle' score of the passed list.  If there is an even
number of scores, the mean of the 2 middle scores is returned.

Usage:   lmedianscore(inlist)
"""

    newlist = copy.deepcopy(inlist)
    newlist.sort()
    if len(newlist) % 2 == 0:   # if even number of scores, average middle 2
	index = len(newlist)/2  # integer division correct
	median = float(newlist[index] + newlist[index-1]) /2
    else:
	index = len(newlist)/2  # int divsion gives mid value when count from 0
	median = newlist[index]
    return median


def lmode(inlist):
    """
Returns a list of the modal (most common) score(s) in the passed
list.  If there is more than one such score, all are returned.  The
bin-count for the mode(s) is also returned.

Usage:   lmode(inlist)
Returns: bin-count for mode(s), a list of modal value(s)
"""

    scores = pstat.unique(inlist)
    scores.sort()
    freq = []
    for item in scores:
	freq.append(inlist.count(item))
    maxfreq = max(freq)
    mode = []
    stillmore = 1
    while stillmore:
	try:
	    indx = freq.index(maxfreq)
	    mode.append(scores[indx])
	    del freq[indx]
	    del scores[indx]
	except ValueError:
	    stillmore=0
    return maxfreq, mode


####################################
############  MOMENTS  #############
####################################

def lmoment(inlist,moment=1):
    """
Calculates the nth moment about the mean for a sample (defaults to
the 1st moment).  Used to calculate coefficients of skewness and kurtosis.

Usage:   lmoment(inlist,moment=1)
Returns: appropriate moment (r) from ... 1/n * SUM((inlist(i)-mean)**r)
"""
    if moment == 1:
	return 0.0
    else:
	mn = mean(inlist)
	n = len(inlist)
	s = 0
	for x in inlist:
	    s = s + (x-mn)**moment
	return s/float(n)


def lvariation(inlist):
    """
Returns the coefficient of variation, as defined in CRC Standard
Probability and Statistics, p.6.

Usage:   lvariation(inlist)
"""
    return 100.0*samplestdev(inlist)/float(mean(inlist))


def lskew(inlist):
    """
Returns the skewness of a distribution, as defined in Numerical
Recipies (alternate defn in CRC Standard Probability and Statistics, p.6.)

Usage:   lskew(inlist)
"""
    return moment(inlist,3)/pow(moment(inlist,2),1.5)


def lkurtosis(inlist):
    """
Returns the kurtosis of a distribution, as defined in Numerical
Recipies (alternate defn in CRC Standard Probability and Statistics, p.6.)

Usage:   lkurtosis(inlist)
"""
    return moment(inlist,4)/pow(moment(inlist,2),2.0)


def ldescribe(inlist):
    """
Returns some descriptive statistics of the passed list (assumed to be 1D).

Usage:   ldescribe(inlist)
Returns: n, mean, standard deviation, skew, kurtosis
"""
    n = len(inlist)
    mm = (min(inlist),max(inlist))
    m = mean(inlist)
    sd = stdev(inlist)
    sk = skew(inlist)
    kurt = kurtosis(inlist)
    return n, mm, m, sd, sk, kurt


####################################
#######  FREQUENCY STATS  ##########
####################################

def litemfreq(inlist):
    """
Returns a list of pairs.  Each pair consists of one of the scores in inlist
and it's frequency count.  Assumes a 1D list is passed.

Usage:   litemfreq(inlist)
Returns: a 2D frequency table (col [0:n-1]=scores, col n=frequencies)
"""
    scores = pstat.unique(inlist)
    scores.sort()
    freq = []
    for item in scores:
	freq.append(inlist.count(item))
    return pstat.abut(scores, freq)


def lscoreatpercentile (inlist, percent):
    """
Returns the score at a given percentile relative to the distribution
given by inlist.

Usage:   lscoreatpercentile(inlist,percent)
"""
    if percent > 1:
	print "\nDividing percent>1 by 100 in lscoreatpercentile().\n"
	percent = percent / 100.0
    targetcf = percent*len(inlist)
    h, lrl, binsize, extras = histogram(inlist)
    cumhist = cumsum(copy.deepcopy(h))
    for i in range(len(cumhist)):
	if cumhist[i] >= targetcf:
	    break
    score = binsize * ((targetcf - cumhist[i-1]) / float(h[i])) + (lrl+binsize*i)
    return score


def lpercentileofscore (inlist, score,histbins=10,defaultlimits=None):
    """
Returns the percentile value of a score relative to the distribution
given by inlist.  Formula depends on the values used to histogram the data(!).

Usage:   lpercentileofscore(inlist,score,histbins=10,defaultlimits=None)
"""

    h, lrl, binsize, extras = histogram(inlist,histbins,defaultlimits)
    cumhist = cumsum(copy.deepcopy(h))
    i = int((score - lrl)/float(binsize))
    pct = (cumhist[i-1]+((score-(lrl+binsize*i))/float(binsize))*h[i])/float(len(inlist)) * 100
    return pct


def lhistogram (inlist,numbins=10,defaultreallimits=None,printextras=0):
    """
Returns (i) a list of histogram bin counts, (ii) the smallest value
of the histogram binning, and (iii) the bin width (the last 2 are not
necessarily integers).  Default number of bins is 10.  If no sequence object
is given for defaultreallimits, the routine picks (usually non-pretty) bins
spanning all the numbers in the inlist.

Usage:   lhistogram (inlist, numbins=10, defaultreallimits=None,suppressoutput=0)
Returns: list of bin values, lowerreallimit, binsize, extrapoints
"""
    if (defaultreallimits <> None):
	if type(defaultreallimits) not in [ListType,TupleType] or len(defaultreallimits)==1: # only one limit given, assumed to be lower one & upper is calc'd
	    lowerreallimit = defaultreallimits
	    upperreallimit = 1.0001 * max(inlist)
	else: # assume both limits given
	    lowerreallimit = defaultreallimits[0]
	    upperreallimit = defaultreallimits[1]
	binsize = (upperreallimit-lowerreallimit)/float(numbins)
    else:     # no limits given for histogram, both must be calc'd
	estbinwidth=(max(inlist)-min(inlist))/float(numbins) + 1 # 1=>cover all
	binsize = ((max(inlist)-min(inlist)+estbinwidth))/float(numbins)
	lowerreallimit = min(inlist) - binsize/2 #lower real limit,1st bin
    bins = [0]*(numbins)
    extrapoints = 0
    for num in inlist:
	try:
	    if (num-lowerreallimit) < 0:
		extrapoints = extrapoints + 1
	    else:
		bintoincrement = int((num-lowerreallimit)/float(binsize))
		bins[bintoincrement] = bins[bintoincrement] + 1
	except:
	    extrapoints = extrapoints + 1
    if (extrapoints > 0 and printextras == 1):
	print '\nPoints outside given histogram range =',extrapoints
    return (bins, lowerreallimit, binsize, extrapoints)


def lcumfreq(inlist,numbins=10,defaultreallimits=None):
    """
Returns a cumulative frequency histogram, using the histogram function.

Usage:   lcumfreq(inlist,numbins=10,defaultreallimits=None)
Returns: list of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(inlist,numbins,defaultreallimits)
    cumhist = cumsum(copy.deepcopy(h))
    return cumhist,l,b,e


def lrelfreq(inlist,numbins=10,defaultreallimits=None):
    """
Returns a relative frequency histogram, using the histogram function.

Usage:   lrelfreq(inlist,numbins=10,defaultreallimits=None)
Returns: list of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(inlist,numbins,defaultreallimits)
    for i in range(len(h)):
	h[i] = h[i]/float(len(inlist))
    return h,l,b,e


####################################
#####  VARIABILITY FUNCTIONS  ######
####################################

def lobrientransform(*args):
    """
Computes a transform on input data (any number of columns).  Used to
test for homogeneity of variance prior to running one-way stats.  From
Maxwell and Delaney, p.112.

Usage:   lobrientransform(*args)
Returns: transformed data for use in an ANOVA
"""
    TINY = 1e-10
    k = len(args)
    n = [0.0]*k
    v = [0.0]*k
    m = [0.0]*k
    nargs = []
    for i in range(k):
	nargs.append(copy.deepcopy(args[i]))
	n[i] = float(len(nargs[i]))
	v[i] = var(nargs[i])
	m[i] = mean(nargs[i])
    for j in range(k):
	for i in range(n[j]):
	    t1 = (n[j]-1.5)*n[j]*(nargs[j][i]-m[j])**2
	    t2 = 0.5*v[j]*(n[j]-1.0)
	    t3 = (n[j]-1.0)*(n[j]-2.0)
	    nargs[j][i] = (t1-t2) / float(t3)
    check = 1
    for j in range(k):
	if v[j] - mean(nargs[j]) > TINY:
	    check = 0
    if check <> 1:
	raise ValueError, 'Problem in obrientransform.'
    else:
	return nargs


def lsamplevar (inlist):
    """
Returns the variance of the values in the passed list using
N for the denominator (i.e., DESCRIBES the sample variance only).

Usage:   lsamplevar(inlist)
"""
    n = len(inlist)
    mn = mean(inlist)
    deviations = []
    for item in inlist:
	deviations.append(item-mn)
    return ss(deviations)/float(n)


def lsamplestdev (inlist):
    """
Returns the standard deviation of the values in the passed list using
N for the denominator (i.e., DESCRIBES the sample stdev only).

Usage:   lsamplestdev(inlist)
"""
    return math.sqrt(samplevar(inlist))


def lvar (inlist):
    """
Returns the variance of the values in the passed list using N-1
for the denominator (i.e., for estimating population variance).

Usage:   lvar(inlist)
"""
    n = len(inlist)
    mn = mean(inlist)
    deviations = [0]*len(inlist)
    for i in range(len(inlist)):
	deviations[i] = inlist[i] - mn
    return ss(deviations)/float(n-1)


def lstdev (inlist):
    """
Returns the standard deviation of the values in the passed list
using N-1 in the denominator (i.e., to estimate population stdev).

Usage:   lstdev(inlist)
"""
    return math.sqrt(var(inlist))


def lsterr(inlist):
    """
Returns the standard error of the values in the passed list using N-1
in the denominator (i.e., to estimate population standard error).

Usage:   lsterr(inlist)
"""
    return stdev(inlist) / float(math.sqrt(len(inlist)))


def lsem (inlist):
    """
Returns the estimated standard error of the mean (sx-bar) of the
values in the passed list.  sem = stdev / sqrt(n)

Usage:   lsem(inlist)
"""
    sd = stdev(inlist)
    n = len(inlist)
    return sd/math.sqrt(n)


def lz (inlist, score):
    """
Returns the z-score for a given input score, given that score and the
list from which that score came.  Not appropriate for population calculations.

Usage:   lz(inlist, score)
"""
    z = (score-mean(inlist))/samplestdev(inlist)
    return z


def lzs (inlist):
    """
Returns a list of z-scores, one for each score in the passed list.

Usage:   lzs(inlist)
"""
    zscores = []
    for item in inlist:
	zscores.append(z(inlist,item))
    return zscores


####################################
#######  TRIMMING FUNCTIONS  #######
####################################

def ltrimboth (l,proportiontocut):
    """
Slices off the passed proportion of items from BOTH ends of the passed
list (i.e., with proportiontocut=0.1, slices 'leftmost' 10% AND 'rightmost'
10% of scores.  Assumes list is sorted by magnitude.  Slices off LESS if
proportion results in a non-integer slice index (i.e., conservatively
slices off proportiontocut).

Usage:   ltrimboth (l,proportiontocut)
Returns: trimmed version of list l
"""
    lowercut = int(proportiontocut*len(l))
    uppercut = len(l) - lowercut
    return l[lowercut:uppercut]


def ltrim1 (l,proportiontocut,tail='right'):
    """
Slices off the passed proportion of items from ONE end of the passed
list (i.e., if proportiontocut=0.1, slices off 'leftmost' or 'rightmost'
10% of scores).  Slices off LESS if proportion results in a non-integer
slice index (i.e., conservatively slices off proportiontocut).

Usage:   ltrim1 (l,proportiontocut,tail='right')  or set tail='left'
Returns: trimmed version of list l
"""
    if tail == 'right':
	lowercut = 0
	uppercut = len(l) - int(proportiontocut*len(l))
    elif tail == 'left':
	lowercut = int(proportiontocut*len(l))
	uppercut = len(l)
    return l[lowercut:uppercut]


####################################
#####  CORRELATION FUNCTIONS  ######
####################################

def lpaired(x,y):
    """
Interactively determines the type of data and then runs the
appropriated statistic for paired group data.

Usage:   lpaired(x,y)
Returns: appropriate statistic name, value, and probability
"""
    samples = ''
    while samples not in ['i','r','I','R','c','C']:
	print '\nIndependent or related samples, or correlation (i,r,c): ',
	samples = raw_input()

    if samples in ['i','I','r','R']:
	print '\nComparing variances ...',
# USE O'BRIEN'S TEST FOR HOMOGENEITY OF VARIANCE, Maxwell & delaney, p.112
	r = obrientransform(x,y)
	f,p = F_oneway(pstat.colex(r,0),pstat.colex(r,1))
	if p<0.05:
	    vartype='unequal, p='+str(round(p,4))
	else:
	    vartype='equal'
	print vartype
	if samples in ['i','I']:
	    if vartype[0]=='e':
		t,p = ttest_ind(x,y,0)
		print '\nIndependent samples t-test:  ', round(t,4),round(p,4)
	    else:
		if len(x)>20 or len(y)>20:
		    z,p = ranksums(x,y)
		    print '\nRank Sums test (NONparametric, n>20):  ', round(z,4),round(p,4)
		else:
		    u,p = mannwhitneyu(x,y)
		    print '\nMann-Whitney U-test (NONparametric, ns<20):  ', round(u,4),round(p,4)

	else:  # RELATED SAMPLES
	    if vartype[0]=='e':
		t,p = ttest_rel(x,y,0)
		print '\nRelated samples t-test:  ', round(t,4),round(p,4)
	    else:
		t,p = ranksums(x,y)
		print '\nWilcoxon T-test (NONparametric):  ', round(t,4),round(p,4)
    else:  # CORRELATION ANALYSIS
	corrtype = ''
	while corrtype not in ['c','C','r','R','d','D']:
	    print '\nIs the data Continuous, Ranked, or Dichotomous (c,r,d): ',
	    corrtype = raw_input()
	if corrtype in ['c','C']:
	    m,b,r,p,see = linregress(x,y)
	    print '\nLinear regression for continuous variables ...'
	    lol = [['Slope','Intercept','r','Prob','SEestimate'],[round(m,4),round(b,4),round(r,4),round(p,4),round(see,4)]]
	    pstat.printcc(lol)
	elif corrtype in ['r','R']:
	    r,p = spearmanr(x,y)
	    print '\nCorrelation for ranked variables ...'
	    print "Spearman's r: ",round(r,4),round(p,4)
	else: # DICHOTOMOUS
	    r,p = pointbiserialr(x,y)
	    print '\nAssuming x contains a dichotomous variable ...'
	    print 'Point Biserial r: ',round(r,4),round(p,4)
    print '\n\n'
    return None


def lpearsonr(x,y):
    """
Calculates a Pearson correlation coefficient and the associated
probability value.  Taken from Heiman's Basic Statistics for the Behav.
Sci (2nd), p.195.

Usage:   lpearsonr(x,y)      where x and y are equal-length lists
Returns: Pearson's r value, two-tailed p-value
"""
    TINY = 1.0e-30
    if len(x) <> len(y):
	raise ValueError, 'Input values not paired in pearsonr.  Aborting.'
    n = len(x)
    x = map(float,x)
    y = map(float,y)
    xmean = mean(x)
    ymean = mean(y)
    r_num = n*(summult(x,y)) - sum(x)*sum(y)
    r_den = math.sqrt((n*ss(x) - square_of_sums(x))*(n*ss(y)-square_of_sums(y)))
    r = (r_num / r_den)  # denominator already a float
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/float(df+t*t))
    return r, prob


def lspearmanr(x,y):
    """
Calculates a Spearman rank-order correlation coefficient.  Taken
from Heiman's Basic Statistics for the Behav. Sci (1st), p.192.

Usage:   lspearmanr(x,y)      where x and y are equal-length lists
Returns: Spearman's r, two-tailed p-value
"""
    TINY = 1e-30
    if len(x) <> len(y):
	raise ValueError, 'Input values not paired in spearmanr.  Aborting.'
    n = len(x)
    rankx = rankdata(x)
    ranky = rankdata(y)
    dsq = sumdiffsquared(rankx,ranky)
    rs = 1 - 6*dsq / float(n*(n**2-1))
    t = rs * math.sqrt((n-2) / ((rs+1.0)*(1.0-rs)))
    df = n-2
    probrs = betai(0.5*df,0.5,df/(df+t*t))  # t already a float
# probability values for rs are from part 2 of the spearman function in
# Numerical Recipies, p.510.  They are close to tables, but not exact. (?)
    return rs, probrs


def lpointbiserialr(x,y):
    """
Calculates a point-biserial correlation coefficient and the associated
probability value.  Taken from Heiman's Basic Statistics for the Behav.
Sci (1st), p.194.

Usage:   lpointbiserialr(x,y)      where x,y are equal-length lists
Returns: Point-biserial r, two-tailed p-value
"""
    TINY = 1e-30
    if len(x) <> len(y):
	raise ValueError, 'INPUT VALUES NOT PAIRED IN pointbiserialr.  ABORTING.'
    data = pstat.abut(x,y)
    categories = pstat.unique(x)
    if len(categories) <> 2:
	raise ValueError, "Exactly 2 categories required for pointbiserialr()."
    else:   # there are 2 categories, continue
	codemap = pstat.abut(categories,range(2))
	recoded = pstat.recode(data,codemap,0)
	x = pstat.linexand(data,0,categories[0])
	y = pstat.linexand(data,0,categories[1])
	xmean = mean(pstat.colex(x,1))
	ymean = mean(pstat.colex(y,1))
	n = len(data)
	adjust = math.sqrt((len(x)/float(n))*(len(y)/float(n)))
	rpb = (ymean - xmean)/samplestdev(pstat.colex(data,1))*adjust
	df = n-2
	t = rpb*math.sqrt(df/((1.0-rpb+TINY)*(1.0+rpb+TINY)))
	prob = betai(0.5*df,0.5,df/(df+t*t))  # t already a float
	return rpb, prob


def lkendalltau(x,y):
    """
Calculates Kendall's tau ... correlation of ordinal data.  Adapted
from function kendl1 in Numerical Recipies.  Needs good test-routine.@@@

Usage:   lkendalltau(x,y)
Returns: Kendall's tau, two-tailed p-value
"""
    n1 = 0
    n2 = 0
    iss = 0
    for j in range(len(x)-1):
	for k in range(j,len(y)):
	    a1 = x[j] - x[k]
	    a2 = y[j] - y[k]
	    aa = a1 * a2
	    if (aa):             # neither list has a tie
		n1 = n1 + 1
		n2 = n2 + 1
		if aa > 0:
		    iss = iss + 1
		else:
		    iss = iss -1
	    else:
		if (a1):
		    n1 = n1 + 1
		else:
		    n2 = n2 + 1
    tau = iss / math.sqrt(n1*n2)
    svar = (4.0*len(x)+10.0) / (9.0*len(x)*(len(x)-1))
    z = tau / math.sqrt(svar)
    prob = erfcc(abs(z)/1.4142136)
    return tau, prob


def llinregress(x,y):
    """
Calculates a regression line on x,y pairs.  

Usage:   llinregress(x,y)      x,y are equal-length lists of x-y coordinates
Returns: slope, intercept, r, two-tailed prob, sterr-of-estimate
"""
    TINY = 1.0e-20
    if len(x) <> len(y):
	raise ValueError, 'Input values not paired in linregress.  Aborting.'
    n = len(x)
    x = map(float,x)
    y = map(float,y)
    xmean = mean(x)
    ymean = mean(y)
    r_num = float(n*(summult(x,y)) - sum(x)*sum(y))
    r_den = math.sqrt((n*ss(x) - square_of_sums(x))*(n*ss(y)-square_of_sums(y)))
    r = r_num / r_den
    z = 0.5*math.log((1.0+r+TINY)/(1.0-r+TINY))
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    slope = r_num / float(n*ss(x) - square_of_sums(x))
    intercept = ymean - slope*xmean
    sterrest = math.sqrt(1-r*r)*samplestdev(y)
    return slope, intercept, r, prob, sterrest


####################################
#####  INFERENTIAL STATISTICS  #####
####################################

def lttest_1samp(a,popmean,printit=0,name='Sample',writemode='a'):
    """
Calculates the t-obtained for the independent samples T-test on ONE group
of scores a, given a population mean.  If printit=1, results are printed
to the screen.  If printit='filename', the results are output to 'filename'
using the given writemode (default=append).  Returns t-value, and prob.

Usage:   lttest_1samp(a,popmean,Name='Sample',printit=0,writemode='a')
Returns: t-value, two-tailed prob
"""
    x = mean(a)
    v = var(a)
    n = len(a)
    df = n-1
    svar = ((n-1)*v)/float(df)
    t = (x-popmean)/math.sqrt(svar*(1.0/n))
    prob = betai(0.5*df,0.5,float(df)/(df+t*t))

    if printit <> 0:
	statname = 'Single-sample T-test.'
	outputpairedstats(printit,writemode,
			  'Population','--',popmean,0,0,0,
			  name,n,x,v,min(a),max(a),
			  statname,t,prob)
    return t,prob


def lttest_ind (a, b, printit=0, name1='Samp1', name2='Samp2', writemode='a'):
    """
Calculates the t-obtained T-test on TWO INDEPENDENT samples of
scores a, and b.  From Numerical Recipies, p.483.  If printit=1, results
are printed to the screen.  If printit='filename', the results are output
to 'filename' using the given writemode (default=append).  Returns t-value,
and prob.

Usage:   lttest_ind(a,b,printit=0,name1='Samp1',name2='Samp2',writemode='a')
Returns: t-value, two-tailed prob
"""
    x1 = mean(a)
    x2 = mean(b)
    v1 = stdev(a)**2
    v2 = stdev(b)**2
    n1 = len(a)
    n2 = len(b)
    df = n1+n2-2
    svar = ((n1-1)*v1+(n2-1)*v2)/float(df)
    t = (x1-x2)/math.sqrt(svar*(1.0/n1 + 1.0/n2))
    prob = betai(0.5*df,0.5,df/(df+t*t))

    if printit <> 0:
	statname = 'Independent samples T-test.'
	outputpairedstats(printit,writemode,
			  name1,n1,x1,v1,min(a),max(a),
			  name2,n2,x2,v2,min(b),max(b),
			  statname,t,prob)
    return t,prob


def lttest_rel (a,b,printit=0,name1='Sample1',name2='Sample2',writemode='a'):
    """
Calculates the t-obtained T-test on TWO RELATED samples of scores,
a and b.  From Numerical Recipies, p.483.  If printit=1, results are
printed to the screen.  If printit='filename', the results are output to
'filename' using the given writemode (default=append).  Returns t-value,
and prob.

Usage:   lttest_rel(a,b,printit=0,name1='Sample1',name2='Sample2',writemode='a')
Returns: t-value, two-tailed prob
"""
    if len(a)<>len(b):
	raise ValueError, 'Unequal length lists in ttest_rel.'
    x1 = mean(a)
    x2 = mean(b)
    v1 = var(a)
    v2 = var(b)
    n = len(a)
    cov = 0
    for i in range(len(a)):
	cov = cov + (a[i]-x1) * (b[i]-x2)
    df = n-1
    cov = cov / float(df)
    sd = math.sqrt((v1+v2 - 2.0*cov)/float(n))
    t = (x1-x2)/sd
    prob = betai(0.5*df,0.5,df/(df+t*t))

    if printit <> 0:
	statname = 'Related samples T-test.'
	outputpairedstats(printit,writemode,
			  name1,n,x1,v1,min(a),max(a),
			  name2,n,x2,v2,min(b),max(b),
			  statname,t,prob)
    return t, prob


def lchisquare(f_obs,f_exp=None):
    """
Calculates a one-way chi square for list of observed frequencies and returns
the result.  If no expected frequencies are given, the total N is assumed to
be equally distributed across all groups.

Usage:   lchisquare(f_obs, f_exp=None)   f_obs = list of observed cell freq.
Returns: chisquare-statistic, associated p-value
"""
    k = len(f_obs)                 # number of groups
    if f_exp == None:
	f_exp = [sum(f_obs)/float(k)] * len(f_obs) # create k bins with = freq.
    chisq = 0
    for i in range(len(f_obs)):
	chisq = chisq + (f_obs[i]-f_exp[i])**2 / float(f_exp[i])
    return chisq, chisqprob(chisq, k-1)


def lks_2samp (data1,data2):
    """
Computes the Kolmogorov-Smirnof statistic on 2 samples.  From
Numerical Recipies in C, page 493.

Usage:   lks_2samp(data1,data2)   data1&2 are lists of values for 2 conditions
Returns: KS D-value, associated p-value
"""
    j1 = 0
    j2 = 0
    fn1 = 0.0
    fn2 = 0.0
    n1 = len(data1)
    n2 = len(data2)
    en1 = n1
    en2 = n2
    d = 0.0
    data1.sort()
    data2.sort()
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
	if math.fabs(dt) > math.fabs(d):
	    d = dt
    try:
        en = math.sqrt(en1*en2/float(en1+en2))
        prob = ksprob((en+0.12+0.11/en)*abs(d))
    except:
        prob = 1.0
    return d, prob


def lmannwhitneyu(x,y):
    """
Calculates a Mann-Whitney U statistic on the provided scores and
returns the result.  Use only when the n in each condition is < 20 and
you have 2 independent samples of ranks.  NOTE: Mann-Whitney U is
significant if the u-obtained is LESS THAN or equal to the critical
value of U found in the tables.  Equivalent to Kruskal-Wallis H with
just 2 groups.

Usage:   lmannwhitneyu(data)
Returns: u-statistic, one-tailed p-value (i.e., p(z(U)))
"""
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(x+y)
    rankx = ranked[0:n1]       # get the x-ranks
    ranky = ranked[n1:]        # the rest are y-ranks
    u1 = n1*n2 + (n1*(n1+1))/2.0 - sum(rankx)  # calc U for x
    u2 = n1*n2 - u1                            # remainder is U for y
    bigu = max(u1,u2)
    smallu = min(u1,u2)
    T = math.sqrt(tiecorrect(ranked))  # correction factor for tied scores
    if T == 0:
	raise ValueError, 'All numbers are identical in lmannwhitneyu'
    sd = math.sqrt(T*n1*n2*(n1+n2+1)/12.0)
    z = abs((bigu-n1*n2/2.0) / sd)  # normal approximation for prob calc
    return smallu, 1.0 - zprob(z)


def ltiecorrect(rankvals):
    """
Corrects for ties in Mann Whitney U and Kruskal Wallis H tests.  See
Siegel, S. (1956) Nonparametric Statistics for the Behavioral Sciences.
New York: McGraw-Hill.  Code adapted from |Stat rankind.c code.

Usage:   ltiecorrect(rankvals)
Returns: T correction factor for U or H
"""
    sorted,posn = shellsort(rankvals)
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


def lranksums(x,y):
    """
Calculates the rank sums statistic on the provided scores and
returns the result.  Use only when the n in each condition is > 20 and you
have 2 independent samples of ranks.

Usage:   lranksums(x,y)
Returns: a z-statistic, two-tailed p-value
"""
    n1 = len(x)
    n2 = len(y)
    alldata = x+y
    ranked = rankdata(alldata)
    x = ranked[:n1]
    y = ranked[n1:]
    s = sum(x)
    expected = n1*(n1+n2+1) / 2.0
    z = (s - expected) / math.sqrt(n1*n2*(n1+n2+1)/12.0)
    prob = 2*(1.0 -zprob(abs(z)))
    return z, prob


def lwilcoxont(x,y):
    """
Calculates the Wilcoxon T-test for related samples and returns the
result.  A non-parametric T-test.

Usage:   lwilcoxont(x,y)
Returns: a t-statistic, two-tail probability estimate
"""
    if len(x) <> len(y):
	raise ValueError, 'Unequal N in wilcoxont.  Aborting.'
    d=[]
    for i in range(len(x)):
	diff = x[i] - y[i]
	if diff <> 0:
	    d.append(diff)
    count = len(d)
    absd = map(abs,d)
    absranked = rankdata(absd)
    r_plus = 0.0
    r_minus = 0.0
    for i in range(len(absd)):
	if d[i] < 0:
	    r_minus = r_minus + absranked[i]
	else:
	    r_plus = r_plus + absranked[i]
    wt = min(r_plus, r_minus)
    mn = count * (count+1) * 0.25
    se =  math.sqrt(count*(count+1)*(2.0*count+1.0)/24.0)
    z = math.fabs(wt-mn) / se
    prob = 2*(1.0 -zprob(abs(z)))
    return wt, prob


def lkruskalwallish(*args):
    """
The Kruskal-Wallis H-test is a non-parametric ANOVA for 3 or more
groups, requiring at least 5 subjects in each group.  This function
calculates the Kruskal-Wallis H-test for 3 or more independent samples
and returns the result.  

Usage:   lkruskalwallish(*args)
Returns: H-statistic (corrected for ties), associated p-value
"""
    args = list(args)
    n = [0]*len(args)
    all = []
    n = map(len,args)
    for i in range(len(args)):
	all = all + args[i]
    ranked = rankdata(all)
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
	raise ValueError, 'All numbers are identical in lkruskalwallish'
    h = h / float(T)
    return h, chisqprob(h,df)


def lfriedmanchisquare(*args):
    """
Friedman Chi-Square is a non-parametric, one-way within-subjects
ANOVA.  This function calculates the Friedman Chi-square test for repeated
measures and returns the result, along with the associated probability
value.  It assumes 3 or more repeated measures.  Only 3 levels requires a
minimum of 10 subjects in the study.  Four levels requires 5 subjects per
level(??).

Usage:   lfriedmanchisquare(*args)
Returns: chi-square statistic, associated p-value
"""
    k = len(args)
    if k < 3:
	raise ValueError, 'Less than 3 levels.  Friedman test not appropriate.'
    n = len(args[0])
    data = apply(pstat.abut,tuple(args))
    for i in range(len(data)):
	data[i] = rankdata(data[i])
    ssbn = 0
    for i in range(k):
	ssbn = ssbn + sum(args[i])**2
    chisq = 12.0 / (k*n*(k+1)) * ssbn - 3*n*(k+1)
    return chisq, chisqprob(chisq,k-1)


####################################
####  PROBABILITY CALCULATIONS  ####
####################################

def lchisqprob(chisq,df):
    """
Returns the (1-tailed) probability value associated with the provided
chi-square value and df.  Adapted from chisq.c in Gary Perlman's |Stat.

Usage:   lchisqprob(chisq,df)
"""
    BIG = 20.0
    def ex(x):
	BIG = 20.0
	if x < -BIG:
	    return 0.0
	else:
	    return math.exp(x)

    if chisq <=0 or df < 1:
	return 1.0
    a = 0.5 * chisq
    if df%2 == 0:
	even = 1
    else:
	even = 0
    if df > 1:
	y = ex(-a)
    if even:
	s = y
    else:
	s = 2.0 * zprob(-math.sqrt(chisq))
    if (df > 2):
	chisq = 0.5 * (df - 1.0)
	if even:
	    z = 1.0
	else:
	    z = 0.5
	if a > BIG:
	    if even:
		e = 0.0
	    else:
		e = math.log(math.sqrt(math.pi))
	    c = math.log(a)
	    while (z <= chisq):
		e = math.log(z) + e
		s = s + ex(c*z-a-e)
		z = z + 1.0
	    return s
	else:
	    if even:
		e = 1.0
	    else:
		e = 1.0 / math.sqrt(math.pi) / math.sqrt(a)
		c = 0.0
		while (z <= chisq):
		    e = e * (a/float(z))
		    c = c + e
		    z = z + 1.0
		return (c*y+s)
    else:
	return s


def lerfcc(x):
    """
Returns the complementary error function erfc(x) with fractional
error everywhere less than 1.2e-7.  Adapted from Numerical Recipies.

Usage:   lerfcc(x)
"""
    z = abs(x)
    t = 1.0 / (1.0+0.5*z)
    ans = t * math.exp(-z*z-1.26551223 + t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
    if x >= 0:
	return ans
    else:
	return 2.0 - ans


def lzprob(z):
    """
Returns the area under the normal curve 'to the left of' the given z value.
Thus, 
    for z<0, zprob(z) = 1-tail probability
    for z>0, 1.0-zprob(z) = 1-tail probability
    for any z, 2.0*(1.0-zprob(abs(z))) = 2-tail probability
Adapted from z.c in Gary Perlman's |Stat.

Usage:   lzprob(z)
"""
    Z_MAX = 6.0    # maximum meaningful z-value
    if z == 0.0:
	x = 0.0
    else:
	y = 0.5 * math.fabs(z)
	if y >= (Z_MAX*0.5):
	    x = 1.0
	elif (y < 1.0):
	    w = y*y
	    x = ((((((((0.000124818987 * w
			-0.001075204047) * w +0.005198775019) * w
		      -0.019198292004) * w +0.059054035642) * w
		    -0.151968751364) * w +0.319152932694) * w
		  -0.531923007300) * w +0.797884560593) * y * 2.0
	else:
	    y = y - 2.0
	    x = (((((((((((((-0.000045255659 * y
			     +0.000152529290) * y -0.000019538132) * y
			   -0.000676904986) * y +0.001390604284) * y
			 -0.000794620820) * y -0.002034254874) * y
		       +0.006549791214) * y -0.010557625006) * y
		     +0.011630447319) * y -0.009279453341) * y
		   +0.005353579108) * y -0.002141268741) * y
		 +0.000535310849) * y +0.999936657524
    if z > 0.0:
	prob = ((x+1.0)*0.5)
    else:
	prob = ((1.0-x)*0.5)
    return prob


def lksprob(alam):
    """
Computes a Kolmolgorov-Smirnov t-test significance level.  Adapted from
Numerical Recipies.

Usage:   lksprob(alam)
"""
    fac = 2.0
    sum = 0.0
    termbf = 0.0
    a2 = -2.0*alam*alam
    for j in range(1,201):
	term = fac*math.exp(a2*j*j)
	sum = sum + term
	if math.fabs(term) <= (0.001*termbf) or math.fabs(term) < (1.0e-8*sum):
	    return sum
	fac = -fac
	termbf = math.fabs(term)
    return 1.0             # Get here only if fails to converge; was 0.0!!


def lfprob (dfnum, dfden, F):
    """
Returns the (1-tailed) significance level (p-value) of an F
statistic given the degrees of freedom for the numerator (dfR-dfF) and
the degrees of freedom for the denominator (dfF).

Usage:   lfprob(dfnum, dfden, F)   where usually dfnum=dfbn, dfden=dfwn
"""
    p = betai(0.5*dfden, 0.5*dfnum, dfden/float(dfden+dfnum*F))
    return p


def lbetacf(a,b,x):
    """
This function evaluates the continued fraction form of the incomplete
Beta function, betai.  (Adapted from: Numerical Recipies in C.)

Usage:   lbetacf(a,b,x)
"""
    ITMAX = 200
    EPS = 3.0e-7

    bm = az = am = 1.0
    qab = a+b
    qap = a+1.0
    qam = a-1.0
    bz = 1.0-qab*x/qap
    for i in range(ITMAX+1):
	em = float(i+1)
	tem = em + em
	d = em*(b-em)*x/((qam+tem)*(a+tem))
	ap = az + d*am
	bp = bz+d*bm
	d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
	app = ap+d*az
	bpp = bp+d*bz
	aold = az
	am = ap/bpp
	bm = bp/bpp
	az = app/bpp
	bz = 1.0
	if (abs(az-aold)<(EPS*abs(az))):
	    return az
    print 'a or b too big, or ITMAX too small in Betacf.'


def lgammln(xx):
    """
Returns the gamma function of xx.
    Gamma(z) = Integral(0,infinity) of t^(z-1)exp(-t) dt.
(Adapted from: Numerical Recipies in C.)

Usage:   lgammln(xx)
"""

    coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516,
	     0.120858003e-2, -0.536382e-5]
    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5)*math.log(tmp)
    ser = 1.0
    for j in range(len(coeff)):
	x = x + 1
	ser = ser + coeff[j]/x
    return -tmp + math.log(2.50662827465*ser)


def lbetai(a,b,x):
    """
Returns the incomplete beta function:

    I-sub-x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)

where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
function of a.  The continued fraction formulation is implemented here,
using the betacf function.  (Adapted from: Numerical Recipies in C.)

Usage:   lbetai(a,b,x)
"""
    if (x<0.0 or x>1.0):
	raise ValueError, 'Bad x in lbetai'
    if (x==0.0 or x==1.0):
	bt = 0.0
    else:
	bt = math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*
		      math.log(1.0-x))
    if (x<(a+1.0)/(a+b+2.0)):
	return bt*betacf(a,b,x)/float(a)
    else:
	return 1.0-bt*betacf(b,a,1.0-x)/float(b)


####################################
#######  ANOVA CALCULATIONS  #######
####################################

def lF_oneway(*lists):
    """
Performs a 1-way ANOVA, returning an F-value and probability given
any number of groups.  From Heiman, pp.394-7.

Usage:   F_oneway(*lists)    where *lists is any number of lists, one per
                                  treatment group
Returns: F value, one-tailed p-value
"""
    a = len(lists)           # ANOVA on 'a' groups, each in it's own list
    means = [0]*a
    vars = [0]*a
    ns = [0]*a
    alldata = []
    tmp = map(N.array,lists)
    means = map(amean,tmp)
    vars = map(avar,tmp)
    ns = map(len,lists)
    for i in range(len(lists)):
	alldata = alldata + lists[i]
    alldata = N.array(alldata)
    bign = len(alldata)
    sstot = ass(alldata)-(asquare_of_sums(alldata)/float(bign))
    ssbn = 0
    for list in lists:
	ssbn = ssbn + asquare_of_sums(N.array(list))/float(len(list))
    ssbn = ssbn - (asquare_of_sums(alldata)/float(bign))
    sswn = sstot-ssbn
    dfbn = a-1
    dfwn = bign - a
    msb = ssbn/float(dfbn)
    msw = sswn/float(dfwn)
    f = msb/msw
    prob = fprob(dfbn,dfwn,f)
    return f, prob


def lF_value (ER,EF,dfnum,dfden):
    """
Returns an F-statistic given the following:
	ER  = error associated with the null hypothesis (the Restricted model)
	EF  = error associated with the alternate hypothesis (the Full model)
	dfR-dfF = degrees of freedom of the numerator
	dfF = degrees of freedom associated with the denominator/Full model

Usage:   lF_value(ER,EF,dfnum,dfden)
"""
    return ((ER-EF)/float(dfnum) / (EF/float(dfden)))


####################################
########  SUPPORT FUNCTIONS  #######
####################################

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
	items = pstat.colex(list2print,col)
	items = map(pstat.makestr,items)
	maxsize[col] = max(map(len,items)) + extra
    for row in listoflists:
	if row == ['\n'] or row == '\n':
	    outfile.write('\n')
	elif row == ['dashes'] or row == 'dashes':
	    dashes = [0]*len(maxsize)
	    for j in range(len(maxsize)):
		dashes[j] = '-'*(maxsize[j]-2)
	    outfile.write(pstat.lineincustcols(dashes,maxsize))
	else:
	    outfile.write(pstat.lineincustcols(row,maxsize))
	outfile.write('\n')
    outfile.close()
    return None


def lincr(l,cap):        # to increment a list up to a max-list of 'cap'
    """
Simulate a counting system from an n-dimensional list.

Usage:   lincr(l,cap)   l=list to increment, cap=max values for each list pos'n
Returns: next set of values for list l, OR -1 (if overflow)
"""
    l[0] = l[0] + 1     # e.g., [0,0,0] --> [2,4,3] (=cap)
    for i in range(len(l)):
	if l[i] > cap[i] and i < len(l)-1: # if carryover AND not done
	    l[i] = 0
	    l[i+1] = l[i+1] + 1
	elif l[i] > cap[i] and i == len(l)-1: # overflow past last column, must be finished
	    l = -1
    return l


def lsum (inlist):
    """
Returns the sum of the items in the passed list.

Usage:   lsum(inlist)
"""
    s = 0
    for item in inlist:
	s = s + item
    return s


def lcumsum (inlist):
    """
Returns a list consisting of the cumulative sum of the items in the
passed list.

Usage:   lcumsum(inlist)
"""
    newlist = copy.deepcopy(inlist)
    for i in range(1,len(newlist)):
	newlist[i] = newlist[i] + newlist[i-1]
    return newlist


def lss(inlist):
    """
Squares each value in the passed list, adds up these squares and
returns the result.

Usage:   lss(inlist)
"""
    ss = 0
    for item in inlist:
	ss = ss + item*item
    return ss


def lsummult (list1,list2):
    """
Multiplies elements in list1 and list2, element by element, and
returns the sum of all resulting multiplications.  Must provide equal
length lists.

Usage:   lsummult(list1,list2)
"""
    if len(list1) <> len(list2):
	raise ValueError, "Lists not equal length in summult."
    s = 0
    for item1,item2 in pstat.abut(list1,list2):
	s = s + item1*item2
    return s


def lsumdiffsquared(x,y):
    """
Takes pairwise differences of the values in lists x and y, squares
these differences, and returns the sum of these squares.

Usage:   lsumdiffsquared(x,y)
Returns: sum[(x[i]-y[i])**2]
"""
    sds = 0
    for i in range(len(x)):
	sds = sds + (x[i]-y[i])**2
    return sds


def lsquare_of_sums(inlist):
    """
Adds the values in the passed list, squares the sum, and returns
the result.

Usage:   lsquare_of_sums(inlist)
Returns: sum(inlist[i])**2
"""
    s = sum(inlist)
    return float(s)*s


def lshellsort(inlist):
    """
Shellsort algorithm.  Sorts a 1D-list.

Usage:   lshellsort(inlist)
Returns: sorted-inlist, sorting-index-vector (for original list)
"""
    n = len(inlist)
    svec = copy.deepcopy(inlist)
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
# svec is now sorted inlist, and ivec has the order svec[i] = vec[ivec[i]]
    return svec, ivec


def lrankdata(inlist):
    """
Ranks the data in inlist, dealing with ties appropritely.  Assumes
a 1D inlist.  Adapted from Gary Perlman's |Stat ranksort.

Usage:   lrankdata(inlist)
Returns: a list of length equal to inlist, containing rank scores
"""
    n = len(inlist)
    svec, ivec = shellsort(inlist)
    sumranks = 0
    dupcount = 0
    newlist = [0]*n
    for i in range(n):
	sumranks = sumranks + i
	dupcount = dupcount + 1
	if i==n-1 or svec[i] <> svec[i+1]:
	    averank = sumranks / float(dupcount) + 1
	    for j in range(i-dupcount+1,i+1):
		newlist[ivec[j]] = averank
	    sumranks = 0
	    dupcount = 0
    return newlist


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
	pstat.printcc(lofl)
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
	file.write(pstat.list2string(['\nTest statistic = ',round(stat,4),'   p = ',round(prob,4),suffix,'\n\n']))
	file.close()
    return None


def lfindwithin (data):
    """
Returns an integer representing a binary vector, where 1=within-
subject factor, 0=between.  Input equals the entire data 2D list (i.e.,
column 0=random factor, column -1=measured values (those two are skipped).
Note: input data is in |Stat format ... a list of lists ("2D list") with 
one row per measured value, first column=subject identifier, last column=
score, one in-between column per factor (these columns contain level
designations on each factor).  See also stats.anova.__doc__.

Usage:   lfindwithin(data)     data in |Stat format
"""

    numfact = len(data[0])-1
    withinvec = 0
    for col in range(1,numfact):
	examplelevel = pstat.unique(pstat.colex(data,col))[0]
	rows = pstat.linexand(data,col,examplelevel)  # get 1 level of this factor
	factsubjs = pstat.unique(pstat.colex(rows,0))
	allsubjs = pstat.unique(pstat.colex(data,0))
	if len(factsubjs) == len(allsubjs):  # fewer Ss than scores on this factor?
	    withinvec = withinvec + (1 << col)
    return withinvec


#########################################################
#########################################################
####### DISPATCH LISTS AND TUPLES TO ABOVE FCNS #########
#########################################################
#########################################################

## CENTRAL TENDENCY:
geometricmean = Dispatch ( (lgeometricmean, (ListType, TupleType)), )
harmonicmean = Dispatch ( (lharmonicmean, (ListType, TupleType)), )
mean = Dispatch ( (lmean, (ListType, TupleType)), )
median = Dispatch ( (lmedian, (ListType, TupleType)), )
medianscore = Dispatch ( (lmedianscore, (ListType, TupleType)), )
mode = Dispatch ( (lmode, (ListType, TupleType)), )

## MOMENTS:
moment = Dispatch ( (lmoment, (ListType, TupleType)), )
variation = Dispatch ( (lvariation, (ListType, TupleType)), )
skew = Dispatch ( (lskew, (ListType, TupleType)), )
kurtosis = Dispatch ( (lkurtosis, (ListType, TupleType)), )
describe = Dispatch ( (ldescribe, (ListType, TupleType)), )

## FREQUENCY STATISTICS:
itemfreq = Dispatch ( (litemfreq, (ListType, TupleType)), )
scoreatpercentile = Dispatch ( (lscoreatpercentile, (ListType, TupleType)), )
percentileofscore = Dispatch ( (lpercentileofscore, (ListType, TupleType)), )
histogram = Dispatch ( (lhistogram, (ListType, TupleType)), )
cumfreq = Dispatch ( (lcumfreq, (ListType, TupleType)), )
relfreq = Dispatch ( (lrelfreq, (ListType, TupleType)), )

## VARIABILITY:
obrientransform = Dispatch ( (lobrientransform, (ListType, TupleType)), )
samplevar = Dispatch ( (lsamplevar, (ListType, TupleType)), )
samplestdev = Dispatch ( (lsamplestdev, (ListType, TupleType)), )
var = Dispatch ( (lvar, (ListType, TupleType)), )
stdev = Dispatch ( (lstdev, (ListType, TupleType)), )
sterr = Dispatch ( (lsterr, (ListType, TupleType)), )
sem = Dispatch ( (lsem, (ListType, TupleType)), )
z = Dispatch ( (lz, (ListType, TupleType)), )
zs = Dispatch ( (lzs, (ListType, TupleType)), )

## TRIMMING FCNS:
trimboth = Dispatch ( (ltrimboth, (ListType, TupleType)), )
trim1 = Dispatch ( (ltrim1, (ListType, TupleType)), )

## CORRELATION FCNS:
paired = Dispatch ( (lpaired, (ListType, TupleType)), )
pearsonr = Dispatch ( (lpearsonr, (ListType, TupleType)), )
spearmanr = Dispatch ( (lspearmanr, (ListType, TupleType)), )
pointbiserialr = Dispatch ( (lpointbiserialr, (ListType, TupleType)), )
kendalltau = Dispatch ( (lkendalltau, (ListType, TupleType)), )
linregress = Dispatch ( (llinregress, (ListType, TupleType)), )

## INFERENTIAL STATS:
ttest_1samp = Dispatch ( (lttest_1samp, (ListType, TupleType)), )
ttest_ind = Dispatch ( (lttest_ind, (ListType, TupleType)), )
ttest_rel = Dispatch ( (lttest_rel, (ListType, TupleType)), )
chisquare = Dispatch ( (lchisquare, (ListType, TupleType)), )
ks_2samp = Dispatch ( (lks_2samp, (ListType, TupleType)), )
mannwhitneyu = Dispatch ( (lmannwhitneyu, (ListType, TupleType)), )
ranksums = Dispatch ( (lranksums, (ListType, TupleType)), )
tiecorrect = Dispatch ( (ltiecorrect, (ListType, TupleType)), )
wilcoxont = Dispatch ( (lwilcoxont, (ListType, TupleType)), )
kruskalwallish = Dispatch ( (lkruskalwallish, (ListType, TupleType)), )
friedmanchisquare = Dispatch ( (lfriedmanchisquare, (ListType, TupleType)), )

## PROBABILITY CALCS:
chisqprob = Dispatch ( (lchisqprob, (IntType, FloatType)), )
zprob = Dispatch ( (lzprob, (IntType, FloatType)), )
ksprob = Dispatch ( (lksprob, (IntType, FloatType)), )
fprob = Dispatch ( (lfprob, (IntType, FloatType)), )
betacf = Dispatch ( (lbetacf, (IntType, FloatType)), )
betai = Dispatch ( (lbetai, (IntType, FloatType)), )
erfcc = Dispatch ( (lerfcc, (IntType, FloatType)), )
gammln = Dispatch ( (lgammln, (IntType, FloatType)), )

## ANOVA FUNCTIONS:
F_oneway = Dispatch ( (lF_oneway, (ListType, TupleType)), )
F_value = Dispatch ( (lF_value, (ListType, TupleType)), )

## SUPPORT FUNCTIONS:
incr = Dispatch ( (lincr, (ListType, TupleType)), )
sum = Dispatch ( (lsum, (ListType, TupleType)), )
cumsum = Dispatch ( (lcumsum, (ListType, TupleType)), )
ss = Dispatch ( (lss, (ListType, TupleType)), )
summult = Dispatch ( (lsummult, (ListType, TupleType)), )
square_of_sums = Dispatch ( (lsquare_of_sums, (ListType, TupleType)), )
sumdiffsquared = Dispatch ( (lsumdiffsquared, (ListType, TupleType)), )
shellsort = Dispatch ( (lshellsort, (ListType, TupleType)), )
rankdata = Dispatch ( (lrankdata, (ListType, TupleType)), )
findwithin = Dispatch ( (lfindwithin, (ListType, TupleType)), )


#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============
#=============  THE ARRAY-VERSION OF THE STATS FUNCTIONS  ===============

try:                         # DEFINE THESE *ONLY* IF NUMERIC IS AVAILABLE
 import scipy.numeric as Numeric
 N = Numeric
 import LinearAlgebra
 LA = LinearAlgebra


#####################################
########  ACENTRAL TENDENCY  ########
#####################################

 def ageometricmean (inarray,dimension=None,keepdims=0):
    """
Calculates the geometric mean of the values in the passed array.
That is:  n-th root of (x1 * x2 * ... * xn).  Defaults to ALL values in
the passed array.  Use dimension=None to flatten array first.  REMEMBER: if
dimension=0, it collapses over dimension 0 ('rows' in a 2D array) only, and
if dimension is a sequence, it collapses over all specified dimensions.  If
keepdims is set to 1, the resulting array will have as many dimensions as
inarray, with only 1 'level' per dim that was collapsed over.

Usage:   ageometricmean(inarray,dimension=None,keepdims=0)
Returns: geometric mean computed over dim(s) listed in dimension
"""
    inarray = N.array(inarray,N.Float)
    if dimension == None:
	inarray = N.ravel(inarray)
	size = len(inarray)
	mult = N.power(inarray,1.0/size)
	mult = N.multiply.reduce(mult)
    elif type(dimension) in [IntType,FloatType]:
	size = inarray.shape[dimension]
	mult = N.power(inarray,1.0/size)
	mult = N.multiply.reduce(mult,dimension)
	if keepdims == 1:
	    shp = list(inarray.shape)
	    shp[dimension] = 1
	    sum = N.reshape(sum,shp)
    else: # must be a SEQUENCE of dims to average over
	dims = list(dimension)
	dims.sort()
	dims.reverse()
	size = N.array(N.multiply.reduce(N.take(inarray.shape,dims)),N.Float)
	mult = N.power(inarray,1.0/size)
	for dim in dims:
	    mult = N.multiply.reduce(mult,dim)
	if keepdims == 1:
	    shp = list(inarray.shape)
	    for dim in dims:
		shp[dim] = 1
	    mult = N.reshape(mult,shp)
    return mult


 def aharmonicmean (inarray,dimension=None,keepdims=0):
    """
Calculates the harmonic mean of the values in the passed array.
That is:  n / (1/x1 + 1/x2 + ... + 1/xn).  Defaults to ALL values in
the passed array.  Use dimension=None to flatten array first.  REMEMBER: if
dimension=0, it collapses over dimension 0 ('rows' in a 2D array) only, and
if dimension is a sequence, it collapses over all specified dimensions.  If
keepdims is set to 1, the resulting array will have as many dimensions as
inarray, with only 1 'level' per dim that was collapsed over.

Usage:   aharmonicmean(inarray,dimension=None,keepdims=0)
Returns: harmonic mean computed over dim(s) in dimension
"""
    inarray = inarray.astype(N.Float)
    if dimension == None:
	inarray = N.ravel(inarray)
	size = len(inarray)
	s = N.add.reduce(1.0 / inarray)
    elif type(dimension) in [IntType,FloatType]:
	size = float(inarray.shape[dimension])
	s = N.add.reduce(1.0/inarray, dimension)
	if keepdims == 1:
	    shp = list(inarray.shape)
	    shp[dimension] = 1
	    s = N.reshape(s,shp)
    else: # must be a SEQUENCE of dims to average over
	dims = list(dimension)
	dims.sort()
	nondims = []
	for i in range(len(inarray.shape)):
	    if i not in dims:
		nondims.append(i)
	tinarray = N.transpose(inarray,nondims+dims) # put keep-dims first
	idx = [0] *len(nondims)
	if idx == []:
	    size = len(N.ravel(inarray))
	    s = asum(1.0 / inarray)
	    if keepdims == 1:
		s = N.reshape([s],N.ones(len(inarray.shape)))
	else:
	    idx[0] = -1
	    loopcap = N.array(tinarray.shape[0:len(nondims)]) -1
	    s = N.zeros(loopcap+1,N.Float)
	    while incr(idx,loopcap) <> -1:
		s[idx] = asum(1.0/tinarray[idx])
	    size = N.multiply.reduce(N.take(inarray.shape,dims))
	    if keepdims == 1:
		shp = list(inarray.shape)
		for dim in dims:
		    shp[dim] = 1
		s = N.reshape(s,shp)
    return size / s


 def amean (inarray,dimension=None,keepdims=0):
    """
Calculates the arithmatic mean of the values in the passed array.
That is:  1/n * (x1 + x2 + ... + xn).  Defaults to ALL values in the
passed array.  Use dimension=None to flatten array first.  REMEMBER: if
dimension=0, it collapses over dimension 0 ('rows' in a 2D array) only, and
if dimension is a sequence, it collapses over all specified dimensions.  If
keepdims is set to 1, the resulting array will have as many dimensions as
inarray, with only 1 'level' per dim that was collapsed over.

Usage:   amean(inarray,dimension=None,keepdims=0)
Returns: arithematic mean calculated over dim(s) in dimension
"""
    if inarray.typecode() in ['l','s','b']:
        inarray = inarray.astype(N.Float)
    if dimension == None:
	inarray = N.ravel(inarray)
	sum = N.add.reduce(inarray)
	denom = float(len(inarray))
    elif type(dimension) in [IntType,FloatType]:
	sum = asum(inarray,dimension)
	denom = float(inarray.shape[dimension])
	if keepdims == 1:
	    shp = list(inarray.shape)
	    shp[dimension] = 1
	    sum = N.reshape(sum,shp)
    else: # must be a TUPLE of dims to average over
	dims = list(dimension)
	dims.sort()
	dims.reverse()
	sum = inarray *1.0
	for dim in dims:
	    sum = N.add.reduce(sum,dim)
	denom = N.array(N.multiply.reduce(N.take(inarray.shape,dims)),N.Float)
	if keepdims == 1:
	    shp = list(inarray.shape)
	    for dim in dims:
		shp[dim] = 1
	    sum = N.reshape(sum,shp)
    return sum/denom


 def amedian (inarray,numbins=1000):
    """
Calculates the COMPUTED median value of an array of numbers, given the
number of bins to use for the histogram (more bins approaches finding the
precise median value of the array; default number of bins = 1000).  From
G.W. Heiman's Basic Stats, or CRC Probability & Statistics.
NOTE:  THIS ROUTINE ALWAYS uses the entire passed array (flattens it first).

Usage:   amedian(inarray,numbins=1000)
Returns: median calculated over ALL values in inarray
"""
    inarray = N.ravel(inarray)
    (hist, smallest, binsize, extras) = ahistogram(inarray,numbins)
    cumhist = N.cumsum(hist)            # make cumulative histogram
    otherbins = N.greater_equal(cumhist,len(inarray)/2.0)
    otherbins = list(otherbins)         # list of 0/1s, 1s start at median bin
    cfbin = otherbins.index(1)		# get 1st(!) index holding 50%ile score
    LRL = smallest + binsize*cfbin	# get lower read limit of that bin
    cfbelow = N.add.reduce(hist[0:cfbin])	# cum. freq. below bin
    freq = hist[cfbin]			# frequency IN the 50%ile bin
    median = LRL + ((len(inarray)/2.0-cfbelow)/float(freq))*binsize # MEDIAN
    return median


 def amedianscore (inarray,dimension=None):
    """
Returns the 'middle' score of the passed array.  If there is an even
number of scores, the mean of the 2 middle scores is returned.  Can function
with 1D arrays, or on the FIRST dimension of 2D arrays (i.e., dimension can
be None, to pre-flatten the array, or else dimension must equal 0).

Usage:   amedianscore(inarray,dimension=None)
Returns: 'middle' score of the array, or the mean of the 2 middle scores
"""
    if dimension == None:
	inarray = N.ravel(inarray)
	dimension = 0
    inarray = N.sort(inarray,dimension)
    if inarray.shape[dimension] % 2 == 0:   # if even number of elements
	indx = inarray.shape[dimension]/2   # integer division correct
	median = N.asarray(inarray[indx]+inarray[indx-1]) / 2.0
    else:
	indx = inarray.shape[dimension] / 2 # integer division correct
	median = N.take(inarray,[indx],dimension)
	if median.shape == (1,):
	    median = median[0]
    return median


 def amode(a, dimension=None):
    """
Returns an array of the modal (most common) score in the passed array.
If there is more than one such score, ONLY THE FIRST is returned.
The bin-count for the modal values is also returned.  Operates on whole
array (dimension=None), or on a given dimension.

Usage:   amode(a, dimension=None)
Returns: array of bin-counts for mode(s), array of corresponding modal values
"""

    if dimension == None:
	a = N.ravel(a)
	dimension = 0
    scores = pstat.aunique(N.ravel(a))       # get ALL unique values
    testshape = list(a.shape)
    testshape[dimension] = 1
    oldmostfreq = N.zeros(testshape)
    oldcounts = N.zeros(testshape)
    for score in scores:
	template = N.equal(a,score)
	counts = asum(template,dimension,1)
	mostfrequent = N.where(N.greater(counts,oldcounts),score,oldmostfreq)
	oldcounts = N.where(N.greater(counts,oldcounts),counts,oldcounts)
	oldmostfreq = mostfrequent
    return oldcounts, mostfrequent


 def atmean(a,limits=None,inclusive=(1,1)):
     """
Returns the arithmetic mean of all values in an array, ignoring values
strictly outside the sequence passed to 'limits'.   Note: either limit
in the sequence, or the value of limits itself, can be set to None.  The
inclusive list/tuple determines whether the lower and upper limiting bounds
(respectively) are open/exclusive (0) or closed/inclusive (1).

Usage:   atmean(a,limits=None,inclusive=(1,1))
"""
     if a.typecode() in ['l','s','b']:
         a = a.astype(N.Float)
     if limits == None:
	 return mean(a)
     assert type(limits) in [ListType,TupleType,N.ArrayType], "Wrong type for limits in atmean"
     if inclusive[0]:	 lowerfcn = N.greater_equal
     else:               lowerfcn = N.greater
     if inclusive[1]:	 upperfcn = N.less_equal
     else:               upperfcn = N.less
     if limits[0] > N.maximum.reduce(N.ravel(a)) or limits[1] < N.minimum.reduce(N.ravel(a)):
         raise ValueError, "No array values within given limits (atmean)."
     elif limits[0]==None and limits[1]<>None:
         mask = upperfcn(a,limits[1])
     elif limits[0]<>None and limits[1]==None:
         mask = lowerfcn(a,limits[0])
     elif limits[0]<>None and limits[1]<>None:
         mask = lowerfcn(a,limits[0])*upperfcn(a,limits[1])
     s = float(N.add.reduce(N.ravel(a*mask)))
     n = float(N.add.reduce(N.ravel(mask)))
     return s/n


 def atvar(a,limits=None,inclusive=(1,1)):
     """
Returns the sample variance of values in an array, (i.e., using N-1),
ignoring values strictly outside the sequence passed to 'limits'.  
Note: either limit in the sequence, or the value of limits itself,
can be set to None.  The inclusive list/tuple determines whether the lower
and upper limiting bounds (respectively) are open/exclusive (0) or
closed/inclusive (1).

Usage:   atvar(a,limits=None,inclusive=(1,1))
"""
     a = a.astype(N.Float)
     if limits == None or limits == [None,None]:
         term1 = N.add.reduce(N.ravel(a*a))
         n = float(len(N.ravel(a))) - 1
         term2 = N.add.reduce(N.ravel(a))**2 / n
	 return (term1 - term2) / n
     assert type(limits) in [ListType,TupleType,N.ArrayType], "Wrong type for limits in atvar"
     if inclusive[0]:	 lowerfcn = N.greater_equal
     else:               lowerfcn = N.greater
     if inclusive[1]:	 upperfcn = N.less_equal
     else:               upperfcn = N.less
     if limits[0] > N.maximum.reduce(N.ravel(a)) or limits[1] < N.minimum.reduce(N.ravel(a)):
         raise ValueError, "No array values within given limits (atvar)."
     elif limits[0]==None and limits[1]<>None:
         mask = upperfcn(a,limits[1])
     elif limits[0]<>None and limits[1]==None:
         mask = lowerfcn(a,limits[0])
     elif limits[0]<>None and limits[1]<>None:
         mask = lowerfcn(a,limits[0])*upperfcn(a,limits[1])
     term1 = N.add.reduce(N.ravel(a*a*mask))
     n = float(N.add.reduce(N.ravel(mask))) - 1
     term2 = N.add.reduce(N.ravel(a*mask))**2 / n
     return (term1 - term2) / n


 def atmin(a,lowerlimit=None,dimension=None,inclusive=1):
     """
Returns the minimum value of a, along dimension, including only values less
than (or equal to, if inclusive=1) lowerlimit.  If the limit is set to None,
all values in the array are used.

Usage:   atmin(a,lowerlimit=None,dimension=None,inclusive=1)
"""
     if inclusive:	 lowerfcn = N.greater
     else:               lowerfcn = N.greater_equal
     if dimension == None:
	 a = N.ravel(a)
	 dimension = 0
     if lowerlimit == None:
	 lowerlimit = N.minimum.reduce(N.ravel(a))-11
     biggest = N.maximum.reduce(N.ravel(a))
     ta = N.where(lowerfcn(a,lowerlimit),a,biggest)
     return N.minimum.reduce(ta,dimension)


 def atmax(a,upperlimit,dimension=None,inclusive=1):
     """
Returns the maximum value of a, along dimension, including only values greater
than (or equal to, if inclusive=1) upperlimit.  If the limit is set to None,
a limit larger than the max value in the array is used.

Usage:   atmax(a,upperlimit,dimension=None,inclusive=1)
"""
     if inclusive:	 upperfcn = N.less
     else:               upperfcn = N.less_equal
     if dimension == None:
	 a = N.ravel(a)
	 dimension = 0
     if upperlimit == None:
	 upperlimit = N.maximum.reduce(N.ravel(a))+1
     smallest = N.minimum.reduce(N.ravel(a))
     ta = N.where(upperfcn(a,upperlimit),a,smallest)
     return N.maximum.reduce(ta,dimension)


 def atstdev(a,limits=None,inclusive=(1,1)):
     """
Returns the standard deviation of all values in an array, ignoring values
strictly outside the sequence passed to 'limits'.   Note: either limit
in the sequence, or the value of limits itself, can be set to None.  The
inclusive list/tuple determines whether the lower and upper limiting bounds
(respectively) are open/exclusive (0) or closed/inclusive (1).

Usage:   atstdev(a,limits=None,inclusive=(1,1))
"""
     return N.sqrt(tvar(a,limits,inclusive))


 def atsem(a,limits=None,inclusive=(1,1)):
     """
Returns the standard error of the mean for the values in an array,
(i.e., using N for the denominator), ignoring values strictly outside
the sequence passed to 'limits'.   Note: either limit in the sequence,
or the value of limits itself, can be set to None.  The inclusive list/tuple
determines whether the lower and upper limiting bounds (respectively) are
open/exclusive (0) or closed/inclusive (1).

Usage:   atsem(a,limits=None,inclusive=(1,1))
"""
     sd = tstdev(a,limits,inclusive)
     if limits == None or limits == [None,None]:
         n = float(len(N.ravel(a)))
     assert type(limits) in [ListType,TupleType,N.ArrayType], "Wrong type for limits in atsem"
     if inclusive[0]:	 lowerfcn = N.greater_equal
     else:               lowerfcn = N.greater
     if inclusive[1]:	 upperfcn = N.less_equal
     else:               upperfcn = N.less
     if limits[0] > N.maximum.reduce(N.ravel(a)) or limits[1] < N.minimum.reduce(N.ravel(a)):
         raise ValueError, "No array values within given limits (atsem)."
     elif limits[0]==None and limits[1]<>None:
         mask = upperfcn(a,limits[1])
     elif limits[0]<>None and limits[1]==None:
         mask = lowerfcn(a,limits[0])
     elif limits[0]<>None and limits[1]<>None:
         mask = lowerfcn(a,limits[0])*upperfcn(a,limits[1])
     term1 = N.add.reduce(N.ravel(a*a*mask))
     n = float(N.add.reduce(N.ravel(mask)))
     return sd/math.sqrt(n)


#####################################
############  AMOMENTS  #############
#####################################

 def amoment(a,moment=1,dimension=None):
    """
Calculates the nth moment about the mean for a sample (defaults to the
1st moment).  Generally used to calculate coefficients of skewness and
kurtosis.  Dimension can equal None (ravel array first), an integer
(the dimension over which to operate), or a sequence (operate over
multiple dimensions).

Usage:   amoment(a,moment=1,dimension=None)
Returns: appropriate moment along given dimension
"""
    if dimension == None:
	a = N.ravel(a)
	dimension = 0
    if moment == 1:
	return 0.0
    else:
	mn = amean(a,dimension,1)  # 1=keepdims
	s = N.power((a-mn),moment)
	return amean(s,dimension)


 def avariation(a,dimension=None):
    """
Returns the coefficient of variation, as defined in CRC Standard
Probability and Statistics, p.6. Dimension can equal None (ravel array
first), an integer (the dimension over which to operate), or a
sequence (operate over multiple dimensions).

Usage:   avariation(a,dimension=None)
"""
    return 100.0*asamplestdev(a,dimension)/amean(a,dimension)


 def askew(a,dimension=None): 
    """ 
Returns the skewness of a distribution (normal ==> 0.0; >0 means extra
weight in left tail).  Use askewtest() to see if it's close enough.
Dimension can equal None (ravel array first), an integer (the
dimension over which to operate), or a sequence (operate over multiple
dimensions).

Usage:   askew(a, dimension=None)
Returns: skew of vals in a along dimension, returning ZERO where all vals equal
"""
    denom = N.power(amoment(a,2,dimension),1.5)
    zero = N.equal(denom,0)
    if type(denom) == N.ArrayType and asum(zero) <> 0:
	print "Number of zeros in askew: ",asum(zero)
    denom = denom + zero  # prevent divide-by-zero
    return N.where(zero, 0, amoment(a,3,dimension)/denom)


 def akurtosis(a,dimension=None):
    """
Returns the kurtosis of a distribution (normal ==> 3.0; >3 means
heavier in the tails, and usually more peaked).  Use akurtosistest()
to see if it's close enough.  Dimension can equal None (ravel array
first), an integer (the dimension over which to operate), or a
sequence (operate over multiple dimensions).

Usage:   akurtosis(a,dimension=None)
Returns: kurtosis of values in a along dimension, and ZERO where all vals equal
"""
    denom = N.power(amoment(a,2,dimension),2)
    zero = N.equal(denom,0)
    if type(denom) == N.ArrayType and asum(zero) <> 0:
	print "Number of zeros in akurtosis: ",asum(zero)
    denom = denom + zero  # prevent divide-by-zero
    return N.where(zero,0,amoment(a,4,dimension)/denom)


 def adescribe(inarray,dimension=None):
     """
Returns several descriptive statistics of the passed array.  Dimension
can equal None (ravel array first), an integer (the dimension over
which to operate), or a sequence (operate over multiple dimensions).

Usage:   adescribe(inarray,dimension=None)
Returns: n, (min,max), mean, standard deviation, skew, kurtosis
"""
     if dimension == None:
	 inarray = N.ravel(inarray)
	 dimension = 0
     n = inarray.shape[dimension]
     mm = (N.minimum.reduce(inarray),N.maximum.reduce(inarray))
     m = amean(inarray,dimension)
     sd = astdev(inarray,dimension)
     skew = askew(inarray,dimension)
     kurt = akurtosis(inarray,dimension)
     return n, mm, m, sd, skew, kurt


#####################################
########  NORMALITY TESTS  ##########
#####################################

 def askewtest(a,dimension=None):
    """
Tests whether the skew is significantly different from a normal
distribution.  Dimension can equal None (ravel array first), an
integer (the dimension over which to operate), or a sequence (operate
over multiple dimensions).

Usage:   askewtest(a,dimension=None)
Returns: z-score and 2-tail z-probability
"""
    if dimension == None:
	a = N.ravel(a)
	dimension = 0
    b2 = askew(a,dimension)
    n = float(a.shape[dimension])
    y = b2 * N.sqrt(((n+1)*(n+3)) / (6.0*(n-2)) )
    beta2 = ( 3.0*(n*n+27*n-70)*(n+1)*(n+3) ) / ( (n-2.0)*(n+5)*(n+7)*(n+9) )
    W2 = -1 + N.sqrt(2*(beta2-1))
    delta = 1/N.sqrt(N.log(N.sqrt(W2)))
    alpha = N.sqrt(2/(W2-1))
    y = N.where(N.equal(y,0),1,y)
    Z = delta*N.log(y/alpha + N.sqrt((y/alpha)**2+1))
    return Z, (1.0-zprob(Z))*2


 def akurtosistest(a,dimension=None):
    """
Tests whether a dataset has normal kurtosis (i.e.,
kurtosis=3(n-1)/(n+1)) Valid only for n>20.  Dimension can equal None
(ravel array first), an integer (the dimension over which to operate),
or a sequence (operate over multiple dimensions).

Usage:   akurtosistest(a,dimension=None)
Returns: z-score and 2-tail z-probability, returns 0 for bad pixels
"""
    if dimension == None:
	a = N.ravel(a)
	dimension = 0
    n = float(a.shape[dimension])
    if n<20:
	print "akurtosistest only valid for n>=20 ... continuing anyway, n=",n
    b2 = akurtosis(a,dimension)
    E = 3.0*(n-1) /(n+1)
    varb2 = 24.0*n*(n-2)*(n-3) / ((n+1)*(n+1)*(n+3)*(n+5))
    x = (b2-E)/N.sqrt(varb2)
    sqrtbeta1 = 6.0*(n*n-5*n+2)/((n+7)*(n+9)) * N.sqrt((6.0*(n+3)*(n+5))/
						       (n*(n-2)*(n-3)))
    A = 6.0 + 8.0/sqrtbeta1 *(2.0/sqrtbeta1 + N.sqrt(1+4.0/(sqrtbeta1**2)))
    term1 = 1 -2/(9.0*A)
    denom = 1 +x*N.sqrt(2/(A-4.0))
    denom = N.where(N.less(denom,0), 99, denom)
    term2 = N.where(N.equal(denom,0), term1, N.power((1-2.0/A)/denom,1/3.0))
    Z = ( term1 - term2 ) / N.sqrt(2/(9.0*A))
    Z = N.where(N.equal(denom,99), 0, Z)
    return Z, (1.0-zprob(Z))*2


 def anormaltest(a,dimension=None):
    """
Tests whether skew and/OR kurtosis of dataset differs from normal
curve.  Can operate over multiple dimensions.  Dimension can equal
None (ravel array first), an integer (the dimension over which to
operate), or a sequence (operate over multiple dimensions).

Usage:   anormaltest(a,dimension=None)
Returns: z-score and 2-tail probability
"""
    if dimension == None:
	a = N.ravel(a)
	dimension = 0
    s,p = askewtest(a,dimension)
    k,p = akurtosistest(a,dimension)
    k2 = N.power(s,2) + N.power(k,2)
    return k2, achisqprob(k2,2)


#####################################
######  AFREQUENCY FUNCTIONS  #######
#####################################

 def aitemfreq(a):
    """
Returns a 2D array of item frequencies.  Column 1 contains item values,
column 2 contains their respective counts.  Assumes a 1D array is passed.

Usage:   aitemfreq(a)
Returns: a 2D frequency table (col [0:n-1]=scores, col n=frequencies)
"""
    scores = pstat.aunique(a)
    scores = N.sort(scores)
    freq = N.zeros(len(scores))
    for i in range(len(scores)):
	freq[i] = N.add.reduce(N.equal(a,scores[i]))
    return N.array(pstat.aabut(scores, freq))


 def ascoreatpercentile (inarray, percent):
    """
Usage:   ascoreatpercentile(inarray,percent)   0<percent<100
Returns: score at given percentile, relative to inarray distribution
"""
    percent = percent / 100.0
    targetcf = percent*len(inarray)
    h, lrl, binsize, extras = histogram(inarray)
    cumhist = cumsum(h*1)
    for i in range(len(cumhist)):
	if cumhist[i] >= targetcf:
	    break
    score = binsize * ((targetcf - cumhist[i-1]) / float(h[i])) + (lrl+binsize*i)
    return score


 def apercentileofscore (inarray,score,histbins=10,defaultlimits=None):
    """
Note: result of this function depends on the values used to histogram
the data(!).

Usage:   apercentileofscore(inarray,score,histbins=10,defaultlimits=None)
Returns: percentile-position of score (0-100) relative to inarray
"""
    h, lrl, binsize, extras = histogram(inarray,histbins,defaultlimits)
    cumhist = cumsum(h*1)
    i = int((score - lrl)/float(binsize))
    pct = (cumhist[i-1]+((score-(lrl+binsize*i))/float(binsize))*h[i])/float(len(inarray)) * 100
    return pct


 def ahistogram (inarray,numbins=10,defaultlimits=None,printextras=1):
    """
Returns (i) an array of histogram bin counts, (ii) the smallest value
of the histogram binning, and (iii) the bin width (the last 2 are not
necessarily integers).  Default number of bins is 10.  Defaultlimits
can be None (the routine picks bins spanning all the numbers in the
inarray) or a 2-sequence (lowerlimit, upperlimit).  Returns all of the
following: array of bin values, lowerreallimit, binsize, extrapoints.

Usage:   ahistogram(inarray,numbins=10,defaultlimits=None,printextras=1)
Returns: (array of bin counts, bin-minimum, min-width, #-points-outside-range)
"""
    inarray = N.ravel(inarray)               # flatten any >1D arrays
    if (defaultlimits <> None):
	lowerreallimit = defaultlimits[0]
	upperreallimit = defaultlimits[1]
	binsize = (upperreallimit-lowerreallimit) / float(numbins)
    else:
	Min = N.minimum.reduce(inarray)
	Max = N.maximum.reduce(inarray)
	estbinwidth = float(Max - Min)/float(numbins) + 1
	binsize = (Max-Min+estbinwidth)/float(numbins)
	lowerreallimit = Min - binsize/2.0  #lower real limit,1st bin
    bins = N.zeros(numbins)
    extrapoints = 0
    for num in inarray:
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


 def acumfreq(a,numbins=10,defaultreallimits=None):
    """
Returns a cumulative frequency histogram, using the histogram function.
Defaultreallimits can be None (use all data), or a 2-sequence containing
lower and upper limits on values to include.

Usage:   acumfreq(a,numbins=10,defaultreallimits=None)
Returns: array of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(a,numbins,defaultreallimits)
    cumhist = cumsum(h*1)
    return cumhist,l,b,e


 def arelfreq(a,numbins=10,defaultreallimits=None):
    """
Returns a relative frequency histogram, using the histogram function.
Defaultreallimits can be None (use all data), or a 2-sequence containing
lower and upper limits on values to include.

Usage:   arelfreq(a,numbins=10,defaultreallimits=None)
Returns: array of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(a,numbins,defaultreallimits)
    h = N.array(h/float(a.shape[0]))
    return h,l,b,e


#####################################
######  AVARIABILITY FUNCTIONS  #####
#####################################

 def aobrientransform(*args):
    """
Computes a transform on input data (any number of columns).  Used to
test for homogeneity of variance prior to running one-way stats.  Each
array in *args is one level of a factor.  If an F_oneway() run on the
transformed data and found significant, variances are unequal.   From
Maxwell and Delaney, p.112.

Usage:   aobrientransform(*args)    *args = 1D arrays, one per level of factor
Returns: transformed data for use in an ANOVA
"""
    TINY = 1e-10
    k = len(args)
    n = N.zeros(k,N.Float)
    v = N.zeros(k,N.Float)
    m = N.zeros(k,N.Float)
    nargs = []
    for i in range(k):
	nargs.append(args[i].astype(N.Float))
	n[i] = float(len(nargs[i]))
	v[i] = var(nargs[i])
	m[i] = mean(nargs[i])
    for j in range(k):
	for i in range(n[j]):
	    t1 = (n[j]-1.5)*n[j]*(nargs[j][i]-m[j])**2
	    t2 = 0.5*v[j]*(n[j]-1.0)
	    t3 = (n[j]-1.0)*(n[j]-2.0)
	    nargs[j][i] = (t1-t2) / float(t3)
    check = 1
    for j in range(k):
	if v[j] - mean(nargs[j]) > TINY:
	    check = 0
    if check <> 1:
	raise ValueError, 'Lack of convergence in obrientransform.'
    else:
	return N.array(nargs)


 def asamplevar (inarray,dimension=None,keepdims=0):
    """
Returns the sample standard deviation of the values in the passed
array (i.e., using N).  Dimension can equal None (ravel array first),
an integer (the dimension over which to operate), or a sequence
(operate over multiple dimensions).  Set keepdims=1 to return an array
with the same number of dimensions as inarray.

Usage:   asamplevar(inarray,dimension=None,keepdims=0)
"""
    if dimension == None:
	inarray = N.ravel(inarray)
	dimension = 0
    if dimension == 1:
	mn = amean(inarray,dimension)[:,N.NewAxis]
    else:
	mn = amean(inarray,dimension,keepdims=1)
    deviations = inarray - mn 
    if type(dimension) == ListType:
	n = 1
	for d in dimension:
	    n = n*inarray.shape[d]
    else:
	n = inarray.shape[dimension]
    svar = ass(deviations,dimension,keepdims) / float(n)
    return svar


 def asamplestdev (inarray, dimension=None, keepdims=0):
    """
Returns the sample standard deviation of the values in the passed
array (i.e., using N).  Dimension can equal None (ravel array first),
an integer (the dimension over which to operate), or a sequence
(operate over multiple dimensions).  Set keepdims=1 to return an array
with the same number of dimensions as inarray.

Usage:   asamplestdev(inarray,dimension=None,keepdims=0)
"""
    return N.sqrt(asamplevar(inarray,dimension,keepdims))


 def asignaltonoise(instack,dimension=0):
    """
Calculates signal-to-noise.  Dimension can equal None (ravel array
first), an integer (the dimension over which to operate), or a
sequence (operate over multiple dimensions).

Usage:   asignaltonoise(instack,dimension=0):
Returns: array containing the value of (mean/stdev) along dimension,
         or 0 when stdev=0
"""
    m = mean(instack,dimension)
    sd = stdev(instack,dimension)
    return N.where(N.equal(sd,0),0,m/sd)


 def avar (inarray, dimension=None,keepdims=0):
    """
Returns the estimated population variance of the values in the passed
array (i.e., N-1).  Dimension can equal None (ravel array first), an
integer (the dimension over which to operate), or a sequence (operate
over multiple dimensions).  Set keepdims=1 to return an array with the
same number of dimensions as inarray.

Usage:   avar(inarray,dimension=None,keepdims=0)
"""
    if dimension == None:
	inarray = N.ravel(inarray)
	dimension = 0
    mn = amean(inarray,dimension,1)
    deviations = inarray - mn
    if type(dimension) == ListType:
	n = 1
	for d in dimension:
	    n = n*inarray.shape[d]
    else:
	n = inarray.shape[dimension]
    var = ass(deviations,dimension,keepdims)/float(n-1)
    return var


 def astdev (inarray, dimension=None, keepdims=0):
    """
Returns the estimated population standard deviation of the values in
the passed array (i.e., N-1).  Dimension can equal None (ravel array
first), an integer (the dimension over which to operate), or a
sequence (operate over multiple dimensions).  Set keepdims=1 to return
an array with the same number of dimensions as inarray.

Usage:   astdev(inarray,dimension=None,keepdims=0)
"""
    return N.sqrt(avar(inarray,dimension,keepdims))


 def asterr (inarray, dimension=None, keepdims=0):
    """
Returns the estimated population standard error of the values in the
passed array (i.e., N-1).  Dimension can equal None (ravel array
first), an integer (the dimension over which to operate), or a
sequence (operate over multiple dimensions).  Set keepdims=1 to return
an array with the same number of dimensions as inarray.

Usage:   asterr(inarray,dimension=None,keepdims=0)
"""
    if dimension == None:
	inarray = N.ravel(inarray)
	dimension = 0
    return astdev(inarray,dimension,keepdims) / float(N.sqrt(inarray.shape[dimension]))


 def asem (inarray, dimension=None, keepdims=0):
    """
Returns the standard error of the mean (i.e., using N) of the values
in the passed array.  Dimension can equal None (ravel array first), an
integer (the dimension over which to operate), or a sequence (operate
over multiple dimensions).  Set keepdims=1 to return an array with the
same number of dimensions as inarray.

Usage:   asem(inarray,dimension=None, keepdims=0)
"""
    if dimension == None:
	inarray = N.ravel(inarray)
	dimension = 0
    if type(dimension) == ListType:
	n = 1
	for d in dimension:
	    n = n*inarray.shape[d]
    else:
	n = inarray.shape[dimension]
    s = asamplestdev(inarray,dimension,keepdims) / N.sqrt(n-1)
    return s


 def az (a, score):
    """
Returns the z-score of a given input score, given thearray from which
that score came.  Not appropriate for population calculations, nor for
arrays > 1D.

Usage:   az(a, score)
"""
    z = (score-amean(a)) / asamplestdev(a)
    return z


 def azs (a):
    """
Returns a 1D array of z-scores, one for each score in the passed array,
computed relative to the passed array.

Usage:   azs(a)
"""
    zscores = []
    for item in a:
	zscores.append(z(a,item))
    return N.array(zscores)


 def azmap (scores, compare, dimension=0):
    """
Returns an array of z-scores the shape of scores (e.g., [x,y]), compared to
array passed to compare (e.g., [time,x,y]).  Assumes collapsing over dim 0
of the compare array.

Usage:   azs(scores, compare, dimension=0)
"""
    mns = amean(compare,dimension)
    sstd = asamplestdev(compare,0)
    return (scores - mns) / sstd


#####################################
#######  ATRIMMING FUNCTIONS  #######
#####################################

 def around(a,digits=1):
     """
Rounds all values in array a to 'digits' decimal places.

Usage:   around(a,digits)
Returns: a, where each value is rounded to 'digits' decimals
"""
     def ar(x,d=digits):
	 return round(x,d)

     if type(a) <> N.ArrayType:
	 try:
	     a = N.array(a)
	 except:
	     a = N.array(a,'O')
     shp = a.shape
     if a.typecode() in ['f','F','d','D']:
	 b = N.ravel(a)
	 b = N.array(map(ar,b))
	 b.shape = shp
     elif a.typecode() in ['o','O']:
	 b = N.ravel(a)*1
	 for i in range(len(b)):
	     if type(b[i]) == FloatType:
		 b[i] = round(b[i],digits)
	 b.shape = shp
     else:  # not a float, double or Object array
	 b = a*1
     return b


 def athreshold(a,threshmin=None,threshmax=None,newval=0):
    """
Like Numeric.clip() except that values <threshmid or >threshmax are replaced
by newval instead of by threshmin/threshmax (respectively).

Usage:   athreshold(a,threshmin=None,threshmax=None,newval=0)
Returns: a, with values <threshmin or >threshmax replaced with newval
"""
    mask = N.zeros(a.shape)
    if threshmin <> None:
	mask = mask + N.where(N.less(a,threshmin),1,0)
    if threshmax <> None:
	mask = mask + N.where(N.greater(a,threshmax),1,0)
    mask = N.clip(mask,0,1)
    return N.where(mask,newval,a)


 def atrimboth (a,proportiontocut):
    """
Slices off the passed proportion of items from BOTH ends of the passed
array (i.e., with proportiontocut=0.1, slices 'leftmost' 10% AND
'rightmost' 10% of scores.  You must pre-sort the array if you want
"proper" trimming.  Slices off LESS if proportion results in a
non-integer slice index (i.e., conservatively slices off
proportiontocut).

Usage:   atrimboth (a,proportiontocut)
Returns: trimmed version of array a
"""
    lowercut = int(proportiontocut*len(a))
    uppercut = len(a) - lowercut
    return a[lowercut:uppercut]


 def atrim1 (a,proportiontocut,tail='right'):
    """
Slices off the passed proportion of items from ONE end of the passed
array (i.e., if proportiontocut=0.1, slices off 'leftmost' or 'rightmost'
10% of scores).  Slices off LESS if proportion results in a non-integer
slice index (i.e., conservatively slices off proportiontocut).

Usage:   atrim1(a,proportiontocut,tail='right')  or set tail='left'
Returns: trimmed version of array a
"""
    if string.lower(tail) == 'right':
	lowercut = 0
	uppercut = len(a) - int(proportiontocut*len(a))
    elif string.lower(tail) == 'left':
	lowercut = int(proportiontocut*len(a))
	uppercut = len(a)
    return a[lowercut:uppercut]


#####################################
#####  ACORRELATION FUNCTIONS  ######
#####################################

 def acovariance(X):
    """
Computes the covariance matrix of a matrix X.  Requires a 2D matrix input.

Usage:   acovariance(X)
Returns: covariance matrix of X
"""
    if len(X.shape) <> 2:
	raise TypeError, "acovariance requires 2D matrices"
    n = X.shape[0]
    mX = amean(X,0)
    return N.dot(N.transpose(X),X) / float(n) - N.multiply.outer(mX,mX)


 def acorrelation(X):
    """
Computes the correlation matrix of a matrix X.  Requires a 2D matrix input.

Usage:   acorrelation(X)
Returns: correlation matrix of X
"""
    C = acovariance(X)
    V = N.diagonal(C)
    return C / N.sqrt(N.multiply.outer(V,V))


 def apaired(x,y):
    """
Interactively determines the type of data in x and y, and then runs the
appropriated statistic for paired group data.

Usage:   apaired(x,y)     x,y = the two arrays of values to be compared
Returns: appropriate statistic name, value, and probability
"""
    samples = ''
    while samples not in ['i','r','I','R','c','C']:
	print '\nIndependent or related samples, or correlation (i,r,c): ',
	samples = raw_input()

    if samples in ['i','I','r','R']:
	print '\nComparing variances ...',
# USE O'BRIEN'S TEST FOR HOMOGENEITY OF VARIANCE, Maxwell & delaney, p.112
	r = obrientransform(x,y)
	f,p = F_oneway(pstat.colex(r,0),pstat.colex(r,1))
	if p<0.05:
	    vartype='unequal, p='+str(round(p,4))
	else:
	    vartype='equal'
	print vartype
	if samples in ['i','I']:
	    if vartype[0]=='e':
		t,p = ttest_ind(x,y,None,0)
		print '\nIndependent samples t-test:  ', round(t,4),round(p,4)
	    else:
		if len(x)>20 or len(y)>20:
		    z,p = ranksums(x,y)
		    print '\nRank Sums test (NONparametric, n>20):  ', round(z,4),round(p,4)
		else:
		    u,p = mannwhitneyu(x,y)
		    print '\nMann-Whitney U-test (NONparametric, ns<20):  ', round(u,4),round(p,4)

	else:  # RELATED SAMPLES
	    if vartype[0]=='e':
		t,p = ttest_rel(x,y,0)
		print '\nRelated samples t-test:  ', round(t,4),round(p,4)
	    else:
		t,p = ranksums(x,y)
		print '\nWilcoxon T-test (NONparametric):  ', round(t,4),round(p,4)
    else:  # CORRELATION ANALYSIS
	corrtype = ''
	while corrtype not in ['c','C','r','R','d','D']:
	    print '\nIs the data Continuous, Ranked, or Dichotomous (c,r,d): ',
	    corrtype = raw_input()
	if corrtype in ['c','C']:
	    m,b,r,p,see = linregress(x,y)
	    print '\nLinear regression for continuous variables ...'
	    lol = [['Slope','Intercept','r','Prob','SEestimate'],[round(m,4),round(b,4),round(r,4),round(p,4),round(see,4)]]
	    pstat.printcc(lol)
	elif corrtype in ['r','R']:
	    r,p = spearmanr(x,y)
	    print '\nCorrelation for ranked variables ...'
	    print "Spearman's r: ",round(r,4),round(p,4)
	else: # DICHOTOMOUS
	    r,p = pointbiserialr(x,y)
	    print '\nAssuming x contains a dichotomous variable ...'
	    print 'Point Biserial r: ',round(r,4),round(p,4)
    print '\n\n'
    return None


 def apearsonr(x,y):
    """
Calculates a Pearson correlation coefficient and returns p.  Taken
from Heiman's Basic Statistics for the Behav. Sci (2nd), p.195.

Usage:   apearsonr(x,y)      where x,y are equal length arrays
Returns: Pearson's r, two-tailed p-value
"""
    TINY = 1.0e-20
    n = len(x)
    xmean = amean(x)
    ymean = amean(y)
    r_num = n*(N.add.reduce(x*y)) - N.add.reduce(x)*N.add.reduce(y)
    r_den = math.sqrt((n*ass(x) - asquare_of_sums(x))*(n*ass(y)-asquare_of_sums(y)))
    r = (r_num / r_den)
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = abetai(0.5*df,0.5,df/(df+t*t))
    return r,prob


 def aspearmanr(x,y):
    """
Calculates a Spearman rank-order correlation coefficient.  Taken
from Heiman's Basic Statistics for the Behav. Sci (1st), p.192.

Usage:   aspearmanr(x,y)      where x,y are equal-length arrays
Returns: Spearman's r, two-tailed p-value
"""
    TINY = 1e-30
    n = len(x)
    rankx = rankdata(x)
    ranky = rankdata(y)
    dsq = N.add.reduce((rankx-ranky)**2)
    rs = 1 - 6*dsq / float(n*(n**2-1))
    t = rs * math.sqrt((n-2) / ((rs+1.0)*(1.0-rs)))
    df = n-2
    probrs = abetai(0.5*df,0.5,df/(df+t*t))
# probability values for rs are from part 2 of the spearman function in
# Numerical Recipies, p.510.  They close to tables, but not exact.(?)
    return rs, probrs


 def apointbiserialr(x,y):
    """
Calculates a point-biserial correlation coefficient and the associated
probability value.  Taken from Heiman's Basic Statistics for the Behav.
Sci (1st), p.194.

Usage:   apointbiserialr(x,y)      where x,y are equal length arrays
Returns: Point-biserial r, two-tailed p-value
"""
    TINY = 1e-30
    categories = pstat.aunique(x)
    data = pstat.aabut(x,y)
    if len(categories) <> 2:
	raise ValueError, "Exactly 2 categories required (in x) for pointbiserialr()."
    else:   # there are 2 categories, continue
	codemap = pstat.aabut(categories,N.arange(2))
	recoded = pstat.arecode(data,codemap,0)
	x = pstat.alinexand(data,0,categories[0])
	y = pstat.alinexand(data,0,categories[1])
	xmean = amean(pstat.acolex(x,1))
	ymean = amean(pstat.acolex(y,1))
	n = len(data)
	adjust = math.sqrt((len(x)/float(n))*(len(y)/float(n)))
	rpb = (ymean - xmean)/asamplestdev(pstat.acolex(data,1))*adjust
	df = n-2
	t = rpb*math.sqrt(df/((1.0-rpb+TINY)*(1.0+rpb+TINY)))
	prob = abetai(0.5*df,0.5,df/(df+t*t))
	return rpb, prob


 def akendalltau(x,y):
    """
Calculates Kendall's tau ... correlation of ordinal data.  Adapted
from function kendl1 in Numerical Recipies.  Needs good test-cases.@@@

Usage:   akendalltau(x,y)
Returns: Kendall's tau, two-tailed p-value
"""
    n1 = 0
    n2 = 0
    iss = 0
    for j in range(len(x)-1):
	for k in range(j,len(y)):
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
		if (a1):
		    n1 = n1 + 1
		else:
		    n2 = n2 + 1
    tau = iss / math.sqrt(n1*n2)
    svar = (4.0*len(x)+10.0) / (9.0*len(x)*(len(x)-1))
    z = tau / math.sqrt(svar)
    prob = erfcc(abs(z)/1.4142136)
    return tau, prob


 def alinregress(*args):
    """
Calculates a regression line on two arrays, x and y, corresponding to x,y
pairs.  If a single 2D array is passed, alinregress finds dim with 2 levels
and splits data into x,y pairs along that dim.

Usage:   alinregress(*args)    args=2 equal-length arrays, or one 2D array
Returns: slope, intercept, r, two-tailed prob, sterr-of-the-estimate
"""
    TINY = 1.0e-20
    if len(args) == 1:  # more than 1D array?
	args = args[0]
	if len(args) == 2:
	    x = args[0]
	    y = args[1]
	else:
	    x = args[:,0]
	    y = args[:,1]
    else:
	x = args[0]
	y = args[1]
    n = len(x)
    xmean = amean(x)
    ymean = amean(y)
    r_num = n*(N.add.reduce(x*y)) - N.add.reduce(x)*N.add.reduce(y)
    r_den = math.sqrt((n*ass(x) - asquare_of_sums(x))*(n*ass(y)-asquare_of_sums(y)))
    r = r_num / r_den
    z = 0.5*math.log((1.0+r+TINY)/(1.0-r+TINY))
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = abetai(0.5*df,0.5,df/(df+t*t))
    slope = r_num / (float(n)*ass(x) - asquare_of_sums(x))
    intercept = ymean - slope*xmean
    sterrest = math.sqrt(1-r*r)*asamplestdev(y)
    return slope, intercept, r, prob, sterrest


#####################################
#####  AINFERENTIAL STATISTICS  #####
#####################################

 def attest_1samp(a,popmean,printit=0,name='Sample',writemode='a'):
    """
Calculates the t-obtained for the independent samples T-test on ONE group
of scores a, given a population mean.  If printit=1, results are printed
to the screen.  If printit='filename', the results are output to 'filename'
using the given writemode (default=append).  Returns t-value, and prob.

Usage:   attest_1samp(a,popmean,Name='Sample',printit=0,writemode='a')
Returns: t-value, two-tailed prob
"""
    if type(a) != N.ArrayType:
	a = N.array(a)
    x = amean(a)
    v = avar(a)
    n = len(a)
    df = n-1
    svar = ((n-1)*v) / float(df)
    t = (x-popmean)/math.sqrt(svar*(1.0/n))
    prob = abetai(0.5*df,0.5,df/(df+t*t))

    if printit <> 0:
	statname = 'Single-sample T-test.'
	outputpairedstats(printit,writemode,
			  'Population','--',popmean,0,0,0,
			  name,n,x,v,N.minimum.reduce(N.ravel(a)),
			  N.maximum.reduce(N.ravel(a)),
			  statname,t,prob)
    return t,prob


 def attest_ind (a, b, dimension=None, printit=0, name1='Samp1', name2='Samp2',writemode='a'):
    """
Calculates the t-obtained T-test on TWO INDEPENDENT samples of scores
a, and b.  From Numerical Recipies, p.483.  If printit=1, results are
printed to the screen.  If printit='filename', the results are output
to 'filename' using the given writemode (default=append).  Dimension
can equal None (ravel array first), or an integer (the dimension over
which to operate on a and b).

Usage:   attest_ind (a,b,dimension=None,printit=0,
		     Name1='Samp1',Name2='Samp2',writemode='a')
Returns: t-value, two-tailed p-value
"""
    if dimension == None:
	a = N.ravel(a)
	b = N.ravel(b)
	dimension = 0
    x1 = amean(a,dimension)
    x2 = amean(b,dimension)
    v1 = avar(a,dimension)
    v2 = avar(b,dimension)
    n1 = a.shape[dimension]
    n2 = b.shape[dimension]
    df = n1+n2-2
    svar = ((n1-1)*v1+(n2-1)*v2) / float(df)
    zerodivproblem = N.equal(svar,0)
    t = (x1-x2)/N.sqrt(svar*(1.0/n1 + 1.0/n2))  # N-D COMPUTATION HERE!!!!!!
    t = N.where(zerodivproblem,1.0,t)           # replace NaN t-values with 1.0
    probs = abetai(0.5*df,0.5,float(df)/(df+t*t))

    if type(t) == N.ArrayType:
	probs = N.reshape(probs,t.shape)
    if len(probs) == 1:
	probs = probs[0]
	
    if printit <> 0:
        if type(t) == N.ArrayType:
            t = t[0]
        if type(probs) == N.ArrayType:
            probs = probs[0]
	statname = 'Independent samples T-test.'
	outputpairedstats(printit,writemode,
			  name1,n1,x1,v1,N.minimum.reduce(N.ravel(a)),
			  N.maximum.reduce(N.ravel(a)),
			  name2,n2,x2,v2,N.minimum.reduce(N.ravel(b)),
			  N.maximum.reduce(N.ravel(b)),
			  statname,t,probs)
	return
    return t, probs


 def attest_rel (a,b,dimension=None,printit=0,name1='Samp1',name2='Samp2',writemode='a'):
    """
Calculates the t-obtained T-test on TWO RELATED samples of scores, a
and b.  From Numerical Recipies, p.483.  If printit=1, results are
printed to the screen.  If printit='filename', the results are output
to 'filename' using the given writemode (default=append).  Dimension
can equal None (ravel array first), or an integer (the dimension over
which to operate on a and b).

Usage:   attest_rel(a,b,dimension=None,printit=0,
		    name1='Samp1',name2='Samp2',writemode='a')
Returns: t-value, two-tailed p-value
"""
    if dimension == None:
	a = N.ravel(a)
	b = N.ravel(b)
	dimension = 0
    if len(a)<>len(b):
	raise ValueError, 'Unequal length arrays.'
    x1 = amean(a,dimension)
    x2 = amean(b,dimension)
    v1 = avar(a,dimension)
    v2 = avar(b,dimension)
    n = a.shape[dimension]
    df = float(n-1)
    d = (a-b).astype('d')

    denom = N.sqrt((n*N.add.reduce(d*d,dimension) - N.add.reduce(d,dimension)**2) /df)
    zerodivproblem = N.equal(denom,0)
    t = N.add.reduce(d,dimension) / denom      # N-D COMPUTATION HERE!!!!!!
    t = N.where(zerodivproblem,1.0,t)          # replace NaN t-values with 1.0
    t = N.where(zerodivproblem,1.0,t)           # replace NaN t-values with 1.0
    probs = abetai(0.5*df,0.5,float(df)/(df+t*t))
    if type(t) == N.ArrayType:
	probs = N.reshape(probs,t.shape)
    if len(probs) == 1:
	probs = probs[0]

    if printit <> 0:
	statname = 'Related samples T-test.'
	outputpairedstats(printit,writemode,
			  name1,n,x1,v1,N.minimum.reduce(N.ravel(a)),
			  N.maximum.reduce(N.ravel(a)),
			  name2,n,x2,v2,N.minimum.reduce(N.ravel(b)),
			  N.maximum.reduce(N.ravel(b)),
			  statname,t,probs)
	return
    return t, probs


 def achisquare(f_obs,f_exp=None):
    """
Calculates a one-way chi square for array of observed frequencies and returns
the result.  If no expected frequencies are given, the total N is assumed to
be equally distributed across all groups.

Usage:   achisquare(f_obs, f_exp=None)   f_obs = array of observed cell freq.
Returns: chisquare-statistic, associated p-value
"""

    k = len(f_obs)
    if f_exp == None:
	f_exp = N.array([sum(f_obs)/float(k)] * len(f_obs),N.Float)
    f_exp = f_exp.astype(N.Float)
    chisq = N.add.reduce((f_obs-f_exp)**2 / f_exp)
    return chisq, chisqprob(chisq, k-1)


 def aks_2samp (data1,data2):
    """
Computes the Kolmogorov-Smirnof statistic on 2 samples.  Modified from
Numerical Recipies in C, page 493.  Returns KS D-value, prob.  Not ufunc-
like.

Usage:   aks_2samp(data1,data2)  where data1 and data2 are 1D arrays
Returns: KS D-value, p-value
"""
    j1 = 0    # N.zeros(data1.shape[1:]) TRIED TO MAKE THIS UFUNC-LIKE
    j2 = 0    # N.zeros(data2.shape[1:])
    fn1 = 0.0 # N.zeros(data1.shape[1:],N.Float)
    fn2 = 0.0 # N.zeros(data2.shape[1:],N.Float)
    n1 = data1.shape[0]
    n2 = data2.shape[0]
    en1 = n1*1
    en2 = n2*1
    d = N.zeros(data1.shape[1:],N.Float)
    data1 = N.sort(data1,0)
    data2 = N.sort(data2,0)
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
        prob = aksprob((en+0.12+0.11/en)*N.fabs(d))
    except:
        prob = 1.0
    return d, prob


 def amannwhitneyu(x,y):
    """
Calculates a Mann-Whitney U statistic on the provided scores and
returns the result.  Use only when the n in each condition is < 20 and
you have 2 independent samples of ranks.  REMEMBER: Mann-Whitney U is
significant if the u-obtained is LESS THAN or equal to the critical
value of U.

Usage:   amannwhitneyu(x,y)     where x,y are arrays of values for 2 conditions
Returns: u-statistic, one-tailed p-value (i.e., p(z(U)))
"""
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(N.concatenate((x,y)))
    rankx = ranked[0:n1]       # get the x-ranks
    ranky = ranked[n1:]        # the rest are y-ranks
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


 def atiecorrect(rankvals):
    """
Tie-corrector for ties in Mann Whitney U and Kruskal Wallis H tests.
See Siegel, S. (1956) Nonparametric Statistics for the Behavioral
Sciences.  New York: McGraw-Hill.  Code adapted from |Stat rankind.c
code.

Usage:   atiecorrect(rankvals)
Returns: T correction factor for U or H
"""
    sorted,posn = ashellsort(N.array(rankvals))
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


 def aranksums(x,y):
    """
Calculates the rank sums statistic on the provided scores and returns
the result.

Usage:   aranksums(x,y)     where x,y are arrays of values for 2 conditions
Returns: z-statistic, two-tailed p-value
"""
    n1 = len(x)
    n2 = len(y)
    alldata = N.concatenate((x,y))
    ranked = arankdata(alldata)
    x = ranked[:n1]
    y = ranked[n1:]
    s = sum(x)
    expected = n1*(n1+n2+1) / 2.0
    z = (s - expected) / math.sqrt(n1*n2*(n1+n2+1)/12.0)
    prob = 2*(1.0 -zprob(abs(z)))
    return z, prob


 def awilcoxont(x,y):
    """
Calculates the Wilcoxon T-test for related samples and returns the
result.  A non-parametric T-test.

Usage:   awilcoxont(x,y)     where x,y are equal-length arrays for 2 conditions
Returns: t-statistic, two-tailed p-value
"""
    if len(x) <> len(y):
	raise ValueError, 'Unequal N in awilcoxont.  Aborting.'
    d = x-y
    d = N.compress(N.not_equal(d,0),d) # Keep all non-zero differences
    count = len(d)
    absd = abs(d)
    absranked = arankdata(absd)
    r_plus = 0.0
    r_minus = 0.0
    for i in range(len(absd)):
	if d[i] < 0:
	    r_minus = r_minus + absranked[i]
	else:
	    r_plus = r_plus + absranked[i]
    wt = min(r_plus, r_minus)
    mn = count * (count+1) * 0.25
    se =  math.sqrt(count*(count+1)*(2.0*count+1.0)/24.0)
    z = math.fabs(wt-mn) / se
    z = math.fabs(wt-mn) / se
    prob = 2*(1.0 -zprob(abs(z)))
    return wt, prob


 def akruskalwallish(*args):
    """
The Kruskal-Wallis H-test is a non-parametric ANOVA for 3 or more
groups, requiring at least 5 subjects in each group.  This function
calculates the Kruskal-Wallis H and associated p-value for 3 or more
independent samples.

Usage:   akruskalwallish(*args)     args are separate arrays for 3+ conditions
Returns: H-statistic (corrected for ties), associated p-value
"""
    assert len(args) == 3, "Need at least 3 groups in stats.akruskalwallish()"
    args = list(args)
    n = [0]*len(args)
    n = map(len,args)
    all = []
    for i in range(len(args)):
	all = all + args[i].tolist()
    ranked = rankdata(all)
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
	raise ValueError, 'All numbers are identical in akruskalwallish'
    h = h / float(T)
    return h, chisqprob(h,df)


 def afriedmanchisquare(*args):
    """
Friedman Chi-Square is a non-parametric, one-way within-subjects
ANOVA.  This function calculates the Friedman Chi-square test for
repeated measures and returns the result, along with the associated
probability value.  It assumes 3 or more repeated measures.  Only 3
levels requires a minimum of 10 subjects in the study.  Four levels
requires 5 subjects per level(??).

Usage:   afriedmanchisquare(*args)   args are separate arrays for 2+ conditions
Returns: chi-square statistic, associated p-value
"""
    k = len(args)
    if k < 3:
	raise ValueError, '\nLess than 3 levels.  Friedman test not appropriate.\n'
    n = len(args[0])
    data = apply(pstat.aabut,args)
    data = data.astype(N.Float)
    for i in range(len(data)):
	data[i] = arankdata(data[i])
    ssbn = asum(asum(args,1)**2)
    chisq = 12.0 / (k*n*(k+1)) * ssbn - 3*n*(k+1)
    return chisq, chisqprob(chisq,k-1)


#####################################
####  APROBABILITY CALCULATIONS  ####
#####################################

 def achisqprob(chisq,df):
    """
Returns the (1-tail) probability value associated with the provided chi-square
value and df.  Heavily modified from chisq.c in Gary Perlman's |Stat.  Can
handle multiple dimensions.

Usage:   achisqprob(chisq,df)    chisq=chisquare stat., df=degrees of freedom
"""
    BIG = 200.0
    def ex(x):
	BIG = 200.0
        exponents = N.where(N.less(x,-BIG),-BIG,x)
	return N.exp(exponents)

    if type(chisq) == N.ArrayType:
	arrayflag = 1
    else:
	arrayflag = 0
	chisq = N.array([chisq])
    if df < 1:
	return N.ones(chisq.shape,N.float)
    probs = N.zeros(chisq.shape,N.Float)
    probs = N.where(N.less_equal(chisq,0),1.0,probs)  # set prob=1 for chisq<0
    a = 0.5 * chisq
    if df > 1:
	y = ex(-a)
    if df%2 == 0:
	even = 1
	s = y*1
	s2 = s*1
    else:
	even = 0
	s = 2.0 * azprob(-N.sqrt(chisq))
	s2 = s*1
    if (df > 2):
	chisq = 0.5 * (df - 1.0)
	if even:
	    z = N.ones(probs.shape,N.Float)
	else:
	    z = 0.5 *N.ones(probs.shape,N.Float)
	if even:
	    e = N.zeros(probs.shape,N.Float)
	else:
	    e = N.log(N.sqrt(N.pi)) *N.ones(probs.shape,N.Float)
	c = N.log(a)
	mask = N.zeros(probs.shape)
	a_big = N.greater(a,BIG)
	a_big_frozen = -1 *N.ones(probs.shape,N.Float)
	totalelements = N.multiply.reduce(N.array(probs.shape))
	while asum(mask)<>totalelements:
	    e = N.log(z) + e
	    s = s + ex(c*z-a-e)
	    z = z + 1.0
	    print z, e, s
	    newmask = N.greater(z,chisq)
	    a_big_frozen = N.where(newmask*N.equal(mask,0)*a_big, s, a_big_frozen)
	    mask = N.clip(newmask+mask,0,1)
	if even:
	    z = N.ones(probs.shape,N.Float)
	    e = N.ones(probs.shape,N.Float)
	else:
	    z = 0.5 *N.ones(probs.shape,N.Float)
	    e = 1.0 / N.sqrt(N.pi) / N.sqrt(a) * N.ones(probs.shape,N.Float)
	c = 0.0
	mask = N.zeros(probs.shape)
	a_notbig_frozen = -1 *N.ones(probs.shape,N.Float)
	while asum(mask)<>totalelements:
	    e = e * (a/z.astype(N.Float))
	    c = c + e
	    z = z + 1.0
	    print '#2', z, e, c, s, c*y+s2
	    newmask = N.greater(z,chisq)
	    a_notbig_frozen = N.where(newmask*N.equal(mask,0)*(1-a_big),
				      c*y+s2, a_notbig_frozen)
	    mask = N.clip(newmask+mask,0,1)
	probs = N.where(N.equal(probs,1),1,
			N.where(N.greater(a,BIG),a_big_frozen,a_notbig_frozen))
	return probs
    else:
	return s


 def aerfcc(x):
    """
Returns the complementary error function erfc(x) with fractional error
everywhere less than 1.2e-7.  Adapted from Numerical Recipies.  Can
handle multiple dimensions.

Usage:   aerfcc(x)
"""
    z = abs(x)
    t = 1.0 / (1.0+0.5*z)
    ans = t * N.exp(-z*z-1.26551223 + t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
    return N.where(N.greater_equal(x,0), ans, 2.0-ans)


 def azprob(z):
    """
Returns the area under the normal curve 'to the left of' the given z value.
Thus, 
    for z<0, zprob(z) = 1-tail probability
    for z>0, 1.0-zprob(z) = 1-tail probability
    for any z, 2.0*(1.0-zprob(abs(z))) = 2-tail probability
Adapted from z.c in Gary Perlman's |Stat.  Can handle multiple dimensions.

Usage:   azprob(z)    where z is a z-value
"""
    def yfunc(y):
	x = (((((((((((((-0.000045255659 * y
			 +0.000152529290) * y -0.000019538132) * y
		       -0.000676904986) * y +0.001390604284) * y
		     -0.000794620820) * y -0.002034254874) * y
		   +0.006549791214) * y -0.010557625006) * y
		 +0.011630447319) * y -0.009279453341) * y
	       +0.005353579108) * y -0.002141268741) * y
	     +0.000535310849) * y +0.999936657524
	return x

    def wfunc(w):
	x = ((((((((0.000124818987 * w
		    -0.001075204047) * w +0.005198775019) * w
		  -0.019198292004) * w +0.059054035642) * w
		-0.151968751364) * w +0.319152932694) * w
	      -0.531923007300) * w +0.797884560593) * N.sqrt(w) * 2.0
	return x

    Z_MAX = 6.0    # maximum meaningful z-value
    x = N.zeros(z.shape,N.Float) # initialize
    y = 0.5 * N.fabs(z)
    x = N.where(N.less(y,1.0),wfunc(y*y),yfunc(y-2.0)) # get x's
    x = N.where(N.greater(y,Z_MAX*0.5),1.0,x)          # kill those with big Z
    prob = N.where(N.greater(z,0),(x+1)*0.5,(1-x)*0.5)
    return prob


 def aksprob(alam):
     """
Returns the probability value for a K-S statistic computed via ks_2samp.
Adapted from Numerical Recipies.  Can handle multiple dimensions.

Usage:   aksprob(alam)
"""
     if type(alam) == N.ArrayType:
	 frozen = -1 *N.ones(alam.shape,N.Float64)
	 alam = alam.astype(N.Float64)
	 arrayflag = 1
     else:
	 frozen = N.array(-1.)
	 alam = N.array(alam,N.Float64)
     mask = N.zeros(alam.shape)
     fac = 2.0 *N.ones(alam.shape,N.Float)
     sum = N.zeros(alam.shape,N.Float)
     termbf = N.zeros(alam.shape,N.Float)
     a2 = N.array(-2.0*alam*alam,N.Float64)
     totalelements = N.multiply.reduce(N.array(mask.shape))
     for j in range(1,201):
	 if asum(mask) == totalelements:
	     break
         exponents = (a2*j*j)
	 overflowmask = N.less(exponents,-746)
	 frozen = N.where(overflowmask,0,frozen)
	 mask = mask+overflowmask
	 term = fac*N.exp(exponents)
	 sum = sum + term
	 newmask = N.where(N.less_equal(abs(term),(0.001*termbf)) +
			   N.less(abs(term),1.0e-8*sum), 1, 0)
	 frozen = N.where(newmask*N.equal(mask,0), sum, frozen)
	 mask = N.clip(mask+newmask,0,1)
	 fac = -fac
	 termbf = abs(term)
     if arrayflag:
	 return N.where(N.equal(frozen,-1), 1.0, frozen)  # 1.0 if doesn't converge
     else:
	 return N.where(N.equal(frozen,-1), 1.0, frozen)[0]  # 1.0 if doesn't converge


 def afprob (dfnum, dfden, F):
    """
Returns the 1-tailed significance level (p-value) of an F statistic
given the degrees of freedom for the numerator (dfR-dfF) and the degrees
of freedom for the denominator (dfF).  Can handle multiple dims for F.

Usage:   afprob(dfnum, dfden, F)   where usually dfnum=dfbn, dfden=dfwn
"""
    if type(F) == N.ArrayType:
	return abetai(0.5*dfden, 0.5*dfnum, dfden/(1.0*dfden+dfnum*F))
    else:
	return abetai(0.5*dfden, 0.5*dfnum, dfden/float(dfden+dfnum*F))


 def abetacf(a,b,x):
    """
Evaluates the continued fraction form of the incomplete Beta function,
betai.  (Adapted from: Numerical Recipies in C.)  Can handle multiple
dimensions for x.

Usage:   abetacf(a,b,x)
"""
    ITMAX = 200
    EPS = 3.0e-7

    arrayflag = 1
    if type(x) == N.ArrayType:
	frozen = N.ones(x.shape,N.Float) *-1  #start out w/ -1s, should replace all
    else:
	arrayflag = 0
	frozen = N.array([-1])
	x = N.array([x])
    mask = N.zeros(x.shape)
    bm = az = am = 1.0
    qab = a+b
    qap = a+1.0
    qam = a-1.0
    bz = 1.0-qab*x/qap
    for i in range(ITMAX+1):
	if N.sum(N.ravel(N.equal(frozen,-1)))==0:
	    break
	em = float(i+1)
	tem = em + em
	d = em*(b-em)*x/((qam+tem)*(a+tem))
	ap = az + d*am
	bp = bz+d*bm
	d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
	app = ap+d*az
	bpp = bp+d*bz
	aold = az*1
	am = ap/bpp
	bm = bp/bpp
	az = app/bpp
	bz = 1.0
	newmask = N.less(abs(az-aold),EPS*abs(az))
	frozen = N.where(newmask*N.equal(mask,0), az, frozen)
	mask = N.clip(mask+newmask,0,1)
    noconverge = asum(N.equal(frozen,-1))
    if noconverge <> 0:
	print 'a or b too big, or ITMAX too small in Betacf for ',noconverge,' elements'
    if arrayflag:
	return frozen
    else:
	return frozen[0]


 def agammln(xx):
    """
Returns the gamma function of xx.
    Gamma(z) = Integral(0,infinity) of t^(z-1)exp(-t) dt.
Adapted from: Numerical Recipies in C.  Can handle multiple dims ... but
probably doesn't normally have to.

Usage:   agammln(xx)
"""
    coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516,
	     0.120858003e-2, -0.536382e-5]
    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5)*N.log(tmp)
    ser = 1.0
    for j in range(len(coeff)):
	x = x + 1
	ser = ser + coeff[j]/x
    return -tmp + N.log(2.50662827465*ser)


 def abetai(a,b,x):
    """
Returns the incomplete beta function:

    I-sub-x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)

where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
function of a.  The continued fraction formulation is implemented
here, using the betacf function.  (Adapted from: Numerical Recipies in
C.)  Can handle multiple dimensions.

Usage:   abetai(a,b,x)
"""
    TINY = 1e-15
    if type(a) == N.ArrayType:
	if asum(N.less(x,0)+N.greater(x,1)) <> 0:
	    raise ValueError, 'Bad x in abetai'
    x = N.where(N.equal(x,0),TINY,x)
    x = N.where(N.equal(x,1.0),1-TINY,x)

    bt = N.where(N.equal(x,0)+N.equal(x,1), 0, -1)
    exponents = ( gammln(a+b)-gammln(a)-gammln(b)+a*N.log(x)+b*
		  N.log(1.0-x) )
    # 746 (below) is the MAX POSSIBLE BEFORE OVERFLOW
    exponents = N.where(N.less(exponents,-740),-740,exponents)
    bt = N.exp(exponents)
    if type(x) == N.ArrayType:
	ans = N.where(N.less(x,(a+1)/(a+b+2.0)),
		      bt*abetacf(a,b,x)/float(a),
		      1.0-bt*abetacf(b,a,1.0-x)/float(b))
    else:
	if x<(a+1)/(a+b+2.0):
            ans = bt*abetacf(a,b,x)/float(a)
	else:
            ans = 1.0-bt*abetacf(b,a,1.0-x)/float(b)
    return ans


#####################################
#######  AANOVA CALCULATIONS  #######
#####################################

 import stats, LinearAlgebra, math, operator
 LA = LinearAlgebra

 def aglm(data,para):
    """
Calculates a linear model fit ... anova/ancova/lin-regress/t-test/etc. Taken
from:
    Peterson et al. Statistical limitations in functional neuroimaging
    I. Non-inferential methods and statistical models.  Phil Trans Royal Soc
    Lond B 354: 1239-1260.

Usage:   aglm(data,para)
Returns: statistic, p-value ???
"""
    if len(para) <> len(data):
	print "data and para must be same length in aglm"
	return
    n = len(para)
    p = pstat.aunique(para)
    x = N.zeros((n,len(p)))  # design matrix
    for l in range(len(p)):
	x[:,l] = N.equal(para,p[l])
    b = N.dot(N.dot(LA.inverse(N.dot(N.transpose(x),x)),  # i.e., b=inv(X'X)X'Y
		    N.transpose(x)),
	      data)
    diffs = (data - N.dot(x,b))
    s_sq = 1./(n-len(p)) * N.dot(N.transpose(diffs), diffs)

    if len(p) == 2:  # ttest_ind
	c = N.array([1,-1])
	df = n-2
	fact = asum(1.0/asum(x,0))  # i.e., 1/n1 + 1/n2 + 1/n3 ...
	t = N.dot(c,b) / N.sqrt(s_sq*fact)
	probs = abetai(0.5*df,0.5,float(df)/(df+t*t))
	return t, probs

 def aanova(data,effects=['A','B','C','D','E','F','G','H','I','J','K']):
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
        Axe mess of 'global' statements (particularly for Drestrict fcns)

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
    global Bscols, Bbetweens, SSlist, SSsources, DM, DN, Bwonly_sources, D
    global Bwithins, alleffects, alleffsources
    outputlist = []
    SSbtw = []
    SSbtwsources = []
    SSwb = []
    SSwbsources = []
    alleffects = []
    alleffsources = []
    SSlist = []
    SSsources = []

    print
    variables = 1       # this function only handles one measured variable
    if type(data)==N.ArrayType:
	data = data.tolist()

## Create a list of all unique values in each column, and a list of these Ns
    alluniqueslist = [0]*(len(data[0])-variables) # all cols but data cols
    Nlevels = [0]*(len(data[0])-variables)        # (as above)
    for column in range(len(Nlevels)): 
	alluniqueslist[column] = pstat.unique(pstat.colex(data,column))
	Nlevels[column] = len(alluniqueslist[column])

    Ncells = N.multiply.reduce(Nlevels[1:]) # total num cells (w/i AND btw)
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
    Nwlevels = N.take(N.array(Nlevels),Wscols) # no.lvls for each w/i subj fact
    Nbtwfactors = len(Bscols) - 1 # WASNfactors - Nwifactors + 1
    Nblevels = N.take(N.array(Nlevels),Bscols)

    Nwsources = 2**Nwifactors - 1 # num within-subject factor-combos
    Nbsources = Nallsources - Nwsources

    #
    # CALC M-VARIABLE (LIST) and Marray/Narray VARIABLES (ARRAY OF CELL MNS/NS)
    #
    # Eliminate replications for the same subject in same condition as well as
    # within-subject repetitions, keep as list
    M = pstat.collapse(data,Bscols,-1,None,None,mean)
    # Create an arrays of Nblevels shape (excl. subj dim)
    Marray = N.zeros(Nblevels[1:],'f')
    Narray = N.zeros(Nblevels[1:],'f')
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
    cdata = pstat.collapse(data,range(Nfactors+1),-1,None,None,mean)

    # Find a value that's not a data score with which to fill the array DA
    dummyval = -1
    datavals = pstat.colex(data,-1)
    while dummyval in datavals:  # find a value that's not a data score
        dummyval = dummyval - 1
    DA = N.ones(Nlevels,'f')*dummyval # create plenty of data-slots to fill

    if len(Bscols) == 1: # ie., if no btw-subj factors
        # 1 (below) needed because we need 2D array even w/ only 1 group of subjects
        subjslots = N.ones((Nsubjects,1)) 
    else: # create array to hold 1s (subj present) and 0s (subj absent)
        subjslots = N.zeros(Nblevels)
    for i in range(len(data)): # for every datapoint given as input
        idx = []
        for j in range(Nfactors+1): # get n-D bin idx for this datapoint
            new = alluniqueslist[j].index(data[i][j])
            idx.append(new)
        DA[idx] = data[i][-1] # put this data point in proper place in DA
        btwidx = N.take(idx,N.array(Bscols))
        subjslots[btwidx] = 1
    # DONE CREATING DATA ARRAY, DA ... #dims = numfactors+1, dim 0=subjects
    # dim -1=measured values, dummyval = values used to fill empty slots in DA

    # PREPARE FOR MAIN LOOP
    dcount = -1     # prepare for pre-increment of D-variable pointer
    Bwsources = []  # binary #s, each=source containing w/i subj factors
    Bwonly_sources = [] # binary #s, each=source of w/i-subj-ONLY factors
    D = N.zeros(Nwsources,N.PyObject) # one slot for each Dx,2**Nwifactors
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

            # Next, make a list of the above.  The list is a list of dimensions in DA
            # because DA has the same number of dimensions as there are factors
            # (including subjects), but with extra dummyval='-1' values the original
            # data array (assuming between-subj vars exist)
            Lwithinnonsource = makelist(Bwithinnonsource,Nfactors+1)

            # Collapse all non-source, w/i subj dims, FROM THE END (otherwise the
            # dim-numbers change as you collapse).  THIS WORKS BECAUSE WE'RE
            # COLLAPSING ACROSS W/I SUBJECT DIMENSIONS, WHICH WILL ALL HAVE THE
            # SAME SUBJ IN THE SAME ARRAY LOCATIONS (i.e., dummyvals will still exist
            # but should remain the same value through the amean() function
            for i in range(len(Lwithinnonsource)-1,-1,-1):
                dwsc = amean(dwsc,Lwithinnonsource[i])
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
            dvarshape = N.array(N.take(mns.shape,Lwithinsourcecol[1:])) -1
            idxarray = N.indices(dvarshape)
            newshape = N.array([idxarray.shape[0],
                                N.multiply.reduce(idxarray.shape[1:])])
            indxlist = N.swapaxes(N.reshape(idxarray,newshape),0,1)

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
                coeffmatrix = N.ones(mns.shape,N.Float) # fewer dims than DA (!!)
                # Make a list of dim #s that are both in source AND w/i subj fact, incl subj
                Wsourcecol = makelist(Bwscols&source,Nfactors+1)
                # Fill coeffmatrix with a complete set of coeffs (1 per w/i-source factor)
                for wfactor in range(len(Lwithinsourcecol[1:])):
                    #put correct coeff. axis as first axis, or "swap it up"
                    coeffmatrix = N.swapaxes(coeffmatrix,0,
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
                    coeffmatrix = N.swapaxes(coeffmatrix,0,
                                             Lwithinsourcecol[wfactor+1])

                # CALCULATE D VARIABLE
                scratch = coeffmatrix * mns
                # Collapse all dimensions EXCEPT subjects dim (dim 0)
                for j in range(len(coeffmatrix.shape[1:])):
                    scratch = N.add.reduce(scratch,1)
                if len(scratch.shape) == 1:
                    scratch.shape = list(scratch.shape)+[1]
                try:
                    # Tack this column onto existing ones
                    tmp = D[dcount].shape
                    D[dcount] = pstat.aabut(D[dcount],scratch)
                except AttributeError: # i.e., D[dcount]=integer/float
                    # If this is the first, plug it in
                    D[dcount] = scratch


            # Big long thing to create DMarray (list of DM variables) for this source
            variables = D[dcount].shape[1] # Num variables for this source
            tidx = range(1,len(subjslots.shape)) + [0] # [0] = Ss dim
            tsubjslots = N.transpose(subjslots,tidx) # put Ss in last dim
            DMarray = N.zeros(list(tsubjslots.shape[0:-1]) +
                              [variables],'f') # btw-subj dims, then vars
            DNarray = N.zeros(list(tsubjslots.shape[0:-1]) +
                              [variables],'f') # btw-subj dims, then vars
            idx = [0] *len(tsubjslots.shape[0:-1])
            idx[0] = -1
            loopcap = N.array(tsubjslots.shape[0:-1]) -1
            while incr(idx,loopcap) <> -1:
                DNarray[idx] = float(asum(tsubjslots[idx]))
                thismean =  (N.add.reduce(tsubjslots[idx] * # 1=subj dim
                                          N.transpose(D[dcount]),1) /
                             DNarray[idx])
                thismean = N.array(thismean,N.PyObject)
                DMarray[idx] = thismean
            DM[dcount] = DMarray
            DN[dcount] = DNarray

        #
        # DONE CREATING M AND D VARIABLES ... TIME FOR SOME SS WORK
        # DONE CREATING M AND D VARIABLES ... TIME FOR SOME SS WORK
        #
        if Bscols[1:] <> []:
            BNs = pstat.colex([Nlevels],Bscols[1:])
        else:
            BNs = [1]
            #
            # FIGURE OUT WHICH VARS TO RESTRICT, see p.680 (Maxwell&Delaney)
            #
            # BETWEEN-SUBJECTS VARIABLES ONLY, use M variable for analysis
            #
        if ((source-1) & Bwithins) == 0:  # btw-subjects vars only?
            sourcecols = makelist(source-1,Nfactors+1)

            # Determine cols (from input list) required for n-way interaction
            Lsource = makelist((Nallsources-1)&Bbetweens,Nfactors+1)
            # NOW convert this list of between-subject column numbers to a list of
            # DIMENSIONS in M, since M has fewer dims than the original data array
            # (assuming within-subj vars exist); Bscols has list of between-subj cols
            # from input list, the indices of which correspond to that var's loc'n in M
            btwcols = map(Bscols.index,Lsource)
            # Obviously-needed loop to get cell means is embedded in the collapse fcn, -1
            # represents last (measured-variable) column, None=std, 1=retain Ns

            hn = aharmonicmean(Narray,-1) # -1=unravel first

            # CALCULATE SSw ... SUBTRACT APPROPRIATE CELL MEAN FROM EACH SUBJ SCORE
            SSw = 0.0
            idxlist = pstat.unique(pstat.colex(M,btwcols))
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
            # DIMENSIONS in M, since M has fewer dims than the original data array
            # (assuming within-subj vars exist); Bscols has list of between-subj cols
            # from input list, the indices of which correspond to that var's loc'n in M
            btwsourcecols = (N.array(map(Bscols.index,Lsource))-1).tolist()

            # Average Marray and get harmonic means of Narray OVER NON-SOURCE DIMS
            Bbtwnonsourcedims = ~source & Bbetweens
            Lbtwnonsourcedims = makelist(Bbtwnonsourcedims,Nfactors+1)
            btwnonsourcedims = (N.array(map(Bscols.index,Lbtwnonsourcedims))-1).tolist()

    ## Average Marray over non-source dimensions (1=keep squashed dims)
            sourceMarray = amean(Marray,btwnonsourcedims,1)

    ## Calculate harmonic means for each level in source
            sourceNarray = aharmonicmean(Narray,btwnonsourcedims,1)

    ## Calc grand average (ga), used for ALL effects
            ga = asum((sourceMarray*sourceNarray)/
                            asum(sourceNarray))
            ga = N.reshape(ga,N.ones(len(Marray.shape)))

    ## If GRAND interaction, use harmonic mean of ALL cell Ns
            if source == Nallsources-1:
                sourceNarray = aharmonicmean(Narray)

    ## Calc all SUBSOURCES to be subtracted from sourceMarray (M&D p.320)
            sub_effects = 1.0 * ga # start with grand mean
            for subsource in range(3,source,2):
        ## Make a list of the non-subsource dimensions
                if subset(subsource-1,source-1):
                    sub_effects = (sub_effects +
                                   alleffects[alleffsources.index(subsource)])
        ## Calc this effect (a(j)'s, b(k)'s, ab(j,k)'s, whatever)
            effect = sourceMarray - sub_effects

        ## Save it so you don't have to calculate it again next time
            alleffects.append(effect)
            alleffsources.append(source)

    ## Calc and save sums of squares for this source
            SS = asum((effect**2 *sourceNarray) *
                      N.multiply.reduce(N.take(Marray.shape,btwnonsourcedims)))
        ## Save it so you don't have to calculate it again next time
            SSlist.append(SS)
            SSsources.append(source)

            collapsed = pstat.collapse(M,btwcols,-1,None,len,mean)
            # Obviously needed for-loop to get source cell-means embedded in collapse fcns
            contrastmns = pstat.collapse(collapsed,btwsourcecols,-2,sterr,len,mean)
            # Collapse again, this time SUMMING instead of averaging (to get cell Ns)
            contrastns = pstat.collapse(collapsed,btwsourcecols,-1,None,None,
					N.sum)
            # Collapse again, this time calculating harmonicmeans (for hns)
            contrasthns = pstat.collapse(collapsed,btwsourcecols,-1,None,None,
					 harmonicmean)
            # CALCULATE *BTW-SUBJ* dfnum, dfden
            sourceNs = pstat.colex([Nlevels],makelist(source-1,Nfactors+1))
            dfnum = N.multiply.reduce(N.ravel(N.array(sourceNs)-1))
            dfden = Nsubjects - N.multiply.reduce(N.ravel(BNs))

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
            workD = D[Bwonly_sources.index(sourcewithins)]  

            # CALCULATE Er, Ef
    ## Set up workD and subjslots for upcoming calcs
            if len(workD.shape)==1:
                workD = workD[:,N.NewAxis]
            if len(subjslots.shape)==1:
                subjslots = subjslots[:,N.NewAxis]

    ## Calculate full-model sums of squares
            ef = Dfull_model(workD,subjslots) # Uses cell-means model

            #
            # **ONLY** WITHIN-SUBJECT VARIABLES TO CONSIDER
            #
            if subset((source-1),Bwithins):
                # restrict grand mean, as per M&D p.680
                er = Drestrict_mean(workD,subjslots) 
        #
        # **BOTH** WITHIN- AND BETWEEN-SUBJECTS VARIABLES TO CONSIDER
        #
            else:
                er = Drestrict_source(workD,subjslots,source) + ef
            SSw = LA.determinant(ef)
            SS = LA.determinant(er) - SSw

        # CALCULATE *W/I-SUBJ* dfnum, dfden
            sourceNs = pstat.colex([Nlevels],makelist(source,Nfactors+1))
            # Calculation of dfnum is straightforward regardless
            dfnum = N.multiply.reduce(N.ravel(N.array(sourceNs)-1)[1:])
            # If only within-subject factors are involved, dfden is straightforward
            if subset(source-1,Bwithins):
                dfden = Nsubjects -N.multiply.reduce(N.ravel(BNs))-dfnum +1
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
                    p = workD.shape[1]
                except IndexError:
                    p = 1
                k = N.multiply.reduce(N.ravel(BNs))
                m = Nsubjects -1 -(p+k)/2.0
                d_en = float(p**2 + (k-1)**2 - 5)
                if d_en == 0.0:
                    s = 1.0
                else:
                    s = math.sqrt(((p*(k-1))**2-4) / d_en)
                dfden = m*s - dfnum/2.0 + 1

                # Given a within-between combined source, Wilk's Lambda is appropriate
                if LA.determinant(er) <> 0:
                    lmbda = LA.determinant(ef) / LA.determinant(er)
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
        adjsourcecols = N.array(makelist(source-1,Nfactors+1)) -1
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
        collapsed = pstat.collapse(cdata,Lsource,-1,sterr,len,mean)

        # First, get the list of level-combos for source cells
        prefixcols = range(len(collapsed[0][:-3]))
        outlist = pstat.colex(collapsed,prefixcols)
        # Start w/ factor names (A,B,C, or ones input to anova())
        eff = []
        for col in Lsource:
            eff.append(effects[col-1])
        # Add in the mean and N labels for printout
        for item in ['MEAN','STERR','N']:
            eff.append(item)
        # To the list of level-combos, abut the corresp. means and Ns
        outlist = pstat.abut(outlist,
                             map(round4,pstat.colex(collapsed,-3)),
                             map(round4,pstat.colex(collapsed,-2)),
                             map(round4,pstat.colex(collapsed,-1)))
        outlist = [eff] + outlist # add titles to the top of the list
        pstat.printcc(outlist)    # print it in customized columns
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
    pstat.printcc(title)
    print

    title = [['Effect','SS','DF','MS','F','p','sig']] + ['dashes']
    outputlist = title + outputlist
    pstat.printcc(outputlist)
    return


 def Dfull_model(workd,subjslots):
     """
     RESTRICTS NOTHING (i.e., FULL MODEL CALCULATION).  Subtracts D-variable
cell-mean for each between-subj group and then calculates the SS array.
     """
     workd = subtr_cellmeans(workd,subjslots)
     sserr = multivar_SScalc(workd)
     return sserr


 def Drestrict_mean(workd,subjslots):
     """
     RESTRICTS GRAND MEAN.  Subtracts D-variable cell-mean for each between-
subj group, and then adds back each D-variable's grand mean.
     """
     # subtract D-variable cell-mean for each (btw-subj) group
     errors = subtr_cellmeans(workd,subjslots)

     # add back in appropriate grand mean from individual scores
     grandDmeans = amean(workd,0,1)
     errors = errors + N.transpose(grandDmeans) # errors has reversed dims!!
     # SS for mean-restricted model is calculated below.  Note: already put
     # subj as last dim because later code expects this code here to leave
     # workd that way
     sserr = multivar_SScalc(errors)
     return sserr


 def Drestrict_source(workd,subjslots,source):
     """
Calculates error for a given model on array workd.  Subjslots is an
array of 1s and 0s corresponding to whether or not the subject is a
member of that between-subjects variable combo.  source is the code
for the type of model to calculate.  source=-1 means no restriction;
source=0 means to restrict workd's grand mean; source>0 means to
restrict the columns of the main data array, DA, specified (in binary)
by the source-value.

Usage:   Derrorcalc(workd,subjslots,source)  source:-1=nothing, 0=mean
Returns: SS array for multivariate F calculation
"""
###
### RESTRICT COLUMNS/DIMENSIONS SPECIFIED IN source (BINARY)
### (i.e., is the value of source not equal to 0 or -1?)
###
     if source > 0:
	 sourcewithins = (source-1) & Bwithins
	 sourcebetweens = (source-1) & Bbetweens
	 dindex = Bwonly_sources.index(sourcewithins)
	 all_cellmeans = N.transpose(DM[dindex],[-1]+range(0,len(DM[dindex].shape)-1))
	 all_cellns = N.transpose(DN[dindex],[-1]+range(0,len(DN[dindex].shape)-1))
	 hn = aharmonicmean(all_cellns)

	 levels = D[dindex].shape[1]  # GENERAL, 'cause each workd is always 2D
	 SSm = N.zeros((levels,levels),'f') #called RCm=SCm in Lindman,p.317-8
	 tworkd = N.transpose(D[dindex])

     ## Calculate SSw, within-subj variance (Lindman approach)
	 RSw = N.zeros((levels,levels),'f')
	 RSinter = N.zeros((levels,levels),N.PyObject)	
	 for i in range(levels):
	     for j in range(i,levels):
		 RSw[i,j] = RSw[j,i] = N.sum(tworkd[i]*tworkd[j])
		 cross = all_cellmeans[i] * all_cellmeans[j]
		 multfirst = asum(cross*all_cellns[i])
		 RSinter[i,j] = RSinter[j,i] = N.asarray(multfirst)
		 SSm[i,j] = SSm[j,i] = (amean(all_cellmeans[i]) *
					amean(all_cellmeans[j]) *
					len(all_cellmeans[i]) *hn)
	 SSw = RSw - RSinter

### HERE BEGINS THE MAXWELL & DELANEY APPROACH TO CALCULATING SS
	 Lsource = makelist(sourcebetweens,Nfactors+1)
	 btwsourcecols = (N.array(map(Bscols.index,Lsource))-1).tolist()
	 Bbtwnonsourcedims = ~source & Bbetweens
	 Lbtwnonsourcedims = makelist(Bbtwnonsourcedims,Nfactors+1)
	 btwnonsourcedims = (N.array(map(Bscols.index,Lbtwnonsourcedims))-1).tolist()

       ## Average Marray over non-source dimensions
	 sourceDMarray = DM[dindex] *1.0
	 for dim in btwnonsourcedims: # collapse all non-source dims
	     if dim == len(DM[dindex].shape)-1:
		 raise ValueError, "Crashing ... shouldn't ever collapse ACROSS variables"
	     sourceDMarray = amean(sourceDMarray,dim,1)

       ## Calculate harmonic means for each level in source
	 sourceDNarray = aharmonicmean(DN[dindex],btwnonsourcedims,1)

       ## Calc grand average (ga), used for ALL effects
	 variableNs = asum(sourceDNarray,
                           range(len(sourceDMarray.shape)-1))
	 ga = asum((sourceDMarray*sourceDNarray) /
                   variableNs,
                   range(len(sourceDMarray.shape)-1),1)

       ## If GRAND interaction, use harmonic mean of ALL cell Ns
	 if source == Nallsources-1:
	     sourceDNarray = aharmonicmean(DN[dindex],
                                           range(len(sourceDMarray.shape)-1))
		
       ## Calc all SUBSOURCES to be subtracted from sourceMarray (M&D p.320)
	 sub_effects = ga *1.0   # start with grand mean
	 for subsource in range(3,source-2,2):
       ## Make a list of the non-subsource dimensions
	     subsourcebtw = (subsource-1) & Bbetweens
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
	 SS = N.zeros((levels,levels),'f')
	 SS = asum((effect**2 *sourceDNarray) *
                   N.multiply.reduce(N.take(DM[dindex].shape,btwnonsourcedims)),
			 range(len(sourceDMarray.shape)-1))
       ## Save it so you don't have to calculate it again next time
	 SSlist.append(SS)
	 SSsources.append(source)

	 return SS


 def multivar_SScalc(workd):
###
### DO SS CALCS ON THE OUTPUT FROM THE SOURCE=0 AND SOURCE=-1 CASES
###
     # this section expects workd to have subj. in LAST dimension!!!!!!
     if len(workd.shape) == 1:
	 levels = 1
     else:
	 levels = workd.shape[0] # works because workd is always 2D

     sserr = N.zeros((levels,levels),'f')
     for i in range(levels):
	 for j in range(i,levels):
	     ssval = N.add.reduce(workd[i]*workd[j])
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
     tsubjslots = N.transpose(subjslots,transidx) # get all Ss for this idx
     tworkd = N.transpose(workd) # swap subj. and variable dims
     errors = 1.0 * tworkd

     if len(sourcedims) == 0:
	 idx = [-1]
	 loopcap = [0]
     if len(sourcedims) <> 0:
	 btwsourcedims = map(Bscols.index,sourcedims)
	 idx = [0] * len(btwsourcedims)
	 idx[0] = -1 # compensate for pre-increment of 1st slot in incr()
	     
	 # Get a list of the maximum values each factor can handle
	 loopcap = N.take(N.array(Nlevels),sourcedims)-1

### WHILE STILL MORE GROUPS, CALCULATE GROUP MEAN FOR EACH D-VAR
     while incr(idx,loopcap) <> -1:  # loop through source btw level-combos
	 mask = tsubjslots[idx]
	 thisgroup = tworkd*mask[N.NewAxis,:]
	 groupmns = amean(N.compress(mask,thisgroup),1)

### THEN SUBTRACT THEM FROM APPROPRIATE SUBJECTS
	 errors = errors - N.multiply.outer(groupmns,mask)
     return errors


 def F_value_wilks_lambda(ER, EF, dfnum, dfden, a, b):
     """
Calculation of Wilks lambda F-statistic for multivarite data, per
Maxwell & Delaney p.657.

Usage:   F_value_wilks_lambda(ER,EF,dfnum,dfden,a,b)
"""
     if type(ER) in [IntType, FloatType]:
	 ER = N.array([[ER]])
     if type(EF) in [IntType, FloatType]:
	 EF = N.array([[EF]])
     lmbda = LA.determinant(EF) / LA.determinant(ER)
     if (a-1)**2 + (b-1)**2 == 5:
	 q = 1
     else:
	 q = math.sqrt( ((a-1)**2*(b-1)**2 - 2) / ((a-1)**2 + (b-1)**2 -5) )
     n_um = (1 - lmbda**(1.0/q))*(a-1)*(b-1)
     d_en = lmbda**(1.0/q) / (m*q - 0.5*(a-1)*(b-1) + 1)
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
	 return round(num,4)
     except:
	 return 'N/A'


 def aF_oneway(*args):
    """
Performs a 1-way ANOVA, returning an F-value and probability given
any number of groups.  From Heiman, pp.394-7.

Usage:   aF_oneway (*args)    where *args is 2 or more arrays, one per
                                  treatment group
Returns: f-value, probability
"""
    na = len(args)            # ANOVA on 'na' groups, each in it's own array
    means = [0]*na
    vars = [0]*na
    ns = [0]*na
    alldata = []
    tmp = map(N.array,args)
    means = map(amean,tmp)
    vars = map(avar,tmp)
    ns = map(len,args)
    alldata = N.concatenate(args)
    bign = len(alldata)
    sstot = ass(alldata)-(asquare_of_sums(alldata)/float(bign))
    ssbn = 0
    for a in args:
	ssbn = ssbn + asquare_of_sums(N.array(a))/float(len(a))
    ssbn = ssbn - (asquare_of_sums(alldata)/float(bign))
    sswn = sstot-ssbn
    dfbn = na-1
    dfwn = bign - na
    msb = ssbn/float(dfbn)
    msw = sswn/float(dfwn)
    f = msb/msw
    prob = fprob(dfbn,dfwn,f)
    return f, prob


 def aF_value (ER,EF,dfR,dfF):
    """
Returns an F-statistic given the following:
	ER  = error associated with the null hypothesis (the Restricted model)
	EF  = error associated with the alternate hypothesis (the Full model)
	dfR = degrees of freedom the Restricted model
	dfF = degrees of freedom associated with the Restricted model
"""
    return ((ER-EF)/float(dfR-dfF) / (EF/float(dfF)))


 def outputfstats(Enum, Eden, dfnum, dfden, f, prob):
     Enum = round(Enum,3)
     Eden = round(Eden,3)
     dfnum = round(Enum,3)
     dfden = round(dfden,3)
     f = round(f,3)
     prob = round(prob,3)
     suffix = ''                       # for *s after the p-value
     if  prob < 0.001:  suffix = '  ***'
     elif prob < 0.01:  suffix = '  **'
     elif prob < 0.05:  suffix = '  *'
     title = [['EF/ER','DF','Mean Square','F-value','prob','']]
     lofl = title+[[Enum, dfnum, round(Enum/float(dfnum),3), f, prob, suffix],
		   [Eden, dfden, round(Eden/float(dfden),3),'','','']]
     pstat.printcc(lofl)
     return


 def F_value_multivariate(ER, EF, dfnum, dfden):
     """
Returns an F-statistic given the following:
	ER  = error associated with the null hypothesis (the Restricted model)
	EF  = error associated with the alternate hypothesis (the Full model)
	dfR = degrees of freedom the Restricted model
	dfF = degrees of freedom associated with the Restricted model
where ER and EF are matrices from a multivariate F calculation.
"""
     if type(ER) in [IntType, FloatType]:
	 ER = N.array([[ER]])
     if type(EF) in [IntType, FloatType]:
	 EF = N.array([[EF]])
     n_um = (LA.determinant(ER) - LA.determinant(EF)) / float(dfnum)
     d_en = LA.determinant(EF) / float(dfden)
     return n_um / d_en


#####################################
#######  ASUPPORT FUNCTIONS  ########
#####################################

 def asign(a):
    """
Usage:   asign(a)
Returns: array shape of a, with -1 where a<0 and +1 where a>=0
"""
    a = N.asarray(a)
    if ((type(a) == type(1.4)) or (type(a) == type(1))):
	return a-a-N.less(a,0)+N.greater(a,0)
    else:
	return N.zeros(N.shape(a))-N.less(a,0)+N.greater(a,0)


 def asum (a, dimension=None,keepdims=0):
     """
An alternative to the Numeric.add.reduce function, which allows one to
(1) collapse over multiple dimensions at once, and/or (2) to retain
all dimensions in the original array (squashing one down to size.
Dimension can equal None (ravel array first), an integer (the
dimension over which to operate), or a sequence (operate over multiple
dimensions).  If keepdims=1, the resulting array will have as many
dimensions as the input array.

Usage:   asum(a, dimension=None, keepdims=0)
Returns: array summed along 'dimension'(s), same _number_ of dims if keepdims=1
"""
     if type(a) == N.ArrayType and a.typecode() in ['l','s','b']:
	 a = a.astype(N.Float)
     if dimension == None:
	 s = N.sum(N.ravel(a))
     elif type(dimension) in [IntType,FloatType]:
	 s = N.add.reduce(a, dimension)
	 if keepdims == 1:
	     shp = list(a.shape)
	     shp[dimension] = 1
	     s = N.reshape(s,shp)
     else: # must be a SEQUENCE of dims to sum over
	dims = list(dimension)
	dims.sort()
	dims.reverse()
	s = a *1.0
	for dim in dims:
	    s = N.add.reduce(s,dim)
	if keepdims == 1:
	    shp = list(a.shape)
	    for dim in dims:
		shp[dim] = 1
	    s = N.reshape(s,shp)
     return s


 def acumsum (a,dimension=None):
    """
Returns an array consisting of the cumulative sum of the items in the
passed array.  Dimension can equal None (ravel array first), an
integer (the dimension over which to operate), or a sequence (operate
over multiple dimensions, but this last one just barely makes sense).

Usage:   acumsum(a,dimension=None)
"""
    if dimension == None:
	a = N.ravel(a)
	dimension = 0
    if type(dimension) in [ListType, TupleType, N.ArrayType]:
	dimension = list(dimension)
	dimension.sort()
	dimension.reverse()
	for d in dimension:
	    a = N.add.accumulate(a,d)
	return a
    else:
	return N.add.accumulate(a,dimension)


 def ass(inarray, dimension=None, keepdims=0):
    """
Squares each value in the passed array, adds these squares & returns
the result.  Unfortunate function name. :-) Defaults to ALL values in
the array.  Dimension can equal None (ravel array first), an integer
(the dimension over which to operate), or a sequence (operate over
multiple dimensions).  Set keepdims=1 to maintain the original number
of dimensions.

Usage:   ass(inarray, dimension=None, keepdims=0)
Returns: sum-along-'dimension' for (inarray*inarray)
"""
    if dimension == None:
	inarray = N.ravel(inarray)
	dimension = 0
    return asum(inarray*inarray,dimension,keepdims)


 def asummult (array1,array2,dimension=None,keepdims=0):
    """
Multiplies elements in array1 and array2, element by element, and
returns the sum (along 'dimension') of all resulting multiplications.
Dimension can equal None (ravel array first), an integer (the
dimension over which to operate), or a sequence (operate over multiple
dimensions).  A trivial function, but included for completeness.

Usage:   asummult(array1,array2,dimension=None,keepdims=0)
"""
    if dimension == None:
	array1 = N.ravel(array1)
	array2 = N.ravel(array2)
	dimension = 0
    return asum(array1*array2,dimension,keepdims)


 def asquare_of_sums(inarray, dimension=None, keepdims=0):
    """
Adds the values in the passed array, squares that sum, and returns the
result.  Dimension can equal None (ravel array first), an integer (the
dimension over which to operate), or a sequence (operate over multiple
dimensions).  If keepdims=1, the returned array will have the same
NUMBER of dimensions as the original.

Usage:   asquare_of_sums(inarray, dimension=None, keepdims=0)
Returns: the square of the sum over dim(s) in dimension
"""
    if dimension == None:
	inarray = N.ravel(inarray)
	dimension = 0
    s = asum(inarray,dimension,keepdims)
    if type(s) == N.ArrayType:
        return s.astype(N.Float)*s
    else:
        return float(s)*s


 def asumdiffsquared(a,b, dimension=None, keepdims=0):
    """
Takes pairwise differences of the values in arrays a and b, squares
these differences, and returns the sum of these squares.  Dimension
can equal None (ravel array first), an integer (the dimension over
which to operate), or a sequence (operate over multiple dimensions).
keepdims=1 means the return shape = len(a.shape) = len(b.shape)

Usage:   asumdiffsquared(a,b)
Returns: sum[ravel(a-b)**2]
"""
    if dimension == None:
	inarray = N.ravel(a)
	dimension = 0
    return asum((a-b)**2,dimension,keepdims)


 def ashellsort(inarray):
    """
Shellsort algorithm.  Sorts a 1D-array.

Usage:   ashellsort(inarray)
Returns: sorted-inarray, sorting-index-vector (for original array)
"""
    n = len(inarray)
    svec = inarray *1.0
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
    return svec, ivec


 def arankdata(inarray):
    """
Ranks the data in inarray, dealing with ties appropritely.  Assumes
a 1D inarray.  Adapted from Gary Perlman's |Stat ranksort.

Usage:   arankdata(inarray)
Returns: array of length equal to inarray, containing rank scores
"""
    n = len(inarray)
    svec, ivec = ashellsort(inarray)
    sumranks = 0
    dupcount = 0
    newarray = N.zeros(n,N.Float)
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


 def afindwithin(data):
    """
Returns a binary vector, 1=within-subject factor, 0=between.  Input
equals the entire data array (i.e., column 1=random factor, last
column = measured values.

Usage:   afindwithin(data)     data in |Stat format
"""
    numfact = len(data[0])-2
    withinvec = [0]*numfact
    for col in range(1,numfact+1):
	rows = pstat.linexand(data,col,pstat.unique(pstat.colex(data,1))[0])  # get 1 level of this factor
	if len(pstat.unique(pstat.colex(rows,0))) < len(rows):   # if fewer subjects than scores on this factor
	    withinvec[col-1] = 1
    return withinvec


 #########################################################
 #########################################################
 ######  RE-DEFINE DISPATCHES TO INCLUDE ARRAYS  #########
 #########################################################
 #########################################################

## CENTRAL TENDENCY:
 geometricmean = Dispatch ( (lgeometricmean, (ListType, TupleType)),
                            (ageometricmean, (N.ArrayType,)) )
 harmonicmean = Dispatch ( (lharmonicmean, (ListType, TupleType)),
                           (aharmonicmean, (N.ArrayType,)) )
 mean = Dispatch ( (lmean, (ListType, TupleType)),
                   (amean, (N.ArrayType,)) )
 median = Dispatch ( (lmedian, (ListType, TupleType)),
                     (amedian, (N.ArrayType,)) )
 medianscore = Dispatch ( (lmedianscore, (ListType, TupleType)),
                          (amedianscore, (N.ArrayType,)) )
 mode = Dispatch ( (lmode, (ListType, TupleType)),
                   (amode, (N.ArrayType,)) )
 tmean = Dispatch ( (atmean, (N.ArrayType,)) )
 tvar = Dispatch ( (atvar, (N.ArrayType,)) )
 tstdev = Dispatch ( (atstdev, (N.ArrayType,)) )
 tsem = Dispatch ( (atsem, (N.ArrayType,)) )

## VARIATION:
 moment = Dispatch ( (lmoment, (ListType, TupleType)),
                     (amoment, (N.ArrayType,)) )
 variation = Dispatch ( (lvariation, (ListType, TupleType)),
                        (avariation, (N.ArrayType,)) )
 skew = Dispatch ( (lskew, (ListType, TupleType)),
                   (askew, (N.ArrayType,)) )
 kurtosis = Dispatch ( (lkurtosis, (ListType, TupleType)),
                       (akurtosis, (N.ArrayType,)) )
 describe = Dispatch ( (ldescribe, (ListType, TupleType)),
                       (adescribe, (N.ArrayType,)) )

## DISTRIBUTION TESTS

 skewtest = Dispatch ( (askewtest, (ListType, TupleType)),
		       (askewtest, (N.ArrayType,)) )
 kurtosistest = Dispatch ( (akurtosistest, (ListType, TupleType)),
			   (akurtosistest, (N.ArrayType,)) )
 normaltest = Dispatch ( (anormaltest, (ListType, TupleType)),
			 (anormaltest, (N.ArrayType,)) )

## FREQUENCY STATS:
 itemfreq = Dispatch ( (litemfreq, (ListType, TupleType)),
		       (aitemfreq, (N.ArrayType,)) )
 scoreatpercentile = Dispatch ( (lscoreatpercentile, (ListType, TupleType)),
                                (ascoreatpercentile, (N.ArrayType,)) )
 percentileofscore = Dispatch ( (lpercentileofscore, (ListType, TupleType)),
                                 (apercentileofscore, (N.ArrayType,)) )
 histogram = Dispatch ( (lhistogram, (ListType, TupleType)),
                        (ahistogram, (N.ArrayType,)) )
 cumfreq = Dispatch ( (lcumfreq, (ListType, TupleType)),
                      (acumfreq, (N.ArrayType,)) )
 relfreq = Dispatch ( (lrelfreq, (ListType, TupleType)),
                      (arelfreq, (N.ArrayType,)) )
 
## VARIABILITY:
 obrientransform = Dispatch ( (lobrientransform, (ListType, TupleType)),
			      (aobrientransform, (N.ArrayType,)) )
 samplevar = Dispatch ( (lsamplevar, (ListType, TupleType)),
                        (asamplevar, (N.ArrayType,)) )
 samplestdev = Dispatch ( (lsamplestdev, (ListType, TupleType)),
                          (asamplestdev, (N.ArrayType,)) )
 signaltonoise = Dispatch( (asignaltonoise, (N.ArrayType,)),)
 var = Dispatch ( (lvar, (ListType, TupleType)),
                  (avar, (N.ArrayType,)) )
 stdev = Dispatch ( (lstdev, (ListType, TupleType)),
                    (astdev, (N.ArrayType,)) )
 sterr = Dispatch ( (lsterr, (ListType, TupleType)),
                    (asterr, (N.ArrayType,)) )
 sem = Dispatch ( (lsem, (ListType, TupleType)),
		  (asem, (N.ArrayType,)) )
 z = Dispatch ( (lz, (ListType, TupleType)),
                (az, (N.ArrayType,)) )
 zs = Dispatch ( (lzs, (ListType, TupleType)),
                 (azs, (N.ArrayType,)) )
 
## TRIMMING FCNS:
 threshold = Dispatch( (athreshold, (N.ArrayType,)),)
 trimboth = Dispatch ( (ltrimboth, (ListType, TupleType)),
                       (atrimboth, (N.ArrayType,)) )
 trim1 = Dispatch ( (ltrim1, (ListType, TupleType)),
                    (atrim1, (N.ArrayType,)) )
 
## CORRELATION FCNS:
 paired = Dispatch ( (lpaired, (ListType, TupleType)),
		     (apaired, (N.ArrayType,)) )
 pearsonr = Dispatch ( (lpearsonr, (ListType, TupleType)),
                       (apearsonr, (N.ArrayType,)) )
 spearmanr = Dispatch ( (lspearmanr, (ListType, TupleType)),
                        (aspearmanr, (N.ArrayType,)) )
 pointbiserialr = Dispatch ( (lpointbiserialr, (ListType, TupleType)),
                             (apointbiserialr, (N.ArrayType,)) )
 kendalltau = Dispatch ( (lkendalltau, (ListType, TupleType)),
                         (akendalltau, (N.ArrayType,)) )
 linregress = Dispatch ( (llinregress, (ListType, TupleType)),
                         (alinregress, (N.ArrayType,)) )
 
## INFERENTIAL STATS:
 ttest_1samp = Dispatch ( (lttest_1samp, (ListType, TupleType)),
                          (attest_1samp, (N.ArrayType,)) )
 ttest_ind = Dispatch ( (lttest_ind, (ListType, TupleType)),
                        (attest_ind, (N.ArrayType,)) )
 ttest_rel = Dispatch ( (lttest_rel, (ListType, TupleType)),
                        (attest_rel, (N.ArrayType,)) )
 chisquare = Dispatch ( (lchisquare, (ListType, TupleType)),
                        (achisquare, (N.ArrayType,)) )
 ks_2samp = Dispatch ( (lks_2samp, (ListType, TupleType)),
		       (aks_2samp, (N.ArrayType,)) )
 mannwhitneyu = Dispatch ( (lmannwhitneyu, (ListType, TupleType)),
                           (amannwhitneyu, (N.ArrayType,)) )
 tiecorrect = Dispatch ( (ltiecorrect, (ListType, TupleType)),
			 (atiecorrect, (N.ArrayType,)) )
 ranksums = Dispatch ( (lranksums, (ListType, TupleType)),
                       (aranksums, (N.ArrayType,)) )
 wilcoxont = Dispatch ( (lwilcoxont, (ListType, TupleType)),
                        (awilcoxont, (N.ArrayType,)) )
 kruskalwallish = Dispatch ( (lkruskalwallish, (ListType, TupleType)),
                             (akruskalwallish, (N.ArrayType,)) )
 friedmanchisquare = Dispatch ( (lfriedmanchisquare, (ListType, TupleType)),
                                (afriedmanchisquare, (N.ArrayType,)) )
 
## PROBABILITY CALCS:
 chisqprob = Dispatch ( (lchisqprob, (IntType, FloatType)),
                        (achisqprob, (N.ArrayType,)) )
 zprob = Dispatch ( (lzprob, (IntType, FloatType)),
                    (azprob, (N.ArrayType,)) )
 ksprob = Dispatch ( (lksprob, (IntType, FloatType)),
		     (aksprob, (N.ArrayType,)) )
 fprob = Dispatch ( (lfprob, (IntType, FloatType)),
                    (afprob, (N.ArrayType,)) )
 betacf = Dispatch ( (lbetacf, (IntType, FloatType)),
                     (abetacf, (N.ArrayType,)) )
 betai = Dispatch ( (lbetai, (IntType, FloatType)),
                    (abetai, (N.ArrayType,)) )
 erfcc = Dispatch ( (lerfcc, (IntType, FloatType)),
                    (aerfcc, (N.ArrayType,)) )
 gammln = Dispatch ( (lgammln, (IntType, FloatType)),
                     (agammln, (N.ArrayType,)) )
 
## ANOVA FUNCTIONS:
 anova = Dispatch ( (aanova, (ListType, TupleType, N.ArrayType)),)
 F_oneway = Dispatch ( (lF_oneway, (ListType, TupleType)),
                       (aF_oneway, (N.ArrayType,)) )
 F_value = Dispatch ( (lF_value, (ListType, TupleType)),
                      (aF_value, (N.ArrayType,)) )

## SUPPORT FUNCTIONS:
 incr = Dispatch ( (lincr, (ListType, TupleType, N.ArrayType)), )
 sum = Dispatch ( (lsum, (ListType, TupleType)),
                  (asum, (N.ArrayType,)) )
 cumsum = Dispatch ( (lcumsum, (ListType, TupleType)),
		     (acumsum, (N.ArrayType,)) )
 ss = Dispatch ( (lss, (ListType, TupleType)),
                 (ass, (N.ArrayType,)) )
 summult = Dispatch ( (lsummult, (ListType, TupleType)),
		      (asummult, (N.ArrayType,)) )
 square_of_sums = Dispatch ( (lsquare_of_sums, (ListType, TupleType)),
			     (asquare_of_sums, (N.ArrayType,)) )
 sumdiffsquared = Dispatch ( (lsumdiffsquared, (ListType, TupleType)),
			     (asumdiffsquared, (N.ArrayType,)) )
 shellsort = Dispatch ( (lshellsort, (ListType, TupleType)),
                        (ashellsort, (N.ArrayType,)) )
 rankdata = Dispatch ( (lrankdata, (ListType, TupleType)),
                       (arankdata, (N.ArrayType,)) )
 findwithin = Dispatch ( (lfindwithin, (ListType, TupleType)),
                         (afindwithin, (N.ArrayType,)) )

######################  END OF NUMERIC FUNCTION BLOCK  #####################

######################  END OF STATISTICAL FUNCTIONS  ######################

except ImportError:
 pass
