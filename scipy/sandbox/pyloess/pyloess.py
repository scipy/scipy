# pylint: disable-msg=E1101
"""
Wrapper to lowess and stl routines.

LOWESS:
Initial Fortran code available at:
http://netlib.bell-labs.com/netlib/go/lowess.f.gz

initial author: W. S. Cleveland, 1979.
Simple to double precision conversion of the Fortran code by Pierre
Gerard-Marchant, 2007/03.

STL:
Initial Fortran code available at:
http://netlib.bell-labs.com/netlib/a/stl.gz
Initial Authors: R. B. Cleveland, W. S. Cleveland, J. E. McRae, and
I. Terpenning, 1990.
Simple-to-double precision conversion of the Fortran code by Pierre
Gerard-Marchant, 2007/03.

LOESS:
Initial C/Fortran package avialable at
http://netlib.bell-labs.com/netlib/a/dloess.gz
Initial authors: W. S. Cleveland, E. Grosse and Shyu
Adaptation to Pyrex/Python by Pierre Gerard-Marchant, 2007/03

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date$
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy
from numpy import bool_, complex_, float_, int_, str_, object_
import numpy.core.numeric as numeric
from numpy.core.records import recarray

from numpy.core import array as narray
from numpy.core import empty as nempty

import _lowess, _stl, _loess


#####---------------------------------------------------------------------------
#--- --- STL ---
#####---------------------------------------------------------------------------
def flowess(x,y,span=0.5,nsteps=2,delta=0):
    """Performs a robust locally weighted regression (lowess).

    Outputs a *3xN* array of fitted values, residuals and fit weights.


:Parameters:
    x : ndarray
        Abscissas of the points on the scatterplot; the values in X must be
        ordered from smallest to largest.
    y : ndarray
        Ordinates of the points on the scatterplot.
    span : Float *[0.5]*
        Fraction of the total number of points used to compute each fitted value.
        As f increases the smoothed values become smoother. Choosing f in the range
        .2 to .8 usually results in a good fit.
    nsteps : Integer *[2]*
        Number of iterations in the robust fit. If nsteps=0, the nonrobust fit
        is returned; setting nsteps=2 should serve most purposes.
    delta : Integer *[0]*
        Nonnegative parameter which may be used to save computations.
        If N (the number of elements in x) is less than 100, set delta=0.0;
        if N is greater than 100 you should find out how delta works by reading
        the additional instructions section.

:Returns:
    A recarray of smoothed values ('smooth'), residuals ('residuals') and local
    robust weights ('weights').


Additional instructions
-----------------------

Fro the original author:

        DELTA can be used to save computations.   Very  roughly  the
        algorithm  is  this:   on the initial fit and on each of the
        NSTEPS iterations locally weighted regression fitted  values
        are computed at points in X which are spaced, roughly, DELTA
        apart; then the fitted values at the  remaining  points  are
        computed  using  linear  interpolation.   The  first locally
        weighted regression (l.w.r.) computation is carried  out  at
        X(1)  and  the  last  is  carried  out at X(N).  Suppose the
        l.w.r. computation is carried out at  X(I).   If  X(I+1)  is
        greater  than  or  equal  to  X(I)+DELTA,  the  next  l.w.r.
        computation is carried out at X(I+1).   If  X(I+1)  is  less
        than X(I)+DELTA, the next l.w.r.  computation is carried out
        at the largest X(J) which is greater than or equal  to  X(I)
        but  is not greater than X(I)+DELTA.  Then the fitted values
        for X(K) between X(I)  and  X(J),  if  there  are  any,  are
        computed  by  linear  interpolation  of the fitted values at
        X(I) and X(J).  If N is less than 100 then DELTA can be  set
        to  0.0  since  the  computation time will not be too great.
        For larger N it is typically not necessary to carry out  the
        l.w.r.  computation for all points, so that much computation
        time can be saved by taking DELTA to be  greater  than  0.0.
        If  DELTA =  Range  (X)/k  then,  if  the  values  in X were
        uniformly  scattered  over  the  range,  the   full   l.w.r.
        computation  would be carried out at approximately k points.
        Taking k to be 50 often works well.

Method
------

        The fitted values are computed by using the nearest neighbor
        routine  and  robust locally weighted regression of degree 1
        with the tricube weight function.  A few additional features
        have  been  added.  Suppose r is FN truncated to an integer.
        Let  h  be  the  distance  to  the  r-th  nearest   neighbor
        from X[i].   All  points within h of X[i] are used.  Thus if
        the r-th nearest neighbor is exactly the  same  distance  as
        other  points,  more  than r points can possibly be used for
        the smooth at  X[i].   There  are  two  cases  where  robust
        locally  weighted regression of degree 0 is actually used at
        X[i].  One case occurs when  h  is  0.0.   The  second  case
        occurs  when  the  weighted  standard error of the X[i] with
        respect to the weights w[j] is  less  than  .001  times  the
        range  of the X[i], where w[j] is the weight assigned to the
        j-th point of X (the tricube  weight  times  the  robustness
        weight)  divided by the sum of all of the weights.  Finally,
        if the w[j] are all zero for the smooth at X[i], the  fitted
        value is taken to be Y[i].

References
----------
    W. S. Cleveland. 1978. Visual and Computational Considerations in
    Smoothing Scatterplots by Locally Weighted Regression. In
    Computer Science and Statistics: Eleventh Annual Symposium on the
    Interface, pages 96-100. Institute of Statistics, North Carolina
    State University, Raleigh, North Carolina, 1978.

    W. S. Cleveland, 1979. Robust Locally Weighted Regression and
    Smoothing Scatterplots. Journal of the American Statistical
    Association, 74:829-836, 1979.

    W. S. Cleveland, 1981. LOWESS: A Program for Smoothing Scatterplots
    by Robust Locally Weighted Regression. The American Statistician,
    35:54.

    """
    x = narray(x, copy=False, subok=True, dtype=float_)
    y = narray(y, copy=False, subok=True, dtype=float_)
    if x.size != y.size:
        raise ValueError("Incompatible size between observations and response!")


    out_dtype = [('smooth',float_), ('weigths', float_), ('residuals', float_)]
    return numeric.fromiter(zip(*_lowess.lowess(x,y,span,nsteps,delta,)),
                            dtype=out_dtype).view(recarray)


class lowess:
    """An object for robust locally weighted regression.

:IVariables:
    inputs : An object storing the inputs.
        x : A (n,) ndarray of observations (sorted by increasing values).
        y : A (n,) ndarray of responses (sorted by increasing x).
    parameters : An object storing the control parameters.
        span : Fraction of the total number of points used in the smooth.
        nsteps : Number of iterations of the robust fit.
        delta : Parameter used to save computation time
    outputs : An object storing the outputs.
        smooth : A (n,) ndarray of fitted values.
        residuals : A (n,) ndarray of fitted residuals.
        weights : A (n,) ndarray of robust weights.


Method
------
    The fitted values are computed by using the nearest neighbor
    routine  and  robust locally weighted regression of degree 1
    with the tricube weight function.  A few additional features
    have  been  added.  Suppose r is FN truncated to an integer.
    Let  h  be  the  distance  to  the  r-th  nearest   neighbor
    from X[i].   All  points within h of X[i] are used.  Thus if
    the r-th nearest neighbor is exactly the  same  distance  as
    other  points,  more  than r points can possibly be used for
    the smooth at  X[i].   There  are  two  cases  where  robust
    locally  weighted regression of degree 0 is actually used at
    X[i].  One case occurs when  h  is  0.0.   The  second  case
    occurs  when  the  weighted  standard error of the X[i] with
    respect to the weights w[j] is  less  than  .001  times  the
    range  of the X[i], where w[j] is the weight assigned to the
    j-th point of X (the tricube  weight  times  the  robustness
    weight)  divided by the sum of all of the weights.  Finally,
    if the w[j] are all zero for the smooth at X[i], the  fitted
    value is taken to be Y[i].

References
----------
    W. S. Cleveland. 1978. Visual and Computational Considerations in
    Smoothing Scatterplots by Locally Weighted Regression. In
    Computer Science and Statistics: Eleventh Annual Symposium on the
    Interface, pages 96-100. Institute of Statistics, North Carolina
    State University, Raleigh, North Carolina, 1978.

    W. S. Cleveland, 1979. Robust Locally Weighted Regression and
    Smoothing Scatterplots. Journal of the American Statistical
    Association, 74:829-836, 1979.

    W. S. Cleveland, 1981. LOWESS: A Program for Smoothing Scatterplots
    by Robust Locally Weighted Regression. The American Statistician,
    35:54.
    """
    #............................................
    class _inputs(object):
        """Inputs of the lowess fit.

:IVariables:
    x : ndarray
        A (n,) float ndarray of observations (sorted by increasing values).
    y : ndarray
        A (n,) float ndarray of responses (sorted by increasing x).
        """
        def __init__(self, x, y):
            x = narray(x, copy=False, subok=True, dtype=float_).ravel()
            y = narray(y, copy=False, subok=True, dtype=float_).ravel()
            if x.size != y.size:
                msg = "Incompatible size between observations (%s) and response (%s)!"
                raise ValueError(msg % (x.size, y.size))
            idx = x.argsort()
            self._x = x[idx]
            self._y = y[idx]
        #.....
        x = property(fget=lambda self:self._x)
        y = property(fget=lambda self:self._y)
    #............................................
    class _parameters(object):
        """Parameters of the lowess fit.

:IVariables:
    span : float *[0.5]*
        Fraction of the total number of points used to compute each fitted value.
        As f increases the smoothed values become smoother. Choosing f in the range
        .2 to .8 usually results in a good fit.
    nsteps : integer *[2]*
        Number of iterations in the robust fit. If nsteps=0, the nonrobust fit
        is returned; setting nsteps=2 should serve most purposes.
    delta : integer *[0]*
        Nonnegative parameter which may be used to save computations.
        If N (the number of observations) is less than 100, set delta=0.0;
        if N is greater than 100 you should find out how delta works by reading
        the additional instructions section.
        """
        def __init__(self, span, nsteps, delta, caller):
            self.activated = False
            self._span = span
            self._nsteps = nsteps
            self._delta = delta
            self._caller = caller
        #.....
        def _get_span(self):
            "Gets the current span."
            return self._span
        def _set_span(self, span):
            "Sets the current span, and refit if needed."
            if span <= 0 or span > 1:
                raise ValueError("span should be between zero and one!")
            self._span = span
            if self.activated:
                self._caller.fit()
        span = property(fget=_get_span, fset=_set_span)
        #.....
        def _get_nsteps(self):
            "Gets the current number of iterations."
            return self._nsteps
        def _set_nsteps(self, nsteps):
            "Sets the current number of iterations, and refit if needed."
            if nsteps < 0:
                raise ValueError("nsteps should be positive!")
            self._nsteps = nsteps
            if self.activated:
                self._caller.fit()
        nsteps = property(fget=_get_nsteps, fset=_set_nsteps)
        #.....
        def _get_delta(self):
            "Gets the current delta."
            return self._delta
        def _set_delta(self, delta):
            "Sets the current delta, and refit if needed."
            if delta < 0:
                raise ValueError("delta should be positive!")
            self._delta = delta
            if self.activated:
                self._caller.fit()
        delta = property(fget=_get_delta, fset=_set_delta)
    #............................................
    class _outputs(object):
        """Outputs of the lowess fit.

:IVariables:
    fitted_values : ndarray
        A (n,) ndarray of fitted values (readonly).
    fitted_residuals : ndarray
        A (n,) ndarray of residuals (readonly).
    weights : ndarray
        A (n,) ndarray of robust weights (readonly).
        """
        def __init__(self, n):
            self._fval = nempty((n,), float_)
            self._rw = nempty((n,), float_)
            self._fres = nempty((n,), float_)
        #.....
        fitted_values = property(fget=lambda self:self._fval)
        robust_weights = property(fget=lambda self:self._rw)
        fitted_residuals = property(fget=lambda self:self._fres)

    #............................................
    def __init__(self, x, y, span=0.5, nsteps=2, delta=0):
        """
:Parameters:
    x : ndarray
        Abscissas of the points on the scatterplot; the values in X must be
        ordered from smallest to largest.
    y : ndarray
        Ordinates of the points on the scatterplot.
    span : Float *[0.5]*
        Fraction of the total number of points used to compute each fitted value.
        As span increases the smoothed values become smoother. Choosing span in
        the range .2 to .8 usually results in a good fit.
    nsteps : Integer *[2]*
        Number of iterations in the robust fit. If nsteps=0, the nonrobust fit
        is returned; setting nsteps=2 should serve most purposes.
    delta : Integer *[0]*
        Nonnegative parameter which may be used to save computations.
        If N (the number of elements in x) is less than 100, set delta=0.0;
        if N is greater than 100 you should find out how delta works by reading
        the additional instructions section.
        """
        # Chek the input data .........
        # Initialize the attributes ...
        self.inputs = lowess._inputs(x,y)
        self.parameters = lowess._parameters(span, nsteps, delta, self)
        self.outputs = lowess._outputs(self.inputs._x.size)
        # Force a fit .................
        self.fit()

    #............................................
    def fit(self):
        """Computes the lowess fit. Returns a lowess.outputs object."""
        (x, y) = (self.inputs._x, self.inputs._y)
        # Get the parameters .....
        self.parameters.activated = True
        f = self.parameters._span
        nsteps = self.parameters._nsteps
        delta = self.parameters._delta
        (tmp_s, tmp_w, tmp_r) = _lowess.lowess(x, y, f, nsteps, delta)
        # Process the outputs .....
        #... set the values
        self.outputs.fitted_values[:] = tmp_s.flat
        self.outputs.robust_weights[:] = tmp_w.flat
        self.outputs.fitted_residuals[:] = tmp_r.flat
        # Clean up the mess .......
        del(tmp_s, tmp_w, tmp_r)
        return self.outputs


#####---------------------------------------------------------------------------
#--- --- STL ---
#####---------------------------------------------------------------------------
def stl(y, np=12, ns=7, nt=None, nl=13, isdeg=1, itdeg=1, ildeg=1,
        nsjump=None,ntjump=None,nljump=None, robust=True, ni=None,no=None):
    """Decomposes a time series into seasonal and trend  components.

:Parameters:
    y : Numerical array
        Time Series to be decomposed.
    np : Integer *[12]*
        Period of the seasonal component.
        For example, if  the  time series is monthly with a yearly cycle, then
        np=12.
    ns : Integer *[7]*
        Length of the seasonal smoother.
        The value of  ns should be an odd integer greater than or equal to 3.
        A value ns>6 is recommended. As ns  increases  the  values  of  the
        seasonal component at a given point in the seasonal cycle (e.g., January
        values of a monthly series with  a  yearly cycle) become smoother.
    nt : Integer *[None]*
        Length of the trend smoother.
        The  value  of  nt should be an odd integer greater than or equal to 3.
        A value of nt between 1.5*np and 2*np is  recommended. As nt increases,
        the values of the trend component become  smoother.
        If nt is None, it is estimated as the smallest odd integer greater
        or equal to (1.5*np)/[1-(1.5/ns)]
    nl : Integer *[None]*
        Length of the low-pass filter.
        The value of nl should  be an odd integer greater than or equal to 3.
        The smallest odd integer greater than or equal to np is used by default.
    isdeg : Integer *[1]*
        Degree of locally-fitted polynomial in seasonal smoothing.
        The value is 0 or 1.
    itdeg : Integer *[1]*
        Degree of locally-fitted polynomial in trend smoothing.
        The value is 0 or 1.
    ildeg : Integer *[1]*
        Degree of locally-fitted polynomial in low-pass smoothing.
        The value is 0 or 1.
    nsjump : Integer *[None]*
        Skipping value for seasonal smoothing.
        The seasonal smoother skips ahead nsjump points and then linearly
        interpolates in between.  The value  of nsjump should be a positive
        integer; if nsjump=1, a seasonal smooth is calculated at all n points.
        To make the procedure run faster, a reasonable choice for nsjump is
        10%-20% of ns. By default, nsjump= 0.1*ns.
    ntjump : Integer *[1]*
        Skipping value for trend smoothing. If None, ntjump= 0.1*nt
    nljump : Integer *[1]*
        Skipping value for low-pass smoothing. If None, nljump= 0.1*nl
    robust : Boolean *[True]*
        Flag indicating whether robust fitting should be performed.
    ni : Integer *[None]*
        Number of loops for updating the seasonal and trend  components.
        The value of ni should be a positive integer.
        See the next argument for advice on the  choice of ni.
        If ni is None, ni is set to 1 for robust fitting, to 5 otherwise.
    no : Integer *[0]*
        Number of iterations of robust fitting. The value of no should
        be a nonnegative integer. If the data are well behaved without
        outliers, then robustness iterations are not needed. In this case
        set no=0, and set ni=2 to 5 depending on how much security
        you want that  the seasonal-trend looping converges.
        If outliers are present then no=3 is a very secure value unless
        the outliers are radical, in which case no=5 or even 10 might
        be better.  If no>0 then set ni to 1 or 2.
        If None, then no is set to 15 for robust fitting, to 0 otherwise.

Returns:
    A recarray of estimated trend values ('trend'), estimated seasonal
    components ('seasonal'), local robust weights ('weights') and fit
    residuals ('residuals').
    The final local robust weights are all 1 if no=0.

Reference
---------

    R. B. Cleveland, W. S. Cleveland, J. E. McRae and I. Terpenning.
    1990. STL: A Seasonal-Trend Decomposition Procedure Based on LOESS
    (with Discussion). Journal of Official Statistics, 6:3-73.


    """
    ns = max(ns, 3)
    if ns%2 == 0:
        ns += 1
    np = max(2, np)
    if nt is None:
        nt = max(int((1.5*np/(1.-1.5/ns))+0.5), 3)
        if not nt%2:
            nt += 1
    if nl is None:
        nl = max(3,np)
        if not nl%2:
            nl += 1
    if nsjump is None:
        nsjump = int(0.1*ns + 0.9)
    if ntjump is None:
        ntjump = int(0.1*nt + 0.9)
    if nljump is None:
        nljump = int(0.1*nl + 0.9)
    if robust:
        if ni is None:
            ni = 1
        if no is None:
            no = 15
    else:
        if ni is None:
            ni = 5
        if no is None:
            no = 0

    if hasattr(y,'_mask') and numpy.any(y._mask):
        raise ValueError,"Missing values should first be filled !"
    y = numeric.array(y, subok=True, copy=False).ravel()
    (rw,szn,trn,work) = _stl.stl(y,np,ns,nt,nl,isdeg,itdeg,ildeg,
                                 nsjump,ntjump,nljump,ni,no,)
    dtyp = [('trend', float_), ('seasonal', float_),
            ('residuals', float_), ('weights', float_)]
    result = numeric.fromiter(zip(trn,szn,y-trn-szn,rw), dtype=dtyp)
    return result.view(recarray)

#####---------------------------------------------------------------------------
#--- --- Loess ---
#####---------------------------------------------------------------------------
loess = _loess.loess
"""
loess : locally weighted estimates. Multi-variate version

:Keywords:
    x : ndarray
        A (n,p) ndarray of independent variables, with n the number of observations
        and p the number of variables.
    y : ndarray
        A (n,) ndarray of observations
    weights : ndarray
        A (n,) ndarray of weights to be given to individual observations in the
        sum of squared residuals that forms the local fitting criterion. If not
        None, the weights should be non negative. If the different observations
        have non-equal variances, the weights should be inversely proportional
        to the variances.
        By default, an unweighted fit is carried out (all the weights are one).
    surface : string ["interpolate"]
        Determines whether the fitted surface is computed directly at all points
        ("direct") or whether an interpolation method is used ("interpolate").
        The default ("interpolate") is what most users should use unless special
        circumstances warrant.
    statistics : string ["approximate"]
        Determines whether the statistical quantities are computed exactly
        ("exact") or approximately ("approximate"). "exact" should only be used
        for testing the approximation in statistical development and is not meant
        for routine usage because computation time can be horrendous.
    trace_hat : string ["wait.to.decide"]
        Determines how the trace of the hat matrix should be computed. The hat
        matrix is used in the computation of the statistical quantities.
        If "exact", an exact computation is done; this could be slow when the
        number of observations n becomes large. If "wait.to.decide" is selected,
        then a default is "exact" for n < 500 and "approximate" otherwise.
        This option is only useful when the fitted surface is interpolated. If
        surface is "exact", an exact computation is always done for the trace.
        Setting trace_hat to "approximate" for large dataset will substantially
        reduce the computation time.
    iterations : integer
        Number of iterations of the robust fitting method. If the family is
        "gaussian", the number of iterations is set to 0.
    cell : integer
        Maximum cell size of the kd-tree. Suppose k = floor(n*cell*span),
        where n is the number of observations, and span the smoothing parameter.
        Then, a cell is further divided if the number of observations within it
        is greater than or equal to k. This option is only used if the surface
        is interpolated.
    span : float [0.75]
        Smoothing factor, as a fraction of the number of points to take into
        account.
    degree : integer [2]
        Overall degree of locally-fitted polynomial. 1 is locally-linear
        fitting and 2 is locally-quadratic fitting.  Degree should be 2 at most.
    normalize : boolean [True]
        Determines whether the independent variables should be normalized.
        If True, the normalization is performed by setting the 10% trimmed
        standard deviation to one. If False, no normalization is carried out.
        This option is only useful for more than one variable. For spatial
        coordinates predictors or variables with a common scale, it should be
        set to False.
    family : string ["gaussian"]
        Determines the assumed distribution of the errors. The values are
        "gaussian" or "symmetric". If "gaussian" is selected, the fit is
        performed with least-squares. If "symmetric" is selected, the fit
        is performed robustly by redescending M-estimators.
    parametric_flags : sequence [ [False]*p ]
        Indicates which independent variables should be conditionally-parametric
       (if there are two or more independent variables). The argument should
       be a sequence of booleans, with the same size as the number of independent
       variables, specified in the order of the predictor group ordered in x.
    drop_square : sequence [ [False]* p]
        When there are two or more independent variables and when a 2nd order
        polynomial is used, "drop_square_flags" specifies those numeric predictors
        whose squares should be dropped from the set of fitting variables.
        The method of specification is the same as for parametric.

:Outputs:
    fitted_values : ndarray
        The (n,) ndarray of fitted values.
    fitted_residuals : ndarray
        The (n,) ndarray of fitted residuals (observations - fitted values).
    enp : float
        Equivalent number of parameters.
    s : float
        Estimate of the scale of residuals.
    one_delta: float
        Statistical parameter used in the computation of standard errors.
    two_delta : float
        Statistical parameter used in the computation of standard errors.
    pseudovalues : ndarray
        The (n,) ndarray of adjusted values of the response when robust estimation
        is used.
    trace_hat : float
        Trace of the operator hat matrix.
    diagonal :
        Diagonal of the operator hat matrix.
    robust : ndarray
        The (n,) ndarray of robustness weights for robust fitting.
    divisor : ndarray
        The (p,) array of normalization divisors for numeric predictors.


    newdata : ndarray
        The (m,p) array of independent variables where the surface must be estimated.
    values : ndarray
        The (m,) ndarray of loess values evaluated at newdata
    stderr : ndarray
        The (m,) ndarray of the estimates of the standard error on the estimated
        values.
    residual_scale : float
        Estimate of the scale of the residuals
    df : integer
        Degrees of freedom of the t-distribution used to compute pointwise
        confidence intervals for the evaluated surface.
    nest : integer
        Number of new observations.
"""

loess_anova = _loess.anova

################################################################################
if __name__ == '__main__':
    from maskedarray.testutils import assert_almost_equal
    from maskedarray import masked_values
