# pylint: disable-msg=E1101
"""
Wrapper to lowess, loess and stl routines, with support for masked arrays.

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
narray = numeric.array
nempty = numeric.empty
nlogical_not = numpy.logical_not

from scipy.sandbox.maskedarray.core import masked, nomask, mask_or, \
     masked_array as marray

import _lowess, _stl, _mloess


#####---------------------------------------------------------------------------
#--- --- STL ---
#####---------------------------------------------------------------------------
class lowess:
    """An object for robust locally weighted regression.

:IVariables:
    inputs :
    parameters :
    outputs :


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
            x = marray(x, copy=False, subok=True, dtype=float_, order='F').ravel()
            y = marray(y, copy=False, subok=True, dtype=float_, order='F').ravel()
            if x.size != y.size:
                msg = "Incompatible size between observations (%s) and response (%s)!"
                raise ValueError(msg % (x.size, y.size))
            idx = x.argsort()
            self._x = x[idx]
            self._y = y[idx]
            self._mask = mask_or(self._x._mask, self._y._mask,
                                 copy=False, small_mask=False)
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
            self._fval = marray(nempty((n,), dtype=float_, order='F'))
            self._rw = marray(nempty((n,), dtype=float_, order='F'))
            self._fres = marray(nempty((n,), dtype=float_, order='F'))
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
        self.outputs = self._outputs(self.inputs._x.size)
        # Force a fit .................
        self.fit()

    #............................................
    def fit(self):
        # Check the mask .........
        mask = self.inputs._mask
        if mask.any():
            unmask = nlogical_not(mask)
            (x, y) = (self.inputs._x[unmask], self.inputs._y[unmask])
        else:
            unmask = slice(None,None)
            (x, y) = (self.inputs._x, self.inputs._y)
        # Get the parameters .....
        self.parameters.activated = True
        f = self.parameters._span
        nsteps = self.parameters._nsteps
        delta = self.parameters._delta
        (tmp_s, tmp_w, tmp_r) = _lowess.lowess(x, y, f, nsteps, delta)
        # Process the outputs .....
        #... set the values
        self.outputs._fval[unmask] = tmp_s[:]
        self.outputs._rw[unmask] = tmp_w[:]
        self.outputs._fres[unmask] = tmp_r[:]
        #... set the masks
        self.outputs._fval._set_mask(mask)
        self.outputs._rw._set_mask(mask)
        self.outputs._fres._set_mask(mask)
        # Clean up the mess .......
        del(tmp_s, tmp_w, tmp_r)
        return self.outputs

#####---------------------------------------------------------------------------
#--- --- STL ---
#####---------------------------------------------------------------------------
class stl:
    class _inputs:
        def __init__(self, y):
            self.y = marray(y, subok=True, copy=False).ravel()
            self._mask = self.y._mask
            if self._mask.any():
                raise ValueError("Masked arrays should be filled first!")
            self.y_eff = self.y.compressed()
    #............................................
    class _model(object):
        """Model parameters of the STL fit.

:IVariables:
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
        """
        def __init__(self,
                     np=12, ns=7, nt=None, nl=13,
                     isdeg=1, itdeg=1, ildeg=1, caller=None):
            self._np = np
            self._ns = ns
            #
            self._nt = nt
            if nt is None:
                self._nt = max(int((1.5*np/(1.-1.5/ns))+0.5), 3)
            else:
                self._nt = nt
            if not self._nt % 2:
                self._nt += 1
            #
            if nl is None:
                self._nl = max(np, 3)
            else:
                self._nl = nl
            if not self._nl % 2:
                self._nl += 1
            #
            self._isdeg = isdeg
            self._itdeg = itdeg
            self._ildeg = ildeg
            self.activated = False
            self.caller = caller
        #.....
        def _get_np(self):
            "Gets the current seasonal period."
            return self._np
        def _set_np(self, np):
            "Sets the current seasonal period."
            self._np = max(np,2)
            if self.activated:
                self.caller.fit()
        np = property(fget=_get_np, fset=_set_np)
        #.....
        def _get_ns(self):
            "Gets the length of the seasonal smoother."
            return self._ns
        def _set_ns(self, ns):
            "Sets the length of the seasonal smoother."
            self._ns = max(ns, 3)
            if self._ns %2 == 0:
                self._ns += 1
            if self.activated:
                self.caller.fit()
        ns = property(fget=_get_ns, fset=_set_ns)
        #.....
        def _get_nt(self):
            "Gets the length of the trend smoother."
            return self._nt
        def _set_nt(self, nt):
            "Sets the length of the trend smoother."
            self._nt = nt
            if self.activated:
                self.caller.fit()
        nt = property(fget=_get_nt, fset=_set_nt)
        #.....
        def _get_nl(self):
            "Gets the length of the trend smoother."
            return self._nl
        def _set_nl(self, nl):
            "Sets the length of the trend smoother."
            self._nl = nl
            if self.activated:
                self.caller.fit()
        nl = property(fget=_get_nl, fset=_set_nl)
        #.....
        def _get_isdeg(self):
            "Gets the degree of the seasonal smoother."
            return self._isdeg
        def _set_isdeg(self, isdeg):
            "Sets the degree of the seasonal smoother."
            if isdeg > 2 or isdeg < 0:
                raise ValueError("The degree of the seasonal smoother should be 1 or 0.")
            self._isdeg = int(isdeg)
            if self.activated:
                self.caller.fit()
        isdeg = property(fget=_get_isdeg, fset=_set_isdeg)
        #.....
        def _get_itdeg(self):
            "Gets the degree of the trend smoother."
            return self._itdeg
        def _set_itdeg(self, itdeg):
            "Sets the degree of the trend smoother."
            if itdeg > 2 or itdeg < 0:
                raise ValueError("The degree of the trend smoother should be 1 or 0.")
            self._itdeg = int(itdeg)
            if self.activated:
                self.caller.fit()
        itdeg = property(fget=_get_itdeg, fset=_set_itdeg)
        #.....
        def _get_ildeg(self):
            "Gets the degree of the low-pass smoother."
            return self._ildeg
        def _set_ildeg(self, ildeg):
            "Sets the degree of the low-pass smoother."
            if ildeg > 2 or ildeg < 0:
                raise ValueError("The degree of the low-pass smoother should be 1 or 0.")
            self._ildeg = int(ildeg)
            if self.activated:
                self.caller.fit()
        ildeg = property(fget=_get_ildeg, fset=_set_ildeg)

    #............................................
    class _control(object):
        """Control parameters of the STL fit.

:IVariables:
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
        """
        def __init__(self,
                     nsjump=None,ntjump=None,nljump=None,
                     robust=True, ni=None,no=None, caller=None):
            (self._nsjump, self._ntjump, self._nljump) = (nsjump, ntjump, nljump)
            #...
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
            (self._robust, self._ni, self._no) = (robust, ni, no)
            #...
            self.activated = False
            self.caller = caller
        #....
        def _get_nsjump(self):
            "Gets the skipping value for seasonal smoothing."
            return self._nsjump
        def _set_nsjump(self, nsjump):
            "Sets the skipping value for seasonal smoothing."
            self._nsjump = nsjump
            if self.activated:
                self.caller.fit()
        nsjump = property(fget=_get_nsjump, fset=_set_nsjump)
        #....
        def _get_ntjump(self):
            "Gets the skipping value for trend smoothing."
            return self._ntjump
        def _set_ntjump(self, ntjump):
            "Sets the skipping value for trend smoothing."
            self._ntjump = ntjump
            if self.activated:
                self.caller.fit()
        ntjump = property(fget=_get_ntjump, fset=_set_ntjump)
        #....
        def _get_nljump(self):
            "Gets the skipping value for low-pass smoothing."
            return self._nljump
        def _set_nljump(self, nljump):
            "Set  the skipping value for low-pass smoothings"
            self._nljump = nljump
            if self.activated:
                self.caller.fit()
        nljump = property(fget=_get_nljump, fset=_set_nljump)
        #....
        def _get_robust(self):
            "Gets whether robust fitting should be performed."
            return self._robust
        def _set_robust(self, robust):
            "Sets whether robust fitting should be performed."
            self._robust = robust
            if self.activated:
                self.caller.fit()
        robust = property(fget=_get_robust, fset=_set_robust)
        #....
        def _get_ni(self):
            "Gets the number of loops."
            return self._ni
        def _set_ni(self, ni):
            "Sets the number of loops."
            if ni < 0:
                raise ValueError("The number of loops should be positive!")
            self._ni = ni
            if self.activated:
                self.caller.fit()
        ni = property(fget=_get_ni, fset=_set_ni)
        #....
        def _get_no(self):
            "Gets the number of iterations for robust fitting."
            return self._no
        def _set_no(self, no):
            "Sets the number of iterations for robust fitting."
            if no < 0 :
                raise ValueError("The number of iterations should be positive!")
            self._no = no
            if self.activated:
                self.caller.fit()
        no = property(fget=_get_no, fset=_set_no)
    #............................................
    class _outputs(object):
        """Outputs of the STL fit.

:IVariables:
    seasonal : ndarray
        Seasonal fitted values.
    trend : ndarray
        Trend fitted values.
    residuals : ndarray
        Fitted residuals.
    weights : ndarray
        Local robust weights. The final local robust weights are all 1 if no=0.
        """
        def __init__(self, n):
            self._seasonal = marray(nempty((n,), float_))
            self._trend = marray(nempty((n,), float_))
            self._weights = marray(nempty((n,), float_))
            self._residuals = marray(nempty((n,), float_))
        #.....
        seasonal = property(fget=lambda self:self._seasonal)
        trend = property(fget=lambda self:self._trend)
        weights = property(fget=lambda self:self._weights)
        residuals = property(fget=lambda self:self._residuals)
    #.............................................
    def __init__(self, y, **options):
        """Decomposes a time series into seasonal and trend  components.

:Parameters:
    y : ndarray
        Time series to be decomposed.
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
        """
        self.inputs = stl._inputs(y)
        self.model = stl._model(**dict(np=options.get('np',12),
                                       ns=options.get('ns',7),
                                       nt=options.get('nt',None),
                                       nl=options.get('nl',13),
                                       isdeg=options.get('isdeg',1),
                                       itdeg=options.get('itdeg',1),
                                       ildeg=options.get('ildeg',1),
                                       caller=self))
        optcontrol = dict(
                     nsjump=options.get('nsjump',int(0.1*self.model.ns+0.9)),
                     ntjump=options.get('ntjump',int(0.1*self.model.nt+0.9)),
                     nljump=options.get('nljump',int(0.1*self.model.nl+0.9)),
                     robust=options.get('robust',True),
                     ni=options.get('ni',None),
                     no=options.get('no',None),)
        self.control = stl._control(**optcontrol)
        self.outputs = stl._outputs(len(self.inputs.y))
        # Force a fit .................
        self.fit()

    #............................................
    def fit(self):
        # Get the input ...............
        y = self.inputs.y_eff
        mask = self.inputs._mask
        if mask is nomask:
            unmask = slice(None,None)
        else:
            unmask = nlogical_not(mask)
        # Get the parameters ..........
        model = self.model
        (np, ns, nt, nl) = (model.np, model.ns, model.nt, model.nl)
        (isdeg, itdeg, ildeg) = (model.isdeg, model.itdeg, model.ildeg)
        control = self.control
        (nsjump, ntjump, nljump) = (control.nsjump, control.ntjump, control.nljump)
        (ni, no) = (control.ni, control.no)
        # Compute the fit .............
        (rw,szn,trn,work) = _stl.stl(y,np,ns,nt,nl,isdeg,itdeg,ildeg,
                                     nsjump,ntjump,nljump,ni,no,)
        # Process the outputs .....
        #... set the values
        self.outputs.trend[unmask] = trn.flat
        self.outputs.seasonal[unmask] = szn.flat
        self.outputs.weights[unmask] = rw.flat
        self.outputs.residuals[unmask] = (y-trn-szn)
        #... set the masks
        self.outputs.trend._set_mask(mask)
        self.outputs.seasonal._set_mask(mask)
        self.outputs.weights._set_mask(mask)
        self.outputs.residuals._set_mask(mask)
        # Clean up the mess .......
        self.model.activated = self.control.activated = True
        del(trn, rw, szn)
        return self.outputs





def fstl(y, np=12, ns=7, nt=None, nl=13, isdeg=1, itdeg=1, ildeg=1,
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
    if ns % 2 == 0:
        ns += 1
    np = max(2, np)
    if nt is None:
        nt = max(int((1.5*np/(1.-1.5/ns))+0.5), 3)
        if not nt % 2:
            nt += 1
    if nl is None:
        nl = max(3,np)
        if not nl % 2:
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
loess = _mloess.loess
loess_anova = _mloess.anova

################################################################################
if __name__ == '__main__':
    from maskedarray.testutils import assert_almost_equal
    from maskedarray import masked_values
    from numpy import fromiter
    import os

    if 1:
        NOx = marray([4.818, 2.849, 3.275, 4.691, 4.255, 5.064, 2.118, 4.602,
                      2.286, 0.970, 3.965, 5.344, 3.834, 1.990, 5.199, 5.283,
                      -9999, -9999, 3.752, 0.537, 1.640, 5.055, 4.937, 1.561])
        NOx = maskedarray.masked_values(NOx, -9999)
        E = marray([0.831, 1.045, 1.021, 0.970, 0.825, 0.891, 0.71, 0.801,
                    1.074, 1.148, 1.000, 0.928, 0.767, 0.701, 0.807, 0.902,
                    -9999, -9999, 0.997, 1.224, 1.089, 0.973, 0.980, 0.665])
        gas_fit_E = numpy.array([0.665, 0.949, 1.224])
        newdata = numpy.array([0.6650000, 0.7581667, 0.8513333, 0.9445000,
                               1.0376667, 1.1308333, 1.2240000])
        coverage = 0.99

        rfile = open(os.path.join('tests','gas_result'), 'r')
        results = []
        for i in range(8):
            rfile.readline()
            z = fromiter((float(v) for v in rfile.readline().rstrip().split()),
                         float_)
            results.append(z)
        #
        gas = loess(E,NOx)
        gas.model.span = 2./3.
        gas.fit()
        assert_almost_equal(gas.outputs.fitted_values.compressed(), results[0], 6)
        assert_almost_equal(gas.outputs.enp, 5.5, 1)
        assert_almost_equal(gas.outputs.s, 0.3404, 4)
