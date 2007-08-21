# -*- Mode: Python -*-  
cimport c_python
cimport c_numpy
from c_numpy cimport ndarray, npy_intp, \
    PyArray_SIZE, PyArray_EMPTY, PyArray_FROMANY, \
    NPY_INT, NPY_DOUBLE, NPY_OWNDATA, NPY_ALIGNED, NPY_FORTRAN, \
    PyArray_SimpleNewFromData
import numpy
narray = numpy.array
float_ = numpy.float_

# NumPy must be initialized
c_numpy.import_array()

cimport c_loess

cdef floatarray_from_data(int rows, int cols, double *data):
    cdef ndarray a_ndr
    cdef npy_intp size
    size = rows*cols
    a_ndr = <object>PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, data)
    if cols > 1:
        a_ndr.shape = (rows, cols)
    return a_ndr

cdef boolarray_from_data(int rows, int cols, int *data):
    cdef ndarray a_ndr
    cdef npy_intp size
    size = rows*cols
    a_ndr = <object>PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, data)
    if cols > 1:
        a_ndr.shape = (rows, cols)
    return a_ndr.astype(numpy.bool)


"""
        

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


#####---------------------------------------------------------------------------
#---- ---- loess model ---
#####---------------------------------------------------------------------------
cdef class loess_inputs:
    """Loess inputs
    
:IVariables:
    x : ndarray
        A (n,p) ndarray of independent variables, with n the number of observations
        and p the number of variables.
    y : ndarray
        A (n,) ndarray of observations
    nobs : integer
        Number of observations: nobs=n.
    npar : integer
        Number of independent variables: npar=p.
    weights : ndarray
        A (n,) ndarray of weights to be given to individual observations in the 
        sum of squared residuals that forms the local fitting criterion. If not
        None, the weights should be non negative. If the different observations
        have non-equal variances, the weights should be inversely proportional 
        to the variances.
        By default, an unweighted fit is carried out (all the weights are one).
    """
    cdef c_loess.c_loess_inputs *_base
    cdef ndarray w_ndr
    cdef readonly ndarray x, y, masked, x_eff, y_eff
    cdef readonly int nobs, npar
    #.........
    def __init__(self, x_data, y_data):
        cdef ndarray unmasked
        self.x = narray(x_data, copy=False, subok=True, dtype=float_, order='C')
        self.y = narray(y_data, copy=False, subok=True, dtype=float_, order='C')
        # Check the dimensions ........
        if self.x.ndim == 1:
            self.npar = 1
        elif self.x.ndim == 2:
            self.npar = self.x.shape[-1]
        else:
            raise ValueError("The array of indepedent varibales should be 2D at most!")
        self.nobs = len(self.x)         
        # Get the effective data ......
        self.x_eff = self.x.ravel()
        self.y_eff = self.y
        self.w_ndr = numpy.ones((self.nobs,), dtype=float_)

    #.........    
    property weights:
        """A (n,) ndarray of weights to be given to individual observations in the 
        sum of squared residuals that forms the local fitting criterion. If not
        None, the weights should be non negative. If the different observations
        have non-equal variances, the weights should be inversely proportional 
        to the variances.
        By default, an unweighted fit is carried out (all the weights are one).
        """
        def __get__(self):
            return self.w_ndr
        
        def __set__(self, w):
            cdef npy_intp *dims
            cdef ndarray w_ndr
            w_ndr = narray(w, copy=False, subok=False)
            if w_ndr.ndim > 1 or w_ndr.size != self.nobs:
                raise ValueError, "Invalid size of the 'weights' vector!"
            self.w_ndr = w_ndr
            self._base.weights = <double *>w_ndr.data
#       
######---------------------------------------------------------------------------
##---- ---- loess control ---
######---------------------------------------------------------------------------
cdef class loess_control:
    """Loess control parameters.
    
:IVariables:
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
    
    """
    cdef c_loess.c_loess_control *_base
    #.........    
    property surface:
        """
    surface : string ["interpolate"]
        Determines whether the fitted surface is computed directly at all points
        ("direct") or whether an interpolation method is used ("interpolate").
        The default ("interpolate") is what most users should use unless special 
        circumstances warrant.
        """
        def __get__(self):
            return self._base.surface
        def __set__(self, surface):
            if surface.lower() not in ('interpolate', 'direct'):
                raise ValueError("Invalid value for the 'surface' argument: "+
                                 "should be in ('interpolate', 'direct').")
            tmpx = surface.lower()
            self._base.surface = tmpx
    #.........
    property statistics:
        """
    statistics : string ["approximate"]
        Determines whether the statistical quantities are computed exactly 
        ("exact") or approximately ("approximate"). "exact" should only be used 
        for testing the approximation in statistical development and is not meant 
        for routine usage because computation time can be horrendous.
        """
        def __get__(self):
            return self._base.statistics
        def __set__(self, statistics):
            if statistics.lower() not in ('approximate', 'exact'):
                raise ValueError("Invalid value for the 'statistics' argument: "\
                                 "should be in ('approximate', 'exact').")
            tmpx = statistics.lower()
            self._base.statistics = tmpx
    #.........
    property trace_hat:
        """
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
        """
        def __get__(self):
            return self._base.trace_hat
        def __set__(self, trace_hat):
            if trace_hat.lower() not in ('approximate', 'exact'):
                raise ValueError("Invalid value for the 'trace_hat' argument: "\
                                 "should be in ('approximate', 'exact').")
            tmpx = trace_hat.lower()
            self._base.trace_hat = tmpx
    #.........
    property iterations:
        """
    iterations : integer
        Number of iterations of the robust fitting method. If the family is 
        "gaussian", the number of iterations is set to 0.
        """
        def __get__(self):
            return self._base.iterations
        def __set__(self, iterations):
            if iterations < 0:
                raise ValueError("Invalid number of iterations: should be positive")
            self._base.iterations = iterations
    #.........
    property cell:
        """
    cell : integer
        Maximum cell size of the kd-tree. Suppose k = floor(n*cell*span),
        where n is the number of observations, and span the smoothing parameter.
        Then, a cell is further divided if the number of observations within it 
        is greater than or equal to k. This option is only used if the surface 
        is interpolated.
        """     
        def __get__(self):
            return self._base.cell
        def __set__(self, cell):
            if cell <= 0:
                raise ValueError("Invalid value for the cell argument: should be positive")
            self._base.cell = cell
    #.........
    def update(self, **cellargs):
        """Updates several parameters at once."""
        surface = cellargs.get('surface', None)
        if surface is not None:
            self.surface = surface
        #
        statistics = cellargs.get('statistics', None)
        if statistics is not None:
            self.statistics = statistics
        #    
        trace_hat = cellargs.get('trace_hat', None)
        if trace_hat is not None:
            self.trace_hat = trace_hat
        #
        iterations = cellargs.get('iterations', None)
        if iterations is not None:
            self.iterations = iterations
        #
        cell = cellargs.get('cell', None)
        if cell is not None:
            self.parametric_flags = cell
        #
    #.........
    def __str__(self):
        strg = ["Control          :",
                "Surface type     : %s" % self.surface,
                "Statistics       : %s" % self.statistics,
                "Trace estimation : %s" % self.trace_hat,
                "Cell size        : %s" % self.cell,
                "Nb iterations    : %s" % self.iterations,]
        return '\n'.join(strg)
        
#    
######---------------------------------------------------------------------------
##---- ---- loess kd_tree ---
######---------------------------------------------------------------------------
cdef class loess_kd_tree:
    cdef c_loess.c_loess_kd_tree *_base

######---------------------------------------------------------------------------
##---- ---- loess model ---
######---------------------------------------------------------------------------
cdef class loess_model:
    """loess_model contains parameters required for a loess fit.
    
:IVariables:
    normalize : boolean [True]
        Determines whether the independent variables should be normalized.  
        If True, the normalization is performed by setting the 10% trimmed 
        standard deviation to one. If False, no normalization is carried out. 
        This option is only useful for more than one variable. For spatial
        coordinates predictors or variables with a common scale, it should be 
        set to False.
    span : float [0.75]
        Smoothing factor, as a fraction of the number of points to take into
        account. 
    degree : integer [2]
        Overall degree of locally-fitted polynomial. 1 is locally-linear 
        fitting and 2 is locally-quadratic fitting.  Degree should be 2 at most.
    family : string ["gaussian"]
        Determines the assumed distribution of the errors. The values are 
        "gaussian" or "symmetric". If "gaussian" is selected, the fit is 
        performed with least-squares. If "symmetric" is selected, the fit
        is performed robustly by redescending M-estimators.
    parametric_flags : sequence [ [False]*p ]
        Indicates which independent variables should be conditionally-parametric
        (if there are two or more independent variables). The argument should be
        a sequence of booleans, with the same size as the number of independent 
        variables, specified in the order of the predictor group ordered in x. 
        Note: elements of the sequence cannot be modified individually: the whole
        sequence must be given.
    drop_square : sequence [ [False]* p]
        When there are two or more independent variables and when a 2nd order
        polynomial is used, "drop_square_flags" specifies those numeric predictors 
        whose squares should be dropped from the set of fitting variables. 
        The method of specification is the same as for parametric.  
        Note: elements of the sequence cannot be modified individually: the whole
        sequence must be given.

    """
        
    cdef c_loess.c_loess_model *_base
    cdef long npar
    #.........
    cdef setup(self, c_loess.c_loess_model *base, long npar):
        self._base = base
        self.npar = npar
        return self    
    #.........
    property normalize:
        def __get__(self):
            return bool(self._base.normalize)
        def __set__(self, normalize):
            self._base.normalize = normalize
    #.........
    property span:
        def __get__(self):
            return self._base.span
        def __set__(self, span):
            if span <= 0. or span > 1.:
                raise ValueError("Span should be between 0 and 1!")
            self._base.span = span
    #.........
    property degree:
        def __get__(self):
            return self._base.degree
        def __set__(self, degree):
            if degree < 0 or degree > 2:
                raise ValueError("Degree should be between 0 and 2!")
            self._base.degree = degree
    #.........
    property family:
        def __get__(self):
            return self._base.family
        def __set__(self, family):
            if family.lower() not in ('symmetric', 'gaussian'):
                raise ValueError("Invalid value for the 'family' argument: "\
                                 "should be in ('symmetric', 'gaussian').")
            self._base.family = family
    #.........
    property parametric_flags:
        def __get__(self):
            return boolarray_from_data(self.npar, 1, self._base.parametric)
        def __set__(self, paramf):
            cdef ndarray p_ndr
            cdef int i
            p_ndr = numpy.atleast_1d(narray(paramf, copy=False, subok=True, 
                                            dtype=numpy.bool))
            for i from 0 <= i < min(self.npar, p_ndr.size):
                self._base.parametric[i] = p_ndr[i]
    #.........
    property drop_square_flags:
        def __get__(self):
            return boolarray_from_data(self.npar, 1, self._base.drop_square)
        def __set__(self, drop_sq):
            cdef ndarray d_ndr
            cdef int i
            d_ndr = numpy.atleast_1d(narray(drop_sq, copy=False, subok=True, 
                                            dtype=numpy.bool))
            for i from 0 <= i < min(self.npar, d_ndr.size):
                self._base.drop_square[i] = d_ndr[i]
    #........
    def update(self, **modelargs):
        family = modelargs.get('family', None)
        if family is not None:
            self.family = family
        #
        span = modelargs.get('span', None)
        if span is not None:
            self.span = span
        #    
        degree = modelargs.get('degree', None)
        if degree is not None:
            self.degree = degree
        #
        normalize = modelargs.get('normalize', None)
        if normalize is not None:
            self.normalize = normalize
        #
        parametric = modelargs.get('parametric', None)
        if parametric is not None:
            self.parametric_flags = parametric
        #
        drop_square = modelargs.get('drop_square', None)
        if drop_square is not None:
            self.drop_square_flags = drop_square
    #.........
    def __repr__(self):
        return "<loess object: model parameters>"
    #.........
    def __str__(self):
        strg = ["Model parameters.....",
                "Family          : %s" % self.family,
                "Span            : %s" % self.span,
                "Degree          : %s" % self.degree,
                "Normalized      : %s" % self.normalize,
                "Parametric      : %s" % self.parametric_flags[:self.npar],
                "Drop_square     : %s" % self.drop_square_flags[:self.npar]
                ]
        return '\n'.join(strg)
        
#####---------------------------------------------------------------------------
#---- ---- loess outputs ---
#####---------------------------------------------------------------------------
cdef class loess_outputs:
    """Outputs of a loess fit. This object is automatically created with empty
values when a new loess object is instantiated. The object gets filled when the 
loess.fit() method is called.
    
:IVariables:
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
    diagonal : ndarray
        Diagonal of the operator hat matrix.
    robust : ndarray
        The (n,) ndarray of robustness weights for robust fitting.
    divisor : ndarray
        The (p,) array of normalization divisors for numeric predictors.
    """
    cdef c_loess.c_loess_outputs *_base
    cdef long nobs, npar
    cdef readonly int activated
    cdef setup(self,c_loess.c_loess_outputs *base, long nobs, long npar):
        self._base = base
        self.nobs = nobs
        self.npar = npar
        self.activated = False
    #........
    property fitted_values:    
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.fitted_values)
    #.........
    property fitted_residuals:
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.fitted_residuals)
    #.........
    property pseudovalues:
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.pseudovalues)
    #.........
    property diagonal:
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.diagonal)
    #.........
    property robust:
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.robust)
    #.........
    property divisor:
        def __get__(self):
            return floatarray_from_data(self.npar, 1, self._base.divisor)
    #.........    
    property enp:
        def __get__(self):
            return self._base.enp
    #.........
    property s:
        def __get__(self):
            return self._base.s
    #.........
    property one_delta:
        def __get__(self):
            return self._base.one_delta 
    #.........
    property two_delta:
        def __get__(self):
            return self._base.two_delta
    #.........
    property trace_hat:
        def __get__(self):
            return self._base.trace_hat
    #.........
    def __str__(self):
        strg = ["Outputs................",
                "Fitted values         : %s\n" % self.fitted_values,
                "Fitted residuals      : %s\n" % self.fitted_residuals,
                "Eqv. nb of parameters : %s" % self.enp,
                "Residual error        : %s" % self.s,
                "Deltas                : %s - %s" % (self.one_delta, self.two_delta),
                "Normalization factors : %s" % self.divisor,]
        return '\n'.join(strg)


        
#####---------------------------------------------------------------------------
#---- ---- loess confidence ---
#####---------------------------------------------------------------------------
cdef class conf_intervals:
    """Pointwise confidence intervals of a loess-predicted object:
    
:IVariables:
    fit : ndarray
        Predicted values.
    lower : ndarray
        Lower bounds of the confidence intervals.
    upper : ndarray
        Upper bounds of the confidence intervals.
    """
    cdef c_loess.c_conf_inv _base
    cdef readonly ndarray lower, fit, upper
    #.........
    def __dealloc__(self):
        c_loess.pw_free_mem(&self._base)
    #.........
    cdef setup(self, c_loess.c_conf_inv base, long nest):
        self._base = base
        self.fit = floatarray_from_data(nest, 1, base.fit)
        self.upper = floatarray_from_data(nest, 1, base.upper)
        self.lower = floatarray_from_data(nest, 1, base.lower)
    #.........
    def __str__(self):
        cdef ndarray tmp_ndr
        tmp_ndr = numpy.r_[[self.lower,self.fit,self.upper]].T
        return "Confidence intervals....\nLower b./ fit / upper b.\n%s" % \
               tmp_ndr 

#####---------------------------------------------------------------------------
#---- ---- loess predictions ---
#####---------------------------------------------------------------------------
cdef class loess_predicted:
    """Predicted values and standard errors of a loess object

:IVariables:
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
        
    cdef c_loess.c_prediction _base
    cdef readonly long nest
    cdef readonly conf_intervals confidence_intervals
    #.........
    def __dealloc__(self):
        c_loess.pred_free_mem(&self._base)      
    #.........
    cdef setup(self, c_loess.c_prediction base, long nest):
        self._base = base
        self.nest = nest
    #.........
    property values:
        def __get__(self):
            return floatarray_from_data(self.nest, 1, self._base.fit)
    #.........
    property stderr:
        def __get__(self):
            return floatarray_from_data(self.nest, 1, self._base.se_fit)
    #.........
    property residual_scale:
        def __get__(self):
            return self._base.residual_scale
    #.........
    property df:
        def __get__(self):
            return self._base.df        
    #.........
    def confidence(self, coverage=0.95):
        """Returns the pointwise confidence intervals for each predicted values,
at the given confidence interval coverage.
        
:Parameters:
    coverage : float
        Confidence level of the confidence intervals limits, as a fraction.
        
:Returns:
    A new conf_intervals object, consisting of:
    fit : ndarray
        Predicted values.
    lower : ndarray
        Lower bounds of the confidence intervals.
    upper : ndarray
        Upper bounds of the confidence intervals.
        """
        cdef c_loess.c_conf_inv _confintv
        if coverage < 0.5:
            coverage = 1. - coverage 
        if coverage > 1. :
            raise ValueError("The coverage precentage should be between 0 and 1!")
        c_loess.c_pointwise(&self._base, self.nest, coverage, &_confintv)
        self.confidence_intervals = conf_intervals()
        self.confidence_intervals.setup(_confintv, self.nest)
        return self.confidence_intervals
    #.........
    def __str__(self):
        strg = ["Outputs................",
                "Predicted values      : %s\n" % self.values,
                "Predicted std error   : %s\n" % self.stderr,
                "Residual scale        : %s" % self.residual_scale,
                "Degrees of freedom    : %s" % self.df,
#                "Confidence intervals  : %s" % self.confidence,
                ]
        return '\n'.join(strg)
    

#####---------------------------------------------------------------------------
#---- ---- loess base class ---
#####---------------------------------------------------------------------------
cdef class loess:
    """
    
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

    """
    cdef c_loess.c_loess _base
    cdef readonly loess_inputs inputs
    cdef readonly loess_model model
    cdef readonly loess_control control
    cdef readonly loess_kd_tree kd_tree
    cdef readonly loess_outputs outputs
    cdef readonly loess_predicted predicted
    
    def __init__(self, object x, object y, object weights=None, **options):
        #
        cdef ndarray x_ndr, y_ndr
        cdef double *x_dat, *y_dat
        cdef int i
        # Initialize the inputs .......
        self.inputs = loess_inputs(x, y)
        x_dat = <double *>self.inputs.x_eff.data
        y_dat = <double *>self.inputs.y_eff.data
        n = self.inputs.nobs
        p = self.inputs.npar
        c_loess.loess_setup(x_dat, y_dat, n, p, &self._base)
        # Sets the _base input .........
        self.inputs._base = &self._base.inputs
        # Initialize the model parameters
        self.model = loess_model()
        self.model.setup(&self._base.model, p)
        # Initialize the control parameters
        self.control = loess_control()
        self.control._base = &self._base.control
        # Initialize the kd tree ......
        self.kd_tree = loess_kd_tree()
        self.kd_tree._base = &self._base.kd_tree
        # Initialize the outputs ......
        self.outputs = loess_outputs()
        self.outputs.setup(&self._base.outputs, n, p,)
        # Process options .............
        modelopt = {}
        controlopt = {}
        for (k,v) in options.iteritems():
            if k in ('family', 'span', 'degree', 'normalize', 
                     'parametric', 'drop_square',):
                modelopt[k] = v
            elif k in ('surface', 'statistics', 'trace_hat', 
                       'iterations', 'cell'):
                controlopt[k] = v
        self.control.update(**controlopt)
        self.model.update(**modelopt)
    #......................................................
    def fit(self):
        """Computes the loess parameters on the current inputs and sets of parameters."""
        c_loess.loess_fit(&self._base)
        self.outputs.activated = True
        if self._base.status.err_status:
            raise ValueError(self._base.status.err_msg)
        return
    #......................................................
    def input_summary(self):
        """Returns some generic information about the loess parameters.
        """
        toprint = [str(self.inputs), str(self.model), str(self.control)]
        return "\n".join(toprint)
        
    def output_summary(self):
        """Returns some generic information about the loess fit."""
        print "Number of Observations         : %d" % self.inputs.nobs
        print "Fit flag                       : %d" % bool(self.outputs.activated)
        print "Equivalent Number of Parameters: %.1f" % self.outputs.enp
        if self.model.family == "gaussian":
            print "Residual Standard Error        : %.4f" % self.outputs.s
        else:
            print "Residual Scale Estimate        : %.4f" % self.outputs.s
    #......................................................
    def predict(self, newdata, stderror=False):
        """Computes loess estimates at the given new data points newdata. Returns
a loess_predicted object, whose attributes are described below.
        
:Parameters:
    newdata : ndarray
        The (m,p) array of independent variables where the surface must be estimated,
        with m the number of new data points, and p the number of independent
        variables.
    stderror : boolean
        Whether the standard error should be computed
        
:Returns:
    A new loess_predicted object, consisting of:
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
        cdef ndarray p_ndr
        cdef double *p_dat
        cdef c_loess.c_prediction _prediction
        cdef int i, m
        # Make sure there's been a fit earlier ...
        if self.outputs.activated == 0:
            c_loess.loess_fit(&self._base)
            self.outputs.activated = True
            if self._base.status.err_status:
                raise ValueError(self._base.status.err_msg)
        # Note : we need a copy as we may have to normalize
        p_ndr = narray(newdata, copy=True, subok=True, order='C').ravel()
        p_dat = <double *>p_ndr.data
        # Test the compatibility of sizes .......
        if p_ndr.size == 0:
            raise ValueError("Can't predict without input data !")
        (m, notOK) = divmod(len(p_ndr), self.inputs.npar)
        if notOK:
            raise ValueError(
                  "Incompatible data size: there should be as many rows as parameters")
        #.....
        c_loess.c_predict(p_dat, m, &self._base, &_prediction, stderror)
        if self._base.status.err_status:
            raise ValueError(self._base.status.err_msg)
        self.predicted = loess_predicted()
        self.predicted._base = _prediction
        self.predicted.nest = m
#        self.predicted.setup(_prediction, m)
        return self.predicted
    #.........
    def __dealloc__(self):
        c_loess.loess_free_mem(&self._base)
    #......................................................
    

#####---------------------------------------------------------------------------
#---- ---- loess anova ---
#####---------------------------------------------------------------------------
cdef class anova:
    cdef readonly double dfn, dfd, F_value, Pr_F
    #
    def __init__(self, loess_one, loess_two):
        cdef double one_d1, one_d2, one_s, two_d1, two_d2, two_s, rssdiff,\
                    d1diff, tmp, df1, df2
        #
        if not isinstance(loess_one, loess) or not isinstance(loess_two, loess):
            raise ValueError("Arguments should be valid loess objects!"\
                             "got '%s' instead" % type(loess_one))
        #
        out_one = loess_one.outputs
        out_two = loess_two.outputs
        #
        one_d1 = out_one.one_delta
        one_d2 = out_one.two_delta
        one_s = out_one.s
        #
        two_d1 = out_two.one_delta
        two_d2 = out_two.two_delta
        two_s = out_two.s
        #
        rssdiff = abs(one_s * one_s * one_d1 - two_s * two_s * two_d1)
        d1diff = abs(one_d1 - two_d1)
        self.dfn = d1diff * d1diff / abs(one_d2 - two_d2)
        df1 = self.dfn
        #
        if out_one.enp > out_two.enp:
            self.dfd = one_d1 * one_d1 / one_d2
            tmp = one_s
        else:
            self.dfd = two_d1 * two_d1 / two_d2
            tmp = two_s
        df2 = self.dfd
        F_value = (rssdiff / d1diff) / (tmp * tmp)
        
        self.Pr_F = 1. - c_loess.ibeta(F_value*df1/(df2+F_value*df1), df1/2, df2/2)
        self.F_value = F_value

        