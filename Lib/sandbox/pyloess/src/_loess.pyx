# -*- Mode: Python -*-  
cimport c_python
cimport c_numpy
from c_numpy cimport ndarray, npy_intp, \
    PyArray_SIZE, PyArray_EMPTY, PyArray_FROMANY, \
    NPY_INT, NPY_DOUBLE, NPY_OWNDATA, NPY_ALIGNED, NPY_FORTRAN, \
    PyArray_SimpleNewFromData
import numpy
narray = numpy.array

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

##cimport modelflags
##import modelflags
#
#cdef list_to_clist(object p_list):
#    cdef int i, imax
#    p_list = list(p_list)
#    imax = min(8, len(p_list))
#    for i from 0 <= i < imax:
#        c_list[i] = p_list[i]
#    return c_list[0]
#cdef object clist_to_list(int c_list[8]):
#    cdef int i, imax
#    p_list = [False] * 8
#    imax = min(8, len(p_list))
#    for i from 0 <= i < imax:
#        p_list[i] = c_list[i]
#    return p_list
#        
#
#class modelflags:
#    def __init__(self):
#        self.str_list = [False] * 8
#        self.c_list = list_to_clist(self.str_list)
#    def __getitem__(self, idx):
#        return self.str_list[idx]
#    def __setitem__(self, idx, val):
#        cdef int tmpval
#        tmpval = val
#        self.c_list[idx] = tmpval
#        self.str_list[idx] = bool(val)
#    def __str__(self):
#        return str(self.str_list)
#
##class modelflags(c_modelflags):
##    def __init__(self):
##        c_modelflags.__init__(self)
##        


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
    cdef c_loess.c_loess_inputs *_base
    #.........
    property x:
        def __get__(self):
            return floatarray_from_data(self._base.n, self._base.p, self._base.x)
    #.........    
    property y:
        def __get__(self):
            return floatarray_from_data(self._base.n, 1, self._base.y)
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
            return floatarray_from_data(self._base.n, 1, self._base.weights)
        
        def __set__(self, w):
            cdef npy_intp *dims
            cdef ndarray w_ndr
            w_ndr = <ndarray>PyArray_FROMANY(w, NPY_DOUBLE, 1, 1, NPY_OWNDATA)
            if w_ndr.ndim > 1 or w_ndr.size != self._base.n:
                raise ValueError, "Invalid size of the 'weights' vector!"
            self._base.weights = <double *>w_ndr.data
    #.........    
    property nobs:
        "Number of observations."
        def __get__(self):
            return self._base.n
    #.........
    property npar:
        "Number of independent variables."
        def __get__(self):
            return self._base.p
#       
######---------------------------------------------------------------------------
##---- ---- loess control ---
######---------------------------------------------------------------------------
cdef class loess_control:
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
    cdef c_loess.c_loess_model *_base
    cdef long npar
#    cdef public double span
#    cdef public int degree
#    cdef public char *family
#    cdef public parametric_mflags, drop_square_mflags
    #.........
    cdef setup(self, c_loess.c_loess_model *base, long npar):
        self._base = base
        self.npar = npar
#        self.parametric_flags = modelflags()
#        self.parametric_flags.c_list[0] = base.parametric[0]
#        self.drop_square_flags = modelflags()
#        self.drop_square_flags.c_list[0] = base.drop_square[0]
#        self.span = self._base.span
#        self.degree = self._base.degree
#        self.family = self._base.family
#        self.parametric_flags = boolarray_from_data(self.npar, 1, self._base.parametric)
#        self.drop_square_flags = boolarray_from_data(self.npar, 1, self._base.drop_square)
        return self    
    #.........
    property normalize:
        """
    normalize : boolean [True]
        Determines whether the independent variables should be normalized.  
        If True, the normalization is performed by setting the 10% trimmed 
        standard deviation to one. If False, no normalization is carried out. 
        This option is only useful for more than one variable. For spatial
        coordinates predictors or variables with a common scale, it should be 
        set to False.
        """
        def __get__(self):
            return bool(self._base.normalize)
        def __set__(self, normalize):
            self._base.normalize = normalize
    #.........
    property span:
        """Smoothing factor, as a fraction of the number of points to take into
    account. By default, span=0.75."""
        def __get__(self):
            return self._base.span
        def __set__(self, span):
            if span <= 0. or span > 1.:
                raise ValueError("Span should be between 0 and 1!")
            self._base.span = span
    #.........
    property degree:
        """
    degree : integer [2]
        Overall degree of locally-fitted polynomial. 1 is locally-linear 
        fitting and 2 is locally-quadratic fitting.  Degree should be 2 at most.
        """
        def __get__(self):
            return self._base.degree
        def __set__(self, degree):
            if degree < 0 or degree > 2:
                raise ValueError("Degree should be between 0 and 2!")
    #.........
    property family:
        """
    family : string ["gaussian"]
        Determines the assumed distribution of the errors. The values are 
        "gaussian" or "symmetric". If "gaussian" is selected, the fit is 
        performed with least-squares. If "symmetric" is selected, the fit
        is performed robustly by redescending M-estimators.
        """    
        def __get__(self):
            return self._base.family
        def __set__(self, family):
            if family.lower() not in ('symmetric', 'gaussian'):
                raise ValueError("Invalid value for the 'family' argument: "\
                                 "should be in ('symmetric', 'gaussian').")
            self._base.family = family
    #.........
    property parametric_flags:
        """
    parametric_flags : sequence [ [False]*p ]
        Indicates which independent variables should be conditionally-parametric
       (if there are two or more independent variables). The argument should 
       be a sequence of booleans, with the same size as the number of independent 
       variables, specified in the order of the predictor group ordered in x. 
        """
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
        """
    drop_square : sequence [ [False]* p]
        When there are two or more independent variables and when a 2nd order
        polynomial is used, "drop_square_flags" specifies those numeric predictors 
        whose squares should be dropped from the set of fitting variables. 
        The method of specification is the same as for parametric.  
        """
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
        return "loess model parameters @%s" % id(self)
    #.........
    def __str__(self):
        strg = ["Model parameters.....",
                "family      : %s" % self.family,
                "span        : %s" % self.span,
                "degree      : %s" % self.degree,
                "normalized  : %s" % self.normalize,
                "parametric  : %s" % self.parametric_flags[:self.npar],
                "drop_square : %s" % self.drop_square_flags[:self.npar]
                ]
        return '\n'.join(strg)
        
#####---------------------------------------------------------------------------
#---- ---- loess outputs ---
#####---------------------------------------------------------------------------
cdef class loess_outputs:
    cdef c_loess.c_loess_outputs *_base
    cdef long nobs, npar
    cdef readonly int activated
    #........
    property fitted_values:    
        """
    fitted_values : ndarray
        The (n,) ndarray of fitted values.
        """
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.fitted_values)
    #.........
    property fitted_residuals:
        """
    fitted_residuals : ndarray
        The (n,) ndarray of fitted residuals (observations - fitted values).
        """
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.fitted_residuals)
    #.........
    property pseudovalues:
        """
    pseudovalues : ndarray
        The (n,) ndarray of adjusted values of the response when robust estimation 
        is used.
        """        
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.pseudovalues)
    #.........
    property diagonal:
        """
    diagonal :
        Diagonal of the operator hat matrix.
        """    
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.diagonal)
    #.........
    property robust:
        """
    robust : ndarray
        The (n,) ndarray of robustness weights for robust fitting.
        """
        def __get__(self):
            return floatarray_from_data(self.nobs, 1, self._base.robust)
    #.........
    property divisor:
        "Equivalent number of parameters."
        def __get__(self):
            return floatarray_from_data(self.npar, 1, self._base.divisor)
    #.........    
    property enp:
        """
    enp : float
        Equivalent number of parameters.
        """
        def __get__(self):
            return self._base.enp
    #.........
    property s:
        """
    s : float
        Estimate of the scale of residuals.
        """
        def __get__(self):
            return self._base.s
    #.........
    property one_delta:
        """
    one_delta: float
        Statistical parameter used in the computation of standard errors.
        """
        def __get__(self):
            return self._base.one_delta 
    #.........
    property two_delta:
        """
    two_delta : float
        Statistical parameter used in the computation of standard errors.
       """
        def __get__(self):
            return self._base.two_delta
    #.........
    property trace_hat:
        """
    trace_hat : float    
        Trace of the operator hat matrix.
        """
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
    cdef c_loess.c_conf_inv _base
    cdef readonly ndarray lower, fit, upper
    #.........
#    def __dealloc__(self):
#        c_loess.pw_free_mem(self._base)
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
    cdef c_loess.c_prediction _base
    cdef readonly long nest
    cdef readonly conf_intervals confidence_intervals
#    cdef readonly ndarray values, stderr
#    cdef readonly double residual_scale, df
    #.........
    def __dealloc__(self):
        c_loess.pred_free_mem(&self._base)      
    #.........
    cdef setup(self, c_loess.c_prediction base, long nest):
        self._base = base
        self.nest = nest
#    cdef setup(self, c_loess.c_loess loess_base, object newdata, stderror):
#        cdef ndarray p_ndr
#        cdef double *p_dat
#        cdef c_loess.c_prediction _prediction
#        cdef int i, m
#        #
#        # Note : we need a copy as we may have to normalize
#        p_ndr = narray(newdata, copy=True, subok=True, order='C').ravel()
#        p_dat = <double *>p_ndr.data
#        # Test the compatibility of sizes .......
#        if p_ndr.size == 0:
#            raise ValueError("Can't predict without input data !")
#        (m, notOK) = divmod(len(p_ndr), loess_base.inputs.p)
#        if notOK:
#            raise ValueError(
#                  "Incompatible data size: there should be as many rows as parameters")
#        #.....
#        c_loess.c_predict(p_dat, m, &loess_base, &_prediction, stderror)
#        if loess_base.status.err_status:
#            raise ValueError(loess_base.status.err_msg)
#        self._base = _prediction
#        self.nest = m
##        self.values = floatarray_from_data(m, 1, _prediction.fit)
##        self.stderr = floatarray_from_data(m, 1, _prediction.se_fit)
##        self.residual_scale = _prediction.residual_scale
##        self.df = _prediction.df
    #.........
    property values:
        """
    values : ndarray
        The (m,) ndarray of loess values evaluated at newdata
        """
        def __get__(self):
            return floatarray_from_data(self.nest, 1, self._base.fit)
    #.........
    property stderr:
        """
    stderr : ndarray
        The (m,) ndarray of the estimates of the standard error on the estimated
        values.
        """
        def __get__(self):
            return floatarray_from_data(self.nest, 1, self._base.se_fit)
    #.........
    property residual_scale:
        """
    residual_scale : float
        Estimate of the scale of the residuals
        """
        def __get__(self):
            return self._base.residual_scale
    #.........
    property df:
        """
    df : integer
        Degrees of freedom of the t-distribution used to compute pointwise 
        confidence intervals for the evaluated surface.
        """
        def __get__(self):
            return self._base.df        
    #.........
    def confidence(self, coverage=0.95):
        """Returns the pointwise confidence intervals for each predicted values,
at the given confidence interval coverage.
        
:Parameters:
    coverage : float
        Confidence level of the confidence intervals limits, as a fraction.
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
    cdef c_loess.c_loess _base
    cdef readonly loess_inputs inputs
    cdef readonly loess_model model
    cdef readonly loess_control control
    cdef readonly loess_kd_tree kd_tree
    cdef readonly loess_outputs outputs
    cdef readonly loess_predicted predicted
    cdef public long nobs, npar
    
    def __init__(self, object x, object y, object weights=None, **options):
        #
        cdef ndarray x_ndr, y_ndr
        cdef double *x_dat, *y_dat
        cdef int i
        # Get the predictor array
        x_ndr = narray(x, copy=True, subok=True, order='C')
        x_dat = <double *>x_ndr.data
        n = len(x_ndr)
        p = x_ndr.size / n
        self.npar = p
        self.nobs = n
        # Ravel the predictor array ...
        if p > 1:
            x_ndr.shape = (n*p,)
        # Get the response array ......
        y_ndr = narray(y, copy=False, subok=True, order='C')
        y_dat = <double *>y_ndr.data
        if y_ndr.size != n:
            raise ValueError("Incompatible size between the response array (%i)"\
                             " and the predictor array (%i)" % (y_ndr,n))
        # Initialization ..............
        c_loess.loess_setup(x_dat, y_dat, n, p, &self._base)
        #
        self.inputs = loess_inputs()
        self.inputs._base = &self._base.inputs
        #
        self.model = loess_model()
        self.model.setup(&self._base.model, p)
#        self.model._base = &self._base.model
#        self.model.npar = p
        #
        self.control = loess_control()
        self.control._base = &self._base.control
        #
        self.kd_tree = loess_kd_tree()
        self.kd_tree._base = &self._base.kd_tree
        #
        self.outputs = loess_outputs()
        self.outputs._base = &self._base.outputs
        self.outputs.activated = False
        self.outputs.nobs = n
        self.outputs.npar = p
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
        c_loess.loess_fit(&self._base)
        self.outputs.activated = True
        if self._base.status.err_status:
            raise ValueError(self._base.status.err_msg)
        return
    #......................................................
    def summary(self):
        print "Number of Observations         : %d" % self.inputs.nobs
        print "Equivalent Number of Parameters: %.1f" % self.outputs.enp
        if self.model.family == "gaussian":
            print "Residual Standard Error        : %.4f" % self.outputs.s
        else:
            print "Residual Scale Estimate        : %.4f" % self.outputs.s
    #......................................................
    def predict(self, newdata, stderror=False):
        """
    newdata: ndarray
        A (m,p) ndarray specifying the values of the predictors at which the 
        evaluation is to be carried out.
    stderr: Boolean
        Logical flag for computing standard errors at newdata.
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
        (m, notOK) = divmod(len(p_ndr), self.npar)
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

        