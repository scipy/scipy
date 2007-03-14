# -*- Mode: Python -*-  
cimport c_python
cimport c_numpy
from c_numpy cimport ndarray, npy_intp, \
    PyArray_SIZE, PyArray_EMPTY, PyArray_FROMANY, \
    NPY_INT, NPY_DOUBLE, NPY_OWNDATA, NPY_ALIGNED, NPY_FORTRAN
import numpy

#narray = _N.array
#ndouble = _N.float
#nbool = _N.bool

# NumPy must be initialized
c_numpy.import_array()

cimport c_loess


cdef floatarray_from_data(rows, cols, double *data):
    cdef ndarray a_ndr
    cdef double *a_dat
    a_ndr = numpy.empty((rows*cols,), dtype=numpy.float)
    a_dat = <double *>a_ndr.data
    for i from 0 <= i < a_ndr.size:
        a_dat[i] = data[i]
    if cols > 1:
        a_ndr.shape = (rows, cols)
    return a_ndr

cdef boolarray_from_data(rows, cols, int *data):
    cdef ndarray a_ndr
    cdef int *a_dat
    a_ndr = numpy.empty((rows*cols,), dtype=numpy.int)
    a_dat = <int *>a_ndr.data
    for i from 0 <= i < a_ndr.size:
        a_dat[i] = data[i]
    if cols > 1:
        a_ndr.shape = (rows, cols)
    return a_ndr.astype(numpy.bool)

    
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
        "Weights"
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
    property nvar:
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
        def __get__(self):
            return self._base.surface
        def __set__(self, surface):
            self._base.surface = surface
    #.........
    property statistics:
        def __get__(self):
            return self._base.statistics
        def __set__(self, statistics):
            self._base.statistics = statistics
    #.........
    property trace_hat:
        def __get__(self):
            return self._base.trace_hat
        def __set__(self, trace_hat):
            self._base.trace_hat = trace_hat
    #.........
    property iterations:
        def __get__(self):
            return self._base.iterations
        def __set__(self, iterations):
            self._base.iterations = iterations
    #.........
    property cell:
        def __get__(self):
            return self._base.cell
        def __set__(self, cell):
            self._base.cell = cell
    #.........
    def update(self, **cellargs):
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
#    #.........
#    property parameter:
#        def __get__(self):
#            return self._base.parameter
#    #.........
#    property a:
#        def __get__(self):
#            return self._base.a
#    #.........
#    property vval:
#        def __get__(self):
#            return self._base.vert
#    #.........
#    property xi:
#        def __get__(self):
#            return self._base.xi
#    #.........
#    property vert:
#        def __get__(self):
#            return self._base.vert
#        return
#    
######---------------------------------------------------------------------------
##---- ---- loess model ---
######---------------------------------------------------------------------------
cdef class loess_model:
    cdef c_loess.c_loess_model *_base
    cdef long npar
    #.........
    property span:
        def __get__(self):
            return self._base.span
        def __set__(self, span):
            self._base.span = span
    #.........
    property degree:
        def __get__(self):
            return self._base.degree
    #.........
    property normalize:
        "Normalize the variables. Only useful if more than one variable..."
        def __get__(self):
            return bool(self._base.normalize)
        def __set__(self, normalize):
            self._base.normalize = normalize
    #.........
    property family:
        def __get__(self):
            return self._base.family
    #.........
    property parametric_flags:
        def __get__(self):
            return boolarray_from_data(8, 1, self._base.parametric)
        def __set__(self, paramf):
            cdef ndarray p_ndr
            cdef long *p_dat
            cdef int i
            p_ndr = <ndarray>PyArray_FROMANY(paramf, NPY_LONG, 1, 1, NPY_OWNDATA)
            p_dat = <long *>p_ndr.data
            for i in 0 <= i < max(8, p_ndr.size):
                self._base.parametric[i] = p_dat[i]
    #.........
    property drop_square_flags:
        def __get__(self):
            return boolarray_from_data(8, 1, self._base.drop_square)
        def __set__(self, drop_sq):
            cdef ndarray d_ndr
            cdef long *d_dat
            cdef int i
            d_ndr = <ndarray>PyArray_FROMANY(drop_sq, NPY_LONG, 1, 1, NPY_OWNDATA)
            d_dat = <long *>d_ndr.data
            for i in 0 <= i < max(8, d_ndr.size):
                self._base.drop_square[i] = d_dat[i]
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
    cdef long nobs
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
            return floatarray_from_data(self.nobs, 1, self._base.divisor)
    #.........    
    property enp:
        "Equivalent number of parameters."
        def __get__(self):
            return self._base.enp
    #.........
    property s:
#        ""
        def __get__(self):
            return self._base.s
    #.........
    property one_delta:
#        ""
        def __get__(self):
            return self._base.one_delta
    #.........
    property two_delta:
#        ""
        def __get__(self):
            return self._base.two_delta
    #.........
    property trace_hat:
#        ""
        def __get__(self):
            return self._base.trace_hat
#    #.........
#    def __str__(self):
#        strg = ["Outputs          :",
#                "enp       : %s" % self.enp,
#                "s : %s" % self.s,
#                "Deltas        : %s/%s" % (self.one_delta, self.two_delta),
#                "Divisor    : %s" % self.divisor,]
#        return '\n'.join(strg)


#####---------------------------------------------------------------------------
#---- ---- loess anova ---
#####---------------------------------------------------------------------------
cdef class loess_anova:
    cdef c_loess.c_anova *_base
    cdef long nest
    #.........
    property dfn:
        def __get__(self):
            return self._base.dfn
    #.........
    property dfd:
        def __get__(self):
            return self._base.dfd
    #.........
    property F_value:
        def __get__(self):
            return self._base.F_value
    #.........
    property Pr_F:
        def __get__(self):
            return self._base.Pr_F
        
#####---------------------------------------------------------------------------
#---- ---- loess confidence ---
#####---------------------------------------------------------------------------
cdef class confidence_interval:
    cdef c_loess.c_conf_inv *_base
    cdef nest
    #.........
    def __dealloc__(self):
        c_loess.pw_free_mem(self)
    #.........
    property fit:
        def __get__(self):
            return floatarray_from_data(nest, 1, self._base.fit)
    #.........
    property upper:
        def __get__(self):
            return floatarray_from_data(nest, 1, self._base.upper)
    #.........
    property lower:
        def __get__(self):
            return floatarray_from_data(nest, 1, self._base.lower)

#####---------------------------------------------------------------------------
#---- ---- loess predictions ---
#####---------------------------------------------------------------------------
cdef class loess_predicted:
    cdef c_loess.c_prediction *_base
    cdef long nest
    cdef confidence_interval conf_interval
    #.........
    def __dealloc__(self):
        c_loess.pred_free_mem(self._base)
    #.........
    property predicted:
        def __get__(self):
            return floatarray_from_data(nest, 1, self._base.fit)
    #.........
    property predicted_stderr:
        def __get__(self):
            return floatarray_from_data(nest, 1, self._base.se_fit)
    #.........
    property residual_scale:
        def __get__(self):
            return self.residual_scale
    #.........
    property df:
        def __get__(self):
            return self.df
    #.........
    def confidence(self, coverage=0.95):
        """
    coverage : float
        Confidence level of the confidence intervals limits as a fraction.
        """
        c_loess.pointwise(self._base, self.nest, coverage,
                          self.conf_interval._base)
    

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
    
    def __init__(self, object x, object y, object weights=None):
        #
        cdef ndarray x_ndr, y_ndr
        cdef double *x_dat, *y_dat
        cdef int i
        #
        x_ndr = <ndarray>PyArray_FROMANY(x, NPY_DOUBLE, 1, 1, NPY_FORTRAN)
        y_ndr = <ndarray>PyArray_FROMANY(y, NPY_DOUBLE, 1, 1, NPY_FORTRAN)
        x_dat = <double *>x_ndr.data
        y_dat = <double *>y_ndr.data
        n = len(x_ndr)
        p = x_ndr.size / n
        c_loess.loess_setup(x_dat, y_dat, n, p, &self._base)
        #
        self.inputs = loess_inputs()
        self.inputs._base = &self._base.inputs
        #
        self.model = loess_model()
        self.model._base = &self._base.model
        self.model.npar = p
        #
        self.control = loess_control()
        self.control._base = &self._base.control
        #
        self.kd_tree = loess_kd_tree()
        self.kd_tree._base = &self._base.kd_tree
        #
        self.outputs = loess_outputs()
        self.outputs._base = &self._base.outputs
        self.outputs.nobs = n
    #......................................................
    def fit(self):
        c_loess.loess_fit(&self._base)
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
    def predict(self, newdata, stderr=False):
        """
    newdata: ndarray
        A (m,p) ndarray specifying the values of the predictors at which the 
        evaluation is to be carried out.
    stderr: Boolean
        Logical flag for computing standard errors at newdata.
        """
        cdef ndarray p_ndr
        cdef double *p_dat
        cdef int i, m
        #
        p_ndr = <ndarray>PyArray_FROMANY(newdata, NPY_DOUBLE, 1, self.nvar, NPY_FORTRAN)
        p_dat = <double *>p_ndr.data
        m = len(p_ndr)
        c_loess.predict(p_dat, m, &self._base, self.predicted._base, stderr)
        self.predicted.nest = m
    #......................................................
#    def pointwisevoid(predicted *pre, int m, double coverage,
#                      struct ci_struct *ci)
#        c_loess.pointwise(predicted *pre, int m, double coverage,
#                      struct ci_struct *ci)
#
#def anova(loess_one, loess_two):
#    cdef c_loess.c_anova result
#    
#    c_loess.anova(loess_one._base, loess_two._base, &result)
#        
        