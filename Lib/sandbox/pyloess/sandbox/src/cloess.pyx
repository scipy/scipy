# -*- Mode: Python -*-  
cimport c_python
cimport c_numpy
from c_numpy cimport ndarray, npy_intp, \
    PyArray_SIZE, PyArray_EMPTY, PyArray_FROMANY, \
    NPY_INT, NPY_DOUBLE, NPY_OWNDATA, NPY_ALIGNED
import numpy as _N
# NumPy must be initialized
c_numpy.import_array()

cimport c_loess
    
#####---------------------------------------------------------------------------
#---- ---- loess model ---
#####---------------------------------------------------------------------------

cdef class loess_inputs:
    cdef c_loess.c_loess_inputs _inputs
    cdef public long nobs, nvar
    cdef readonly ndarray x, y, weights
    
    def __init__(self, object x, object y, object weights=None):
        cdef ndarray x_ndr, y_ndr
        cdef double *x_dat, *y_in_dat, *w_dat
        cdef npy_intp n, p
        cdef npy_intp *dims
        print "DEBUG: Initializing loess_inputs...",
        self.x = <ndarray>PyArray_FROMANY(x, NPY_DOUBLE, 1, 1, NPY_OWNDATA)
        self.y = <ndarray>PyArray_FROMANY(y, NPY_DOUBLE, 1, 1, NPY_OWNDATA)
        #
        if self.x.ndim > 2:
            raise ValueError,"Argument 'x' should be 2D at most!"
        n = len(self.x)
        p = self.x.size / n
        # Initialize the python side ............
        self.nobs = <long>n
        self.nvar = <long>p
        dims[0] = n
        # ... input weights ...
        if weights is None:
            self.weights = <ndarray>PyArray_EMPTY(1, dims, NPY_DOUBLE, NPY_ALIGNED)
            w_dat = <double *>self.weights.data
            for i from 0 <= i < dims[0]:
                w_dat[i] = 1
        else:
            self.weights = <ndarray>PyArray_FROMANY(weights, NPY_DOUBLE, 1, 1, NPY_OWNDATA)
            if self.weights.ndim > 1 or self.weights.size != n:
                raise ValueError, "Invalid size of the 'weights' vector!"
            w_dat = <double *>self.weights.data
        # Initialize the underlying C object ....
        self._inputs.n = self.nobs
        self._inputs.p = self.nvar
        self._inputs.x = <double *>self.x.data
        self._inputs.y = <double *>self.y.data
        self._inputs.weights = <double *>self.weights.data
        
#        for i from 0 <= i < (self.nobs*self.nvar):
#            self._inputs.x[i] = self.x.data[i]
#        for i from 0 <= i < self.nobs:
#            self._inputs.y[i] = self.y.data[i]
#            self._inputs.weights[i] = w_dat[i]
        print " OK."
        return

       
#####---------------------------------------------------------------------------
#---- ---- loess control ---
#####---------------------------------------------------------------------------
cdef class loess_control:
    cdef c_loess.c_loess_control _control
    cdef public char *surface, *statistics, *trace_hat
    cdef public double cell
    cdef public int iterations
    
    def __init__(self):
        print "DEBUG: Initializing loess_control...",
        self.surface = self._control.surface = "interpolate"
        self.statistics = self._control.statistics = "approximate"
        self.cell = self._control.cell = 0.2
        self.trace_hat = self._control.trace_hat = "wait.to.decide"
        self.iterations = self._control.iterations = 4
        print "OK."
        return 
        
    def __str__(self):
        strg = ["Control          :",
                "Surface type     : %s" % self.surface,
                "Statistics       : %s" % self.statistics,
                "Trace estimation : %s" % self.trace_hat,
                "Cell size        : %s" % self.cell,
                "Nb iterations    : %s" % self.iterations,]
        return '\n'.join(strg)
        
    
#####---------------------------------------------------------------------------
#---- ---- loess outputs ---
#####---------------------------------------------------------------------------
cdef class loess_outputs:
    cdef c_loess.c_loess_outputs _outputs
    cdef ndarray fitted_values, fitted_residuals, pseudovalues, diagonal, robust, divisor
    cdef double enp, s, one_delta, two_delta, trace_hat
    
    def __init__(self, n, p):
        cdef npy_intp *rows, *cols
        #cdef double *fv_dat, *fr_dat, *pv_dat, *diag_dat, *rw_dat, *div_dat
        rows[0] = <npy_intp>n
        cols[0] = <npy_intp>p
        print "DEBUG: Initializing loess_outputs...",
        # Initialize the python side ............
        self.fitted_values = <ndarray>PyArray_EMPTY(1, rows, NPY_DOUBLE, NPY_ALIGNED)
        self.fitted_residuals = <ndarray>PyArray_EMPTY(1, rows, NPY_DOUBLE, NPY_ALIGNED)
        self.pesudovalues = <ndarray>PyArray_EMPTY(1, rows, NPY_DOUBLE, NPY_ALIGNED)
        self.diagonal = <ndarray>PyArray_EMPTY(1, rows, NPY_DOUBLE, NPY_ALIGNED)
        self.robust = <ndarray>PyArray_EMPTY(1, rows, NPY_DOUBLE, NPY_ALIGNED)
        self.divisor = <ndarray>PyArray_EMPTY(1, cols, NPY_DOUBLE, NPY_ALIGNED)
        # Initialize the C side .................
        self._outputs.fitted_values = <double *>self.fitted_values.data
        self._outputs.fitted_residuals = <double *>self.fitted_residuals.data
        self._outputs.pseudovalues = <double *>self.pseudovalues.data
        self._outputs.diagonal = <double *>self.diagonal.data
        self._outputs.robust = <double *>self.robust.data
        self._outputs.divisor = <double *>self.divisor.data
        # Common initialization .................
        self.enp = self._outputs.enp = 0
        self.s = self._outputs.s = 0
        self.one_delta = self._outputs.one_delta = 0
        self.two_delta = self._outputs.two_delta = 0
        self.trace_hat = self._outputs.trace_hat = 0 
        print "OK."   
        
#    def __str__(self):
#        strg = ["Outputs          :",
#                "enp       : %s" % self.enp,
#                "s : %s" % self.s,
#                "Deltas        : %s/%s" % (self.one_delta, self.two_delta),
#                "Divisor    : %s" % self.divisor,]
#        return '\n'.join(strg)
    
#####---------------------------------------------------------------------------
#---- ---- loess kd_tree ---
#####---------------------------------------------------------------------------
cdef class loess_kd_tree:
    cdef c_loess.c_loess_kd_tree _kdtree
    cdef ndarray parameter, a, xi, vert, vval
    
    def __init__(self, long n, long p):
        cdef long maxkd, nval, nvert
        cdef npy_intp *nmaxkd, *nnval, *nnvert, *npars
        #
        print "DEBUG: Initializing loess_kdtree...",
        maxkd = max(n, 200)
        nvert = p * 2
        nval = (p+1) * maxkd
        # Initialize the python side ............
        print "(python side)",
        nmaxkd[0] = <npy_intp>maxkd
        nnvert[0] = <npy_intp>nvert
        nnval[0] = <npy_intp>nval
        npars[0] = <npy_intp>8
        self.parameter = <ndarray>PyArray_EMPTY(1, npars, NPY_LONG, NPY_ALIGNED)
        self.a = <ndarray>PyArray_EMPTY(1, nmaxkd, NPY_LONG, NPY_ALIGNED)
        self.xi = <ndarray>PyArray_EMPTY(1, nmaxkd, NPY_DOUBLE, NPY_ALIGNED)
        self.vert = <ndarray>PyArray_EMPTY(1, nnvert, NPY_DOUBLE, NPY_ALIGNED)
        self.vval = <ndarray>PyArray_EMPTY(1, nnval, NPY_DOUBLE, NPY_ALIGNED)
#        self.parameter = <ndarray>_N.empty((8,), _N.long, order='F')
#        self.a = <ndarray>_N.empty((maxkd,), _N.long, order='F')
#        self.xi = <ndarray>_N.empty((maxkd,), _N.float, order='F')
#        self.vert = <ndarray>_N.empty((p*2,), _N.float, order='F')
#        self.vval = <ndarray>_N.empty((nnval,), _N.float, order='F')
        # Initialize the C side .................
        print "(C side)",
        self._kdtree.parameter = <long *>self.parameter.data
        self._kdtree.a = <long *>self.a.data
        self._kdtree.xi = <double *>self.xi.data
        self._kdtree.vert = <double *>self.vert.data
        self._kdtree.vval = <double *>self.vval.data
#        # FIXME : Do we need to fill the arrays ?
#        for i from 0 <= i < 8:
#            self._kdtree.parameter[i] = 0
#        for i from 0 <= i < len(self.a):
#            self._kdtree.a[i] = 0
#            self._kdtree.xi[i] = 0
#        for i from 0 <= i < (p*2):
#            self._kdtree.vert[i] = 0
#        for i from 0 <= i < nnval:
#            self._kdtree.vval[i] = 0
        print "OK."
        return
    
#####---------------------------------------------------------------------------
#---- ---- loess model ---
#####---------------------------------------------------------------------------

cdef class loess_model:
    cdef c_loess.c_loess_model _model
    cdef double span
    cdef int degree, normalize
    cdef char *family
#    cdef ndarray parametric_flags, drop_square_flags
    cdef object parametric_flags, drop_square_flags
    #
    def __init__(self, double span=0.75, int degree=2, int normalize=1, 
                 object parametric_in=False, object drop_square_in=False, 
                 object family="gaussian"):
        cdef int i
        cdef int *parmf_dat, *dropsq_dat
        cdef npy_intp *npars
        print "DEBUG: Initializing loess_model...",
        self.span = self._model.span = span
        self.degree = self._model.degree = degree
        self.normalize = self._model.normalize = normalize
        self.family = self._model.family = family
        # FIXME : trying to use a ndarray crashes around here....
        self.parametric_flags = [None]*8
        if hasattr(parametric_in, '__len__'):
            for i from 0 <= i < len(parametric_in):
                self.parametric_flags[i] = parametric_in[i]
                self._model.parametric[i] = parametric_in[i]
        else:
            for i from 0 <= i <=7:
                self.parametric_flags[i] = parametric_in
                self._model.parametric[i] = parametric_in
        #....
        self.drop_square_flags = [None]*8
        if hasattr(drop_square_in, '__len__'):
            for i from 0 <= i < len(drop_square_in):
                self.drop_square_flags[i] = drop_square_in[i]
                self._model.drop_square[i] = drop_square_in[i]
        else:
            for i from 0 <= i < 8:
                self.drop_square_flags[i] = drop_square_in
                self._model.drop_square[i] = drop_square_in
        print "OK."
#        npars[0] = 8
#        self.parametric = <ndarray>PyArray_EMPTY(1, npars, NPY_INT, NPY_ALIGNED)
#        self.drop_square = <ndarray>PyArray_EMPTY(1, npars, NPY_INT, NPY_ALIGNED)
#        parmf_dat = <int *>self.parametric_flags.data
#        dropsq_dat = <int *>self.drop_square_flags.data
#        for i from 0 <= i < 8:
#            print "i:%i" % i
#            parmf_dat[i] = 0
#            self._model.parametric[i] = parmf_dat[i]
#            dropsq_dat[i] = 0
#            self._model.drop_square[i] = dropsq_dat[i]
#        print "DEBUG: loess_model: initialized"
        return
    #
    def __repr__(self):
        return "loess model parameters @%s" % id(self)
    def __str__(self):
        strg = ["Object      : %s" % self.__name__,
                "family      : %s" % self._model.family,
                "span        : %s" % self._model.span,
                "degree      : %s" % self._model.degree,
                "normalized  : %s" % self._model.normalize,
                "parametric  : %s" % self.parametric,
                "drop_square : %s" % self.drop_square]
        return '\n'.join(strg)
    
#####---------------------------------------------------------------------------
#---- ---- loess base class ---
#####---------------------------------------------------------------------------

cdef class loess:
#    cdef c_loess.c_loess *_base # If we try the pure C way
    cdef c_loess.c_loess _base
    cdef loess_inputs inputs
    cdef loess_model model
    cdef loess_control control
    cdef loess_kd_tree kd_tree
    cdef loess_outputs outputs
    
    def __init__(self, object x, object y, object weights=None):
        #
        cdef ndarray x_ndr, y_ndr
        cdef double *x_dat, *y_dat
        cdef long n, p
        cdef int i
        #
#        print "Try setup"
#        # CHECK : The following function segfaults :(
#        c_loess.loess_setup(x_dat, y_dat, n, d, self._base)
#        # CHECK : So we gonna try the hard way
        # Initialize the python side ............
        self.inputs = loess_inputs(x, y, weights)
        n = self.inputs.nobs
        p = self.inputs.nvar
        self.model = loess_model()
        self.control = loess_control()
        self.outputs = loess_outputs(n,p)
        self.kd_tree = loess_kd_tree(n,p)
        # Initialize the C side .................
        print "DEBUG:Initializing loess_cside"
        self._base.inputs = self.inputs._inputs
        self._base.model = self.model._model
        self._base.control = self.control._control
        self._base.kd_tree = self.kd_tree._kdtree
        self._base.outputs = self.outputs._outputs



#        self.inputs = base.in
#        self.model = base.model
#        self.kd_tree = base.kd_tree
#        self.outputs = base.out
        
    
    #......................................................
    def summary(self):
        print "Number of Observations         : %d" % self.inputs.n
        print "Equivalent Number of Parameters: %.1f" % self.outputs.enp
        if self.model.family == "gaussian":
            print "Residual Standard Error        : %.4f" % self.outputs.s
        else:
            print "Residual Scale Estimate        : %.4f" % self.outputs.s