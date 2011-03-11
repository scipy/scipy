# maxentropy.py: Routines for fitting maximum entropy models.

# Copyright: Ed Schofield, 2003-2006
# License: BSD-style (see LICENSE.txt in main source directory)

# Future imports must come before any code in 2.5
from __future__ import division

__author__ = "Ed Schofield"
__version__ = '2.1'
__changelog__ = """
This module is an adaptation of "ftwmaxent" by Ed Schofield, first posted
on SourceForge as part of the "textmodeller" project in 2002.  The
official repository is now SciPy (since Nov 2005); the SourceForge
ftwmaxent code will not be developed further.

------------

Change log:

Since 2.0:
* Code simplification.  Removed dualapprox(), gradapprox() and other
  alias methods for bigmodel objects.  Use dual(), grad() etc. instead.
* Added support for testing on an external sample during optimization.
* Removed incomplete support for the (slow) GIS algorithm

Since 2.0-alpha4:
* Name change maxent -> maxentropy
* Removed online (sequential) estimation of feature expectations and
  variances.

Since v2.0-alpha3:
(1) Name change ftwmaxent -> scipy/maxent
(2) Modification for inclusion in scipy tree.  Broke one big class into
    two smaller classes, one for small models, the other for large models.
    Here a 'small' model is one defined on a sample space small enough to sum
    over in practice, whereas a 'large' model is on a sample space that is
    high-dimensional and continuous or discrete but too large to sum over,
    and requires Monte Carlo simulation.
(3) Refactoring:
    self.Eapprox -> self.mu
    p_0 -> aux_dist
    p0 -> aux_dist
    p_dot -> aux_dist_dot
    qdot -> p_dot
    q_dot -> p_dot
    q_theta -> p_theta
    E_p -> E_p_tilde
    E_q -> E_p

Since v2.0-alpha2:
Using multiple static feature matrices is now supported.  The generator
function supplied to generate feature matrices is called matrixtrials'
times each iteration.  This is useful for variance estimation of the E
and log Z estimators across the trials, without drawing another sample
each iteration (when staticsample = True).

Since v2.0-alpha1:
Sample feature matrices, if used, are sampled on the fly with a supplied
generator function, optionally multiple times to estimate the sample
variance of the feature expectation estimates.  An alternative is the
online estimation alg.

Since v0.8.5:
Added code for online (sequential) estimation of feature expectations and
variances.


"""


import math, types, cPickle
import numpy as np
from numpy import exp, asarray
from scipy import optimize
from scipy.linalg import norm
from scipy.misc import logsumexp
from scipy.maxentropy.maxentutils import arrayexp, \
        innerprod, innerprodtranspose, columnmeans, columnvariances, \
        flatten, DivergenceError, sparsefeaturematrix


class basemodel(object):
    """A base class providing generic functionality for both small and
    large maximum entropy models.  Cannot be instantiated.
    """

    def __init__(self):
        self.format = self.__class__.__name__[:4]
        if self.format == 'base':
            raise ValueError("this class cannot be instantiated directly")
        self.verbose = False

        self.maxgtol = 1e-5
        # Required tolerance of gradient on average (closeness to zero,axis=0) for
        # CG optimization:
        self.avegtol = 1e-3
        # Default tolerance for the other optimization algorithms:
        self.tol = 1e-4
        # Default tolerance for stochastic approximation: stop if
        # ||params_k - params_{k-1}|| < paramstol:
        self.paramstol = 1e-5

        self.maxiter = 1000
        self.maxfun = 1500
        self.mindual = -100.    # The entropy dual must actually be
                                # non-negative, but the estimate may be
                                # slightly out with bigmodel instances
                                # without implying divergence to -inf
        self.callingback = False
        self.iters = 0          # the number of iterations so far of the
                                # optimization algorithm
        self.fnevals = 0
        self.gradevals = 0

        # Variances for a Gaussian prior on the parameters for smoothing
        self.sigma2 = None

        # Store the duals for each fn evaluation during fitting?
        self.storeduals = False
        self.duals = {}
        self.storegradnorms = False
        self.gradnorms = {}

        # Do we seek to minimize the KL divergence between the model and a
        # prior density p_0?  If not, set this to None; then we maximize the
        # entropy.  If so, set this to an array of the log probability densities
        # p_0(x) for each x in the sample space.  For bigmodel objects, set this
        # to an array of the log probability densities p_0(x) for each x in the
        # random sample from the auxiliary distribution.
        self.priorlogprobs = None

        # By default, use the sample matrix sampleF to estimate the
        # entropy dual and its gradient.  Otherwise, set self.external to
        # the index of the sample feature matrix in the list self.externalFs.
        # This applies to 'bigmodel' objects only, but setting this here
        # simplifies the code in dual() and grad().
        self.external = None
        self.externalpriorlogprobs = None


    def fit(self, K, algorithm='CG'):
        """Fit the maxent model p whose feature expectations are given
        by the vector K.

        Model expectations are computed either exactly or using Monte
        Carlo simulation, depending on the 'func' and 'grad' parameters
        passed to this function.

        For 'model' instances, expectations are computed exactly, by summing
        over the given sample space.  If the sample space is continuous or too
        large to iterate over, use the 'bigmodel' class instead.

        For 'bigmodel' instances, the model expectations are not computed
        exactly (by summing or integrating over a sample space) but
        approximately (by Monte Carlo simulation).  Simulation is necessary
        when the sample space is too large to sum or integrate over in
        practice, like a continuous sample space in more than about 4
        dimensions or a large discrete space like all possible sentences in a
        natural language.

        Approximating the expectations by sampling requires an instrumental
        distribution that should be close to the model for fast convergence.
        The tails should be fatter than the model.  This instrumental
        distribution is specified by calling setsampleFgen() with a
        user-supplied generator function that yields a matrix of features of a
        random sample and its log pdf values.

        The algorithm can be 'CG', 'BFGS', 'LBFGSB', 'Powell', or
        'Nelder-Mead'.

        The CG (conjugate gradients) method is the default; it is quite fast
        and requires only linear space in the number of parameters, (not
        quadratic, like Newton-based methods).

        The BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm is a
        variable metric Newton method.  It is perhaps faster than the CG
        method but requires O(N^2) instead of O(N) memory, so it is
        infeasible for more than about 10^3 parameters.

        The Powell algorithm doesn't require gradients.  For small models
        it is slow but robust.  For big models (where func and grad are
        simulated) with large variance in the function estimates, this
        may be less robust than the gradient-based algorithms.
        """
        dual = self.dual
        grad = self.grad

        if isinstance(self, bigmodel):
            # Ensure the sample matrix has been set
            if not hasattr(self, 'sampleF') and hasattr(self, 'samplelogprobs'):
                raise AttributeError("first specify a sample feature matrix"
                                      " using sampleFgen()")
        else:
            # Ensure the feature matrix for the sample space has been set
            if not hasattr(self, 'F'):
                raise AttributeError("first specify a feature matrix"
                                      " using setfeaturesandsamplespace()")

        # First convert K to a numpy array if necessary
        K = np.asarray(K, float)

        # Store the desired feature expectations as a member variable
        self.K = K

        # Sanity checks
        try:
            self.params
        except AttributeError:
            self.reset(len(K))
        else:
            assert len(self.params) == len(K)

        # Don't reset the number of function and gradient evaluations to zero
        # self.fnevals = 0
        # self.gradevals = 0

        # Make a copy of the parameters
        oldparams = np.array(self.params)

        callback = self.log

        if algorithm == 'CG':
            retval = optimize.fmin_cg(dual, oldparams, grad, (), self.avegtol, \
                                      maxiter=self.maxiter, full_output=1, \
                                      disp=self.verbose, retall=0,
                                      callback=callback)

            (newparams, fopt, func_calls, grad_calls, warnflag) = retval

        elif algorithm == 'LBFGSB':
            if callback is not None:
                raise NotImplementedError("L-BFGS-B optimization algorithm"
                        " does not yet support callback functions for"
                        " testing with an external sample")
            retval = optimize.fmin_l_bfgs_b(dual, oldparams, \
                        grad, args=(), bounds=self.bounds, pgtol=self.maxgtol,
                        maxfun=self.maxfun)
            (newparams, fopt, d) = retval
            warnflag, func_calls = d['warnflag'], d['funcalls']
            if self.verbose:
                print algorithm + " optimization terminated successfully."
                print "\tFunction calls: " + str(func_calls)
                # We don't have info on how many gradient calls the LBFGSB
                # algorithm makes

        elif algorithm == 'BFGS':
            retval = optimize.fmin_bfgs(dual, oldparams, \
                                        grad, (), self.tol, \
                                        maxiter=self.maxiter, full_output=1, \
                                        disp=self.verbose, retall=0, \
                                        callback=callback)

            (newparams, fopt, gopt, Lopt, func_calls, grad_calls, warnflag) = retval

        elif algorithm == 'Powell':
            retval = optimize.fmin_powell(dual, oldparams, args=(), \
                                   xtol=self.tol, ftol = self.tol, \
                                   maxiter=self.maxiter, full_output=1, \
                                   disp=self.verbose, retall=0, \
                                   callback=callback)

            (newparams, fopt, direc, numiter, func_calls, warnflag) = retval

        elif algorithm == 'Nelder-Mead':
            retval = optimize.fmin(dual, oldparams, args=(), \
                                   xtol=self.tol, ftol = self.tol, \
                                   maxiter=self.maxiter, full_output=1, \
                                   disp=self.verbose, retall=0, \
                                   callback=callback)

            (newparams, fopt, numiter, func_calls, warnflag) = retval

        else:
            raise AttributeError("the specified algorithm '" + str(algorithm)
                    + "' is unsupported.  Options are 'CG', 'LBFGSB', "
                    "'Nelder-Mead', 'Powell', and 'BFGS'")

        if np.any(self.params != newparams):
            self.setparams(newparams)
        self.func_calls = func_calls



    def dual(self, params=None, ignorepenalty=False, ignoretest=False):
        """Computes the Lagrangian dual L(theta) of the entropy of the
        model, for the given vector theta=params.  Minimizing this
        function (without constraints) should fit the maximum entropy
        model subject to the given constraints.  These constraints are
        specified as the desired (target) values self.K for the
        expectations of the feature statistic.

        This function is computed as:
            L(theta) = log(Z) - theta^T . K

        For 'bigmodel' objects, it estimates the entropy dual without
        actually computing p_theta.  This is important if the sample
        space is continuous or innumerable in practice.  We approximate
        the norm constant Z using importance sampling as in
        [Rosenfeld01whole].  This estimator is deterministic for any
        given sample.  Note that the gradient of this estimator is equal
        to the importance sampling *ratio estimator* of the gradient of
        the entropy dual [see my thesis], justifying the use of this
        estimator in conjunction with grad() in optimization methods that
        use both the function and gradient. Note, however, that
        convergence guarantees break down for most optimization
        algorithms in the presence of stochastic error.

        Note that, for 'bigmodel' objects, the dual estimate is
        deterministic for any given sample.  It is given as:

            L_est = log Z_est - sum_i{theta_i K_i}

        where
            Z_est = 1/m sum_{x in sample S_0} p_dot(x) / aux_dist(x),

        and m = # observations in sample S_0, and K_i = the empirical
        expectation E_p_tilde f_i (X) = sum_x {p(x) f_i(x)}.
        """

        if self.external is None and not self.callingback:
            if self.verbose:
                print "Function eval #", self.fnevals

        if params is not None:
            self.setparams(params)

        # Subsumes both small and large cases:
        L = self.lognormconst() - np.dot(self.params, self.K)

        if self.verbose and self.external is None:
            print "  dual is ", L

        # Use a Gaussian prior for smoothing if requested.
        # This adds the penalty term \sum_{i=1}^m \params_i^2 / {2 \sigma_i^2}.
        # Define 0 / 0 = 0 here; this allows a variance term of
        # sigma_i^2==0 to indicate that feature i should be ignored.
        if self.sigma2 is not None and ignorepenalty==False:
            ratios = np.nan_to_num(self.params**2 / self.sigma2)
            # Why does the above convert inf to 1.79769e+308?

            L += 0.5 * ratios.sum()
            if self.verbose and self.external is None:
                print "  regularized dual is ", L

        if not self.callingback and self.external is None:
            if hasattr(self, 'callback_dual') \
                               and self.callback_dual is not None:
                # Prevent infinite recursion if the callback function
                # calls dual():
                self.callingback = True
                self.callback_dual(self)
                self.callingback = False

        if self.external is None and not self.callingback:
            self.fnevals += 1

        # (We don't reset self.params to its prior value.)
        return L


    # An alias for the dual function:
    entropydual = dual

    def log(self, params):
        """This method is called every iteration during the optimization
        process.  It calls the user-supplied callback function (if any),
        logs the evolution of the entropy dual and gradient norm, and
        checks whether the process appears to be diverging, which would
        indicate inconsistent constraints (or, for bigmodel instances,
        too large a variance in the estimates).
        """

        if self.external is None and not self.callingback:
            if self.verbose:
                print "Iteration #", self.iters

        # Store new dual and/or gradient norm
        if not self.callingback:
            if self.storeduals:
                self.duals[self.iters] = self.dual()
            if self.storegradnorms:
                self.gradnorms[self.iters] = norm(self.grad())

        if not self.callingback and self.external is None:
            if hasattr(self, 'callback'):
                # Prevent infinite recursion if the callback function
                # calls dual():
                self.callingback = True
                self.callback(self)
                self.callingback = False

        # Do we perform a test on external sample(s) every iteration?
        # Only relevant to bigmodel objects
        if hasattr(self, 'testevery') and self.testevery > 0:
            if (self.iters + 1) % self.testevery != 0:
                if self.verbose:
                    print "Skipping test on external sample(s) ..."
            else:
                self.test()

        if not self.callingback and self.external is None:
            if self.mindual > -np.inf and self.dual() < self.mindual:
                raise DivergenceError("dual is below the threshold 'mindual'"
                        " and may be diverging to -inf.  Fix the constraints"
                        " or lower the threshold!")

        self.iters += 1


    def grad(self, params=None, ignorepenalty=False):
        """Computes or estimates the gradient of the entropy dual.
        """

        if self.verbose and self.external is None and not self.callingback:
            print "Grad eval #" + str(self.gradevals)

        if params is not None:
            self.setparams(params)

        G = self.expectations() - self.K

        if self.verbose and self.external is None:
            print "  norm of gradient =",  norm(G)

        # (We don't reset params to its prior value.)

        # Use a Gaussian prior for smoothing if requested.  The ith
        # partial derivative of the penalty term is \params_i /
        # \sigma_i^2.  Define 0 / 0 = 0 here; this allows a variance term
        # of sigma_i^2==0 to indicate that feature i should be ignored.
        if self.sigma2 is not None and ignorepenalty==False:
            penalty = self.params / self.sigma2
            G += penalty
            features_to_kill = np.where(np.isnan(penalty))[0]
            G[features_to_kill] = 0.0
            if self.verbose and self.external is None:
                normG = norm(G)
                print "  norm of regularized gradient =", normG

        if not self.callingback and self.external is None:
            if hasattr(self, 'callback_grad') \
                               and self.callback_grad is not None:
                # Prevent infinite recursion if the callback function
                # calls grad():
                self.callingback = True
                self.callback_grad(self)
                self.callingback = False

        if self.external is None and not self.callingback:
            self.gradevals += 1

        return G


    def crossentropy(self, fx, log_prior_x=None, base=np.e):
        """Returns the cross entropy H(q, p) of the empirical
        distribution q of the data (with the given feature matrix fx)
        with respect to the model p.  For discrete distributions this is
        defined as:

            H(q, p) = - n^{-1} \sum_{j=1}^n log p(x_j)

        where x_j are the data elements assumed drawn from q whose
        features are given by the matrix fx = {f(x_j)}, j=1,...,n.

        The 'base' argument specifies the base of the logarithm, which
        defaults to e.

        For continuous distributions this makes no sense!
        """
        H = -self.logpdf(fx, log_prior_x).mean()
        if base != np.e:
            # H' = H * log_{base} (e)
            return H / np.log(base)
        else:
            return H


    def normconst(self):
        """Returns the normalization constant, or partition function, for
        the current model.  Warning -- this may be too large to represent;
        if so, this will result in numerical overflow.  In this case use
        lognormconst() instead.

        For 'bigmodel' instances, estimates the normalization term as
        Z = E_aux_dist [{exp (params.f(X))} / aux_dist(X)] using a sample
        from aux_dist.
        """
        return np.exp(self.lognormconst())


    def setsmooth(self, sigma):
        """Specifies that the entropy dual and gradient should be
        computed with a quadratic penalty term on magnitude of the
        parameters.  This 'smooths' the model to account for noise in the
        target expectation values or to improve robustness when using
        simulation to fit models and when the sampling distribution has
        high variance.  The smoothing mechanism is described in Chen and
        Rosenfeld, 'A Gaussian prior for smoothing maximum entropy
        models' (1999).

        The parameter 'sigma' will be squared and stored as self.sigma2.
        """
        self.sigma2 = sigma**2


    def setparams(self, params):
        """Set the parameter vector to params, replacing the existing
        parameters.  params must be a list or numpy array of the same
        length as the model's feature vector f.
        """

        self.params = np.array(params, float)        # make a copy

        # Log the new params to disk
        self.logparams()

        # Delete params-specific stuff
        self.clearcache()


    def clearcache(self):
        """Clears the interim results of computations depending on the
        parameters and the sample.
        """
        for var in ['mu', 'logZ', 'logZapprox', 'logv']:
            if hasattr(self, var):
                exec('del self.' + var)

    def reset(self, numfeatures=None):
        """Resets the parameters self.params to zero, clearing the cache
        variables dependent on them.  Also resets the number of function
        and gradient evaluations to zero.
        """

        if numfeatures:
            m = numfeatures
        else:
            # Try to infer the number of parameters from existing state
            if hasattr(self, 'params'):
                m = len(self.params)
            elif hasattr(self, 'F'):
                m = self.F.shape[0]
            elif hasattr(self, 'sampleF'):
                m = self.sampleF.shape[0]
            elif hasattr(self, 'K'):
                m = len(self.K)
            else:
                raise ValueError("specify the number of features / parameters")

        # Set parameters, clearing cache variables
        self.setparams(np.zeros(m, float))

        # These bounds on the param values are only effective for the
        # L-BFGS-B optimizer:
        self.bounds = [(-100., 100.)]*len(self.params)

        self.fnevals = 0
        self.gradevals = 0
        self.iters = 0
        self.callingback = False

        # Clear the stored duals and gradient norms
        self.duals = {}
        self.gradnorms = {}
        if hasattr(self, 'external_duals'):
            self.external_duals = {}
        if hasattr(self, 'external_gradnorms'):
            self.external_gradnorms = {}
        if hasattr(self, 'external'):
            self.external = None


    def setcallback(self, callback=None, callback_dual=None, \
                    callback_grad=None):
        """Sets callback functions to be called every iteration, every
        function evaluation, or every gradient evaluation. All callback
        functions are passed one argument, the current model object.

        Note that line search algorithms in e.g. CG make potentially
        several function and gradient evaluations per iteration, some of
        which we expect to be poor.
        """
        self.callback = callback
        self.callback_dual = callback_dual
        self.callback_grad = callback_grad

    def logparams(self):
        """Saves the model parameters if logging has been
        enabled and the # of iterations since the last save has reached
        self.paramslogfreq.
        """
        if not hasattr(self, 'paramslogcounter'):
            # Assume beginlogging() was never called
            return
        self.paramslogcounter += 1
        if not (self.paramslogcounter % self.paramslogfreq == 0):
            return

        # Check whether the params are NaN
        if not np.all(self.params == self.params):
            raise FloatingPointError("some of the parameters are NaN")

        if self.verbose:
            print "Saving parameters ..."
        paramsfile = open(self.paramslogfilename + '.' + \
                          str(self.paramslogcounter) + '.pickle', 'wb')
        cPickle.dump(self.params, paramsfile, cPickle.HIGHEST_PROTOCOL)
        paramsfile.close()
        #self.paramslog += 1
        #self.paramslogcounter = 0
        if self.verbose:
            print "Done."

    def beginlogging(self, filename, freq=10):
        """Enable logging params for each fn evaluation to files named
        'filename.freq.pickle', 'filename.(2*freq).pickle', ... each
        'freq' iterations.
        """
        if self.verbose:
            print "Logging to files " + filename + "*"
        self.paramslogcounter = 0
        self.paramslogfilename = filename
        self.paramslogfreq = freq
        #self.paramslog = 1

    def endlogging(self):
        """Stop logging param values whenever setparams() is called.
        """
        del self.paramslogcounter
        del self.paramslogfilename
        del self.paramslogfreq





class model(basemodel):
    """A maximum-entropy (exponential-form) model on a discrete sample
    space.
    """
    def __init__(self, f=None, samplespace=None):
        super(model, self).__init__()

        if f is not None and samplespace is not None:
            self.setfeaturesandsamplespace(f, samplespace)
        elif f is not None and samplespace is None:
            raise ValueError("not supported: specify both features and"
                    " sample space or neither")


    def setfeaturesandsamplespace(self, f, samplespace):
        """Creates a new matrix self.F of features f of all points in the
        sample space. f is a list of feature functions f_i mapping the
        sample space to real values.  The parameter vector self.params is
        initialized to zero.

        We also compute f(x) for each x in the sample space and store
        them as self.F.  This uses lots of memory but is much faster.

        This is only appropriate when the sample space is finite.
        """
        self.f = f
        self.reset(numfeatures=len(f))
        self.samplespace = samplespace
        self.F = sparsefeaturematrix(f, samplespace, 'csr_matrix')


    def lognormconst(self):
        """Compute the log of the normalization constant (partition
        function) Z=sum_{x \in samplespace} p_0(x) exp(params . f(x)).
        The sample space must be discrete and finite.
        """
        # See if it's been precomputed
        if hasattr(self, 'logZ'):
            return self.logZ

        # Has F = {f_i(x_j)} been precomputed?
        if not hasattr(self, 'F'):
            raise AttributeError("first create a feature matrix F")

        # Good, assume the feature matrix exists
        log_p_dot = innerprodtranspose(self.F, self.params)

        # Are we minimizing KL divergence?
        if self.priorlogprobs is not None:
            log_p_dot += self.priorlogprobs

        self.logZ = logsumexp(log_p_dot)
        return self.logZ


    def expectations(self):
        """The vector E_p[f(X)] under the model p_params of the vector of
        feature functions f_i over the sample space.
        """
        # For discrete models, use the representation E_p[f(X)] = p . F
        if not hasattr(self, 'F'):
            raise AttributeError("first set the feature matrix F")

        # A pre-computed matrix of features exists
        p = self.pmf()
        return innerprod(self.F, p)

    def logpmf(self):
        """Returns an array indexed by integers representing the
        logarithms of the probability mass function (pmf) at each point
        in the sample space under the current model (with the current
        parameter vector self.params).
        """
        # Have the features already been computed and stored?
        if not hasattr(self, 'F'):
            raise AttributeError("first set the feature matrix F")

        # Yes:
        # p(x) = exp(params.f(x)) / sum_y[exp params.f(y)]
        #      = exp[log p_dot(x) - logsumexp{log(p_dot(y))}]

        log_p_dot = innerprodtranspose(self.F, self.params)
        # Do we have a prior distribution p_0?
        if self.priorlogprobs is not None:
            log_p_dot += self.priorlogprobs
        if not hasattr(self, 'logZ'):
            # Compute the norm constant (quickly!)
            self.logZ = logsumexp(log_p_dot)
        return log_p_dot - self.logZ


    def pmf(self):
        """Returns an array indexed by integers representing the values
        of the probability mass function (pmf) at each point in the
        sample space under the current model (with the current parameter
        vector self.params).

        Equivalent to exp(self.logpmf())
        """
        return arrayexp(self.logpmf())

    # An alias for pmf
    probdist = pmf

    def pmf_function(self, f=None):
        """Returns the pmf p_theta(x) as a function taking values on the
        model's sample space.  The returned pmf is defined as:

            p_theta(x) = exp(theta.f(x) - log Z)

        where theta is the current parameter vector self.params.  The
        returned function p_theta also satisfies
            all([p(x) for x in self.samplespace] == pmf()).

        The feature statistic f should be a list of functions
        [f1(),...,fn(x)].  This must be passed unless the model already
        contains an equivalent attribute 'model.f'.

        Requires that the sample space be discrete and finite, and stored
        as self.samplespace as a list or array.
        """

        if hasattr(self, 'logZ'):
            logZ = self.logZ
        else:
            logZ = self.lognormconst()

        if f is None:
            try:
                f = self.f
            except AttributeError:
                raise AttributeError("either pass a list f of feature"
                           " functions or set this as a member variable self.f")

        # Do we have a prior distribution p_0?
        priorlogpmf = None
        if self.priorlogprobs is not None:
            try:
                priorlogpmf = self.priorlogpmf
            except AttributeError:
                raise AttributeError("prior probability mass function not set")

        def p(x):
            f_x = np.array([f[i](x) for i in range(len(f))], float)

            # Do we have a prior distribution p_0?
            if priorlogpmf is not None:
                priorlogprob_x = priorlogpmf(x)
                return math.exp(np.dot(self.params, f_x) + priorlogprob_x \
                                - logZ)
            else:
                return math.exp(np.dot(self.params, f_x) - logZ)
        return p


class conditionalmodel(model):
    """
    A conditional maximum-entropy (exponential-form) model p(x|w) on a
    discrete sample space.

    This is useful for classification problems:
    given the context w, what is the probability of each class x?

    The form of such a model is::

        p(x | w) = exp(theta . f(w, x)) / Z(w; theta)

    where Z(w; theta) is a normalization term equal to::

        Z(w; theta) = sum_x exp(theta . f(w, x)).

    The sum is over all classes x in the set Y, which must be supplied to
    the constructor as the parameter 'samplespace'.

    Such a model form arises from maximizing the entropy of a conditional
    model p(x | w) subject to the constraints::

        K_i = E f_i(W, X)

    where the expectation is with respect to the distribution::

        q(w) p(x | w)

    where q(w) is the empirical probability mass function derived from
    observations of the context w in a training set.  Normally the vector
    K = {K_i} of expectations is set equal to the expectation of f_i(w,
    x) with respect to the empirical distribution.

    This method minimizes the Lagrangian dual L of the entropy, which is
    defined for conditional models as::

        L(theta) = sum_w q(w) log Z(w; theta)
                   - sum_{w,x} q(w,x) [theta . f(w,x)]

    Note that both sums are only over the training set {w,x}, not the
    entire sample space, since q(w,x) = 0 for all w,x not in the training
    set.

    The partial derivatives of L are::

        dL / dtheta_i = K_i - E f_i(X, Y)

    where the expectation is as defined above.

    """
    def __init__(self, F, counts, numcontexts):
        """The F parameter should be a (sparse) m x size matrix, where m
        is the number of features and size is |W| * |X|, where |W| is the
        number of contexts and |X| is the number of elements X in the
        sample space.

        The 'counts' parameter should be a row vector stored as a (1 x
        |W|*|X|) sparse matrix, whose element i*|W|+j is the number of
        occurrences of x_j in context w_i in the training set.

        This storage format allows efficient multiplication over all
        contexts in one operation.
        """
        # Ideally, the 'counts' parameter could be represented as a sparse
        # matrix of size C x X, whose ith row # vector contains all points x_j
        # in the sample space X in context c_i.  For example:
        #     N = sparse.lil_matrix((len(contexts), len(samplespace)))
        #     for (c, x) in corpus:
        #         N[c, x] += 1

        # This would be a nicer input format, but computations are more
        # efficient internally with one long row vector.  What we really need is
        # for sparse matrices to offer a .reshape method so this conversion
        # could be done internally and transparently.  Then the numcontexts
        # argument to the conditionalmodel constructor could also be inferred
        # from the matrix dimensions.

        super(conditionalmodel, self).__init__()
        self.F = F
        self.numcontexts = numcontexts

        S = F.shape[1] // numcontexts          # number of sample point
        assert isinstance(S, int)

        # Set the empirical pmf:  p_tilde(w, x) = N(w, x) / \sum_c \sum_y N(c, y).
        # This is always a rank-2 beast with only one row (to support either
        # arrays or dense/sparse matrices.
        if not hasattr(counts, 'shape'):
            # Not an array or dense/sparse matrix
            p_tilde = asarray(counts).reshape(1, len(counts))
        else:
            if counts.ndim == 1:
                p_tilde = counts.reshape(1, len(counts))
            elif counts.ndim == 2:
                # It needs to be flat (a row vector)
                if counts.shape[0] > 1:
                    try:
                        # Try converting to a row vector
                        p_tilde = counts.reshape((1, counts.size))
                    except AttributeError:
                        raise ValueError("the 'counts' object needs to be a"
                            " row vector (1 x n) rank-2 array/matrix) or have"
                            " a .reshape method to convert it into one")
                else:
                    p_tilde = counts
        # Make a copy -- don't modify 'counts'
        self.p_tilde = p_tilde / p_tilde.sum()

        # As an optimization, p_tilde need not be copied or stored at all, since
        # it is only used by this function.

        self.p_tilde_context = np.empty(numcontexts, float)
        for w in xrange(numcontexts):
            self.p_tilde_context[w] = self.p_tilde[0, w*S : (w+1)*S].sum()

        # Now compute the vector K = (K_i) of expectations of the
        # features with respect to the empirical distribution p_tilde(w, x).
        # This is given by:
        #
        #     K_i = \sum_{w, x} q(w, x) f_i(w, x)
        #
        # This is independent of the model parameters.
        self.K = flatten(innerprod(self.F, self.p_tilde.transpose()))
        self.numsamplepoints = S


    def lognormconst(self):
        """Compute the elementwise log of the normalization constant
        (partition function) Z(w)=sum_{y \in Y(w)} exp(theta . f(w, y)).
        The sample space must be discrete and finite.  This is a vector
        with one element for each context w.
        """
        # See if it's been precomputed
        if hasattr(self, 'logZ'):
            return self.logZ

        numcontexts = self.numcontexts
        S = self.numsamplepoints
        # Has F = {f_i(x_j)} been precomputed?
        if not hasattr(self, 'F'):
            raise AttributeError("first create a feature matrix F")

        # Good, assume F has been precomputed

        log_p_dot = innerprodtranspose(self.F, self.params)

        # Are we minimizing KL divergence?
        if self.priorlogprobs is not None:
            log_p_dot += self.priorlogprobs

        self.logZ = np.zeros(numcontexts, float)
        for w in xrange(numcontexts):
            self.logZ[w] = logsumexp(log_p_dot[w*S: (w+1)*S])
        return self.logZ


    def dual(self, params=None, ignorepenalty=False):
        """The entropy dual function is defined for conditional models as

            L(theta) = sum_w q(w) log Z(w; theta)
                       - sum_{w,x} q(w,x) [theta . f(w,x)]

        or equivalently as

            L(theta) = sum_w q(w) log Z(w; theta) - (theta . k)

        where K_i = \sum_{w, x} q(w, x) f_i(w, x), and where q(w) is the
        empirical probability mass function derived from observations of the
        context w in a training set.  Normally q(w, x) will be 1, unless the
        same class label is assigned to the same context more than once.

        Note that both sums are only over the training set {w,x}, not the
        entire sample space, since q(w,x) = 0 for all w,x not in the training
        set.

        The entropy dual function is proportional to the negative log
        likelihood.

        Compare to the entropy dual of an unconditional model:
            L(theta) = log(Z) - theta^T . K
        """
        if not self.callingback:
            if self.verbose:
                print "Function eval #", self.fnevals

            if params is not None:
                self.setparams(params)

        logZs = self.lognormconst()

        L = np.dot(self.p_tilde_context, logZs) - np.dot(self.params, self.K)

        if self.verbose and self.external is None:
            print "  dual is ", L

        # Use a Gaussian prior for smoothing if requested.
        # This adds the penalty term \sum_{i=1}^m \theta_i^2 / {2 \sigma_i^2}
        if self.sigma2 is not None and ignorepenalty==False:
            penalty = 0.5 * (self.params**2 / self.sigma2).sum()
            L += penalty
            if self.verbose and self.external is None:
                print "  regularized dual is ", L

        if not self.callingback:
            if hasattr(self, 'callback_dual'):
                # Prevent infinite recursion if the callback function calls
                # dual():
                self.callingback = True
                self.callback_dual(self)
                self.callingback = False
            self.fnevals += 1

        # (We don't reset params to its prior value.)
        return L


    # These do not need to be overridden:
    #     grad
    #     pmf
    #     probdist


    def fit(self, algorithm='CG'):
        """Fits the conditional maximum entropy model subject to the
        constraints

            sum_{w, x} p_tilde(w) p(x | w) f_i(w, x) = k_i

        for i=1,...,m, where k_i is the empirical expectation
            k_i = sum_{w, x} p_tilde(w, x) f_i(w, x).
        """

        # Call base class method
        return model.fit(self, self.K, algorithm)


    def expectations(self):
        """The vector of expectations of the features with respect to the
        distribution p_tilde(w) p(x | w), where p_tilde(w) is the
        empirical probability mass function value stored as
        self.p_tilde_context[w].
        """
        if not hasattr(self, 'F'):
            raise AttributeError("need a pre-computed feature matrix F")

        # A pre-computed matrix of features exists

        numcontexts = self.numcontexts
        S = self.numsamplepoints
        p = self.pmf()
        # p is now an array representing p(x | w) for each class w.  Now we
        # multiply the appropriate elements by p_tilde(w) to get the hybrid pmf
        # required for conditional modelling:
        for w in xrange(numcontexts):
            p[w*S : (w+1)*S] *= self.p_tilde_context[w]

        # Use the representation E_p[f(X)] = p . F
        return flatten(innerprod(self.F, p))

        # # We only override to modify the documentation string.  The code
        # # is the same as for the model class.
        # return model.expectations(self)


    def logpmf(self):
        """Returns a (sparse) row vector of logarithms of the conditional
        probability mass function (pmf) values p(x | c) for all pairs (c,
        x), where c are contexts and x are points in the sample space.
        The order of these is log p(x | c) = logpmf()[c * numsamplepoints
        + x].
        """
        # Have the features already been computed and stored?
        if not hasattr(self, 'F'):
            raise AttributeError("first set the feature matrix F")

        # p(x | c) = exp(theta.f(x, c)) / sum_c[exp theta.f(x, c)]
        #      = exp[log p_dot(x) - logsumexp{log(p_dot(y))}]

        numcontexts = self.numcontexts
        S = self.numsamplepoints
        log_p_dot = flatten(innerprodtranspose(self.F, self.params))
        # Do we have a prior distribution p_0?
        if self.priorlogprobs is not None:
            log_p_dot += self.priorlogprobs
        if not hasattr(self, 'logZ'):
            # Compute the norm constant (quickly!)
            self.logZ = np.zeros(numcontexts, float)
            for w in xrange(numcontexts):
                self.logZ[w] = logsumexp(log_p_dot[w*S : (w+1)*S])
        # Renormalize
        for w in xrange(numcontexts):
            log_p_dot[w*S : (w+1)*S] -= self.logZ[w]
        return log_p_dot


class bigmodel(basemodel):
    """A maximum-entropy (exponential-form) model on a large sample
    space.

    The model expectations are not computed exactly (by summing or
    integrating over a sample space) but approximately (by Monte Carlo
    estimation).  Approximation is necessary when the sample space is too
    large to sum or integrate over in practice, like a continuous sample
    space in more than about 4 dimensions or a large discrete space like
    all possible sentences in a natural language.

    Approximating the expectations by sampling requires an instrumental
    distribution that should be close to the model for fast convergence.
    The tails should be fatter than the model.
    """

    def __init__(self):
        super(bigmodel, self).__init__()

        # Number of sample matrices to generate and use to estimate E and logZ
        self.matrixtrials = 1

        # Store the lowest dual estimate observed so far in the fitting process
        self.bestdual = float('inf')

        # Most of the attributes below affect only the stochastic
        # approximation procedure.  They should perhaps be removed, and made
        # arguments of stochapprox() instead.

        # Use Kersten-Deylon accelerated convergence fo stoch approx
        self.deylon = False

        # By default, use a stepsize decreasing as k^(-3/4)
        self.stepdecreaserate = 0.75

        # If true, check convergence using the exact model.  Only useful for
        # testing small problems (e.g. with different parameters) when
        # simulation is unnecessary.
        self.exacttest = False

        # By default use Ruppert-Polyak averaging for stochastic approximation
        self.ruppertaverage = True

        # Use the stoch approx scaling modification of Andradottir (1996)
        self.andradottir = False

        # Number of iterations to hold the stochastic approximation stepsize
        # a_k at a_0 for before decreasing it
        self.a_0_hold = 0

        # Whether or not to use the same sample for all iterations
        self.staticsample = True

        # How many iterations of stochastic approximation between testing for
        # convergence
        self.testconvergefreq = 0

        # How many sample matrices to average over when testing for convergence
        # in stochastic approx
        self.testconvergematrices = 10

        # Test for convergence every 'testevery' iterations, using one or
        # more external samples. If None, don't test.
        self.testevery = None
        # self.printevery = 1000


    def resample(self):
        """(Re)samples the matrix F of sample features.
        """

        if self.verbose >= 3:
            print "(sampling)"

        # First delete the existing sample matrix to save memory
        # This matters, since these can be very large
        for var in ['sampleF, samplelogprobs, sample']:
            if hasattr(self, var):
                exec('del self.' + var)

        # Now generate a new sample
        output = self.sampleFgen.next()
        try:
            len(output)
        except TypeError:
            raise ValueError("output of sampleFgen.next() not recognized")
        if len(output) == 2:
            # Assume the format is (F, lp)
            (self.sampleF, self.samplelogprobs) = output
        elif len(output) == 3:
            # Assume the format is (F, lp, sample)
            (self.sampleF, self.samplelogprobs, self.sample) = output
        else:
            raise ValueError("output of sampleFgen.next() not recognized")

        # Check whether the number m of features is correct
        try:
            # The number of features is defined as the length of
            # self.params, so first check if it exists:
            self.params
            m = len(self.params)
        except AttributeError:
            (m, n) = self.sampleF.shape
            self.reset(m)
        else:
            if self.sampleF.shape[0] != m:
                raise ValueError("the sample feature generator returned"
                                  " a feature matrix of incorrect dimensions")
        if self.verbose >= 3:
            print "(done)"

        # Now clear the temporary variables that are no longer correct for this
        # sample
        self.clearcache()


    def lognormconst(self):
        """Estimate the normalization constant (partition function) using
        the current sample matrix F.
        """
        # First see whether logZ has been precomputed
        if hasattr(self, 'logZapprox'):
            return self.logZapprox

        # Compute log v = log [p_dot(s_j)/aux_dist(s_j)]   for
        # j=1,...,n=|sample| using a precomputed matrix of sample
        # features.
        logv = self._logv()

        # Good, we have our logv.  Now:
        n = len(logv)
        self.logZapprox = logsumexp(logv) - math.log(n)
        return self.logZapprox


    def expectations(self):
        """Estimates the feature expectations E_p[f(X)] under the current
        model p = p_theta using the given sample feature matrix.  If
        self.staticsample is True, uses the current feature matrix
        self.sampleF.  If self.staticsample is False or self.matrixtrials
        is > 1, draw one or more sample feature matrices F afresh using
        the generator function supplied to sampleFgen().
        """
        # See if already computed
        if hasattr(self, 'mu'):
            return self.mu
        self.estimate()
        return self.mu

    def _logv(self):
        """This function helps with caching of interim computational
        results.  It is designed to be called internally, not by a user.

        This is defined as the array of unnormalized importance sampling
        weights corresponding to the sample x_j whose features are
        represented as the columns of self.sampleF.
            logv_j = p_dot(x_j) / q(x_j),
        where p_dot(x_j) = p_0(x_j) exp(theta . f(x_j)) is the
        unnormalized pdf value of the point x_j under the current model.
        """
        # First see whether logv has been precomputed
        if hasattr(self, 'logv'):
            return self.logv

        # Compute log v = log [p_dot(s_j)/aux_dist(s_j)]   for
        # j=1,...,n=|sample| using a precomputed matrix of sample
        # features.
        if self.external is None:
            paramsdotF = innerprodtranspose(self.sampleF, self.params)
            logv = paramsdotF - self.samplelogprobs
            # Are we minimizing KL divergence between the model and a prior
            # density p_0?
            if self.priorlogprobs is not None:
                logv += self.priorlogprobs
        else:
            e = self.external
            paramsdotF = innerprodtranspose(self.externalFs[e], self.params)
            logv = paramsdotF - self.externallogprobs[e]
            # Are we minimizing KL divergence between the model and a prior
            # density p_0?
            if self.externalpriorlogprobs is not None:
                logv += self.externalpriorlogprobs[e]

        # Good, we have our logv.  Now:
        self.logv = logv
        return logv


    def estimate(self):
        """This function approximates both the feature expectation vector
        E_p f(X) and the log of the normalization term Z with importance
        sampling.

        It also computes the sample variance of the component estimates
        of the feature expectations as: varE = var(E_1, ..., E_T) where T
        is self.matrixtrials and E_t is the estimate of E_p f(X)
        approximated using the 't'th auxiliary feature matrix.

        It doesn't return anything, but stores the member variables
        logZapprox, mu and varE.  (This is done because some optimization
        algorithms retrieve the dual fn and gradient fn in separate
        function calls, but we can compute them more efficiently
        together.)

        It uses a supplied generator sampleFgen whose .next() method
        returns features of random observations s_j generated according
        to an auxiliary distribution aux_dist.  It uses these either in a
        matrix (with multiple runs) or with a sequential procedure, with
        more updating overhead but potentially stopping earlier (needing
        fewer samples).  In the matrix case, the features F={f_i(s_j)}
        and vector [log_aux_dist(s_j)] of log probabilities are generated
        by calling resample().

        We use [Rosenfeld01Wholesentence]'s estimate of E_p[f_i] as:
            {sum_j  p(s_j)/aux_dist(s_j) f_i(s_j) }
              / {sum_j p(s_j) / aux_dist(s_j)}.

        Note that this is consistent but biased.

        This equals:
            {sum_j  p_dot(s_j)/aux_dist(s_j) f_i(s_j) }
              / {sum_j p_dot(s_j) / aux_dist(s_j)}

        Compute the estimator E_p f_i(X) in log space as:
            num_i / denom,
        where
            num_i = exp(logsumexp(theta.f(s_j) - log aux_dist(s_j)
                        + log f_i(s_j)))
        and
            denom = [n * Zapprox]

        where Zapprox = exp(self.lognormconst()).

        We can compute the denominator n*Zapprox directly as:
            exp(logsumexp(log p_dot(s_j) - log aux_dist(s_j)))
          = exp(logsumexp(theta.f(s_j) - log aux_dist(s_j)))
        """

        if self.verbose >= 3:
            print "(estimating dual and gradient ...)"

        # Hereafter is the matrix code

        mus = []
        logZs = []

        for trial in range(self.matrixtrials):
            if self.verbose >= 2 and self.matrixtrials > 1:
                print "(trial " + str(trial) + " ...)"

            # Resample if necessary
            if (not self.staticsample) or self.matrixtrials > 1:
                self.resample()

            logv = self._logv()
            n = len(logv)
            logZ = self.lognormconst()
            logZs.append(logZ)

            # We don't need to handle negative values separately,
            # because we don't need to take the log of the feature
            # matrix sampleF. See my thesis, Section 4.4

            logu = logv - logZ
            if self.external is None:
                averages = innerprod(self.sampleF, arrayexp(logu))
            else:
                averages = innerprod(self.externalFs[self.external], \
                                     arrayexp(logu))
            averages /= n
            mus.append(averages)

        # Now we have T=trials vectors of the sample means.  If trials > 1,
        # estimate st dev of means and confidence intervals
        ttrials = len(mus)   # total number of trials performed
        if ttrials == 1:
            self.mu = mus[0]
            self.logZapprox = logZs[0]
            try:
                del self.varE       # make explicit that this has no meaning
            except AttributeError:
                pass
            return
        else:
            # The log of the variance of logZ is:
            #     -log(n-1) + logsumexp(2*log|Z_k - meanZ|)

            self.logZapprox = logsumexp(logZs) - math.log(ttrials)
            stdevlogZ = np.array(logZs).std()
            mus = np.array(mus)
            self.varE = columnvariances(mus)
            self.mu = columnmeans(mus)
            return


    def setsampleFgen(self, sampler, staticsample=True):
        """Initializes the Monte Carlo sampler to use the supplied
        generator of samples' features and log probabilities.  This is an
        alternative to defining a sampler in terms of a (fixed size)
        feature matrix sampleF and accompanying vector samplelogprobs of
        log probabilities.

        Calling sampler.next() should generate tuples (F, lp), where F is
        an (m x n) matrix of features of the n sample points x_1,...,x_n,
        and lp is an array of length n containing the (natural) log
        probability density (pdf or pmf) of each point under the
        auxiliary sampling distribution.

        The output of sampler.next() can optionally be a 3-tuple (F, lp,
        sample) instead of a 2-tuple (F, lp).  In this case the value
        'sample' is then stored as a class variable self.sample.  This is
        useful for inspecting the output and understanding the model
        characteristics.

        If matrixtrials > 1 and staticsample = True, (which is useful for
        estimating variance between the different feature estimates),
        sampler.next() will be called once for each trial
        (0,...,matrixtrials) for each iteration.  This allows using a set
        of feature matrices, each of which stays constant over all
        iterations.

        We now insist that sampleFgen.next() return the entire sample
        feature matrix to be used each iteration to avoid overhead in
        extra function calls and memory copying (and extra code).

        An alternative was to supply a list of samplers,
        sampler=[sampler0, sampler1, ..., sampler_{m-1}, samplerZ], one
        for each feature and one for estimating the normalization
        constant Z. But this code was unmaintained, and has now been
        removed (but it's in Ed's CVS repository :).

        Example use:
        >>> import spmatrix
        >>> model = bigmodel()
        >>> def sampler():
        ...     n = 0
        ...     while True:
        ...         f = spmatrix.ll_mat(1,3)
        ...         f[0,0] = n+1; f[0,1] = n+1; f[0,2] = n+1
        ...         yield f, 1.0
        ...         n += 1
        ...
        >>> model.setsampleFgen(sampler())
        >>> type(model.sampleFgen)
        <type 'generator'>
        >>> [model.sampleF[0,i] for i in range(3)]
        [1.0, 1.0, 1.0]

        We now set matrixtrials as a class property instead, rather than
        passing it as an argument to this function, where it can be
        written over (perhaps with the default function argument by
        accident) when we re-call this func (e.g. to change the matrix
        size.)
        """

        # if not sequential:
        assert type(sampler) is types.GeneratorType
        self.sampleFgen = sampler
        self.staticsample = staticsample
        if staticsample:
            self.resample()


    def pdf(self, fx):
        """Returns the estimated density p_theta(x) at the point x with
        feature statistic fx = f(x).  This is defined as
            p_theta(x) = exp(theta.f(x)) / Z(theta),
        where Z is the estimated value self.normconst() of the partition
        function.
        """
        return exp(self.logpdf(fx))

    def pdf_function(self):
        """Returns the estimated density p_theta(x) as a function p(f)
        taking a vector f = f(x) of feature statistics at any point x.
        This is defined as:
            p_theta(x) = exp(theta.f(x)) / Z
        """
        log_Z_est = self.lognormconst()

        def p(fx):
            return np.exp(innerprodtranspose(fx, self.params) - log_Z_est)
        return p


    def logpdf(self, fx, log_prior_x=None):
        """Returns the log of the estimated density p(x) = p_theta(x) at
        the point x.  If log_prior_x is None, this is defined as:
            log p(x) = theta.f(x) - log Z
        where f(x) is given by the (m x 1) array fx.

        If, instead, fx is a 2-d (m x n) array, this function interprets
        each of its rows j=0,...,n-1 as a feature vector f(x_j), and
        returns an array containing the log pdf value of each point x_j
        under the current model.

        log Z is estimated using the sample provided with
        setsampleFgen().

        The optional argument log_prior_x is the log of the prior density
        p_0 at the point x (or at each point x_j if fx is 2-dimensional).
        The log pdf of the model is then defined as
            log p(x) = log p0(x) + theta.f(x) - log Z
        and p then represents the model of minimum KL divergence D(p||p0)
        instead of maximum entropy.
        """
        log_Z_est = self.lognormconst()
        if len(fx.shape) == 1:
            logpdf = np.dot(self.params, fx) - log_Z_est
        else:
            logpdf = innerprodtranspose(fx, self.params) - log_Z_est
        if log_prior_x is not None:
            logpdf += log_prior_x
        return logpdf


    def stochapprox(self, K):
        """Tries to fit the model to the feature expectations K using
        stochastic approximation, with the Robbins-Monro stochastic
        approximation algorithm: theta_{k+1} = theta_k + a_k g_k - a_k
        e_k where g_k is the gradient vector (= feature expectations E -
        K) evaluated at the point theta_k, a_k is the sequence a_k = a_0
        / k, where a_0 is some step size parameter defined as self.a_0 in
        the model, and e_k is an unknown error term representing the
        uncertainty of the estimate of g_k.  We assume e_k has nice
        enough properties for the algorithm to converge.
        """
        if self.verbose:
            print "Starting stochastic approximation..."

        # If we have resumed fitting, adopt the previous parameter k
        try:
            k = self.paramslogcounter
            #k = (self.paramslog-1)*self.paramslogfreq
        except:
            k = 0

        try:
            a_k = self.a_0
        except AttributeError:
            raise AttributeError("first define the initial step size a_0")

        avgparams = self.params
        if self.exacttest:
            # store exact error each testconvergefreq iterations
            self.SAerror = []
        while True:
            k += 1
            if k > self.a_0_hold:
                if not self.deylon:
                    n = k - self.a_0_hold
                elif k <= 2 + self.a_0_hold:   # why <= 2?
                    # Initialize n for the first non-held iteration
                    n = k - self.a_0_hold
                else:
                    # Use Kersten-Deylon accelerated SA, based on the rate of
                    # changes of sign of the gradient.  (If frequent swaps, the
                    # stepsize is too large.)
                    #n += (np.dot(y_k, y_kminus1) < 0)   # an indicator fn
                    if np.dot(y_k, y_kminus1) < 0:
                        n += 1
                    else:
                        # Store iterations of sign switches (for plotting
                        # purposes)
                        try:
                            self.nosignswitch.append(k)
                        except AttributeError:
                            self.nosignswitch = [k]
                        print "No sign switch at iteration " + str(k)
                    if self.verbose >= 2:
                        print "(using Deylon acceleration.  n is " + str(n) + " instead of " + str(k - self.a_0_hold) + "...)"
                if self.ruppertaverage:
                    if self.stepdecreaserate is None:
                        # Use log n / n as the default.  Note: this requires a
                        # different scaling of a_0 than a stepsize decreasing
                        # as, e.g., n^(-1/2).
                        a_k = 1.0 * self.a_0 * math.log(n) / n
                    else:
                        # I think that with Ruppert averaging, we need a
                        # stepsize decreasing as n^(-p), where p is in the open
                        # interval (0.5, 1) for almost sure convergence.
                        a_k = 1.0 * self.a_0 / (n ** self.stepdecreaserate)
                else:
                    # I think we need a stepsize decreasing as n^-1 for almost
                    # sure convergence
                    a_k = 1.0 * self.a_0 / (n ** self.stepdecreaserate)
            # otherwise leave step size unchanged
            if self.verbose:
                print "  step size is: " + str(a_k)

            self.matrixtrials = 1
            self.staticsample = False
            if self.andradottir:    # use Andradottir (1996)'s scaling?
                self.estimate()   # resample and reestimate
                y_k_1 = self.mu - K
                self.estimate()   # resample and reestimate
                y_k_2 = self.mu - K
                y_k = y_k_1 / max(1.0, norm(y_k_2)) + \
                      y_k_2 / max(1.0, norm(y_k_1))
            else:
                # Standard Robbins-Monro estimator
                if not self.staticsample:
                    self.estimate()   # resample and reestimate
                try:
                    y_kminus1 = y_k    # store this for the Deylon acceleration
                except NameError:
                    pass               # if we're on iteration k=1, ignore this
                y_k = self.mu - K
            norm_y_k = norm(y_k)
            if self.verbose:
                print "SA: after iteration " + str(k)
                print "  approx dual fn is: " + str(self.logZapprox \
                            - np.dot(self.params, K))
                print "  norm(mu_est - k) = " + str(norm_y_k)

            # Update params (after the convergence tests too ... don't waste the
            # computation.)
            if self.ruppertaverage:
                # Use a simple average of all estimates so far, which
                # Ruppert and Polyak show can converge more rapidly
                newparams = self.params - a_k*y_k
                avgparams = (k-1.0)/k*avgparams + 1.0/k * newparams
                if self.verbose:
                    print "  new params[0:5] are: " + str(avgparams[0:5])
                self.setparams(avgparams)
            else:
                # Use the standard Robbins-Monro estimator
                self.setparams(self.params - a_k*y_k)

            if k >= self.maxiter:
                print "Reached maximum # iterations during stochastic" \
                        " approximation without convergence."
                break


    def settestsamples(self, F_list, logprob_list, testevery=1, priorlogprob_list=None):
        """Requests that the model be tested every 'testevery' iterations
        during fitting using the provided list F_list of feature
        matrices, each representing a sample {x_j} from an auxiliary
        distribution q, together with the corresponding log probabiltiy
        mass or density values log {q(x_j)} in logprob_list.  This is
        useful as an external check on the fitting process with sample
        path optimization, which could otherwise reflect the vagaries of
        the single sample being used for optimization, rather than the
        population as a whole.

        If self.testevery > 1, only perform the test every self.testevery
        calls.

        If priorlogprob_list is not None, it should be a list of arrays
        of log(p0(x_j)) values, j = 0,. ..., n - 1, specifying the prior
        distribution p0 for the sample points x_j for each of the test
        samples.
        """
        # Sanity check
        assert len(F_list) == len(logprob_list)

        self.testevery = testevery
        self.externalFs = F_list
        self.externallogprobs = logprob_list
        self.externalpriorlogprobs = priorlogprob_list

        # Store the dual and mean square error based on the internal and
        # external (test) samples.  (The internal sample is used
        # statically for sample path optimization; the test samples are
        # used as a control for the process.)  The hash keys are the
        # number of function or gradient evaluations that have been made
        # before now.

        # The mean entropy dual and mean square error estimates among the
        # t external (test) samples, where t = len(F_list) =
        # len(logprob_list).
        self.external_duals = {}
        self.external_gradnorms = {}


    def test(self):
        """Estimate the dual and gradient on the external samples,
        keeping track of the parameters that yield the minimum such dual.
        The vector of desired (target) feature expectations is stored as
        self.K.
        """
        if self.verbose:
            print "  max(params**2)    = " + str((self.params**2).max())

        if self.verbose:
            print "Now testing model on external sample(s) ..."

        # Estimate the entropy dual and gradient for each sample.  These
        # are not regularized (smoothed).
        dualapprox = []
        gradnorms = []
        for e in xrange(len(self.externalFs)):
            self.external = e
            self.clearcache()
            if self.verbose >= 2:
                print "(testing with sample %d)" % e
            dualapprox.append(self.dual(ignorepenalty=True, ignoretest=True))
            gradnorms.append(norm(self.grad(ignorepenalty=True)))

        # Reset to using the normal sample matrix sampleF
        self.external = None
        self.clearcache()

        meandual = np.average(dualapprox,axis=0)
        self.external_duals[self.iters] = dualapprox
        self.external_gradnorms[self.iters] = gradnorms

        if self.verbose:
            print "** Mean (unregularized) dual estimate from the %d" \
                  " external samples is %f" % \
                 (len(self.externalFs), meandual)
            print "** Mean mean square error of the (unregularized) feature" \
                    " expectation estimates from the external samples =" \
                    " mean(|| \hat{\mu_e} - k ||,axis=0) =", np.average(gradnorms,axis=0)
        # Track the parameter vector params with the lowest mean dual estimate
        # so far:
        if meandual < self.bestdual:
            self.bestdual = meandual
            self.bestparams = self.params
            if self.verbose:
                print "\n\t\t\tStored new minimum entropy dual: %f\n" % meandual


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
