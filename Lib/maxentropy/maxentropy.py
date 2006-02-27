# maxentropy.py: Routines for fitting maximum entropy models.

# Copyright: Ed Schofield, 2003-2006
# License: BSD-style (see LICENSE.txt in main source directory)

__author__ = "Ed Schofield"
__version__ = '2.0'
__changelog__ = """ 
This module is an adaptation of "ftwmaxent" by Ed Schofield, first posted
on SourceForge as part of the "textmodeller" project in 2002.  The
official repository is now SciPy (since Nov 2005); the SourceForge
ftwmaxent code will not be developed further.

------------

Change log:

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


from __future__ import division
import math, types, cPickle
import numpy
from scipy import optimize, sparse
from scipy.linalg import norm
from scipy.maxentropy.maxentutils import *


class basemodel(object):
    """A base class providing generic functionality for both small and
    large maximum entropy models.  Cannot be instantiated.
    """
    
    def __init__(self):
        self.format = self.__class__.__name__[:4]
        if self.format == 'base':
            raise ValueError, "this class cannot be instantiated directly"
        self.verbose = False
        
        self.maxgtol = 1e-5
        # Required tolerance of gradient on average (closeness to zero) for
        # CG optimization:
        self.avegtol = 1e-3
        # Default tolerance for the other optimization algorithms:
        self.tol = 1e-4       
        # Default tolerance for stochastic approximation: stop if
        # ||theta_k - theta_{k-1}|| < thetatol:
        self.thetatol = 1e-5         
        
        self.maxiter = 1000
        self.maxfun = 1500
        self.mindual = 0        # was -100.  But the entropy dual must be
                                # non-negative, and it seems any negative
                                # estimates indicate over-fitting to a
                                # particular sample
        #self.usecomplexlogs = False
        self.fnevals = 0
        self.gradevals = 0
        self.debug = 0
        
        # Variances for a Gaussian prior for smoothing
        self.sigma2 = None

        # Store the duals for each fn evaluation during fitting?
        self.storeduals = False
        self.duals = {}
        self.storegradnorms = False
        self.gradnorms = {}

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
            try:
                self.sampleF
                self.samplelogprobs
            except AttributeError:
                raise AttributeError, "first specify a feature matrix" \
                                      " using sampleFgen()"
        else:
            # Ensure the feature matrix for the sample space has been set
            try:
                self.F
            except AttributeError:
                raise AttributeError, "first specify a feature matrix" \
                                      " using setfeaturesandsamplespace()"
        
        # First convert K to a numpy array if necessary
        K = numpy.asarray(K, float)
            
        # Store the desired feature expectations as a member variable
        self.K = K
        
        # Sanity checks
        try:
            self.theta
        except AttributeError:
            self.resetparams(len(K))

        # Reset the number of function and gradient evaluations to zero
        self.fnevals = 0
        self.gradevals = 0
        
        # Make a copy of the parameters
        oldtheta = numpy.array(self.theta)
            
        assert len(K) == self.numconstraints()
        
        if algorithm == 'CG':
            retval = optimize.fmin_cg(dual, oldtheta, grad, (), self.avegtol, \
                                      maxiter=self.maxiter, full_output=1, \
                                      disp=self.verbose, retall=0)
            
            (newtheta, fopt, func_calls, grad_calls, warnflag) = retval

        elif algorithm == 'LBFGSB':
            retval = optimize.fmin_l_bfgs_b(dual, oldtheta, \
                        grad, args=(), bounds=self.bounds, pgtol=self.maxgtol,
                        maxfun=self.maxfun)
            (newtheta, fopt, d) = retval
            warnflag, func_calls = d['warnflag'], d['funcalls']
            if self.verbose:
                print algorithm + " optimization terminated successfully."
                print "\tFunction calls: " + str(func_calls)
                # We don't have info on how many gradient calls the LBFGSB
                # algorithm makes

        elif algorithm == 'BFGS':
            retval = optimize.fmin_bfgs(dual, oldtheta, \
                                        grad, (), self.tol, \
                                        maxiter=self.maxiter, full_output=1, \
                                        disp=self.verbose, retall=0)
            
            (newtheta, fopt, gopt, Lopt, func_calls, grad_calls, warnflag) = retval

        elif algorithm == 'Powell':
            retval = optimize.fmin_powell(dual, oldtheta, args=(), \
                                   xtol=self.tol, ftol = self.tol, \
                                   maxiter=self.maxiter, full_output=1, \
                                   disp=self.verbose, retall=0)
            
            (newtheta, fopt, direc, numiter, func_calls, warnflag) = retval
        
        elif algorithm == 'Nelder-Mead':
            retval = optimize.fmin(dual, oldtheta, args=(), \
                                   xtol=self.tol, ftol = self.tol, \
                                   maxiter=self.maxiter, full_output=1, \
                                   disp=self.verbose, retall=0)
            
            (newtheta, fopt, numiter, func_calls, warnflag) = retval

        else:
            raise AttributeError, "the specified algorithm '" + str(algorithm) \
                    + "' is unsupported.  Options are 'CG', 'LBFGSB', " \
                    "'Nelder-Mead', 'Powell', and 'BFGS'"
        
        if numpy.any(self.theta != newtheta):
            self.setparams(newtheta)
        self.func_calls = func_calls
    
  
    def fit_GIS(self, K=None, theta0=None, gtol=1e-6, maxiter=None):
        r"""A version of Darroch and Ratcliff's Generalized Iterative Scaling
        (GIS) algorithm, by default with the step-size approximation of
        Rosenfeld et al.  in "Whole-sentence exponential language models",
        Comp. Ling., 2002.
        
        Pass in 'K' as the desired (target) expectations of the features f(X),
        or set 'K' as a member variable of the model.
    
        By default the step sizes a_i used are those suggested by Rosenfeld et
        al., whose reciprocals are:
        
            F_i = \sum_s [p(s) f_i(s) \sum_k f_k(s)] / Ef_i
        
        The update rule is:
            x_{i+1} = x_i - a_i \bullet d_i,
        where \bullet is the Hadamard (componentwise) vector product, and
        where d_i is a difference of logs:
            d_i = log( E f_i(X) / K_i ) = log E f_i(X) - log K_i.
        Wherever K_i is 0, we define d_i as +inf.
        
        Note that the gradient is
            g_i = E f_i(X) - K_i
        so this differs from a gradient descent algorithm.
        """
        if maxiter is None:
            maxiter = self.maxiter 
        else:
            self.maxiter = maxiter
        if theta0 is None:
            # Initialize parameters with the zero vector by default
            theta = numpy.zeros(self.sampleF.shape[0], dtype=float)
        else:
            theta = theta0
        
        if K is None:
            try:
                K = self.K
            except AttributeError:
                raise AttributeError, "first set the member variable 'K' to" \
                        " the target feature expectations or pass them as an"\
                        " argument"
        assert len(K) == len(theta)
        iteration = 0
        # Cache variables
        (verbose, stepsizes) = (self.verbose, self._GISstepsizes)
        while True:
            self.setparams(theta)
            Ef = expectations()
            err = norm(Ef-K)
            if verbose:
                print "Iteration =", iteration
                print "Error =", err
                print "Entropy dual =", self.dual()
            if err < gtol:
                if verbose:
                    print "Desired tolerance reached"
                break
            # Compute a_i as the reciprocal of F, where
            #     F = \sum_s [p(s) f_i(s) \sum_k f_k(s)] / Ef_i
            a = stepsizes()
            if verbose >= 2:
                print "Step sizes =", a
                print "Elementwise product a * g =", a * g
            theta -= a * numpy.nan_to_num(log(Ef / K))   # elementwise product,
                                                    # replacing NaNs with zeros
            iteration += 1
            self.fnevals += 1   # actually, gradevals, but callback functions
                                # are only trigged after function evals
            if iteration >= self.maxiter:
                print "Warning: maximum number of iterations exceeded"
                break
        return theta
    
    
    def entropydual(self):
        """An alias for dual()
        """
        return self.dual()
    
    
    def dual(self, theta=None, ignorepenalty=False, ignoretest=False):
        """Computes the Lagrangian dual L(theta) of the entropy of the
        model.  Minimizing this function (without constraints) should
        fit the constrained model, where the desired (target) feature
        expectations need to be available as self.K.
                
        This function is computed as:
            L(theta) = log(Z) - theta^T . K

        This function approximates the entropy of the model p_theta
        with exponential form exp(theta^T . f_vec(x)) / Z(theta).
        
        For 'bigmodel' objects, it does this without actually computing
        p_theta.  This is important if the event space is practically
        innumerable.  We approximate the norm constant Z using importance
        sampling as in [Rosenfeld02whole].  Note that the gradient of
        this estimator is equal to the importance-sampling *ratio
        estimator* of the gradient of the entropy dual [see my paper,
        ICSLP 2004], so using this estimator in conjunction with the
        function gradapprox() in gradient-based optimization methods
        should be stable.  Also note that this estimator is deterministic
        for any given sample.
        """
        
        if self.external is None:
            if self.verbose:
                print "Function eval #", self.fnevals
            
            if theta is not None:
                self.setparams(theta)
        
        # Subsumes both small and large cases:
        L = self.lognormconst() - numpy.dot(self.theta, self.K)
            
        if self.verbose and self.external is None:
            print "  dual is ", L
        
        # Use a Gaussian prior for smoothing if requested.
        # This adds the penalty term \sum_{i=1}^m \theta_i^2 / {2 \sigma_i^2}
        if self.sigma2 is not None and ignorepenalty==False:
            penalty = 0.5 * (self.theta**2 / self.sigma2).sum()
            L += penalty
            if self.verbose and self.external is None:
                print "  regularized dual is ", L
        
        # Store new dual
        if self.external is None:
            if self.storeduals:
                self.duals[self.fnevals] = L
            
        # Test the new parameters with the test samples if they exist.
        if not ignoretest and isinstance(self, bigmodel) and \
                self.testevery is not None:
            M = self.testevery
            if self.fnevals % M == M - 1:
                self.test()
        
        # Call the callback function.  The complex-looking attribute-test
        # sequence is necessary to prevent an AttributeError in the
        # callback function from being caught and silenced
        # inappropriately.
        try:
            self.callback
        except AttributeError:
            pass
        else:
            self.callback(self)

        if self.external is None:
            if L < self.mindual:
                raise DivergenceError, \
                  "the dual is below the threshold 'mindual' and may be diverging" \
                  " to -inf.  Fix the constraints or lower the"\
                  " threshold."
            self.fnevals += 1
        
        # (We don't reset theta to its prior value.)
        return L
 
   
    def grad(self, theta=None, ignorepenalty=False):
        """Computes or estimates the gradient of the entropy dual.
        """
        
        if self.verbose:
            print "Grad eval #" + str(self.gradevals)

        # Has theta changed?  If so, use the new one, clearing the
        # speed-enhancing temporary variables that are functions of theta.
        #if numpy.any(self.theta != theta):
        # ** Always recompute
        if theta is not None:
            self.setparams(theta)
        
        G = self.expectations() - self.K
        
        if self.verbose and self.external is None:
            print "  norm of gradient =",  norm(G)

        # (We don't reset theta to its prior value.)
        
        # Use a Gaussian prior for smoothing if requested.
        # The ith partial derivative of the penalty term is 
        # \theta_i / \sigma_i^2
        if self.sigma2 is not None:
            penalty = self.theta / self.sigma2
            G += penalty
            if self.verbose and self.external is None:
                normG = norm(G)
                print "  norm of regularized gradient =", normG
        
        # Store gradient norm
        if self.external is None:
            if self.storegradnorms:
                # Re-use normG if it's already been computed
                try:
                    normG
                except UnboundLocalError:
                    self.gradnorms[self.gradevals] = norm(G)
                else:
                    self.gradnorms[self.gradevals] = normG
        
        if self.external is None:
            self.gradevals += 1
        
        return G
            
        
    def expectations(self):
        """For a 'model' object, the vector E_p[f(X)] under the model
        p_theta of the vector of feature functions f_i over the sample
        space.
        
        For a 'bigmodel' object, estimates of the feature expectations
        E_p[f(X)] under the current model p = p_theta using the given
        sample feature matrix.  If self.staticsample is True, use the
        current feature matrix F.  If self.staticsample is False or
        self.matrixtrials is > 1, draw one or more sample feature
        matrices F afresh using the supplied sampleFgen() generator
        function.
        """
        if isinstance(self, model):
            # If discrete, use the representation E_p[f(X)] = p . F
            try:
                self.F
            except AttributeError:
                raise AttributeError, "need a pre-computed feature matrix F"
 
            # A pre-computed matrix of features exists
            p = self.pmf()
            return innerprod(self.F, p)
        elif isinstance(self, bigmodel):
            try:
                # See if already computed
                self.mu
                if self.verbose >= 3:
                    print "(returning pre-computed expectations)"
                return self.mu
            except AttributeError:
                self.estimate()
                return self.mu
    
    
    def lognormconst(self):
        """Compute the log of the normalization constant (partition
        function) Z=sum_{x \in samplespace} exp(theta . f(x)).  The
        sample space must be discrete and finite.

        For 'bigmodel' instances, this function estimates the
        normalization constant using the current sample matrix F.
        """
        
        if isinstance(self, model):
            # See if it's been precomputed
            try:
                return self.logZ
            except AttributeError:
                pass
            
            # Has F = {f_i(x_j)} been precomputed?
            try:
                self.F
                
                # Good, assume F has been precomputed
                log_p_dot = innerprodtranspose(self.F, self.theta)
                self.logZ = logsumexp(log_p_dot)
                return self.logZ
            
            except AttributeError:
                raise AttributeError, "first create a feature matrix F"
        
        elif isinstance(self, bigmodel):
            # First see whether logZ has been precomputed
            try:
                return self.logZapprox
            except AttributeError:
                pass
 
            # Compute log w = log [p_dot(s_j)/aux_dist(s_j)]   for
            # j=1,...,n=|sample| using a precomputed matrix of sample
            # features.
            logw = self._logw()
            
            # Good, we have our logw.  Now:
            n = len(logw)
            self.logZapprox = logsumexp(logw) - math.log(n)
            return self.logZapprox
    
    
    def normconst(self):
        """Returns the normalization constant, or partition function, for
        the current model.  Warning -- this may be too large to represent,
        resulting in numerical overflow!  In this case use lognormconst()
        instead.
        """
        return math.exp(self.lognormconst())
    
    
    def _GISstepsizes(self):
        """A function that computes step sizes for the GIS iterative scaling
        algorithm.  This choice of step sizes was proposed by
        Rosenfeld et al., "Whole-sentence exponential language models", Comp.
        Ling. (2002).

        It returns a vector a of step sizes, one for each parameter params[i].
        These are given by:
            a_i = 1 / {\sum_x [p(x) f_i(x) \sum_k f_k(x)] /
                       [\sum_x [p(x) f_i(x)]]}
        """
        if isinstance(self, bigmodel):
            F = self.sampleF
        elif isinstance(self, model):
            F = self.F
        
        ## Overwrite the parameters
        #oldtheta = self.theta
        #if numpy.any(oldtheta != theta):
        #    self.setparams(theta)

        # Store the matrix S with elements s_ij = f_i(x_j) \sum_k f_k(x_j)
        try:
            self.S
        except AttributeError:
            F_column_sum = F.sum(axis=0).A.squeeze()
            if isinstance(F, numpy.ndarray):
                F = numpy.asarray(F)    # convert to array if it's a matrix
                                        # to get the elementwise product
                self.S = F * F_column_sum
            else:
                # Assume it's a sparse matrix
                n = F.shape[1]
                # Construct a sparse matrix D whose diagonals are F_column_sum
                # D = sparse.spdiags(F_column_sum, 0, n, n)
                D = sparse.dok_matrix((n, n))
                D.setdiag(F_column_sum)
                D = D.tocsc()
                # Now S is given by the matrix product F D.
                self.S = F * D

        # Compute b as the vector with components
        # b_j = \sum_k f_k(x_j) * p_dot(x_j)} = F_column_sum_j * p_dot(x_j)
        p_dot = flatten(arrayexp(innerprodtranspose(F, self.theta)))
        
        # Now the step sizes are given as follows:
        a = innerprod(F, p_dot) / innerprod(self.S, p_dot)     # elementwise division
        a = numpy.nan_to_num(a) # zero out all NaNs -- stepsize zero.
                                # This should be sufficient, since this occurs
                                # only when both num and denom are zero (so the
                                # expectation is correct already)
        return a 
    
    
    def numconstraints(self):
        return len(self.theta)


    def setsmooth(sigma):
        """Speficies that the entropy dual and gradient should be
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
    
    
    def setparams(self, theta):
        """Set the parameter vector theta, replacing the existing
        parameters.  theta must be a list or numpy array of the same
        length as the model's feature vector f.
        """

        # assert len(theta) == self.numconstraints()
        self.theta = numpy.array(theta, float)        # make a copy
        
        # Log the new params to disk
        self.logparams()
            
        # Delete theta-specific stuff
        self.clearcache()
        
       
    def clearcache(self):
        """Clears the interim results of computations depending on the
        parameters and the sample.
        """
        try:
            del self.mu
        except:
            pass
        try:
            del self.logZ               # used by 'model' objects
        except:
            pass
        try:
            del self.logZapprox         # used by 'bigmodel' objects
        except:
            pass
        try:
            del self.logw               # used by 'bigmodel' objects
        except:
            pass

    def resetparams(self, numfeatures=None):
        """Resets the parameters theta to zero, clearing variables
        dependent on them.
        """

        if numfeatures:
            m = numfeatures
        else:
            m = self.numconstraints()            
            
        # Set parameters, clearing cache variables
        self.setparams(numpy.zeros(m, float))  
        
        # These bounds on the param values are only effective for the
        # L-BFGS-B optimizer:
        self.bounds = [(-100., 100.)]*len(self.theta)  
        
   
    def setcallback(self, callback):
        """Sets a callback function to be called before the parameters
        are changed with setparams(theta).  The function is passed one
        argument, the current model before the parameters are changed for
        the next fn/grad evaluation.

        Note that line search algorithms in e.g. CG make potentially
        several evals per iteration, some of which we expect to be poor.
        """
        self.callback = callback
    
    def logparams(self):
        """Saves the model parameters if logging has been
        enabled and the # of iterations since the last save has reached
        self.thetalogfreq.
        """
        try:
            self.thetalogcounter += 1
            #if not (self.thetalog > 0 and self.thetalogcounter % self.thetalogfreq == 0):
            if not (self.thetalogcounter % self.thetalogfreq == 0):
                return
        except AttributeError:
            # Assume beginlogging() was never called
            return
        
        # Check whether the params are NaN
        if not numpy.all(self.theta == self.theta):
            raise FloatingPointError, "some of the parameters are NaN"

        if self.verbose:                    
            print "Saving parameters ..."
        paramsfile = open(self.thetalogfilename + '.' + \
                          str(self.thetalogcounter) + '.pickle', 'wb')
        cPickle.dump(self.theta, paramsfile, cPickle.HIGHEST_PROTOCOL)
        paramsfile.close()
        #self.thetalog += 1
        #self.thetalogcounter = 0
        if self.verbose:
            print "Done."
    
    def beginlogging(self, filename, freq=10):
        """Enable logging params for each fn evaluation to files named
        'filename.freq.pickle', 'filename.(2*freq).pickle', ... each
        'freq' iterations.
        """
        if self.verbose:
            print "Logging to files " + filename + "*"
        self.thetalogcounter = 0
        self.thetalogfilename = filename
        self.thetalogfreq = freq
        #self.thetalog = 1

    def endlogging(self):
        """Stop logging param values whenever setparams() is called.
        """
        del self.thetalogcounter
        del self.thetalogfilename 
        del self.thetalogfreq
 




class model(basemodel):
    """A maximum-entropy (exponential-form) model on a discrete sample
    space.
    """
    def __init__(self, f=None, samplespace=None):
        super(model, self).__init__()
        
        if f != None and samplespace != None:
            self.setfeaturesandsamplespace(f, samplespace)
        elif f != None and samplespace is None:
            raise ValueError, "not supported: specify both features and" \
                    " sample space or neither"
        
    def setfeaturesandsamplespace(self, f, samplespace):
        """Creates a new exponential model, where f is a list of feature
        functions f_i mapping the sample space to real values.  The
        parameter vector theta is initialized to the zero vector, of the
        same length as the list of feature functions f_i.
        
        We also compute f(x) for each x in the sample space and store
        them as self.F.  This uses lots of memory but is much faster.
        
        This is only appropriate when the sample space is finite.
        """
        self.f = f
        self.resetparams(numfeatures=len(f))
        self.samplespace = samplespace
        self.F = sparsefeaturematrix(f, samplespace, 'csr_matrix')
    
    
    def logpmf(self):
        """Returns an array indexed by integers representing the
        logarithms of the probability mass function (pmf) at each point
        in the sample space under the current model (with the current
        parameter vector self.theta).
        """
        # Have the features already been computed and stored?
        try:
            self.F
                
            # Yes:
            # p(x) = exp(theta.f(x)) / sum_y[exp theta.f(y)]
            #      = exp[log p_dot(x) - logsumexp{log(p_dot(y))}]
            
            log_p_dot = innerprodtranspose(self.F, self.theta)
            try:
                self.logZ
            except AttributeError:
                # Compute the norm constant (quickly!)
                self.logZ = logsumexp(log_p_dot)
            return log_p_dot - self.logZ
        except AttributeError:
            raise AttributeError, "first set the feature matrix F"
    
    
    def pmf(self):
        """Returns an array indexed by integers representing the values
        of the probability mass function (pmf) at each point in the
        sample space under the current model (with the current parameter
        vector self.theta).

        Equivalent to exp(self.logpmf())
        """
        return arrayexp(self.logpmf())
    
    def probdist(self):
        """An alias for pmf()
        """
        return self.pmf()
    
    def pmf_function(self, f=None):
        """Returns the pmf p_theta(x) as a function taking values on the
        model's sample space.  The returned pmf is defined as:
            
            p_theta(x) = exp(theta.f(x) - log Z)
        
        The returned function p_theta also satisfies
            all([p(x) for x in self.samplespace] == pmf()).

        The feature statistic f should be a list of functions
        [f1(),...,fn(x)].  This must be passed unless the model already
        contains an equivalent attribute 'model.f'.

        Requires that the sample space be discrete and finite, and stored
        as self.samplespace as a list or array.
        """
        
        try:
            logZ = self.logZ
        except AttributeError:
            logZ = self.lognormconst()
        
        if f is None:
            try:
                f = self.f
            except AttributeError:
                raise AttributeError, "either pass a list f of feature" \
                           " functions or set this as a member variable self.f"
        
        def p(x):
            f_x = numpy.array([f[i](x) for i in range(len(f))], float)
            #f_x = numpy.empty(len(f), float)
            #for i in range(len(f)):
            #    f_x[i] = f[i](x)
            return math.exp(numpy.dot(self.theta, f_x) - logZ)
        
        return p
    
  
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

        # By default, use the sample matrix sampleF to estimate the
        # entropy dual and its gradient.  Otherwise, set self.external to
        # the index of the sample feature matrix in the list self.externalFs.
        self.external = None
        
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
        
        # For comparing sampling methods and opt algorithms -- specifies that
        # we can compute the exact expectations at any iteration with
        # self.expectations() as a convergence criterion
        self.testconvergecheat = False   
        
        # Number of stdevs either side of the mean for Z and E confidence
        # intervals
        self.z = 3.0         
        
        # Desired precision with expectation estimates
        self.Etol = 5e-5     
        
        # Desired precision with logZ estimates
        self.Ztol = 5e-5        
        
        # Using relative precision for the sampling stopping criterion is
        # disabled by default:
        self.Ertol = 0.0
        
        # Number of samples to compute before tracking the variance
        self.init_samples = 10000     
        
        # Min number of samples to compute the variance of after starting
        # tracking
        self.min_samples = 10000      
        
        # Test for convergence every 'testevery' iterations when using the
        # sequential procedure
        self.testevery = 1000         
        self.printevery = 1000
    
    
    def lognormconstapprox(self):
        """Estimates the normalization constant (partition function)
        using the current sample matrix F.
        """
        return self.lognormconst()
    
    
    def normconstapprox(self):
        """Estimates the normalization constant (partition function) Z =
        E_aux_dist [{exp (theta.f(X))} / aux_dist(X)] using a sample from
        aux_dist.
        
        Warning -- this will cause numerical overflow if it is too large
        to represent. In this case use lognormconstapprox() instead!
        """
        
        return math.exp(self.lognormconst())
        
        
    def resample(self):
        """(Re)samples the matrix F of sample features.  
        """

        if self.verbose >= 3:
            print "(sampling)"
        
        # First delete the existing sample matrix to save memory
        # This matters, since these can be very large
        try:
            del self.sampleF
            del self.samplelogprobs
            del self.sample
        except AttributeError:
            pass
        # Now generate a new sample
        output = self.sampleFgen.next()
        try:
            len(output)
        except TypeError:
            raise ValueError, "output of sampleFgen.next() not recognized"
        if len(output) == 2:
            # Assume the format is (F, lp)
            (self.sampleF, self.samplelogprobs) = output
        elif len(output) == 3:
            # Assume the format is (F, lp, sample)
            (self.sampleF, self.samplelogprobs, self.sample) = output
        else:
            raise ValueError, "output of sampleFgen.next() not recognized"
        
        # Check whether the number m of features is correct
        try:
            # The number of features is defined as the length of
            # self.theta, so first check if it exists:
            self.theta
            m = self.numconstraints()
        except AttributeError:
            (m, n) = self.sampleF.shape
            self.resetparams(m)
        else:
            if self.sampleF.shape[0] != m:
                raise ValueError, "the sample feature generator returned" \
                                  " a feature matrix of incorrect dimensions"
        if self.verbose >= 3:
            print "(done)"

        # Now clear the temporary variables that are no longer correct for this
        # sample
        self.clearcache()

       
    def _logw(self):
        """This function helps with caching of interim computational
        results.  It is designed to be called internally, not by a user.

        This is defined as the array of unnormalized importance sampling
        weights corresponding to the sample x_j whose features are
        represented as the columns of self.sampleF.
            logw_j = p_dot(x_j) / q(x_j),
        where p_dot(x_j) is the unnormalized pdf value of the point x_j
        under the current model.
        """
        # First see whether logw has been precomputed
        try:
            return self.logw
        except AttributeError:
            pass
 
        # Compute log w = log [p_dot(s_j)/aux_dist(s_j)]   for
        # j=1,...,n=|sample| using a precomputed matrix of sample
        # features.
        if self.external is None:
            thetadotF = innerprodtranspose(self.sampleF, self.theta)
            logw = thetadotF - self.samplelogprobs
        else:
            e = self.external
            thetadotF = innerprodtranspose(self.externalFs[e], self.theta)
            logw = thetadotF - self.externallogprobs[e]
        
        # Good, we have our logw.  Now:
        self.logw = logw
        return logw
    
    
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
        
        Compute the estimator E_p f_i(X) using logs as:
            num_i / denom,
        where
            num_i = exp(logsumexp(theta.f(s_j) - log aux_dist(s_j) 
                        + log f_i(s_j)))
        and
            denom = [n * Zapprox]
        
        where Zapprox = exp(self.lognormconstapprox()).
        
        We can compute the denominator n*Zapprox directly as:
            exp(logsumexp(log p_dot(s_j) - log aux_dist(s_j)))
          = exp(logsumexp(theta.f(s_j) - log aux_dist(s_j)))
        """
        
        if self.verbose >= 3:
            print "(estimating dual and gradient ...)"

        if not self.usesamplematrix:
            raise NotImplementedError, "sequential estimation is not" \
                    " supported. Either use sample path optimization (set" \
                    " usesamplematrix=True) or stochastic approximation"
    
        # Hereafter is the matrix code

        Es = []
        logZs = []

        for trial in range(self.matrixtrials):
            if self.verbose >= 2 and self.matrixtrials > 1:
                print "(trial " + str(trial) + " ...)"
            
            # Resample if necessary
            if (not self.staticsample) or self.matrixtrials > 1:
                self.resample()
            
            # 1. Compute log w = log [p_dot(s_j)/aux_dist(s_j)]   for
            #    j=1,...,n=|sample| using a precomputed matrix of sample
            #    features.
            logw = self._logw()
            
            # 2. Good, we have our logw.  Now:
            n = len(logw)
            lse = logsumexp(logw)
            logZs.append(lse - math.log(n))
           
            # We don't need to handle negative values separately,
            # because we don't need to take the log of the feature
            # matrix sampleF.
            
            logwminuslogZ = logw - logZs[-1]
            if self.external is None:
                averages = innerprod(self.sampleF, arrayexp(logwminuslogZ)) 
            else:
                averages = innerprod(self.externalFs[self.external], \
                                     arrayexp(logwminuslogZ)) 
            averages /= n
            Es.append(averages)
                        
        # Now we have T=trials vectors of the sample means.  
        # If trials > 1, estimate st dev of means and confidence intervals
        ttrials = len(Es)   # total number of trials performed
        if ttrials == 1:
            #self.mu = Es[0,:]
            self.mu = Es[0]
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
            self.logZsapprox = logZs
            #logstdevZ = 0.5*(-math.log(n-1) + logsumexp([2.*logdiffexp(logZ_k, self.logZapprox) for logZ_k in logZs]))
            stdevlogZ = numpy.array(logZs).std()
            Etemp = numpy.array(Es)
            self.varE = columnvariances(Etemp)
            self.mu = columnmeans(Etemp)
            
            return
  
    
    def entropydualapprox(self):
        """An alias for dualapprox()
        """
        return self.dualapprox()
   
    
    def dualapprox(self):
        """This function approximates the entropy of the model p_theta
        with exponential form exp(theta^T . f_vec(x)) / Z(theta) without
        actually computing p_theta.  This is important if the event space
        is practically innumerable.  We approximate the norm constant Z
        using importance sampling as in [Rosenfeld02whole].  Note that
        the gradient of this estimator is equal to the
        importance-sampling *ratio estimator* of the gradient of the
        entropy dual [see my paper, ICSLP 2004], so using this estimator
        in conjunction with the function gradapprox() in gradient-based
        optimization methods should be stable.
        
        Note that this estimator is deterministic for any given sample.
        
        We use the following estimator:
        
            L_est = log Z_est - sum_i{theta_i K_i}
        
        where
        
        Z_est(theta) = 1/m sum_{x in sample S_0} p_dot(x) / aux_dist(x),
        
        and m = # observations in sample S_0, and K_i = the empirical
        expectation E_p_tilde f_i (X) = sum_x {p(x) f_i(x)}.
        """
        return self.dual()
    
    
    #def _entropydualapproxgradient(self, theta, K):
    #    """Deprecated.  Now an alias for _gradapprox(...).  See my
    #    paper in Interspeech 04 for a proof that these are equal.
    #    """
    #    return self._gradapprox(theta, K)
   
    
       
    def expectationsapprox(self):
        """Estimates the feature expectations E_p f(X) under the current
        model p = p_theta.  If self.staticsample is True, use the current
        feature matrix F.  If self.staticsample is False or
        self.matrixtrials is > 1, draw one or more sample feature
        matrices F afresh using the supplied sampleFgen() generator
        function.
        """
        return self.expectations()        
   
    
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
        self.usesamplematrix = True
        self.logZsapprox = []
        self.Esapprox = []

        assert type(sampler) is types.GeneratorType
        self.sampleFgen = sampler
        self.staticsample = staticsample
        if staticsample:
            self.resample()
                
       
    def pdfapprox(self, fx):
        """Returns the approximate density p_theta(x) at the point x with
        feature statistic fx = f(x).  This is defined as:
            p_theta(x) = exp(theta.f(x)) / Z
        """
        return exp(self.logpdfapprox(fx))
    
    def pdf_function(self):
        """Returns the approximate density p_theta(x) as a function p(fx)
        taking a vector fx = f(x) of feature statistics at any
        point x.  This is defined as:
            p_theta(x) = exp(theta.f(x)) / Z
        """
        logZ = self.lognormconstapprox()
        
        def p(fx):
            return numpy.exp(innerprodtranspose(fx, self.theta) - logZ)
        return p
     
    
    def logpdfapprox(self, fx):
        """Returns the log of the approximate density p_theta(x) at the
        point x.  This is defined as:
            logp(x) = theta.f(x) - log Z
        where f(x) is given by the (m x 1) array fx.

        If, instead, fx is a 2-d (m x n) array, this function interprets
        each of its rows j=0,...,n-1 as a feature fector f(x_j), and
        returns an array containing the log pdf value of each point x_j
        under the current model .
        
        log Z is estimated using the sample provided with
        setsampleFgen().
        """
        logZ = self.lognormconstapprox()
        if len(fx.shape) == 1:
            return numpy.dot(self.theta, fx) - logZ
        else:
            return innerprodtranspose(fx, self.theta) - logZ


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
            k = self.thetalogcounter
            #k = (self.thetalog-1)*self.thetalogfreq
        except:
            k = 0
        
        try:
            a_k = self.a_0
        except AttributeError:
            raise AttributeError, "first define the initial step size a_0"

        avgtheta = self.theta
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
                    #n += (numpy.dot(y_k, y_kminus1) < 0)   # an indicator fn
                    if numpy.dot(y_k, y_kminus1) < 0:
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
                    if self.stepdecreaserate == None:
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
                            - numpy.dot(self.theta, K))
                print "  norm(mu_est - k) = " + str(norm_y_k)
            
            # Update theta (after the convergence tests too ... don't waste the
            # computation.)
            if self.ruppertaverage:
                # Use a simple average of all estimates so far, which
                # Ruppert and Polyak show can converge more rapidly
                newtheta = self.theta - a_k*y_k
                avgtheta = (k-1.0)/k*avgtheta + 1.0/k * newtheta
                if self.verbose:
                    print "  new params[0:5] are: " + str(avgtheta[0:5])
                self.setparams(avgtheta)
            else:
                # Use the standard Robbins-Monro estimator
                self.setparams(self.theta - a_k*y_k)
            
            if k >= self.maxiter:
                print "Reached maximum # iterations during stochastic" \
                        " approximation without convergence."
                break

                      
    def settestsamples(self, F_list, logprob_list, M):
        """Sets the model to test itself every M iterations during
        fitting using the provided list of feature matrices, each
        representing a sample {x_j} from an auxiliary distribution q,
        together with the corresponding log probabiltiy mass or density
        values log {q(x_j)} in logprob_list.  This is useful as an
        external check on the fitting process with sample path
        optimization, which could otherwise reflect the vagaries of the
        single sample being used for optimization, rather than the
        population as a whole.
        """
        # Sanity check 
        assert len(F_list) == len(logprob_list)
        
        self.externalFs = F_list
        self.externallogprobs = logprob_list
        self.testevery = M
        
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

        This should only be called every self.testevery iterations, e.g.
        with:
            if fnevals % M == M - 1:
                self.testexternal()
        """

        # Reduce clutter
        fnevals, K = self.fnevals, self.K
        
        if self.verbose:
            print "  max(theta**2)    = " + str((self.theta**2).max())
            print "Now testing model on external samples ..."

        # Estimate the entropy dual and gradient for each sample
        dualapprox = []
        dualapproxpenalty = []
        mu_ests = []
        for e in xrange(len(self.externalFs)):
            self.external = e
            self.clearcache()
            if self.verbose >= 2:
                print "(testing with sample %d)" % e
            dualapprox.append(self.dual(ignorepenalty=True, ignoretest=True))
            dualapproxpenalty.append(self.dual(ignoretest=True))
            mu_ests.append(self.expectations())
        
        # Reset to using the normal sample matrix sampleF
        self.external = None
        self.clearcache()
        
        meandual = numpy.average(dualapprox)
        meandualpenalty = numpy.average(dualapproxpenalty)
        self.external_duals[fnevals] = dualapprox
        MSEs = [norm(mu - K) for mu in mu_ests]
        self.external_gradnorms[fnevals] = MSEs
                
        if self.verbose:
            print "** Mean dual estimate from the %d external samples"\
                  " is %f" % \
                 (len(self.externalFs), meandual)
            print "** Mean regularized dual estimate from the %d external" \
                  " samples is %f" % \
                 (len(self.externalFs), meandualpenalty)
            print "** Mean mean square error of the feature expectation" \
                    " estimates from the external samples =" \
                    " mean(|| \hat{\mu_e} - k ||) =", numpy.average(MSEs)
        # Track the parameter vector theta with the lowest mean dual estimate so
        # far:
        if meandual < self.bestdual:
            self.bestdual = meandual
            self.bestparams = self.theta
            if self.verbose:
                print "\n\t\t\tStored new minimum entropy dual: %f\n" % meandual


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()


