# maxent.py: Routines for fitting maximum entropy models.

# Copyright: Ed Schofield, 2003-2006
# License: BSD-style (see LICENSE.txt in main source directory)

__author__ = "Ed Schofield"
__version__ = '2.0-alpha4'
__changelog__ = """ 
This module is an adaptation of "ftwmaxent" by Ed Schofield, first posted
on SourceForge as part of the "textmodeller" project in 2002.  The
official repository is now SciPy (since Nov 2005); the SourceForge maxent
code will not be developed further.

------------

Change log:

Since v2.0-alpha3:
(1) Name change ftwmaxent -> scipy/maxent
(2) Modification for inclusion in scipy tree.  Broke one big class into
two smaller classes, one for small models, one for large models.  Here a
'small' sample space is defined as finite and small enough to iterate
over in practice, and a 'large' model is one that requires simulation.
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
from scipy import optimize
from scipy.maxent.maxentutils import logsumexp,sparsefeaturematrix
from scipy.maxent.maxentutils import innerprod,innerprodtranspose,arrayexp


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
        # Stop stoch approx if ||theta_k - theta_{k-1}|| < thetatol:
        self.thetatol = 1e-5         
        # Number of CPUs to use in each computation of the dual and gradient:
        #self.numprocs = numprocs      
        # Number of this cpu (from 0 to numprocs-1)
        #self.myid = myid
        
        self.maxiter = 800
        self.maxfun = 1500
        self.minL = -100
        #self.usecomplexlogs = False
        self.fnevals = 0
        self.gradevals = 0
        self.debug = 0

    def _fit(self, K, func, grad=None, algorithm='CG'):
        """Fit the maxent model q whose feature expectations are given
        by the vector K.

        Model expectations are computed either exactly or using Monte
        Carlo simulation, depending on the 'func' and 'grad' parameters
        passed to this function.
        
        The algorithm can be 'CG', 'BFGS', 'LBFGSB', 'Powell', or
        'Nelder-Mead'.
        
        The CG (conjugate gradients) method is the default; it is quite
        fast and requires only O(N) memory in the number of parameters.
        
        The BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm is a
        variable metric Newton method.  It is perhaps faster than the CG
        method but requires O(N^2) instead of O(N) memory, so it is
        infeasible for more than about 10^3 parameters.

        The Powell algorithm doesn't require gradients.  For small models
        it is slow but robust.  For big models (where func and grad are
        simulated) with large variance in the function estimates, this
        may be less robust than the gradient-based algorithms.
        """
        
        # Sanity checks
        assert type(func) in (types.FunctionType, types.MethodType)
        if grad != None:
            assert type(grad) in (types.FunctionType, types.MethodType)

        # First convert K to a numpy array if necessary
        K = numpy.asarray(K, float)
            
        # Sanity checks
        try:
            self.theta
        except AttributeError:
            raise AttributeError, "first specify a feature matrix or feature "\
                    "functions"
        # Make a copy of the parameters
        oldtheta = numpy.array(self.theta)
            
        assert len(K) == self.numconstraints()
        
        if algorithm == 'CG':
            retval = optimize.fmin_cg(func, oldtheta, \
                                      grad, (K,), \
                                      self.avegtol, maxiter=self.maxiter, \
                                      full_output=1, disp=self.verbose, retall=0)
            
            (newtheta, fopt, func_calls, grad_calls, warnflag) = retval

        elif algorithm == 'LBFGSB':
            retval = optimize.fmin_l_bfgs_b(func, oldtheta, \
                        grad, args=(K,), bounds=self.bounds, pgtol=self.maxgtol,
                        maxfun=self.maxfun)
            (newtheta, fopt, d) = retval
            warnflag, func_calls = d['warnflag'], d['funcalls']
            if self.verbose:
                print algorithm + " optimization terminated successfully."
                print "\tFunction calls: " + str(func_calls)
                # We don't have info on how many gradient calls the LBFGSB
                # algorithm makes

        elif algorithm == 'BFGS':
            retval = optimize.fmin_bfgs(func, oldtheta, \
                                        grad, (K,), self.tol, \
                                        maxiter=self.maxiter, full_output=1, \
                                        disp=self.verbose, retall=0)
            
            (newtheta, fopt, gopt, Lopt, func_calls, grad_calls, warnflag) = retval

        elif algorithm == 'Powell':
            retval = optimize.fmin_powell(func, oldtheta, args=(K,), \
                                   xtol=self.tol, ftol = self.tol, \
                                   maxiter=self.maxiter, full_output=1, \
                                   disp=self.verbose, retall=0)
            
            (newtheta, fopt, direc, numiter, func_calls, warnflag) = retval
        
        elif algorithm == 'Nelder-Mead':
            retval = optimize.fmin(func, oldtheta, args=(K,), \
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
    
  
    def numconstraints(self):
        return len(self.theta)

    def setparams(self, theta, newiter=True):
        """Set the parameter vector theta, replacing the existing
        parameters.  theta must be a list or numpy array of the same
        length as the model's feature vector f.
        """
        # First call the callback function passing the old model, if necessary,
        # before updating it with its new paramaeters.  
        
        # The complex-looking attribute-test sequence is necessary to prevent
        # an AttributeError in the callback function from being caught and
        # silenced inappropriately.

        # If this is a new iteration, call the callbackiter function:
        if newiter:
            try:
                self.callbackiter
            except AttributeError:
                pass
            else:
                self.callbackiter(self)

        # Call the standard callback function, regardless of whether this is a
        # new iteration or a line search:
        try:
            self.callback
        except AttributeError:
            pass
        else:
            self.callback(self)


        # assert len(theta) == self.numconstraints()
        self.theta = numpy.array(theta, float)        # make a copy
        
        # Log the new params to disk
        self.logparams()
            
        # Delete theta-specific stuff
        try:
            del self.mu
        except:
            pass
        #try:
        #    del self.logw
        #except:
        #    pass
        #try:
        #    del self.log_num
        #except:
        #    pass
        try:
            del self.logZ
        except:
            pass
        try:
            del self.logZapprox
        except:
            pass
        
        # Reset the sampler cache, so we begin drawing samples from the
        # stored matrix
        # self.resetsampler = True
        
    def resetparams(self, numfeatures=None):
        """Resets the parameters theta to zero, clearing variables
        dependent on them.
        """

        self.fnevals = 0
        self.gradevals = 0

        if numfeatures:
            m = numfeatures
        else:
            m = self.numconstraints()            
            
        self.theta = numpy.zeros(m, float)
        
        # These bounds on the param values are only effective for the
        # L-BFGS-B optimizer:
        self.bounds = [(-10., 10.)]*len(self.theta)  
        
        # Delete theta-specific stuff
        try:
            del self.mu
        except:
            pass
        #try:
        #    del self.logw
        #except:
        #    pass
        #try:
        #    del self.log_num
        #except:
        #    pass
        try:
            del self.logZ
        except:
            pass
        try:
            del self.logZapprox
        except:
            pass

   
    def setcallbackiter(self, callback):
        """Sets a callback function to be called each new iteration but
        NOT each function / gradient evaluation in line searches, which
        will be wilder.  The function is passed one argument, the current
        model before the parameters are changed for the next iteration.
        """
        self.callbackiter = callback
     
    def setcallback(self, callback):
        """Sets a callback function to be called before the parameters
        are changed with setparams(theta).  The function is passed one
        argument, the current model before the parameters are changed for
        the next fn/grad evaluation.

        Note that line search algorithms in e.g. CG make potentially
        several evals per iter, some of which we expect to be poor.  If
        this is a problem, use callbackiter() instead.
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
        paramsfile = open(self.thetalogfilename + '.' + str(self.thetalogcounter) + '.pickle', 'wb')
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
        try:
            if self.verbose:
                print "Logging to files " + filename + "*"
        except AttributeError:
            pass
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
        
    
    def fit(self, K, algorithm='CG'):
        """Fit the maxent model q whose feature expectations are given
        by the vector K.

        Model expectations are computed exactly, by summing over the
        given sample space.  If the sample space is continuous or too
        large to iterate over, use the 'bigmodel' class instead.
        
        The algorithm can be 'CG', 'BFGS', 'LBFGSB', 'Powell', or
        'Nelder-Mead'.
        
        The CG (conjugate gradients) method is the default; it is quite
        fast and requires only O(N) memory in the number of parameters.
        
        The BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm is a
        variable metric Newton method.  It is perhaps faster than the CG
        method but requires O(N^2) instead of O(N) memory, so it is
        infeasible for more than about 10^3 parameters.

        The Powell algorithm doesn't require gradients.  It is slow but
        robust.
        """
        
        func = self._entropydual
        grad = self._grad
        self._fit(K, func, grad, algorithm=algorithm)

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
    
    
    def expectations(self):
        """The vector E_q[f(X)] under the model p_theta of the vector of
        feature functions f_i over the sample space.
        """
        
        # If discrete, use the representation E_q[f(X)] = q . F
        try:
            self.F

            # A pre-computed matrix of features exists
            q = self.probdistarray()
            innerproduct = innerprod(self.F, q)
            
        except AttributeError:
            raise AttributeError, "need a pre-computed feature matrix F"
            
        return innerproduct
    

    def lognormconst(self):
        """If the sample space is discrete, compute the log of the
        normalization constant Z=sum_{x \in samplespace} exp(theta .
        f(x))

        This assumes it is finite!  If infinite, it will iterate
        endlessly; in this case, use bigmodel.lognormconstapprox() instead!
        """
        
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
        

    def normconst(self):
        """Returns the normalization constant, or partition function, for
        the current model.  Warning -- this may be too large to represent,
        resulting in numerical overflow!  In this case use lognormconst()
        instead.
        """
        return math.exp(self.lognormconst())

    def probdistarray(self):
        """Like the probdist function but returns an array indexed
        by integers instead of a sparse vector indexed by hashable
        objects.  Faster and more robust.  Requires that samplespace
        be stored as a list or an array.
        """
        samplespace = self.samplespace
        
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
            else:
                logsumexp_p_dot = self.logZ
            
            p = arrayexp(log_p_dot - self.logZ)
            
        except AttributeError:
            raise AttributeError, "first store the feature matrix F"
        
        return p
        
    def pdf(self, f=None):
        """Returns the pmf p_theta(x) as a function taking values on the
        model's sample space.  The returned pmf is:
            
            p_theta(x) = exp(theta.f(x) - log Z)
        
        The returned function p_theta also satisfies
            [p(x) for x in self.samplespace] == probdistarray().

        The feature statistic f should be a list of functions
        [f1(),...,fn(x)].  This must be passed unless the model already
        contains an equivalent attribute 'model.f'.
        """
        
        try:
            logZ = self.logZ
        except AttributeError:
            logZ = self.lognormconst()
        
        if f==None:
            f = self.f
        
        def p(x):
            f_x = numpy.array([f[i](x) for i in range(len(f))], float)
            #f_x = numpy.empty(len(f), float)
            #for i in range(len(f)):
            #    f_x[i] = f[i](x)
            return math.exp(numpy.dot(self.theta, f_x) - logZ)
        
        return p

 
    def entropydual(self, K):
        """Computes the Lagrangian dual L(theta) of the entropy of the
        model.  Minimising this function (without constraints) should
        fit the constrained model.  Given by:
            L(theta) = log(Z) - theta^T . K
        """
        
        logZ = self.lognormconst()
        return logZ - numpy.dot(K, self.theta)


    def _entropydual(self, theta, K):
        """Computes the Lagrangian dual of the entropy function L(p)
        with specified Lagriangian multipliers (parameters theta) and
        constraints E_p(f(X)) = K.  
        
        Minimize this to maximize the entropy among all constrained
        models.  The minimum value equals the model entropy and the model
        with the corresponding parameters satisfies the given
        constraints.
        """

        oldtheta = self.theta
        if numpy.any(oldtheta != theta):
            self.setparams(theta)

        # If constraints E_q f(X) = K have been imposed, compute the entropy
        # of the model p_theta with exponential form
        # exp(theta^T . f_vec(x)) / Z(theta) as:
        #        L(p) = log Z(theta) - sum_i { theta_i K_i }
        # where K_i = empirical expectation E_p_tilde f_i (X) = sum_x {p(x)f_i(x)}.
        
        logZ = self.lognormconst()
        L = logZ - numpy.dot(K, self.theta)
        
        # (We don't reset theta to its prior value.)
        
        if self.verbose:
            print "Entropy dual is: " + str(L)
        
        if L < self.minL:
            raise ValueError, "the dual is below the threshold minL and may " \
                    "be diverging to -inf.  Fix the constraints or lower the "\
                    "threshold."

        return L
     
    def _grad(self, theta, K):
        """Returns the gradient vector of the dual of the entropy
        function:
            G(theta) = E_p_theta[f(X)] - K.

        This function is called by the optimization routines.
        """
        
        oldtheta = self.theta
        if numpy.any(oldtheta != theta):
            self.setparams(theta)
        
        negG = self.expectations() - K
        
        # (We don't reset theta to its prior value.)
        
        if self.verbose:
            print "Computed the gradient."
        
        return negG
 




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
    
    def fit(self, K, algorithm='CG'):
        """Fit the maxent model q whose feature expectations are given
        by the vector K.
       
        The algorithm can be 'CG', 'BFGS', 'LBFGSB', 'Powell', or
        'Nelder-Mead'.
        
        The CG (conjugate gradients) method is the default; it is quite
        fast and requires only O(N) memory in the number of parameters.
        
        The BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm is a
        variable metric Newton method.  It is perhaps faster than the CG
        method but requires O(N^2) instead of O(N) memory, so it is
        infeasible for more than about 10^3 parameters.

        The Powell algorithm doesn't require gradients.  It is slow but
        robust.
        """
        
        func = self._entropydualapprox
        grad = self._gradapprox
        self._fit(K, func, grad, algorithm=algorithm)

    def __init__(self):
        super(bigmodel, self).__init__()
        
        # Number of sample matrices to generate and use to estimate E and logZ
        self.matrixtrials = 1      
        
        self.tol = 1e-4
        
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

        # First see whether logZ has been precomputed
        try:
            return self.logZapprox
        except AttributeError:
            pass
 
        # Compute log w = log [p_dot(s_j)/aux_dist(s_j)]   for
        # j=1,...,n=|sample| using a precomputed matrix of sample
        # features.
        thetadotF = innerprodtranspose(self.sampleF, self.theta)
        logw = thetadotF - self.samplelogprobs
        
        # Good, we have our logw.  Now:
        n = len(logw)
        self.logZapprox = logsumexp(logw) - math.log(n)
        return self.logZapprox
    
    
    def normconstapprox(self):
        """Returns the approximate normalization constant Z = E_aux_dist
        [{exp (theta.f(X))} / aux_dist(X)] for computing the distribution
        q_dist etc., using a sample from aux_dist and the law of large
        numbers.  Warning -- this may cause numerical overflow if it is
        too large to represent!  In this case use lognormconstapprox().
        """
        
        return arrayexp(self.lognormconstapprox())
    
        
    def resample(self):
        """(Re)samples the matrix F of sample features.  
        """

        if self.verbose >= 2:
            print "(sampling " + str(n) + " observations ...)"

        (self.sampleF, self.samplelogprobs) = self.sampleFgen.next()

        try:
            self.theta
            m = self.numconstraints()
        except AttributeError:
            (m, n) = self.sampleF.shape
            if self.numsamples:
                assert n == self.numsamples
            else:
                self.numsamples = n
            self.resetparams(m)
        else:
            if (m, self.numsamples) != self.sampleF.shape:
                raise ValueError, "the sample feature generator returned" \
                        " a feature matrix of incorrect dimensions"
        if self.verbose >= 2:
            print "(done.)"

        # Now clear the temporary variables that are no longer correct for this sample
        try:
            del self.mu
        except:
            pass
        try:
            del self.logZapprox
        except:
            pass

    
    def estimate(self):
        """This function approximates both the feature expectation vector
        E_q f(X) and the log of the normalization term Z with importance
        sampling.
        
        It also computes the sample variance of the component estimates
        of the feature expectations as: varE = var(E_1, ..., E_T) where T
        is self.matrixtrials and E_t is the estimate of E_q f(X)
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

        We use [Rosenfeld01Wholesentence]'s estimate of E_q[f_i] as:
            {sum_j  q(s_j)/aux_dist(s_j) f_i(s_j) } 
              / {sum_j q(s_j) / aux_dist(s_j)}.
        
        Note that this is consistent but biased.

        This equals:
            {sum_j  p_dot(s_j)/aux_dist(s_j) f_i(s_j) } 
              / {sum_j p_dot(s_j) / aux_dist(s_j)}
        
        Compute the estimator E_q f_i(X) using logs as:
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
        
        if self.verbose >= 2:
            print "(estimating dual and gradient ...)"

        if not self.usesamplematrix:
            if self.verbose >= 2:
                print "(estimating sequentially ...)"
            (self.mu, self.varapprox, self.logZapprox, n) = \
                    self.onlineestimates()
            return
    
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
            thetadotF = innerprodtranspose(self.sampleF, self.theta)
            logw = thetadotF - self.samplelogprobs
            
            # 2. Good, we have our logw.  Now:
            n = len(logw)
            logZs.append(logsumexp(logw) - math.log(n))

            # The following code is obsolete.  I'll leave it in until
            # I've documented why the above is better:

            # 3. Now compute the numerator self.num, a vector with ith component
            #        num_i = sum_x {p_dot(x)/aux_dist(x) * f_i(x)}.
            #    In matrix notation this is   w . F, where w_j = Q_dot_j / aux_dist_j
            
            # The following is unsafe -- sometimes leads to NaN instead of 
            # an exception.
            #try:
            #    w = arrayexp(logw)
            #    num = innerprodtranspose(self.sampleF, w)
            #
            #    try:
            #        Es.append(num / numpy.sum(w))
            #        continue
            #    except OverflowError:
            #        # A slightly safer version
            #        lognum = robustarraylog(num)
            #        Es.append(arrayexp(lognum - logsumexp(logw)))
            #        continue
                
            #except OverflowError:
            # For large values like x>=800, and for safety with smaller values
            # close to the overflow limit, we'll use our logsumexp function
            # again.  In the past I did this:
            # num_i = exp[logsumexp[log p_dot(x) - log aux_dist(x) + log f_i(x)]]
            # assuming f_i(x) is positive for all x, or using special handling
            # otherwise (like the complex log if negative).
                 
            #if self.usecomplexlogs:   # a little flag to switch on the 
            #                          # complex log routines
            #    if self.verbose:
            #        print "Using complex logarithm routines ..."
            #        
            #    lognum = emptyarray(m, numpy.Complex64)

            #    # Compute lognum_i as:
            #    #     logsumexp(log p_dot(x) - log aux_dist(x) + log f_i(x))
            #    #   = logsumexp(logw(x) + log f_i(x))

            #    for i in range(m):   # loop over features
            #        # Let f_i be the ith column of the sparse matrix sampleF
            #        # f_i = self.sampleF[:,i] works -- but we want to be able to iterate over this
            #        f_i = column(self.sampleF, i)
            #        
            #        # Was:
            #        # log_f_i = fftwutils.robustarraylog(f_i)
            #        # lognum[i] = fftwutils.logsumexpcomplex(logw + log_f_i)
            #        # Replaced by this SLOWER version to assist profiling:
            #        log_f_i = robustarraylog(f_i)
            #        lognum[i] = logsumexpcomplex_wrapper(logw + log_f_i)
            #    # -----------
            #    # We have computed lognum successfully.
            #    
            #    # Could do this:
            #    #logsum_w = fftwutils.logsumexp(logw)
            #    # But this is faster, since we've already computed logZ = logsumexp(logw) - log(n)
            #    logsum_w = logZs[-1] + math.log(n)
            #    logE = lognum - logsum_w
            #    #Es[trial,:] = arrayexp(logE)
            #    Es.append(arrayexp(logE))

            #else:

            # We don't need to handle negative values separately,
            # because we don't need to take the log of the feature
            # matrix sampleF.
            
            logwminuslogZ = logw - logZs[-1]
            averages = innerprod(self.sampleF, arrayexp(logwminuslogZ)) 
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
            stdevlogZ = var(logZs)**0.5
            Etemp = numpy.array(Es)
            self.varE = columnvariances(Etemp)
            self.mu = columnmeans(Etemp)
            
            return
  
    
    def _entropydualapprox(self, theta, K, newiter=True):
        """This function is not designed to be called by a user, but by
        an optimization routine.  (The parameters theta are stored in the
        model object anyway.)

        If newiter is False, this evaluation corresponds to a line search
        routine, so we should exclude the resulting dual evaluation from
        plots.
        """
        self.fnevals += 1
        
        if self.verbose:
            print "Function eval #" + str(self.fnevals)
            if self.verbose >= 2 and not newiter:
                print "(line search)"

        # Has theta changed?  If so, set it anew, clearing the cache
        # variables that are functions of theta.
        if numpy.any(self.theta != theta):
            self.setparams(theta, newiter)
        
        L = self.entropydualapprox(K)
        
        if L < self.minL:
            raise DivergenceError, \
              "the dual is below the threshold minL and may be diverging" \
              " to -inf.  Fix the constraints or lower the"\
              " threshold."
        
        if self.verbose:
            print "Approx dual is " + str(L)

        # (We don't reset theta to its prior value.)
        return L
    
    
    def entropydualapprox(self, K):
        """This function approximates the entropy of the model p_theta
        with exponential form exp(theta^T . f_vec(x)) / Z(theta) without
        actually computing p_theta.  This is important if the event space
        is practically innumerable.  We approximate the norm constant Z
        using importance sampling as in [Rosenfeld02whole].  Note that
        the gradient of this estimator is equal to the
        importance-sampling estimator of the gradient E_q f(X) of the
        entropy [see my paper, ICSLP 2004], so using this estimator in
        conjunction with the function E_q_IS in gradient-based
        optimization methods should be stable.
        
        Note that this estimator is deterministic for any given sample.
        
        We use the following estimator:
        
            L_est = log Z_est - sum_i{theta_i K_i}
        
        where
        
        Z_est(theta) = 1/m sum_{x in sample S_0} p_dot(x) / aux_dist(x),
        
        and m = # observations in sample S_0, and K_i = the empirical
        expectation E_p_tilde f_i (X) = sum_x {p(x) f_i(x)}.
        """
        
        L = self.lognormconstapprox() - numpy.dot(K, self.theta)

        return L    
    
    
    #def _entropydualapproxgradient(self, theta, K):
    #    """Deprecated.  Now an alias for _gradapprox(...).  See my
    #    paper in Interspeech 04 for a proof that these are equal.
    #    """
    #    return self._gradapprox(theta, K)
   
    
    def _gradapprox(self, theta, K, newiter=True):
        """Estimate the gradient of the entropy dual.
        """
        
        self.gradevals += 1
        
        if self.verbose:
            print "Grad eval #" + str(self.gradevals)
            #if not newiter:
            #    print "(line search)"

        # Has theta changed?  If so, clear the speed-enhancing temporary
        # variables that are functions of theta.
        if numpy.any(self.theta != theta):
            self.setparams(theta, newiter)
        
        G = self.expectationsapprox() - K
        
        # (We don't reset theta to its prior value.)
        
        return G
        
        
    def expectationsapprox(self):
        """Estimates the feature expectations E_p f(X) under the current
        model p = p_theta.  If 'staticsample' is True, use the current
        feature matrix F.  If 'staticsample' is False or 'matrixtrials'
        is > 1, draw one or more sample feature matrices F afresh using
        the supplied sampleFgen() generator function.
        """
        
        try:
            # See if already computed
            self.mu
            #if self.verbose >= 3:
            #    print "(returning pre-computed expectations)"
            return self.mu
        
        except AttributeError:
            self.estimate()
            return self.mu
     
   
    def logprobapprox(self, x):
        """Returns the approximate log prob of the observation x under
        the model.  This assumes the model has the attribute f for
        feature functions and that the model has already been fit
        using one of the fitapprox* routines.
        
        We compute this as theta . f(x) - self.logZapprox
        """
        m, n = self.sampleF.shape
        assert m == len(self.theta)
        
        # 1. First compute log(p_dot(x)) = theta . f(x)
        
        T = numpy.asarray(self.theta)
        fx = sparsefeatures(self.f, x)
        log_p_dot = dotprod(fx, T)
        
        # 2. Now compute log q(x) = log(p_dot(x)/Z) = log p_dot(x) - log(Z)
        try:
            return log_p_dot - self.logZapprox
        except AttributeError:
            return log_p_dot - self.lognormconstapprox()
     
    
    def setsampleFgen(self, sampler, sequential=False, matrixsize=None, sparse=True, staticsample=True):  
        """Initializes the Monte Carlo sampler to use the supplied
        generator of samples' features and log probabilities.  This is an
        alternative to defining a sampler in terms of a (fixed size)
        feature matrix sampleF and accompanying vector samplelogprobs of
        log probabilities.
            
        Calling sampler() should generate tuples (f(x), log_aux_dist_x)
        where x is an observation sampled from the auxiliary distribution
        aux_dist, f(x) is a vector of features of x, and log_aux_dist_x
        is its (natural) log probability.

        If matrixtrials > 1 and staticsample = True, (which is useful for
        estimating variance between the different feature estimates),
        sampler.next() will be called once for each trial
        (0,...,matrixtrials) for each iteration.  This allows using a set
        of feature matrices, each of which stays constant over all
        iterations.
               
        We now insist that sampleFgen.next() return batches of size
        matrixsize to avoid overhead in extra function calls and memory
        copying (and extra code).
        
        An alternative was to supply a list of samplers,
        sampler=[sampler0, sampler1, ..., sampler_{m-1}, samplerZ], one
        for each feature and one for estimating the normalization
        constant Z. But this code was unmaintained, and has now been
        removed (but it's in CVS :).

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
        
        if not sequential:
            self.usesamplematrix = True
            self.numsamples = matrixsize
           
            self.logZsapprox = []
            self.Esapprox = []

        # First check whether we have a separate sampling dist for each feature
        if type(sampler) is list:
            raise NotImplementedError, "code for using multiple samplers is" \
                " out of date"
            # We require one sampler for each feature and one for the 
            # normalization term Z:
            #assert len(sampler) == 1+self.numconstraints()
        else:
            assert type(sampler) is types.GeneratorType
        self.sampleFgen = sampler
        self.staticsample = staticsample
        if staticsample:
            self.resample()
                
       
    def pdfapprox(self):
        """Returns the approximate density p_theta(x) as a function on
        the model's sample space.  The function p_theta is defined as:
            p_theta(x) = exp(theta.f(x) - log Z)
        """
        
        try:
            logZ = self.logZapprox
        except AttributeError:
            logZ = self.lognormconstapprox()
        
        f = self.f
        def p(x):
            f_x = numpy.empty(len(f), float)
            for i in range(len(f)):
                f_x[i] = f[i](x)
            return math.exp(numpy.dot(self.theta, f_x) - logZ)
        
        return p


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
                        # Store iterations of sign switches (for plotting purposes)
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
                    # I think we need a stepsize decreasing as n^-1 for almost sure convergence
                    a_k = 1.0 * self.a_0 / (n ** self.stepdecreaserate)
            # otherwise leave unchanged
           
            self.matrixtrials = 1
            self.staticsample = False
            if self.andradottir:    # use Andradottir (1996)'s scaling?
                self.estimate()   # resample and reestimate
                y_k_1 = self.mu - K
                self.estimate()   # resample and reestimate
                y_k_2 = self.mu - K
                y_k = y_k_1 / max(1.0, norm(y_k_2)) + y_k_2 / max(1.0, norm(y_k_1))
            else:
                # Standard Robbins-Monro estimator
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
                if self.verbose >= 2:
                    print "  params[0:5] are: " + str(self.theta[0:5])
                    print "  norm(mu - k) = " + str(norm_y_k)
            #if self.storenormerrors:
            #    self.normerrors.append(norm_y_k)
            #if self.storeduals:
            #    self.duals.append(self.logZapprox - numpy.dot(self.theta, K))
               
            # Update theta (after the convergence tests too ... don't waste the
            # computation.)
            if self.ruppertaverage:
                # Use a simple average of all estimates so far, which
                # Ruppert and Polyak show can converge more rapidly
                newtheta = self.theta - a_k*y_k
                avgtheta = (k-1.0)/k*avgtheta + 1.0/k * newtheta
                if self.verbose >= 2:
                    print "  new params[0:5] are: " + str(avgtheta[0:5])
                self.setparams(avgtheta)
            else:
                # Use the standard Robbins-Monro estimator
                self.setparams(self.theta - a_k*y_k)
            
            if k >= self.maxiter:
                if self.verbose:
                    print "Reached maximum # iterations during stochastic" \
                            " approximation without convergence."
                break

                      

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()


