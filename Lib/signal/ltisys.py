#
# Author: Travis Oliphant 2001
#

from filter_design import tf2zpk, zpk2tf, normalize
from Numeric import product, zeros, concatenate, \
     array, dot, transpose, arange, ones, Float
import Numeric
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import scipy.linalg as linalg
from scipy_base import r_, c_, eye, real, atleast_1d, atleast_2d, poly, \
     squeeze, diag, asarray
from Matrix import Matrix as Mat

def tf2ss(num, den):
    """Transfer function to state-space representation.

    Inputs:

      num, den -- sequences representing the numerator and denominator polynomials.

    Outputs:

      A, B, C, D -- state space representation of the system.
    """
    # Controller canonical state-space representation.
    #  if M+1 = len(num) and K+1 = len(den) then we must have M <= K 
    #  states are found by asserting that X(s) = U(s) / D(s)
    #  then Y(s) = N(s) * X(s) 
    #
    #   A, B, C, and D follow quite naturally.
    #
    num, den = normalize(num, den)   # Strips zeros, checks arrays
    nn = len(num.shape)
    if nn == 1:
        num = asarray([num],num.typecode())
    M = num.shape[1] 
    K = len(den) 
    if (M > K):
        raise ValueError, "Improper transfer function."
    if (M == 0 or K == 0):  # Null system
        return array([],Float), array([], Float), array([], Float), \
               array([], Float)

    # pad numerator to have same number of columns has denominator
    num = c_[zeros((num.shape[0],K-M),num.typecode()), num]

    if num.shape[-1] > 0:
        D = num[:,0]
    else:
        D = array([],Float)

    if K == 1:
        return array([], Float), array([], Float), array([], Float), D
    
    frow = -array([den[1:]])
    A = r_[frow, eye(K-2, K-1)]
    B = eye(K-1, 1)
    C = num[:,1:] - num[:,0] * den[1:]
    return A, B, C, D

def none_to_empty(arg):
    if arg is None:
        return []
    else:
        return arg

def abcd_normalize(A=None, B=None, C=None, D=None):
    """Check state-space matrices and ensure they are rank-2.
    """
    A, B, C, D = map(none_to_empty, (A, B, C, D))
    A, B, C, D = map(atleast_2d, (A, B, C, D))

    if ((len(A.shape) > 2) or (len(B.shape) > 2) or \
        (len(C.shape) > 2) or (len(D.shape) > 2)):
        raise ValueError, "A, B, C, D arrays can be no larger than rank-2."

    MA, NA = A.shape
    MB, NB = B.shape
    MC, NC = C.shape
    MD, ND = D.shape

    if (MC == 0) and (NC == 0) and (MD != 0) and (NA != 0):
        MC, NC = MD, NA
        C = zeros((MC, NC))
    if (MB == 0) and (NB == 0) and (MA != 0) and (ND != 0):
        MB, NB = MA, ND
        B = zeros(MB, NB)
    if (MD == 0) and (ND == 0) and (MC != 0) and (NB != 0):
        MD, ND = MC, NB
        D = zeros(MD, ND)
    if (MA == 0) and (NA == 0) and (MB != 0) and (NC != 0):
        MA, NA = MB, NC
        A = zeros(MA, NA)

    if MA != NA:
        raise ValueError, "A must be square."
    if MA != MB:
        raise ValueError, "A and B must have the same number of rows."
    if NA != NC:
        raise ValueError, "A and C must have the same number of columns."
    if MD != MC:
        raise ValueError, "C and D must have the same number of rows."
    if ND != NB:
        raise ValueErrro, "B and D must have the same number of columns."

    return A, B, C, D
    
def ss2tf(A, B, C, D, input=0):
    """State-space to transfer function.

    Inputs:

      A, B, C, D -- state-space representation of linear system.
      input -- For multiple-input systems, the input to use.

    Outputs:

      num, den -- Numerator and denominator polynomials (as sequences)
                  respectively.
    """
    # transfer function is C (sI - A)**(-1) B + D
    A, B, C, D = map(asarray, (A, B, C, D))
    # Check consistency and 
    #     make them all rank-2 arrays
    A, B, C, D = abcd_normalize(A, B, C, D)

    nout, nin = D.shape
    if input >= nin:
        raise ValueError, "System does not have the input specified."

    # make MOSI from possibly MOMI system.
    if B.shape[-1] != 0:
        B = B[:,input]
    B.shape = (B.shape[0],1)
    if D.shape[-1] != 0:
        D = D[:,input]

    den = poly(A)

    if (product(B.shape) == 0) and (product(C.shape) == 0):
        num = Numeric.ravel(D)
        if (product(D.shape) == 0) and (product(A.shape) == 0):
            den = []
        end
        return num, den

    num_states = A.shape[0] 
    type_test = A[:,0] + B[:,0] + C[0,:] + D
    num = Numeric.zeros((nout, num_states+1),type_test.typecode())
    for k in range(nout):
        Ck = atleast_2d(C[k,:])
        num[k] = poly(A - dot(B,Ck)) + (D[k]-1)*den

    return num, den

def zpk2ss(z,p,k):
    """Zero-pole-gain representation to state-space representation

    Inputs:

      z, p, k -- zeros, poles (sequences), and gain of system

    Outputs:

      A, B, C, D -- state-space matrices.
    """
    return tf2ss(*zpk2tf(z,p,k))

def ss2zpk(A,B,C,D,input=0):
    """State-space representation to zero-pole-gain representation.

    Inputs:

      A, B, C, D -- state-space matrices.
      input -- for multiple-input systems, the input to use.

    Outputs:

      z, p, k -- zeros and poles in sequences and gain constant.
    """    
    return tf2zpk(*ss2tf(A,B,C,D,input=input))

class lti:
    """Linear Time Invariant class which simplifies representation.
    """
    def __init__(self,*args,**kwords):
        """Initialize the LTI system using either:
           (numerator, denominator)
           (zeros, poles, gain)
           (A, B, C, D) -- state-space.
        """
        N = len(args)
        if N == 2:  # Numerator denominator transfer function input
            self.__dict__['num'], self.__dict__['den'] = normalize(*args)
            self.__dict__['zeros'], self.__dict__['poles'], \
            self.__dict__['gain'] = tf2zpk(*args)
            self.__dict__['A'], self.__dict__['B'], \
                                self.__dict__['C'], \
                                self.__dict__['D'] = tf2ss(*args)
            self.inputs = 1
            if len(self.num.shape) > 1:
                self.outputs = self.num.shape[0]
            else:
                self.outputs = 1
        elif N == 3:      # Zero-pole-gain form
            self.__dict__['zeros'], self.__dict__['poles'], \
                                    self.__dict__['gain'] = args
            self.__dict__['num'], self.__dict__['den'] = zpk2tf(*args)
            self.__dict__['A'], self.__dict__['B'], \
                                self.__dict__['C'], \
                                self.__dict__['D'] = zpk2ss(*args)
            self.inputs = 1
            if len(self.zeros.shape) > 1:
                self.outputs = self.zeros.shape[0]
            else:
                self.outputs = 1
        elif N == 4:       # State-space form
            self.__dict__['A'], self.__dict__['B'], \
                                self.__dict__['C'], \
                                self.__dict__['D'] = abcd_normalize(*args)
            self.__dict__['zeros'], self.__dict__['poles'], \
                                    self.__dict__['gain'] = ss2zpk(*args)
            self.__dict__['num'], self.__dict__['den'] = ss2tf(*args)
            self.inputs = B.shape[-1]
            self.outputs = C.shape[0]
        else:
            raise ValueError, "Needs 2, 3, or 4 arguments."

    def __setattr__(self, attr, val):
        if attr in ['num','den']:
            self.__dict__[attr] = val
            self.__dict__['zeros'], self.__dict__['poles'], \
                                    self.__dict__['gain'] = \
                                    tf2zpk(self.num, self.den)
            self.__dict__['A'], self.__dict__['B'], \
                                self.__dict__['C'], \
                                self.__dict__['D'] = \
                                tf2ss(self.num, self.den)
        elif attr in ['zeros', 'poles', 'gain']:
            self.__dict__[attr] = val
            self.__dict__['num'], self.__dict__['den'] = \
                                  zpk2tf(self.zeros,
                                         self.poles, self.gain)
            self.__dict__['A'], self.__dict__['B'], \
                                self.__dict__['C'], \
                                self.__dict__['D'] = \
                                zpk2ss(self.zeros,
                                       self.poles, self.gain)
        elif attr in ['A', 'B', 'C', 'D']:
            self.__dict__[attr] = val
            self.__dict__['zeros'], self.__dict__['poles'], \
                                    self.__dict__['gain'] = \
                                    ss2zpk(self.A, self.B,
                                           self.C, self.D)
            self.__dict__['num'], self.__dict__['den'] = \
                                  ss2tf(self.A, self.B,
                                        self.C, self.D)            
        else:
            self.__dict__[attr] = val

    def impulse(self, X0=None, T=None, N=None):
        return impulse(self, X0=X0, T=T, N=N)

    def step(self, X0=None, T=None, N=None):
        return step(self, X0=X0, T=T, N=N)

    def output(self, U, T, X0=None):
        return lsim(self, U, T, X0=X0)


def lsim2(system, U, T, X0=None):
    """Simulate output of a continuous-time linear system, using ODE solver.

    Inputs:

      system -- an instance of the LTI class or a tuple describing the
                system.  The following gives the number of elements in
                the tuple and the interpretation.
                  2 (num, den)
                  3 (zeros, poles, gain)
                  4 (A, B, C, D)
      U -- an input array describing the input at each time T
           (linear interpolation is assumed between given times).
           If there are multiple inputs, then each column of the 
           rank-2 array represents an input.
      T -- the time steps at which the input is defined and at which
           the output is desired.
      X0 -- (optional, default=0) the initial conditions on the state vector.

    Outputs: (T, yout, xout)

      T -- the time values for the output.
      yout -- the response of the system.
      xout -- the time-evolution of the state-vector.
    """
    # system is an lti system or a sequence
    #  with 2 (num, den)
    #       3 (zeros, poles, gain)
    #       4 (A, B, C, D)
    #  describing the system
    #  U is an input vector at times T
    #   if system describes multiple outputs
    #   then U can be a rank-2 array with the number of columns
    #   being the number of inputs

    # rather than use lsim, use direct integration and matrix-exponential.
    if isinstance(system, lti):        
        sys = system
    else:
        sys = lti(*system)
    U = atleast_1d(U)
    T = atleast_1d(T)
    if len(U.shape) == 1:
        U.shape = (U.shape[0],1)
    sU = U.shape        
    if len(T.shape) != 1:
        raise ValueError, "T must be a rank-1 array."
    if sU[0] != len(T):
        raise ValueError, "U must have the same number of rows as elements in T."
    if sU[1] != sys.inputs:
        raise ValueError, "System does not define that many inputs."

    if X0 is None:
        X0 = zeros(sys.B.shape[0],sys.A.typecode())

    # for each output point directly integrate assume zero-order hold
    #   or linear interpolation.

    ufunc = interpolate.linear_1d(T, U, axis=0, bounds_error=0, fill_value=0)

    def fprime(x, t, sys, ufunc):
        return dot(sys.A,x) + squeeze(dot(sys.B,ufunc([t])))

    xout = integrate.odeint(fprime, X0, T, args=(sys, ufunc))
    yout = dot(sys.C,transpose(xout)) + dot(sys.D,transpose(U))
    return T, squeeze(transpose(yout)), xout


def lsim(system, U, T, X0=None, interp=1):
    """Simulate output of a continuous-time linear system.

    Inputs:

      system -- an instance of the LTI class or a tuple describing the
                system.  The following gives the number of elements in
                the tuple and the interpretation.
                  2 (num, den)
                  3 (zeros, poles, gain)
                  4 (A, B, C, D)
      U -- an input array describing the input at each time T
           (interpolation is assumed between given times).
           If there are multiple inputs, then each column of the 
           rank-2 array represents an input.
      T -- the time steps at which the input is defined and at which
           the output is desired.
      X0 -- (optional, default=0) the initial conditions on the state vector.
      interp -- linear (1) or zero-order hold (0) interpolation

    Outputs: (T, yout, xout)

      T -- the time values for the output.
      yout -- the response of the system.
      xout -- the time-evolution of the state-vector.
    """
    # system is an lti system or a sequence
    #  with 2 (num, den)
    #       3 (zeros, poles, gain)
    #       4 (A, B, C, D)
    #  describing the system
    #  U is an input vector at times T
    #   if system describes multiple inputs
    #   then U can be a rank-2 array with the number of columns
    #   being the number of inputs
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    U = atleast_1d(U)
    T = atleast_1d(T)
    if len(U.shape) == 1:
        U.shape = (U.shape[0],1)
    sU = U.shape        
    if len(T.shape) != 1:
        raise ValueError, "T must be a rank-1 array."
    if sU[0] != len(T):
        raise ValueError, "U must have the same number of rows as elements in T."
    if sU[1] != sys.inputs:
        raise ValueError, "System does not define that many inputs."

    if X0 is None:
        X0 = zeros(sys.B.shape[0],sys.A.typecode())

    xout = zeros((len(T),sys.B.shape[0]),sys.A.typecode())
    xout[0] = X0
    A = sys.A
    AT, BT = transpose(sys.A), transpose(sys.B)
    dt = T[1]-T[0]
    lam, v = linalg.eig(A)
    vt = transpose(v)
    vti = linalg.inv(vt)
    GT = dot(dot(vti,diag(Numeric.exp(dt*lam))),vt).astype(xout.typecode())
    ATm1 = linalg.inv(AT)
    ATm2 = dot(ATm1,ATm1)
    I = eye(A.shape[0],typecode=A.typecode())
    GTmI = GT-I
    F1T = dot(dot(BT,GTmI),ATm1)
    if interp:
        F2T = dot(BT,dot(GTmI,ATm2)/dt - ATm1)

    for k in xrange(1,len(T)):
        dt1 = T[k] - T[k-1]
        if dt1 != dt:
            dt = dt1
            GT = dot(dot(vti,diag(Numeric.exp(dt*lam))),vt).astype(xout.typecode())
            GTmI = GT-I
            F1T = dot(dot(BT,GTmI),ATm1)
            if interp:
                F2T = dot(BT,dot(GTmI,ATm2)/dt - ATm1)

        xout[k] = dot(xout[k-1],GT) + dot(U[k-1],F1T)
        if interp:
            xout[k] = xout[k] + dot((U[k]-U[k-1]),F2T)

    yout = squeeze(dot(U,transpose(sys.D))) + squeeze(dot(xout,transpose(sys.C)))
    return T, squeeze(yout), squeeze(xout)


def impulse(system, X0=None, T=None, N=None):
    """Impulse response of continuous-time system.

    Inputs:

      system -- an instance of the LTI class or a tuple with 2, 3, or 4
                elements representing (num, den), (zero, pole, gain), or
                (A, B, C, D) representation of the system.
      X0 -- (optional, default = 0) inital state-vector.
      T -- (optional) time points (autocomputed if not given).
      N -- (optional) number of time points to autocompute (100 if not given).

    Ouptuts: (T, yout)

      T -- output time points,
      yout -- impulse response of system (except possible singularities at 0).
    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    if X0 is None:
        B = sys.B
    else:
        B = sys.B + X0
    if N is None:
        N = 100
    if T is None:
        vals = linalg.eigvals(sys.A)
        tc = 1.0/min(abs(real(vals)))
        T = arange(0,7*tc,7*tc / float(N))
    h = zeros(T.shape, sys.A.typecode())
    s,v = linalg.eig(sys.A)
    vi = linalg.inv(v)
    C = sys.C
    for k in range(len(h)):
        es = diag(Numeric.exp(s*T[k]))
        eA = (dot(dot(v,es),vi)).astype(h.typecode())
        h[k] = squeeze(dot(dot(C,eA),B))
    return T, h

def step(system, X0=None, T=None, N=None):
    """Step response of continuous-time system.

    Inputs:

      system -- an instance of the LTI class or a tuple with 2, 3, or 4
                elements representing (num, den), (zero, pole, gain), or
                (A, B, C, D) representation of the system.
      X0 -- (optional, default = 0) inital state-vector.
      T -- (optional) time points (autocomputed if not given).
      N -- (optional) number of time points to autocompute (100 if not given).

    Ouptuts: (T, yout)

      T -- output time points,
      yout -- step response of system.
    """
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    if N is None:
        N = 100
    if T is None:
        vals = linalg.eigvals(sys.A)
        tc = 1.0/min(abs(real(vals)))
        T = arange(0,7*tc,7*tc / float(N))
    U = ones(T.shape, sys.A.typecode())
    vals = lsim(sys, U, T, X0=X0)
    return vals[0], vals[1]
    
