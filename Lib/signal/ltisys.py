from filter_design import tf2zpk, zpk2tf, normalize
from Numeric import product, zeros, asarray, concatenate, \
     array, dot, transpose
import Numeric
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import scipy.linalg as linalg
from scipy import r_, c_, eye
from scipy import r1array, r2array
from scipy import poly, squeeze, Mat

def tf2ss(num, den):
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
    A, B, C, D = map(none_to_empty, (A, B, C, D))
    A, B, C, D = map(r2array, (A, B, C, D))

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
        Ck = r2array(C[k,:])
        num[k] = poly(A - dot(B,Ck)) + (D[k]-1)*den

    return num, den

def zpk2ss(z,p,k):
    return tf2ss(*zpk2tf(z,p,k))

def ss2zpk(A,B,C,D,input=0):
    return tf2zpk(*ss2tf(A,B,C,D,input=input))

class lti:
    def __init__(self,*args,**kwords):
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


def lsim(system, U, T, X0=None):
    # system is an lti system or a sequence
    #  with 2 (num, den)
    #       3 (zeros, poles, gain)
    #       4 (A, B, C, D)
    #  describing the system
    #  U is an input vector at times T
    #   if system describes multiple outputs
    #   then U can be a rank-2 array with the number of columns
    #   being the number of inputs
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    U = r1array(U)
    T = r1array(T)
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

    ufunc = interpolate.linear_1d(T, U, axis=0, bounds_error=0, fill_value=0)

    def fprime(x, t, sys, ufunc):
        return dot(sys.A,x) + squeeze(dot(sys.B,ufunc([t])))

    xout = integrate.odeint(fprime, X0, T, args=(sys, ufunc))
    yout = dot(sys.C,transpose(xout)) + dot(sys.D,transpose(U))
    return squeeze(transpose(yout)), T, xout


def impulse(system, U, T, X0=None):
    if isinstance(system, lti):
        sys = system
    else:
        sys = lti(*system)
    eA,B,C = map(Mat,(linalg.expm(sys.A), sys.B, sys.C))
    

def step(system, U, T, X0=None):
    pass
