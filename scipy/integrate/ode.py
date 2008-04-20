# Authors: Pearu Peterson, Pauli Virtanen
"""
First-order ODE integrators

User-friendly interface to various numerical integrators for solving a
system of first order ODEs with prescribed initial conditions::

    d y(t)[i]
    ---------  = f(t,y(t))[i],
       d t

    y(t=0)[i] = y0[i],

where::

    i = 0, ..., len(y0) - 1

class ode
---------

A generic interface class to numeric integrators. It has the following
methods::

    integrator = ode(f,jac=None)
    integrator = integrator.set_integrator(name,**params)
    integrator = integrator.set_initial_value(y0,t0=0.0)
    integrator = integrator.set_f_params(*args)
    integrator = integrator.set_jac_params(*args)
    y1 = integrator.integrate(t1,step=0,relax=0)
    flag = integrator.successful()

"""

integrator_info = \
"""
Available integrators
---------------------

vode
~~~~

Real-valued Variable-coefficient Ordinary Differential Equation
solver, with fixed-leading-coefficient implementation. It provides
implicit Adams method (for non-stiff problems) and a method based on
backward differentiation formulas (BDF) (for stiff problems).

Source: http://www.netlib.org/ode/vode.f

This integrator accepts the following parameters in set_integrator()
method of the ode class:

- atol : float or sequence
  absolute tolerance for solution
- rtol : float or sequence
  relative tolerance for solution
- lband : None or int
- rband : None or int
  Jacobian band width, jac[i,j] != 0 for i-lband <= j <= i+rband.
  Setting these requires your jac routine to return the jacobian
  in packed format, jac_packed[i-j+lband, j] = jac[i,j].
- method: 'adams' or 'bdf'
  Which solver to use, Adams (non-stiff) or BDF (stiff)
- with_jacobian : bool
  Whether to use the jacobian
- nsteps : int
  Maximum number of (internally defined) steps allowed during one
  call to the solver.
- first_step : float
- min_step : float
- max_step : float
  Limits for the step sizes used by the integrator.
- order : int
  Maximum order used by the integrator,
  order <= 12 for Adams, <= 5 for BDF.

zvode
~~~~~

Complex-valued Variable-coefficient Ordinary Differential Equation
solver, with fixed-leading-coefficient implementation.  It provides
implicit Adams method (for non-stiff problems) and a method based on
backward differentiation formulas (BDF) (for stiff problems).

Source: http://www.netlib.org/ode/zvode.f

This integrator accepts the same parameters in set_integrator()
as the "vode" solver.

:Note:
    When using ZVODE for a stiff system, it should only be used for
    the case in which the function f is analytic, that is, when each f(i)
    is an analytic function of each y(j).  Analyticity means that the
    partial derivative df(i)/dy(j) is a unique complex number, and this
    fact is critical in the way ZVODE solves the dense or banded linear
    systems that arise in the stiff case.  For a complex stiff ODE system
    in which f is not analytic, ZVODE is likely to have convergence
    failures, and for this problem one should instead use DVODE on the
    equivalent real system (in the real and imaginary parts of y).

"""

__doc__ += integrator_info

# XXX: Integrators must have:
# ===========================
# cvode - C version of vode and vodpk with many improvements.
#   Get it from http://www.netlib.org/ode/cvode.tar.gz
#   To wrap cvode to Python, one must write extension module by
#   hand. Its interface is too much 'advanced C' that using f2py
#   would be too complicated (or impossible).
#
# How to define a new integrator:
# ===============================
#
# class myodeint(IntegratorBase):
#
#     runner = <odeint function> or None
#
#     def __init__(self,...):                           # required
#         <initialize>
#
#     def reset(self,n,has_jac):                        # optional
#         # n - the size of the problem (number of equations)
#         # has_jac - whether user has supplied its own routine for Jacobian
#         <allocate memory,initialize further>
#
#     def run(self,f,jac,y0,t0,t1,f_params,jac_params): # required
#         # this method is called to integrate from t=t0 to t=t1
#         # with initial condition y0. f and jac are user-supplied functions
#         # that define the problem. f_params,jac_params are additional
#         # arguments
#         # to these functions.
#         <calculate y1>
#         if <calculation was unsuccesful>:
#             self.success = 0
#         return t1,y1
#
#     # In addition, one can define step() and run_relax() methods (they
#     # take the same arguments as run()) if the integrator can support
#     # these features (see IntegratorBase doc strings).
#
# if myodeint.runner:
#     IntegratorBase.integrator_classes.append(myodeint)

__all__ = ['ode']
__version__ = "$Id$"
__docformat__ = "restructuredtext en"

from numpy import asarray, array, zeros, sin, int32, isscalar
import re, sys

#------------------------------------------------------------------------------
# User interface
#------------------------------------------------------------------------------

class ode(object):
    """\
A generic interface class to numeric integrators.

See also
--------
odeint : an integrator with a simpler interface based on lsoda from ODEPACK
quad : for finding the area under a curve

Examples
--------
A problem to integrate and the corresponding jacobian:

>>> from scipy import eye
>>> from scipy.integrate import ode
>>>
>>> y0, t0 = [1.0j, 2.0], 0
>>>
>>> def f(t, y, arg1):
>>>     return [1j*arg1*y[0] + y[1], -arg1*y[1]**2]
>>> def jac(t, y, arg1):
>>>     return [[1j*arg1, 1], [0, -arg1*2*y[1]]]

The integration:

>>> r = ode(f, jac).set_integrator('zvode', method='bdf', with_jacobian=True)
>>> r.set_initial_value(y0, t0).set_f_params(2.0).set_jac_params(2.0)
>>> t1 = 10
>>> dt = 1
>>> while r.successful() and r.t < t1:
>>>     r.integrate(r.t+dt)
>>>     print r.t, r.y

"""

    __doc__ += integrator_info

    def __init__(self, f, jac=None):
        """
        Define equation y' = f(y,t) where (optional) jac = df/dy.

        Parameters
        ----------
        f : f(t, y, *f_args)
            Rhs of the equation. t is a scalar, y.shape == (n,).
            f_args is set by calling set_f_params(*args)
        jac : jac(t, y, *jac_args)
            Jacobian of the rhs, jac[i,j] = d f[i] / d y[j]
            jac_args is set by calling set_f_params(*args)
        """
        self.stiff = 0
        self.f = f
        self.jac  = jac
        self.f_params = ()
        self.jac_params = ()
        self.y = []

    def set_initial_value(self, y, t=0.0):
        """Set initial conditions y(t) = y."""
        if isscalar(y):
            y = [y]
        n_prev = len(self.y)
        if not n_prev:
            self.set_integrator('') # find first available integrator
        self.y = asarray(y, self._integrator.scalar)
        self.t = t
        self._integrator.reset(len(self.y),self.jac is not None)
        return self

    def set_integrator(self, name, **integrator_params):
        """
        Set integrator by name.

        Parameters
        ----------
        name : str
            Name of the integrator
        integrator_params
            Additional parameters for the integrator.
        """
        integrator = find_integrator(name)
        if integrator is None:
            print 'No integrator name match with %s or is not available.'\
                  %(`name`)
        else:
            self._integrator = integrator(**integrator_params)
            if not len(self.y):
                self.t = 0.0
                self.y = array([0.0], self._integrator.scalar)
            self._integrator.reset(len(self.y),self.jac is not None)
        return self

    def integrate(self, t, step=0, relax=0):
        """Find y=y(t), set y as an initial condition, and return y."""
        if step and self._integrator.supports_step:
            mth = self._integrator.step
        elif relax and self._integrator.supports_run_relax:
            mth = self._integrator.run_relax
        else:
            mth = self._integrator.run
        self.y,self.t = mth(self.f,self.jac or (lambda :None),
                            self.y,self.t,t,
                            self.f_params,self.jac_params)
        return self.y

    def successful(self):
        """Check if integration was successful."""
        try: self._integrator
        except AttributeError: self.set_integrator('')
        return self._integrator.success==1

    def set_f_params(self,*args):
        """Set extra parameters for user-supplied function f."""
        self.f_params = args
        return self

    def set_jac_params(self,*args):
        """Set extra parameters for user-supplied function jac."""
        self.jac_params = args
        return self

#------------------------------------------------------------------------------
# ODE integrators
#------------------------------------------------------------------------------

def find_integrator(name):
    for cl in IntegratorBase.integrator_classes:
        if re.match(name,cl.__name__,re.I):
            return cl
    return

class IntegratorBase(object):

    runner = None            # runner is None => integrator is not available
    success = None           # success==1 if integrator was called successfully
    supports_run_relax = None
    supports_step = None
    integrator_classes = []
    scalar = float

    def reset(self,n,has_jac):
        """Prepare integrator for call: allocate memory, set flags, etc.
        n - number of equations.
        has_jac - if user has supplied function for evaluating Jacobian.
        """

    def run(self,f,jac,y0,t0,t1,f_params,jac_params):
        """Integrate from t=t0 to t=t1 using y0 as an initial condition.
        Return 2-tuple (y1,t1) where y1 is the result and t=t1
        defines the stoppage coordinate of the result.
        """
        raise NotImplementedError,\
        'all integrators must define run(f,jac,t0,t1,y0,f_params,jac_params)'

    def step(self,f,jac,y0,t0,t1,f_params,jac_params):
        """Make one integration step and return (y1,t1)."""
        raise NotImplementedError,'%s does not support step() method' %\
              (self.__class__.__name__)

    def run_relax(self,f,jac,y0,t0,t1,f_params,jac_params):
        """Integrate from t=t0 to t>=t1 and return (y1,t)."""
        raise NotImplementedError,'%s does not support run_relax() method' %\
              (self.__class__.__name__)

    #XXX: __str__ method for getting visual state of the integrator

class vode(IntegratorBase):
    try:
        import vode as _vode
    except ImportError:
        print sys.exc_value
        _vode = None
    runner = getattr(_vode,'dvode',None)

    messages = {-1:'Excess work done on this call. (Perhaps wrong MF.)',
                -2:'Excess accuracy requested. (Tolerances too small.)',
                -3:'Illegal input detected. (See printed message.)',
                -4:'Repeated error test failures. (Check all input.)',
                -5:'Repeated convergence failures. (Perhaps bad'
                ' Jacobian supplied or wrong choice of MF or tolerances.)',
                -6:'Error weight became zero during problem. (Solution'
                ' component i vanished, and ATOL or ATOL(i) = 0.)'
                }
    supports_run_relax = 1
    supports_step = 1

    def __init__(self,
                 method = 'adams',
                 with_jacobian = 0,
                 rtol=1e-6,atol=1e-12,
                 lband=None,uband=None,
                 order = 12,
                 nsteps = 500,
                 max_step = 0.0, # corresponds to infinite
                 min_step = 0.0,
                 first_step = 0.0, # determined by solver
                 ):

        if re.match(method,r'adams',re.I): self.meth = 1
        elif re.match(method,r'bdf',re.I): self.meth = 2
        else: raise ValueError,'Unknown integration method %s'%(method)
        self.with_jacobian = with_jacobian
        self.rtol = rtol
        self.atol = atol
        self.mu = uband
        self.ml = lband

        self.order = order
        self.nsteps = nsteps
        self.max_step = max_step
        self.min_step = min_step
        self.first_step = first_step
        self.success = 1

    def reset(self,n,has_jac):
        # Calculate parameters for Fortran subroutine dvode.
        if has_jac:
            if self.mu is None and self.ml is None:
                miter = 1
            else:
                if self.mu is None: self.mu = 0
                if self.ml is None: self.ml = 0
                miter = 4
        else:
            if self.mu is None and self.ml is None:
                if self.with_jacobian:
                    miter = 2
                else:
                    miter = 0
            else:
                if self.mu is None: self.mu = 0
                if self.ml is None: self.ml = 0
                if self.ml==self.mu==0:
                    miter = 3
                else:
                    miter = 5
        mf = 10*self.meth + miter
        if mf==10:
            lrw = 20 + 16*n
        elif mf in [11,12]:
            lrw = 22 + 16*n + 2*n*n
        elif mf == 13:
            lrw = 22 + 17*n
        elif mf in [14,15]:
            lrw = 22 + 18*n + (3*self.ml+2*self.mu)*n
        elif mf == 20:
            lrw =  20 +  9*n
        elif mf in [21,22]:
            lrw = 22 + 9*n + 2*n*n
        elif mf == 23:
            lrw = 22 + 10*n
        elif mf in [24,25]:
            lrw = 22 + 11*n + (3*self.ml+2*self.mu)*n
        else:
            raise ValueError,'Unexpected mf=%s'%(mf)
        if miter in [0,3]:
            liw = 30
        else:
            liw = 30 + n
        rwork = zeros((lrw,), float)
        rwork[4] = self.first_step
        rwork[5] = self.max_step
        rwork[6] = self.min_step
        self.rwork = rwork
        iwork = zeros((liw,), int32)
        if self.ml is not None:
            iwork[0] = self.ml
        if self.mu is not None:
            iwork[1] = self.mu
        iwork[4] = self.order
        iwork[5] = self.nsteps
        iwork[6] = 2           # mxhnil
        self.iwork = iwork
        self.call_args = [self.rtol,self.atol,1,1,self.rwork,self.iwork,mf]
        self.success = 1

    def run(self,*args):
        y1,t,istate = self.runner(*(args[:5]+tuple(self.call_args)+args[5:]))
        if istate <0:
            print 'vode:',self.messages.get(istate,'Unexpected istate=%s'%istate)
            self.success = 0
        else:
            self.call_args[3] = 2 # upgrade istate from 1 to 2
        return y1,t

    def step(self,*args):
        itask = self.call_args[2]
        self.call_args[2] = 2
        r = self.run(*args)
        self.call_args[2] = itask
        return r

    def run_relax(self,*args):
        itask = self.call_args[2]
        self.call_args[2] = 3
        r = self.run(*args)
        self.call_args[2] = itask
        return r

if vode.runner:
    IntegratorBase.integrator_classes.append(vode)


class zvode(vode):
    try:
        import vode as _vode
    except ImportError:
        print sys.exc_value
        _vode = None
    runner = getattr(_vode,'zvode',None)

    supports_run_relax = 1
    supports_step = 1
    scalar = complex

    def reset(self, n, has_jac):
        # Calculate parameters for Fortran subroutine dvode.
        if has_jac:
            if self.mu is None and self.ml is None:
                miter = 1
            else:
                if self.mu is None: self.mu = 0
                if self.ml is None: self.ml = 0
                miter = 4
        else:
            if self.mu is None and self.ml is None:
                if self.with_jacobian:
                    miter = 2
                else:
                    miter = 0
            else:
                if self.mu is None: self.mu = 0
                if self.ml is None: self.ml = 0
                if self.ml==self.mu==0:
                    miter = 3
                else:
                    miter = 5

        mf = 10*self.meth + miter

        if mf in (10,):
            lzw = 15*n
        elif mf in (11, 12):
            lzw = 15*n + 2*n**2
        elif mf in (-11, -12):
            lzw = 15*n + n**2
        elif mf in (13,):
            lzw = 16*n
        elif mf in (14,15):
            lzw = 17*n + (3*self.ml + 2*self.mu)*n
        elif mf in (-14,-15):
            lzw = 16*n + (2*self.ml + self.mu)*n
        elif mf in (20,):
            lzw = 8*n
        elif mf in (21, 22):
            lzw = 8*n + 2*n**2
        elif mf in (-21,-22):
            lzw = 8*n + n**2
        elif mf in (23,):
            lzw = 9*n
        elif mf in (24, 25):
            lzw = 10*n + (3*self.ml + 2*self.mu)*n
        elif mf in (-24, -25):
            lzw = 9*n + (2*self.ml + self.mu)*n

        lrw = 20 + n

        if miter in (0, 3):
            liw = 30
        else:
            liw = 30 + n

        zwork = zeros((lzw,), complex)
        self.zwork = zwork

        rwork = zeros((lrw,), float)
        rwork[4] = self.first_step
        rwork[5] = self.max_step
        rwork[6] = self.min_step
        self.rwork = rwork

        iwork = zeros((liw,), int32)
        if self.ml is not None:
            iwork[0] = self.ml
        if self.mu is not None:
            iwork[1] = self.mu
        iwork[4] = self.order
        iwork[5] = self.nsteps
        iwork[6] = 2           # mxhnil
        self.iwork = iwork

        self.call_args = [self.rtol,self.atol,1,1,
                          self.zwork,self.rwork,self.iwork,mf]
        self.success = 1

    def run(self,*args):
        y1,t,istate = self.runner(*(args[:5]+tuple(self.call_args)+args[5:]))
        if istate < 0:
            print 'zvode:', self.messages.get(istate,
                                              'Unexpected istate=%s'%istate)
            self.success = 0
        else:
            self.call_args[3] = 2 # upgrade istate from 1 to 2
        return y1, t

if zvode.runner:
    IntegratorBase.integrator_classes.append(zvode)
