# Original Author: Travis Oliphant 2002
# Bug-fixes in 2006 by Tim Leslie


import numpy
from numpy import asarray, tan, exp, ones, squeeze, sign, \
     all, log, sqrt, pi, shape, array, minimum, where
from numpy import random

__all__ = ['anneal']

_double_min = numpy.finfo(float).min
_double_max = numpy.finfo(float).max
class base_schedule(object):
    def __init__(self):
        self.dwell = 20
        self.learn_rate = 0.5
        self.lower = -10
        self.upper = 10
        self.Ninit = 50
        self.accepted = 0
        self.tests = 0
        self.feval = 0
        self.k = 0
        self.T = None

    def init(self, **options):
        self.__dict__.update(options)
        self.lower = asarray(self.lower)
        self.lower = where(self.lower == numpy.NINF, -_double_max, self.lower)
        self.upper = asarray(self.upper)
        self.upper = where(self.upper == numpy.PINF, _double_max, self.upper)
        self.k = 0
        self.accepted = 0
        self.feval = 0
        self.tests = 0

    def getstart_temp(self, best_state):
        """ Find a matching starting temperature and starting parameters vector
        i.e. find x0 such that func(x0) = T0.

        Parameters
        ----------
        best_state : _state
            A _state object to store the function value and x0 found.

        Returns
        -------
        x0 : array
            The starting parameters vector.
        """

        assert(not self.dims is None)
        lrange = self.lower
        urange = self.upper
        fmax = _double_min
        fmin = _double_max
        for _ in range(self.Ninit):
            x0 = random.uniform(size=self.dims)*(urange-lrange) + lrange
            fval = self.func(x0, *self.args)
            self.feval += 1
            if fval > fmax:
                fmax = fval
            if fval < fmin:
                fmin = fval
                best_state.cost = fval
                best_state.x = array(x0)

        self.T0 = (fmax-fmin)*1.5
        return best_state.x

    def accept_test(self, dE):
        T = self.T
        self.tests += 1
        if dE < 0:
            self.accepted += 1
            return 1
        p = exp(-dE*1.0/self.boltzmann/T)
        if (p > random.uniform(0.0, 1.0)):
            self.accepted += 1
            return 1
        return 0

    def update_guess(self, x0):
        pass

    def update_temp(self, x0):
        pass


#  A schedule due to Lester Ingber
class fast_sa(base_schedule):
    def init(self, **options):
        self.__dict__.update(options)
        if self.m is None:
            self.m = 1.0
        if self.n is None:
            self.n = 1.0
        self.c = self.m * exp(-self.n * self.quench)

    def update_guess(self, x0):
        x0 = asarray(x0)
        u = squeeze(random.uniform(0.0, 1.0, size=self.dims))
        T = self.T
        y = sign(u-0.5)*T*((1+1.0/T)**abs(2*u-1)-1.0)
        xc = y*(self.upper - self.lower)
        xnew = x0 + xc
        return xnew

    def update_temp(self):
        self.T = self.T0*exp(-self.c * self.k**(self.quench))
        self.k += 1
        return

class cauchy_sa(base_schedule):
    def update_guess(self, x0):
        x0 = asarray(x0)
        numbers = squeeze(random.uniform(-pi/2, pi/2, size=self.dims))
        xc = self.learn_rate * self.T * tan(numbers)
        xnew = x0 + xc
        return xnew

    def update_temp(self):
        self.T = self.T0/(1+self.k)
        self.k += 1
        return

class boltzmann_sa(base_schedule):
    def update_guess(self, x0):
        std = minimum(sqrt(self.T)*ones(self.dims), (self.upper-self.lower)/3.0/self.learn_rate)
        x0 = asarray(x0)
        xc = squeeze(random.normal(0, 1.0, size=self.dims))

        xnew = x0 + xc*std*self.learn_rate
        return xnew

    def update_temp(self):
        self.k += 1
        self.T = self.T0 / log(self.k+1.0)
        return

class _state(object):
    def __init__(self):
        self.x = None
        self.cost = None

# TODO:
#     allow for general annealing temperature profile
#     in that case use update given by alpha and omega and
#     variation of all previous updates and temperature?

# Simulated annealing

def anneal(func, x0, args=(), schedule='fast', full_output=0,
           T0=None, Tf=1e-12, maxeval=None, maxaccept=None, maxiter=400,
           boltzmann=1.0, learn_rate=0.5, feps=1e-6, quench=1.0, m=1.0, n=1.0,
           lower=-100, upper=100, dwell=50):
    """Minimize a function using simulated annealing.

    Schedule is a schedule class implementing the annealing schedule.
    Available ones are 'fast', 'cauchy', 'boltzmann'

    Parameters
    ----------
    func : callable f(x, *args)
        Function to be optimized.
    x0 : ndarray
        Initial guess.
    args : tuple
        Extra parameters to `func`.
    schedule : base_schedule
        Annealing schedule to use (a class).
    full_output : bool
        Whether to return optional outputs.
    T0 : float
        Initial Temperature (estimated as 1.2 times the largest
        cost-function deviation over random points in the range).
    Tf : float
        Final goal temperature.
    maxeval : int
        Maximum function evaluations.
    maxaccept : int
        Maximum changes to accept.
    maxiter : int
        Maximum cooling iterations.
    learn_rate : float
        Scale constant for adjusting guesses.
    boltzmann : float
        Boltzmann constant in acceptance test
        (increase for less stringent test at each temperature).
    feps : float
        Stopping relative error tolerance for the function value in
        last four coolings.
    quench, m, n : float
        Parameters to alter fast_sa schedule.
    lower, upper : float or ndarray
        Lower and upper bounds on `x`.
    dwell : int
        The number of times to search the space at each temperature.

    Outputs
    -------
    xmin : ndarray
        Point giving smallest value found.
    retval : int
        Flag indicating stopping condition::

                0 : Cooled to global optimum
                1 : Cooled to final temperature
                2 : Maximum function evaluations
                3 : Maximum cooling iterations reached
                4 : Maximum accepted query locations reached

    Jmin : float
        Minimum value of function found.
    T : float
        Final temperature.
    feval : int
        Number of function evaluations.
    iters : int
        Number of cooling iterations.
    accept : int
        Number of tests accepted.

    """
    x0 = asarray(x0)
    lower = asarray(lower)
    upper = asarray(upper)

    schedule = eval(schedule+'_sa()')
    #   initialize the schedule
    schedule.init(dims=shape(x0),func=func,args=args,boltzmann=boltzmann,T0=T0,
                  learn_rate=learn_rate, lower=lower, upper=upper,
                  m=m, n=n, quench=quench, dwell=dwell)

    current_state, last_state, best_state = _state(), _state(), _state()
    if T0 is None:
        x0 = schedule.getstart_temp(best_state)
    else:
        best_state.x = None
        best_state.cost = 300e8

    last_state.x = asarray(x0).copy()
    fval = func(x0,*args)
    schedule.feval += 1
    last_state.cost = fval
    if last_state.cost < best_state.cost:
        best_state.cost = fval
        best_state.x = asarray(x0).copy()
    schedule.T = schedule.T0
    fqueue = [100, 300, 500, 700]
    iters = 0
    while 1:
        for n in range(dwell):
            current_state.x = schedule.update_guess(last_state.x)
            current_state.cost = func(current_state.x,*args)
            schedule.feval += 1

            dE = current_state.cost - last_state.cost
            if schedule.accept_test(dE):
                last_state.x = current_state.x.copy()
                last_state.cost = current_state.cost
                if last_state.cost < best_state.cost:
                    best_state.x = last_state.x.copy()
                    best_state.cost = last_state.cost
        schedule.update_temp()
        iters += 1
        # Stopping conditions
        # 0) last saved values of f from each cooling step
        #     are all very similar (effectively cooled)
        # 1) Tf is set and we are below it
        # 2) maxeval is set and we are past it
        # 3) maxiter is set and we are past it
        # 4) maxaccept is set and we are past it

        fqueue.append(squeeze(last_state.cost))
        fqueue.pop(0)
        af = asarray(fqueue)*1.0
        if all(abs((af-af[0])/af[0]) < feps):
            retval = 0
            if abs(af[-1]-best_state.cost) > feps*10:
                retval = 5
                print "Warning: Cooled to %f at %s but this is not" \
                      % (squeeze(last_state.cost), str(squeeze(last_state.x))) \
                      + " the smallest point found."
            break
        if (Tf is not None) and (schedule.T < Tf):
            retval = 1
            break
        if (maxeval is not None) and (schedule.feval > maxeval):
            retval = 2
            break
        if (iters > maxiter):
            print "Warning: Maximum number of iterations exceeded."
            retval = 3
            break
        if (maxaccept is not None) and (schedule.accepted > maxaccept):
            retval = 4
            break

    if full_output:
        return best_state.x, best_state.cost, schedule.T, \
               schedule.feval, iters, schedule.accepted, retval
    else:
        return best_state.x, retval



if __name__ == "__main__":
    from numpy import cos
    # minimum expected at ~-0.195
    func = lambda x: cos(14.5*x-0.3) + (x+0.2)*x
    print anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='cauchy')
    print anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='fast')
    print anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='boltzmann')

    # minimum expected at ~[-0.195, -0.1]
    func = lambda x: cos(14.5*x[0]-0.3) + (x[1]+0.2)*x[1] + (x[0]+0.2)*x[0]
    print anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='cauchy')
    print anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='fast')
    print anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='boltzmann')
