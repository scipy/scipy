"""
particle_swarm_optimization: The particle swarm optimization algorithm

Author: David Martínez Rodríguez - damarro3@upv.es
Insituto Universitario de Matematica Multidisciplinar - IMM
Universitat Politècnica de València
"""

import numpy as np
from scipy.optimize import OptimizeResult
from scipy.optimize._constraints import old_bound_to_new


__all__ = ['particle_swarm_optimization']

_MACHEPS = np.finfo(np.float64).eps


def particle_swarm_optimization(func, bounds, initial_bounds=None, args=(),
                                n_particles=30, w=0.9, c1=2, c2=2,
                                velocity_strategy='constant',
                                velocity_bounds=None, mutation=0.0,
                                maxiter=1000, seed=None):
    """
    Parameters
    ----------
    func : callable
        Objective function to be minimized.
        Must be in the form ``f(x, *args)``, where ``x`` is the argument in the
        form of a 1-D array and ``args`` is a tuple of any additional fixed
        parameters needed to completely specify the function.
    bounds : sequence or `Bounds`
        Bounds for variables. There are two ways to specify the bounds:
        1. Instance of `Bounds` class.
        2. ``(min, max)`` pairs for each element in ``x``, defining the finite
        lower and upper bounds for the optimizing argument of `func`. It is
        required to have ``len(bounds) == len(x)``. ``len(bounds)`` is used
        to determine the number of parameters in ``x``.
    initial_bounds : sequence or `Bounds`, optional
        Initial bounds for variables, inside the global search bounds. If it is
        not specified, the initial bounds are settle as the global bounds.
        There are two ways to specify the initial bounds:
        1. Instance of `Bounds` class.
        2. ``(min, max)`` pairs for each element in ``x``, defining the finite
        lower and upper initial bounds for the optimizing argument of `func`. It
        is required to have ``len(initial_bounds) == len(x)``.
    args : tuple, optional
        Any additional fixed parameters needed to completely specify the
        objective function.
    n_particles : integer, optional
        Number of particles that form the swarm.
        The default value is 30.
    w : float or tuple(float, float), optional
        Inertia weight.
        The default value is 0.9.
        If specified as a tuple ``(min, max)``, linearly decreasing velocity
        update strategy must be employed.
    c1 : float, optional
        Individual weight.
        The default value is 2.0.
    c2 : float, optional
        Social weight.
        The default value is 2.0.
    velocity_strategy: str, optional
        The velocity update strategy to use. Should be one of:
            - 'constant'
            - 'linearly decreasing'
            - 'random'
        The default value is 'constant'.
    velocity_bounds: sequence, optional
        ``(min, max)`` pairs for each element in ``x``, defining the finite
        lower and upper velocity bounds. It is required to have
        ``len(velocity) == len(x)``.
        The default value is a sequence
        ``(-0.2 * (bound max - bound min), 0.2 * (bound max - bound min))``
    mutation : float, optional
        Mutation probability of a particle, with the aim of avoiding the stuck
        of the algorithm in a local minima. The value must be in the
        interval [0, 1).
        The default value is 0.0.
    maxiter : integer, optional
        Maximum number of swarm updates during the optimization process.
        Notice that the maximum number of function evaluation will be product between
        maxiter and n_particles.
        The default value is 1000.
    seed : int, optional
        Seed for repeatable minimizations of random numbers.
        If not provided, then fresh, unpredictable entropy will be pulled
        from the OS.

    Returns
        -------
        res : OptimizeResult
            The optimization result represented as a ``OptimizeResult`` object.
            Attributes are: ``x`` the solution array, ``fun`` the
            minimum obtained, ``nfev`` the number of function evaluations performed,
            ``nit`` the number of swarm updates and ``message`` the reason of the
            algorithm ending.

        Notes
        -----
        PSO is a bio-inspired meta-heuristic algorithm developed by James Kennedy
        and Russell Eberhart in 1995 [1]_. It is based on the social behaviour of
        bird flocks, which try to make their movements in the most optimal possible
        way to find food. Since its proposal, the PSO algorithm has been
        successfully applied to solve a large amount of optimization problems.

        Mathematically, given a function :math:`f`, the optimization problem
        consists of finding:

        .. math::

            x_{\text{opt}}|f(x_{opt}) \leq f(x), \qquad \forall x \in \mathbb{R}^{D}


        The :math:`D`-dimensional domain of function :math:`f` in
        :math:`\mathbb{R}^{D}` is called the search space. Every point characterized
        by the :math:`D`-dimensional vector represents a candidate solution to the
        problem. These vectors are called particles[2]_.

        Let :math:`X` be a {\em swarm}, formed by :math:`N` particles:

        .. math::
            X = [x_{1}, x_{2}, \ldots, x_{i}, \ldots, x_{N}],

        where each particle is a D-dimensional vector in the domain of :math:`f`,

        .. math::
            x_{i} = [x_{i1}, x_{i2}, \ldots, x_{iD}].

        During the search for the optimum of the function :math:`f`, a maximum
        number of iterations :math:`T` is set. For each iteration :math:`t`,
        particles update their position according to the following formula:

        .. math::
            x_{i}(t+1) = x_{i}(t) + v_{i}(t+1),

        where the velocity vector $v_{i}(t+1)$ is updated according to:

        .. math::
            \begin{aligned}
            v_{i}(t+1) = &w \cdot v_{i}(t) +
                         &c_{1} \cdot (p_{i}-x_{i}(t)) +
                         &c_{2} \cdot (g-x_{i}(t)).
            \end{aligned}

        :math:`w` is the inertia weight, :math:`c_{1}$` is the individual weight
        and :math:`c_{2}` is the social weight. This three parameters are scalars
        and their values are settled independently from each other. :math:`w`
        values are usually in the range [0.5, 1.5], :math:`c1` and :math:`c2`
        are usually in the range [1, 3], and these values depend on the problem
        to be solved. Also, :math:`p_{i}` stands for the current best position of
        :math:`x_{i}` (local best position) and :math:`g` stands for the position
        with the best value among all the particles which have formed the swarm
        (the global best position).

        Examples
        --------
        Let us consider the problem of minimizing the Rosenbrock function. This
        function is implemented in `rosen` in `scipy.optimize`.

        >>> from scipy.optimize import rosen, particle_swarm_optimization
        >>> bounds = [(0,2), (0, 2)]
        >>> result = particle_swarm_optimization(rosen, bounds,
                                                 velocity_strategy='random',
                                                 mutation=0.05,
                                                 maxiter=100000)
        >>> result.x, result.fun
        (array([0.99999724, 0.99999448]), 7.612345716842813e-12)

        References
        ----------
        .. [1]  J. Kennedy and R. Eberhart. Particle swarm optimization.
                InProceedings of ICNN 95 - InternationalConference on Neural
                Networks. IEEE, 1995.
        .. [2] Federico Marini and Beata Walczak.
               Particle swarm optimization (PSO). A tutorial.
               Chemometrics and Intelligent Laboratory Systems, 149:153 – 165, 2015.
    """

    minimization = Swarm(func, bounds, initial_bounds=initial_bounds, args=args, 
                         n_particles=n_particles, w=w, c1=c1, c2=c2,
                         velocity_strategy=velocity_strategy, velocity_bounds=velocity_bounds,
                         mutation=mutation, maxiter=maxiter, seed=seed)
    minimization.solve()

    return minimization.solution


class Particle:

    def __init__(self, position, velocity, bestpos, bestof, func):
        self.position = position
        self.velocity = velocity
        self.bestpos = bestpos
        self.bestof = bestof
        self.func = func

    def updatevel(self, w, c1, c2, bestgpos, vmin, vmax):
        self.velocity = (w * self.velocity +
                         c1 * (self.bestpos - self.position) +
                         c2 * (bestgpos - self.position))

        self.velocity = np.maximum(self.velocity, vmin)
        self.velocity = np.minimum(self.velocity, vmax)

    def updatepos(self, pmin, pmax):
        self.position = self.position + self.velocity

        self.position = np.maximum(self.position, pmin)
        self.position = np.minimum(self.position, pmax)

        pos_value = self.func(self.position)
        if pos_value <= self.bestof:
            self.bestof = pos_value
            self.bestpos = self.position


class Swarm:
    """
    This class implements the particle swarm optimization solver

    Parameters
    ----------
    func : callable
        Objective function to be minimized.
        Must be in the form ``f(x, *args)``, where ``x`` is the argument in the 
        form of a 1-D array and ``args`` is a tuple of any additional fixed 
        parameters needed to completely specify the function.
    bounds : sequence or `Bounds`
        Bounds for variables. There are two ways to specify the bounds:
        1. Instance of `Bounds` class.
        2. ``(min, max)`` pairs for each element in ``x``, defining the finite
        lower and upper bounds for the optimizing argument of `func`. It is
        required to have ``len(bounds) == len(x)``. ``len(bounds)`` is used
        to determine the number of parameters in ``x``.
    initial_bounds : sequence or `Bounds`, optional
        Initial bounds for variables, inside the global search bounds. If it is 
        not specified, the initial bounds are settle as the global bounds. 
        There are two ways to specify the initial bounds:
        1. Instance of `Bounds` class.
        2. ``(min, max)`` pairs for each element in ``x``, defining the finite
        lower and upper initial bounds for the optimizing argument of `func`. It
        is required to have ``len(initial_bounds) == len(x)``.
    args : tuple, optional
        Any additional fixed parameters needed to completely specify the 
        objective function.
    n_particles : integer, optional
        Number of particles that form the swarm. 
        The default value is 30.
    w : float or tuple(float, float), optional
        Inertia weight. 
        The default value is 0.9.
        If specified as a tuple ``(min, max)``, linearly decreasing velocity 
        update strategy must be employed.
    c1 : float, optional
        Individual weight. 
        The default value is 2.0.
    c2 : float, optional
        Social weight.
        The default value is 2.0.
    velocity_strategy: str, optional
        The velocity update strategy to use. Should be one of:
            - 'constant'
            - 'linearly decreasing'
            - 'random'
        The default value is 'constant'.
    velocity_bounds: sequence, optional
        ``(min, max)`` pairs for each element in ``x``, defining the finite
        lower and upper velocity bounds. It is required to have 
        ``len(velocity) == len(x)``.
        The default value is a sequence 
        ``(-0.2 * (bound max - bound min), 0.2 * (bound max - bound min))``
    mutation : float, optional
        Mutation probability of a particle, with the aim of avoiding the stuck
        of the algorithm in a local minima. The value must be in the 
        interval [0, 1).
        The default value is 0.0.
    maxiter : integer, optional
        Maximum number of swarm updates during the optimization process.
        Notice that the maximum number of function evaluation will be product between
        maxiter and n_particles.
        The default value is 1000.
    seed : int, optional
        Seed for repeatable minimizations of random numbers.
        If not provided, then fresh, unpredictable entropy will be pulled 
        from the OS.
    """

    def __init__(self, func, bounds, initial_bounds=None, args=(),
                 n_particles=30, w=0.9, c1=2, c2=2,
                 velocity_strategy='constant', velocity_bounds=None, 
                 mutation=0.0, maxiter=1000, workers=1, seed=None):

        self.bounds = bounds
        self.n_particles = n_particles
        self.w = w
        self.c1 = c1
        self.c2 = c2
        self.velocity_strategy = velocity_strategy
        self.velocity_bounds = velocity_bounds
        self.mutation = mutation
        self.maxiter = maxiter
        self.workers = workers
        self.seed = seed

        if self.workers == -1:
            self.workers = multiprocessing.cpu_count()

        self.n_params = len(self.bounds)

        self.bestgpos = None
        self.bestgof = np.inf
        self.particles = [None for _ in range(self.n_particles)]

        # set position limits
        if initial_bounds is not None:
            self.initial_bounds = initial_bounds
        else:
            self.initial_bounds = bounds

        self.initial_pos_min, self.initial_pos_max = old_bound_to_new(self.initial_bounds)
        self.pos_min, self.pos_max = old_bound_to_new(self.bounds)

        # set velocity limits
        if velocity_bounds is None:
            self.vel_min = - 0.2 * (self.pos_max - self.pos_min)
            self.vel_max = 0.2 * (self.pos_max - self.pos_min)
        else:
            self.vel_min, self.vel_max = old_bound_to_new(velocity_bounds)

        # we create a wrapped function to allow the use of map (and Pool.map
        # in the future)
        self.func = _FunctionWrapper(func, args)

        # random number generator. There is a different one for each particle
        self.rng = np.random.default_rng(seed)

        self.it = 0
        self.fevals = 0

    def _test(self):
        """ Parameters tests"""

        # test objective function
        if not callable(self.func):
            raise ValueError('The objective function must be a callable')

        # test global bounds
        if (np.size(self.bounds, 1) != 2 or not 
            np.all(np.isfinite(self.bounds))):
            raise ValueError('bounds should be a sequence containing '
                             'real valued (min, max) pairs for each value'
                             ' in x')
        if np.any(self.pos_max < self.pos_min):
            raise ValueError('all bound (min) values must be smaller than'
                             ' (max) values')

        # test initial_bounds
        if np.any(self.initial_pos_max < self.initial_pos_min):
            raise ValueError('all initial_bounds (min) values must be smaller'
                             ' than (max) values')
        if (np.any(self.initial_pos_max > self.pos_max) or 
            np.any(self.initial_pos_min < self.pos_min)):
            raise ValueError('initial_bounds must be inside bounds')

        # test num_particles
        if self.n_particles <= 0:
            raise ValueError('The number of particles must be greater'
                             ' than 0')

        # test w
        if (((self.velocity_strategy == 'constant') or (self.velocity_strategy == 'random'))
                and (type(self.w) is not float)):
            raise ValueError('The inertia weight must be float if'
                             ' ``constant`` or ``random`` velocity strategy is'
                             ' selected')

        if (type(self.w) is float) and (self.w <= 0):
            raise ValueError('The inertia weight must be greater'
                             ' than 0')

        if (self.velocity_strategy == 'linearly decreasing') and (type(self.w) is not tuple):
            raise ValueError('The inertia weight must be a tuple if'
                             ' ``linearly decreasing`` velocity strategy is'
                             ' selected')
        if (self.velocity_strategy == 'linearly decreasing') and (self.w[0] >= self.w[1]):
            raise ValueError('The inertia weight must be in the form (w min, w max)')

        # test c1
        if self.c1 <= 0:
            raise ValueError('The individual weight must be greater'
                             ' than 0')

        # test c2
        if self.c2 <= 0:
            raise ValueError('The social weight must be greater'
                             ' than 0')

        # test velocity_strategy
        if self.velocity_strategy not in ['constant', 'linearly decreasing', 'random']:
            raise ValueError('Velocity strategy must be'
                             ' ``constant``, ``linearly decreasing`` or ``random``')

        # test velocity_bounds
        if self.velocity_bounds is not None:
            if np.size(self.velocity_bounds, 1) != 2 or not np.all(np.isfinite(self.velocity_bounds)):
                raise ValueError('velocity_bounds should be a sequence'
                                 ' containing real valued (min, max) pairs for'
                                 ' each value in x')

        # test mutation
        if (self.mutation >= 1) or (self.mutation < 0):
            raise ValueError('The mutation value must be in the'
                             ' interval [0, 1)')

        # test seed
        if (type(self.seed) is not int) or (type(self.seed) is not None):
            raise ValueError('seed must be int')

    def _initialization(self):
        """ Initialization of the swarm particles """
        for i in range(self.n_particles):
            self._particle_initialization(i)
            if self.particles[i].bestof <= self.bestgof:
                self.bestgpos = self.particles[i].bestpos
                self.bestgof = self.particles[i].bestof

        self.fevals += self.n_particles
        self.it += 1

    def _particle_initialization(self, i):
        """ Initialization of particle values """
        init_position = self._particle_random_generation(self.bounds)
        self.particles[i] = Particle(position=init_position,
                                     velocity=0,
                                     bestpos=init_position,
                                     bestof=self.func(init_position),
                                     func=self.func)

    def _update(self):
        """ Update of the swarm particles """
        for i in range(self.n_particles):
            self._particle_update(i)
            if self.particles[i].bestof <= self.bestgof:
                self.bestgpos = self.particles[i].bestpos
                self.bestgof = self.particles[i].bestof

        self.fevals += self.n_particles
        self.it += 1

    def _particle_update(self, i):
        """ Update of particle values """
        if self.rng.random() <= self.mutation:
            init_position = self._particle_random_generation(self.bounds)
            self.particles[i] = Particle(position=init_position,
                                         velocity=0,
                                         bestpos=init_position,
                                         bestof=self.func(init_position),
                                         func=self.func)
            
        else:
            if self.velocity_strategy == 'constant':
                self.particles[i].updatevel(self.w, self.c1, self.c2,
                                            self.bestgpos,
                                            self.vel_min, self.vel_max)
            if self.velocity_strategy == 'linearly decreasing':
                self.particles[i].updatevel(self.w[1] - (self.w[1] - self.w[0]) * self.it / self.maxiter,
                                            self.c1, self.c2,
                                            self.bestgpos,
                                            self.vel_min, self.vel_max)
            if self.velocity_strategy == 'random':
                self.particles[i].updatevel(self.rng.random(self.n_params) * self.w,
                                            self.rng.random(self.n_params) * self.c1,
                                            self.rng.random(self.n_params) * self.c2,
                                            self.bestgpos,
                                            self.vel_min, self.vel_max)

            self.particles[i].updatepos(self.pos_min, self.pos_max)

    def solve(self):
        self._initialization()
        while self.it < self.maxiter:
            self._update()

        if self.it == self.maxiter:
            message = 'Maximum number of iterations reached'
        else:
            message = 'Process not finished'

        self.solution = OptimizeResult(x=self.bestgpos,
                                       fun=self.bestgof,
                                       nfev=self.fevals,
                                       nit=self.it,
                                       message=message)

    def _particle_random_generation(self, bounds):
        """ Generates a random particle inside the given boundaries

        Parameters
        ----------
        bounds : sequence or `Bounds`
            Bounds for variables. There are two ways to specify the bounds:
            1. Instance of `Bounds` class.
            2. ``(min, max)`` pairs for each element in ``x``, defining the 
            finite lower and upper bounds for the optimizing argument of `func`.
            It is required to have ``len(bounds) == len(x)``. ``len(bounds)`` 
            is used to determine the number of parameters in ``x``.
        """
        # * deletes the tuple and transforms it into two numbers
        return np.array([self.rng.uniform(*element) for element in bounds])


class _FunctionWrapper:
    """
    Object to wrap user cost function, allowing picklability
    """
    def __init__(self, f, args):
        self.f = f
        self.args = [] if args is None else args

    def __call__(self, x):
        return self.f(x, *self.args)
