"""
particle_swarm_optimization: The particle swarm optimization algorithm

David Martínez Rodríguez - damarro3@upv.es
Insituto Universitario de Matematica Multidisciplinar - IMM
Universitat Politècnica de València
"""

import numpy as np
from scipy.optimize import OptimizeResult
from scipy.optimize._constraints import old_bound_to_new


__all__ = ['particle_swarm_optimization']

_MACHEPS = np.finfo(np.float64).eps


def particle_swarm_optimization(func, bounds, initial_bounds=None, args=None,
                                num_particles=30, w=0.9, c1=2, c2=2,
                                velocity_strategy='constant',
                                velocity_bounds=None, mutation=0.0, 
                                max_iter=1000):
    """ Algorithm to find the global minimum of a multivariate function.

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
        ``len(initial_bounds)`` is used to determine the number of parameters 
        in ``x``.
    args : tuple, optional
        Any additional fixed parameters needed to completely specify the 
        objective function.
    num_particles : integer, optional
        Number of particles that form the swarm. 
        The default value is 30.
    w : float, optional
        Inertia weight. 
        The default value is 0.9.
    c1 : float, optional
        Individual weight. 
        The default value is 2.0.
    c2 : float, optional
        Social weight.
        The default value is 2.0.
    velocity_strategy: str, optional
        The velocity update strategy to use. Should be one of:
            - 'constant'
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
        interval [0, 1[.
        The default value is 0.0.
    max_iter : integer, optional
        Maximum number of function evaluations of the optimization algorithm. 
        The default value is 1000.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Attributes are: ``x`` the solution array and ``fun`` the
        minimum obtained.

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
                                             max_iter=100000)
    >>> result.x, result.fun
    (array([1.03943528, 1.07281126]), 0.007353124231320812)

    References
    ----------
    .. [1]  J. Kennedy and R. Eberhart. Particle swarm optimization.  
    		InProceedings of ICNN 95 - InternationalConference on Neural 
    		Networks. IEEE, 1995.
    .. [2] Federico Marini and Beata Walczak. 
    	   Particle swarm optimization (PSO). A tutorial. 
    	   Chemometrics and Intelligent Laboratory Systems, 149:153 – 165, 2015.

    """

    solver = PSO(func, bounds, initial_bounds=None, args=None,
                 num_particles=30, w=0.9, c1=2, c2=2,
                 velocity_strategy='constant',
                 velocity_bounds=None, mutation=0.0, 
                 max_iter=1000)
    solution = solver.solve()

    return solution



class PSO:
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
    num_particles : integer, optional
        Number of particles that form the swarm. 
        The default value is 30.
    w : float, optional
        Inertia weight. 
        The default value is 0.9.
    c1 : float, optional
        Individual weight. 
        The default value is 2.0.
    c2 : float, optional
        Social weight.
        The default value is 2.0.
    velocity_strategy: str, optional
        The velocity update strategy to use. Should be one of:
            - 'constant'
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
    max_iter : integer, optional
        Maximum number of function evaluations during the optimization process. 
        The default value is 1000.
    """

    def __init__(self, func, bounds, initial_bounds=None, args=None, 
                 num_particles=30, w=0.9, c1=2, c2=2,
                 velocity_strategy='constant', velocity_bounds=None, 
                 mutation=0.0, max_iter=1000):

        self.func = _FunctionWrapper(func, args)
        self.bounds = bounds
        self.num_particles = num_particles
        self.w = w
        self.c1 = c1
        self.c2 = c2
        self.velocity_strategy = velocity_strategy
        self.velocity_bounds = velocity_bounds
        self.mutation = mutation
        self.max_iter = max_iter

        if initial_bounds is not None:
            self.initial_bounds = initial_bounds
        else:
            self.initial_bounds = bounds

        self.pos_min, self.pos_max = old_bound_to_new(self.bounds)
        self.initial_pos_min, self.initial_pos_max = old_bound_to_new(
        	                                                self.initial_bounds)

        self.num_param = len(self.pos_min)

        # set velocity limits
        if velocity_bounds is None:
            self.vel_min = -0.2 * (self.pos_max - self.pos_min)
            self.vel_max = -self.vel_min
        else:
            self.vel_min, self.vel_max = old_bound_to_new(velocity_bounds)

        self.best_global_function_eval = np.Inf

        self.it = 0

        # test correct form of class parameters
        self._test()

    def _test(self):
        """ Test for parameters"""

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
        if self.num_particles <= 0:
            raise ValueError('The number of particles must be greater'
                             ' than 0')
        if self.num_particles > self.max_iter:
            raise ValueError('The number of particles must be smaller'
                             ' than max_iter')

        # test w
        if self.w <= 0:
            raise ValueError('The inertia weight must be greater'
                             ' than 0')

        # test c1
        if self.c1 <= 0:
            raise ValueError('The individual weight must be greater'
                             ' than 0')

        # test c2
        if self.c2 <= 0:
            raise ValueError('The social weight must be greater'
                             ' than 0')

        # test velocity_strategy
        if self.velocity_strategy not in ['constant', 'random']:
            raise ValueError('Velocity strategy must be'
                             ' ``constant`` or ``random``')

        # test velocity_bounds
        if self.velocity_bounds is not None:
	        if (np.size(self.velocity_bounds, 1) != 2 or not 
	            np.all(np.isfinite(self.velocity_bounds))):
	            raise ValueError('velocity_bounds should be a sequence'
	            	             ' containing real valued (min, max) pairs for'
	            	             ' each value in x')

        # test mutation
        if (self.mutation >= 1) or (self.mutation < 0):
            raise ValueError('The mutation value must be in the'
                             ' interval [0, 1)')

    def _particle_initialization(self):
        """ Initialization of particle values """

        particle_vector = self.swarm

        for part in particle_vector:

            # Random position in the swarm
            part.position = self._random_generation(self.initial_bounds)

            # Velocity initialization
            part.velocity = np.zeros(self.num_param)

            # Evaluation
            part.function_eval = self.func(part.position)
            self.it += 1

            # Individual best update
            part.best_position = part.position
            part.best_function_eval = part.function_eval

            # Global best update
            if part.best_function_eval < self.best_global_function_eval:
                self.best_global_function_eval = part.best_function_eval
                self.best_global_position = part.best_position

    def _particle_update(self):
        """ Update of particle values """

        particle_vector = self.swarm

        for part in particle_vector:

            # Mutation of explorer particles
            if np.random.rand() < self.mutation:
                part.position = self._random_generation(self.bounds)
                part.velocity = np.zeros(self.num_param)
                part.best_function_eval = np.Inf

            # Velocity update
            part.velocity = self._velocity_update(self.velocity_strategy,
                                                  part)


            # Velocity limits
            part.velocity = np.maximum(part.velocity, self.vel_min)
            part.velocity = np.minimum(part.velocity, self.vel_max)

            # Position update
            part.position = part.position + part.velocity

            # Boundary limits
            part.position = np.maximum(part.position, self.pos_min)
            part.position = np.minimum(part.position, self.pos_max)

            # Evaluation
            part.function_eval = self.func(part.position)
            self.it += 1

            # Individual best update
            if part.function_eval < part.best_function_eval:
                part.best_function_eval = part.function_eval
                part.best_position = part.position

                # Global best update
                if part.best_function_eval < self.best_global_function_eval:
                    self.best_global_function_eval = part.best_function_eval
                    self.best_global_position = part.best_position

    def _random_generation(self, bounds):
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
        return np.array([np.random.uniform(*element) for element in bounds])

    def _velocity_update(self, strategy, p):
        """ Updates the velocity of a particle

        Parameters
        ----------
        strategy : str
            The velocity update strategy to use. Should be one of:
            - 'constant'
            - 'random'
        p : Particle class instance
            Particle to update its velocity
        """
        if strategy == 'constant':
            v = ((self.w * p.velocity) + 
                 (self.c1 * (p.best_position - p.position)) +
                 (self.c2 * (self.best_global_position - p.position)))
            return v
        if strategy == 'random':
            v = ((np.random.rand(self.num_param) * self.w * p.velocity) + 
                 (np.random.rand(self.num_param) * self.c1 * 
                                               (p.best_position - p.position)) +
                 (np.random.rand(self.num_param) * self.c2 * 
                                      (self.best_global_position - p.position)))
            return v


    def solve(self):
        """ Runs the PSO algorithm 
        
        Returns
        -------
        res : OptimizeResult
            The optimization result represented as a ``OptimizeResult`` object.
            Important attributes are: ``x`` the solution array and ``fun`` the
            minimum obtained.
        """

        # list of instances of Particle class
        self.swarm = [_Particle() for _ in range(self.num_particles)]

        # initialization of the swarm
        self._particle_initialization()

        # main loop of the algorithm 
        while (self.it < self.max_iter):
            self._particle_update()

        res = OptimizeResult(x=self.best_global_position, 
        					 fun=self.best_global_function_eval)

        return res


class _Particle:
    """
    Particle class definition
    """

    def __init__(self,
                 position=None,
                 velocity=None,
                 function_eval=np.Inf,
                 best_position=None,
                 best_function_eval=np.Inf):
        self.position = position
        self.velocity = velocity
        self.function_eval = function_eval
        self.best_position = best_position
        self.best_function_eval = best_function_eval


class _FunctionWrapper:
    """
    Object to wrap user cost function, allowing picklability
    """
    def __init__(self, f, args):
        self.f = f
        self.args = [] if args is None else args

    def __call__(self, x):
        return self.f(x, *self.args)
