"""Differential Evolution Optimization

:Author: Robert Kern

Copyright 2005 by Robert Kern.
"""

import numpy as np


# Notes: for future modifications:
# Ali, M. M., and A. Toern. Topographical differential evolution using
# pre-calculated differentials. _Stochastic and Global Optimization_. 1--17.
#
#  A good scale value:
#    F = max(l_min, 1-min(abs(f_min/f_max), abs(f_max/f_min)))
#      ~ 0.3 <= l_min <= 0.4
#      ~ f_min and f_max are the minimum and maximum values in the initial
#        population.
#
#  Pre-calculated differentials:
#    Keep a set of differentials A.
#    For each x_i of the population S:
#      Every M steps, randomly select 3 points x_r1, x_r2, x_r3 from S (not x_i).
#        Compute new x_i using x_r1 + F*(x_r2-x_r3).
#        Store differential in A.
#      Each other step:
#        Randomly select x_r1 from S and a differential vector from A.
#      Crossover.
#
#  Convergence criterion:
#    f_max - f_min < eps
#
#  Topographical DEPD:
#    Two populations S and Sa (auxiliary).
#    Phase counter t = 0 and array shift[:] = False.
#    Stopping condition: e.g. t >= 4.
#    g << N, number of nearest neighbors to search for graph minima.
#    Ng << N, number of points for graph.
#    For each x_i in S, do DEPD as described above to get y_i.
#    if f(y_i) < f(x_i):
#      if shift[i] is False:
#        shift[i] = True
#        S[i] = y_i
#      else:
#        Sa[i] = y_i
#      if alltrue(shift,axis=0):
#        Find graph minima of f(x) using the Ng best points in S.
#        Do local search from each minimum.
#        Replace worst Ng points in S with best Ng points in Sa.
#        If best from this phase is better than previous best, t=0.
#        Else: t += 1.
#        shift[:] = False
#    Next generation.

class DiffEvolver(object):
    """Minimize a function using differential evolution.

    Constructors
    ------------
    DiffEvolver(func, pop0, args=(), crossover_rate=0.5, scale=None,
        strategy=('rand', 2, 'bin'), eps=1e-6)
      func -- function to minimize
      pop0 -- sequence of initial vectors
      args -- additional arguments to apply to func
      crossover_rate -- crossover probability [0..1] usually 0.5 or so
      scale -- scaling factor to apply to differences [0..1] usually > 0.5
        if None, then calculated from pop0 using a heuristic
      strategy -- tuple specifying the differencing/crossover strategy
        The first element is one of 'rand', 'best', 'rand-to-best' to specify
        how to obtain an initial trial vector.
        The second element is either 1 or 2 (or only 1 for 'rand-to-best') to
        specify the number of difference vectors to add to the initial trial.
        The third element is (currently) 'bin' to specify binomial crossover.
      eps -- if the maximum and minimum function values of a given generation are
        with eps of each other, convergence has been achieved.
      prng -- a RandomState instance. By default, this is the global
        numpy.random instance.

    DiffEvolver.frombounds(func, lbound, ubound, npop, crossover_rate=0.5,
        scale=None, strategy=('rand', 2, 'bin'), eps=1e-6)
      Randomly initialize the population within given rectangular bounds.
      lbound -- lower bound vector
      ubound -- upper bound vector
      npop -- size of population

    Public Methods
    --------------
    solve(newgens=100)
      Run the minimizer for newgens more generations. Return the best parameter
      vector from the whole run.

    Public Members
    --------------
    best_value -- lowest function value in the history
    best_vector -- minimizing vector
    best_val_history -- list of best_value's for each generation
    best_vec_history -- list of best_vector's for each generation
    population -- current population
    pop_values -- respective function values for each of the current population
    generations -- number of generations already computed
    func, args, crossover_rate, scale, strategy, eps -- from constructor
    """
    def __init__(self, func, pop0, args=(), crossover_rate=0.5, scale=None,
            strategy=('rand', 2, 'bin'), eps=1e-6, prng=np.random):
        self.func = func
        self.population = np.array(pop0)
        self.npop, self.ndim = self.population.shape
        self.args = args
        self.crossover_rate = crossover_rate
        self.strategy = strategy
        self.eps = eps

        self.pop_values = [self.func(m, *args) for m in self.population]
        bestidx = np.argmin(self.pop_values)
        self.best_vector = self.population[bestidx]
        self.best_value = self.pop_values[bestidx]

        if scale is None:
            self.scale = self.calculate_scale()
        else:
            self.scale = scale

        self.generations = 0
        self.best_val_history = []
        self.best_vec_history = []

        self.jump_table = {
            ('rand', 1, 'bin'): (self.choose_rand, self.diff1, self.bin_crossover),
            ('rand', 2, 'bin'): (self.choose_rand, self.diff2, self.bin_crossover),
            ('best', 1, 'bin'): (self.choose_best, self.diff1, self.bin_crossover),
            ('best', 2, 'bin'): (self.choose_best, self.diff2, self.bin_crossover),
            ('rand-to-best', 1, 'bin'):
                (self.choose_rand_to_best, self.diff1, self.bin_crossover),
            }

    def clear(self):
        self.best_val_history = []
        self.best_vec_history = []
        self.generations = 0
        self.pop_values = [self.func(m, *self.args) for m in self.population]

    def frombounds(cls, func, lbound, ubound, npop, crossover_rate=0.5,
            scale=None, strategy=('rand', 2, 'bin'), eps=1e-6, prng=np.random):
        lbound = np.asarray(lbound)
        ubound = np.asarray(ubound)
        pop0 = prng.uniform(lbound, ubound, size=(npop, len(lbound)))
        return cls(func, pop0, crossover_rate=crossover_rate, scale=scale,
            strategy=strategy, eps=eps)
    frombounds = classmethod(frombounds)

    def calculate_scale(self):
        rat = abs(max(self.pop_values)/self.best_value)
        rat = min(rat, 1./rat)
        return max(0.3, 1.-rat)

    def bin_crossover(self, oldgene, newgene):
        mask = self.prng.rand(self.ndim) < self.crossover_rate
        return np.where(mask, newgene, oldgene)

    def select_samples(self, candidate, nsamples):
        possibilities = range(self.npop)
        possibilities.remove(candidate)
        return self.prng.permutation(possibilities)[:nsamples]

    def diff1(self, candidate):
        i1, i2 = self.select_samples(candidate, 2)
        return self.scale * (self.population[i1] - self.population[i2])

    def diff2(self, candidate):
        i1, i2, i3, i4 = self.select_samples(candidate, 4)
        return self.scale * (self.population[i1] - self.population[i2] +
                             self.population[i3] - self.population[i4])

    def choose_best(self, candidate):
        return self.best_vector

    def choose_rand(self, candidate):
        i = self.select_samples(candidate, 1)[0]
        return self.population[i]

    def choose_rand_to_best(self, candidate):
        return ((1-self.scale) * self.population[candidate] +
                self.scale * self.best_vector)

    def get_trial(self, candidate):
        chooser, differ, crosser = self.jump_table[self.strategy]
        trial = crosser(self.population[candidate],
            chooser(candidate) + differ(candidate))
        return trial

    def converged(self):
        return max(self.pop_values) - min(self.pop_values) <= self.eps

    def solve(self, newgens=100):
        """Run for newgens more generations.

        Return best parameter vector from the entire run.
        """
        for gen in xrange(self.generations+1, self.generations+newgens+1):
            for candidate in range(self.npop):
                trial = self.get_trial(candidate)
                trial_value = self.func(trial, *self.args)
                if trial_value < self.pop_values[candidate]:
                    self.population[candidate] = trial
                    self.pop_values[candidate] = trial_value
                    if trial_value < self.best_value:
                        self.best_vector = trial
                        self.best_value = trial_value
            self.best_val_history.append(self.best_value)
            self.best_vec_history.append(self.best_vector)
            if self.converged():
                break
        self.generations = gen
        return self.best_vector
