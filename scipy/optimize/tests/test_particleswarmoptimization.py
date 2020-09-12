"""
Unit test for the particle swarm optimization algorithm
"""

from scipy.optimize import particle_swarm_optimization
from scipy.optimize._particleswarmoptimization import PSO

import numpy as np
from scipy.optimize._constraints import old_bound_to_new
from numpy.testing import assert_equal, assert_allclose, assert_array_equal
from pytest import raises as assert_raises


class TestParticleSwarmOptimization:

    def setup_method(self):
        self.ld_bounds = ((-5, 5),) * 2
        self.hd_bounds = ((-5, 5),) * 10
        self.in_bounds = ((2, 4),) * 2
        self.fun_args = (10,)
        self.particles = 60
        self.inertia_weight = 1.2
        self.individual_weight = 2.5
        self.social_weight = 3.0
        self.vel_up = ['constant', 'random']
        self.vel_bounds = ((-1, 1),) * 2
        self.mut = 0.999

    def squaresum(self, x, args=()):
        # Using Multidimensional Sphere Shifted function for performing tests
        if args:
            shift = args
        else:
            shift = 0

        return sum([el ** 2 for el in x]) + shift

    def test_low_dim_bounds(self):
        """ Test the algorithm solved a simple problem """

        # random velocity strategy gives more accuracy with less iterations
        res = particle_swarm_optimization(self.squaresum, 
                                          self.ld_bounds,
                                          velocity_strategy='random',
                                          max_iter=10000)

        assert_allclose(res.x, np.array([0] * 2), atol=1e-3)
        assert_allclose(res.fun, 0, atol=1e-3)

    def test_high_dim_bounds(self):
        """ Test that the algorithm solved a simple high dimensional problem """

        # random velocity strategy gives more accuracy with less iterations
        res = particle_swarm_optimization(self.squaresum, 
                                          self.hd_bounds,
                                          velocity_strategy='random',
                                          max_iter=10000)

        assert_allclose(res.x, np.array([0] * 10), atol=1e-3)
        assert_allclose(res.fun, 0, atol=1e-3)

    def test_num_particles(self):
        """ Test the number of particles of the swarm is correctly assigned """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds,
                        num_particles=10)
        algorithm.solve()

        assert_equal(10, len(algorithm.swarm))

    def test_random_generation(self):
        """ Test random generation of particles  """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds)

        a = algorithm._random_generation(self.ld_bounds)
        b = algorithm._random_generation(self.ld_bounds)

        assert_raises(AssertionError, assert_array_equal, a, b)


    def test_bounds(self):
        """ Test the algorithm initializes the particles inside 
        the global boundaries """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds,
                        num_particles=100, 
                        max_iter=100)
        # the initial evaluations are taken in account in the maximum number
        # of iterations
        algorithm.solve()

        low_ld_bounds, up_ld_bounds = old_bound_to_new(self.ld_bounds)

        out_bounds=False
        for particle in algorithm.swarm:
            if np.any(particle.position < low_ld_bounds):
               out_bounds=True
               break
            if np.any(particle.position > up_ld_bounds):
               out_bounds=True
               break

        assert out_bounds == False


    def test_initial_bounds(self):
        """ Test the algorithm initializes the particles inside 
        the initial boundaries """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds,
                        initial_bounds=self.in_bounds,
                        num_particles=100, 
                        max_iter=100)
        # the initial evaluations are taken in account in the maximum number
        # of iterations
        algorithm.solve()

        low_in_bounds, up_in_bounds = old_bound_to_new(self.in_bounds)

        out_bounds=False
        for particle in algorithm.swarm:
            if np.any(particle.position < low_in_bounds):
               out_bounds=True
               break
            if np.any(particle.position > up_in_bounds):
               out_bounds=True
               break

        assert out_bounds == False

    def test_function_arguments(self):
        """ Test the algorithm assigns the objective function arguments """

        # random velocity strategy gives more accuracy with less iterations
        res = particle_swarm_optimization(self.squaresum, 
                                          self.ld_bounds, 
                                          args=self.fun_args,
                                          velocity_strategy='random',
                                          max_iter=10000)

        assert_allclose(res.fun, self.fun_args[0], atol=1e-3)

    def test_inertia_weight(self):
        """ Test the correct assignment of the inertia weight """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds, 
                        w=self.inertia_weight)
        algorithm.solve()

        assert_equal(self.inertia_weight, algorithm.w)

    def test_individual_weight(self):
        """ Test the correct assignment of the individual weight """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds, 
                        c1=self.individual_weight)
        algorithm.solve()

        assert_equal(self.individual_weight, algorithm.c1)

    def test_social_weight(self):
        """ Test the correct assignment of the social weight """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds, 
                        c2=self.social_weight)
        algorithm.solve()

        assert_equal(self.social_weight, algorithm.c2)

    def test_velocity_update(self):
        """ Test the velocity strategy """
        algorithm = PSO(self.squaresum, 
                         self.ld_bounds)
        algorithm.solve()

        c1 = algorithm._velocity_update(self.vel_up[0], algorithm.swarm[0])
        c2 = algorithm._velocity_update(self.vel_up[0], algorithm.swarm[0])
        assert_array_equal(c1, c2)

        r1 = algorithm._velocity_update(self.vel_up[1], algorithm.swarm[0])
        r2 = algorithm._velocity_update(self.vel_up[1], algorithm.swarm[0])
        assert_raises(AssertionError, assert_array_equal, r1, r2)

        assert_raises(AssertionError, assert_array_equal, c1, r1)

    def test_velocity_bounds(self):
        """ Test the algorithm keeps the velocity inside 
        the velocity boundaries """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds,
                        velocity_bounds = self.vel_bounds,
                        num_particles=100, 
                        max_iter=200)
        algorithm.solve()

        low_vel_bounds, up_vel_bounds = old_bound_to_new(self.vel_bounds)

        out_bounds=False
        for particle in algorithm.swarm:
            if np.any(particle.velocity < low_vel_bounds):
               out_bounds=True
               break
            if np.any(particle.velocity > up_vel_bounds):
               out_bounds=True
               break

        assert out_bounds == False

    def test_mutation(self):
        """ Test the algorithm mutation method """
        algorithm = PSO(self.squaresum, 
                        self.ld_bounds,
                        mutation = self.mut,
                        max_iter=10000)
        algorithm.solve()

        different_minimum = False
        for particle in algorithm.swarm:
            if not np.allclose(particle.position, np.array([0, 0]), atol=1e-3):
                different_minimum = True
                break

        assert different_minimum == True
