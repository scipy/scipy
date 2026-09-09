"""RBF backend with JAX JIT compiled evaluations.
"""
from ._rbfinterp_xp import *

# bring the _undersored names, too
from ._rbfinterp_xp import (
    _build_and_solve_system, _build_evaluation_coefficients, _build_system,
    kernel_matrix, _monomial_powers, _monomial_powers_impl, polynomial_matrix
)

import jax

compute_interpolation = jax.jit(compute_interpolation, static_argnames=["kernel", "xp"])

