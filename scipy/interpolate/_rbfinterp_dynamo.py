"""RBF backend with torch dynamo JIT compiled evaluations.
"""
from ._rbfinterp_xp import *

# bring the _undersored names, too
from ._rbfinterp_xp import (
    _build_and_solve_system, _build_evaluation_coefficients, _build_system, 
    _kernel_matrix, _monomial_powers, _monomial_powers_impl, _polynomial_matrix
)

import torch

compute_interpolation = torch.compile(fullgraph=True, dynamic=True)(compute_interpolation)

# needed for tests
torch._dynamo.config.cache_size_limit = 160
