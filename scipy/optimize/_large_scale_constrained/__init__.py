"""This module contains the equality constrained SQP solver."""


from .equality_constrained_sqp import equality_constrained_sqp
from .tr_interior_point import tr_interior_point

__all__ = ['equality_constrained_sqp',
           'tr_interior_point']
