"""
Functions for identifying crossings in signals.
"""
from __future__ import division, print_function, absolute_import

import numpy as np


__all__ = []


def _select_cross_comparator(cross_type):
    """
    Select comparators used in finding crossings.
    
    Parameters
    ----------
    cross_type : str
        Accepts 'up' or 'down'.
    
    Return
    ------
    comparator_1 : callable
        Function to use to compare two data points.
        Should take two arrays as arguments.
    comparator_2 : callable
        Function to use to compare two data points.
        Should take two arrays as arguments.

    """
    if cross_type == 'up':
        return np.less_equal, np.greater
    elif cross_type == 'down':
        return np.greater, np.less_equal
    else:
        raise ValueError('cross_type must be "up" or "down"')
