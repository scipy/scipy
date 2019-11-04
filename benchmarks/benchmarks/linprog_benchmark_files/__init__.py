# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

"""
==============================================================================
`` --  Problems for testing linear programming routines
==============================================================================

This module provides a comprehensive set of problems for benchmarking linear 
programming routines, that is, scipy.optimize.linprog with method =
'interior-point' or 'simplex'.

"""

"""
All problems are from the Netlib LP Test Problem Set, courtesy of CUTEr
ftp://ftp.numerical.rl.ac.uk/pub/cutest/netlib/netlib.html

Converted from SIF (MPS) format by Matt Haberland
"""

__all__ = [s for s in dir() if not s.startswith('_')]
