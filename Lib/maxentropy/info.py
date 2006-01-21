"""
Routines for fitting maximum entropy models
===========================================

Contains two classes. One class, "model", is for small discrete sample
spaces, using explicit summation. The other, "bigmodel", is for sample
spaces that are either continuous (and perhaps high-dimensional) or
discrete but too large to sum over, and uses importance sampling or
conditional Monte Carlo methods.

See the file bergerexample.py for a walk-through of how to use these
routines when the sample space can be enumerated.

See bergerexamplesimulated.py for a a similar walk-through using
simulation.

Copyright: Ed Schofield, 2003-2006
License: BSD-style (see LICENSE.txt in main source directory)
"""

postpone_import = 1
depends = ['optimize']
