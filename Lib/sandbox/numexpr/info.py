"""Fast numerical expression evaluator

Usage:

Build an expression up by using E.<variable> for the variables:

>>> ex = 2.0 * E.a + 3 * E.b * E.c

then compile it to a function:

>>> func = numexpr(ex, input_names=('a', 'b', 'c'))

(if input_names is not given to compile, the variables are sort lexically.)

Then, you can use that function for elementwise operations:

>>> func(array([1., 2, 3]), array([4., 5, 6]), array([7., 8, 9]))
array([  86.,  124.,  168.])

Currently, this is only implemented for arrays of float64, and only
for the simple operations +, -, *, and /.

Copyright 2006 David M. Cooke <cookedm@physics.mcmaster.ca>
Licenced under a BSD-style license. See LICENSE.txt in the scipy source
directory.
"""

__all__ = ['E', 'numexpr', 'evaluate']
depends = ['core', 'testing']
