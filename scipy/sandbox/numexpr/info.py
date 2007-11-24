"""Fast numerical expression evaluator

Usage:

The easiest way is to use the ``evaluate`` function:

>>> a = numpy.array([1., 2, 3])
>>> b = numpy.array([4, 5, 6])
>>> c = numpy.array([7., 8, 9])
>>> numexpr.evaluate("2.0 * a + 3 * b * c")
array([  86.,  124.,  168.])

This works for every datatype defined in NumPy and the next operators
are supported:

    Logical operators: &, |, ~
    Comparison operators: <, <=, ==, !=, >=, >
    Unary arithmetic operators: -
    Binary arithmetic operators: +, -, *, /, **, %

Note: Types do not support all operators. Boolean values only support
    logical and strict (in)equality comparison operators, while
    strings only support comparisons, numbers do not work with logical
    operators, and complex comparisons can only check for strict
    (in)equality. Unsupported operations (including invalid castings)
    raise NotImplementedError exceptions.

The next functions are also supported:

    where(bool, number1, number2): number - number1 if the bool
    condition is true, number2 otherwise.

    {sin,cos,tan}(float|complex): float|complex - trigonometric sinus,
    cosinus or tangent.

    {arcsin,arccos,arctan}(float|complex): float|complex -
    trigonometric inverse sinus, cosinus or tangent.

    arctan2(float1, float2): float - trigonometric inverse tangent of
    float1/float2.

    {sinh,cosh,tanh}(float|complex): float|complex - hyperbolic sinus,
    cosinus or tangent.

    sqrt(float|complex): float|complex - square root.

    {real,imag}(complex): float - real or imaginary part of complex.

    complex(float, float): complex - complex from real and imaginary
    parts.


Copyright 2006,2007 David M. Cooke <cookedm@physics.mcmaster.ca>
Licenced under a BSD-style license. See LICENSE.txt in the scipy source
directory.
"""

depends = ['core', 'testing']
