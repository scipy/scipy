"""Collection of Model instances for use with the odrpack fitting package.
"""

from scipy.sandbox.odr.odrpack import Model
import numpy as sb
from types import *

def _lin_fcn(B, x, sum=sb.sum):
    a, b = B[0], B[1:]
    b.shape = (b.shape[0], 1)
    
    return a + sum(x*b)

def _lin_fjb(B, x, concatenate=sb.concatenate, Float=sb.Float, 
             ones=sb.ones, ravel=sb.ravel):
    a = ones((x.shape[-1],), Float)
    res = concatenate((a, ravel(x)))
    res.shape = (B.shape[-1], x.shape[-1])
    return res

def _lin_fjd(B, x, repeat=sb.repeat):
    b = B[1:]
    b = repeat(b, (x.shape[-1],)*b.shape[-1])
    b.shape = x.shape
    return b

def _lin_est(data):
    # Eh. The answer is analytical, so just return all ones.
    # Don't return zeros since that will interfere with 
    # ODRPACK's auto-scaling procedures.
    
    if len(data.x.shape) == 2:
        m = data.x.shape[0]
    else:
        m = 1

    return ones((m + 1,), Float)

def _poly_fcn(B, x, powers, power=sb.power, sum=sb.sum):
    a, b = B[0], B[1:]
    b.shape = (b.shape[0], 1)

    return a + sum(b * power(x, powers))

def _poly_fjacb(B, x, powers, power=sb.power,
                concatenate=sb.concatenate, Float=sb.Float, ones=sb.ones):
    res = concatenate((ones((x.shape[-1],), Float), power(x, powers).flat))
    res.shape = (B.shape[-1], x.shape[-1])
    return res

def _poly_fjacd(B, x, powers, power=sb.power, sum=sb.sum):
    b = B[1:]
    b.shape = (b.shape[0], 1)
    
    b = b * powers
    
    return sum(b * power(x, powers-1))

def _exp_fcn(B, x, exp=sb.exp):
    return B[0] + exp(B[1] * x)

def _exp_fjd(B, x, exp=sb.exp):
    return B[1] * exp(B[1] * x)

def _exp_fjb(B, x, exp=sb.exp, concatenate=sb.concatenate, ones=sb.ones, 
             Float=sb.Float):
    res = concatenate((ones((x.shape[-1],), sb.Float), x * exp(B[1] * x)))
    res.shape = (2, x.shape[-1])
    return res

def _exp_est(data):
    # Eh.
    return array([1., 1.])

multilinear = Model(_lin_fcn, fjacb=_lin_fjb, 
               fjacd=_lin_fjd, estimate=_lin_est,
               meta={'name': 'Arbitrary-dimensional Linear',
                     'equ':'y = B_0 + Sum[i=1..m, B_i * x_i]',
                     'TeXequ':'$y=\\beta_0 + \sum_{i=1}^m \\beta_i x_i$'})

def polynomial(order):
    """Factory function for a general polynomial model.

The argument "order" can be either an integer, where it becomes the 
order of the polynomial to fit, or it can be a sequence of numbers to 
explicitly determine the powers in the polynomial.

Oh yeah, a constant is always included, so don't include 0.

Thus, polynomial(n) is equivalent to polynomial(range(1, n+1)).

  polynomial(order)"""

    if type(order) is int:
        order = range(1, order+1)
    
    powers = sb.asarray(order)
    powers.shape = (len(powers), 1)

    len_beta = len(powers) + 1

    def _poly_est(data, len_beta=len_beta):
        # Eh. Ignore data and return all ones.

        return sb.ones((len_beta,), sb.Float)

    return Model(_poly_fcn, fjacd=_poly_fjd, fjacb=_poly_fjb, 
                 estimate=_poly_est, extra_args=(powers,),
                 meta={'name': 'Sorta-general Polynomial',
                 'equ':'y = B_0 + Sum[i=1..%s, B_i * (x**i)]' % (len_beta-1),
                 'TeXequ':'$y=\\beta_0 + \sum_{i=1}^{%s} \\beta_i x^i$' %\
                        (len_beta-1)})

exponential = Model(_exp_fcn, fjacd=_exp_fjd, fjacb=_exp_fjb,
                    estimate=_exp_est, meta={'name':'Exponential',
                    'equ':'y= B_0 + exp(B_1 * x)',
                    'TeXequ':'$y=\\beta_0 + e^{\\beta_1 x}$'})

def _unilin(B, x):
    return x*B[0] + B[1]

def _unilin_fjd(B, x):
    return sb.ones(x.shape, sb.Float) * B[0]

def _unilin_fjb(B, x, cat=sb.concatenate):
    _ret = cat((x,sb.ones(x.shape, sb.Float)))
    _ret.shape = (2,) + x.shape

    return _ret

def _unilin_est(data):
    return (1., 1.)

def _quadratic(B, x):
    return x*(x*B[0] + B[1]) + B[2]

def _quad_fjd(B, x):
    return 2*x*B[0] + B[1]

def _quad_fjb(B,x,cat=sb.concatenate):
    _ret = cat((x*x, x, sb.ones(x.shape, Float)))
    _ret.shape = (3,) + x.shape

    return _ret

def _quad_est(data):
    return (1.,1.,1.)

unilinear = Model(_unilin, fjacd=_unilin_fjd, fjacb=_unilin_fjb,
                  estimate=_unilin_est, meta={'name': 'Univariate Linear',
                  'equ': 'y = B_0 * x + B_1',
                  'TeXequ': '$y = \\beta_0 x + \\beta_1$'})

quadratic = Model(_quadratic, fjacd=_quad_fjd, fjacb=_quad_fjb,
                  estimate=_quad_est, meta={'name': 'Quadratic',
                  'equ': 'y = B_0*x**2 + B_1*x + B_2',
                  'TeXequ': '$y = \\beta_0 x^2 + \\beta_1 x + \\beta_2'})
