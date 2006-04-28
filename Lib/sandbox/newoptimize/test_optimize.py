import optimize as opt
from base import *

def f( x ):
    return nm.sum( x ** 4 )

def fg( x ):
    return 4 * (x ** 3)

conf = Struct( epsRD = 1e-3,
               epsOF = 1e-4,
               epsOFG = 1e-4,
               iMax = 20,
               norm = nm.Inf,
               ls = True, # Linesearch.
               log = True,
               check = 0,
               delta = 1e-6 )

x0 = nm.arange( 100, dtype = nm.float64 )
x, log = opt.fmin_sd( conf, x0, f, fg )
print x
print log
