# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from slice3 import plane3
from Numeric import *
from scipy_base.fastumath import *
from graftypes import *

class Plane :

   """Plane (normal, point) simply embodies a plane through the
   given point with the given normal.
   """

   def type (self) :
       return PlaneType

   def __init__ (self, normal = array ([1., 0., 0.]),
                 point = array ([0., 0., 0.])) :
       self.coeffs = plane3 (normal, point)

   def __repr__ (self) :
       return `self.coeffs`

   def __str__ (self) :
       return `self.coeffs`

   def __neg__ (self) :
       p = Plane ()
       p.coeffs = - self.coeffs
       return p

   def rep (self) :
       return self.coeffs

   def astype (self, newtype) :
       return (array (self.coeffs, newtype))
