# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.
# I've felt the need for such a test for a long time;
# this tells you whether an item is a scalar or not.

from types import *
from Numeric import *

def is_scalar (x) :
   if type (x) == StringType : return 1
   try :
      y = len (x)
   except (TypeError, AttributeError) :
      return 1
   return 0

# This routine should be able to tell you the size of any object:
def no_of_dims (x) :
   if x == None : return 0
   if (type (x) == ArrayType) : return len (x.shape)
   if (type (x) == ListType or type (x) == TupleType) : return 1
   # I don't know if there are any other possibilities.
   for i in range (10) :
      if is_scalar (x) : return i
      try :
         x = x[0]
      except IndexError, TypeError :
         return i
   return -1

# This routine should define rshape for odd objects:
# It returns a list rather than a tuple; and if it is
# applied to a list, it returns a list each entry of which
# tells the shape of the corresponding list entry, recursively.
def rshape (x) :
   ndim = no_of_dims (x)
   if ndim == 0 : return 0
   s = []
   if type (x) == ListType or type (x) == TupleType :
      for i in range (len (x)) :
         s.append (rshape (x [0]))
         x = x [1:]
   else : # an array or ???
      for i in range ( ndim ) :
         s.append (len (x))
         x = x[0]
   return s
 
