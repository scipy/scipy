# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.
from Numeric import *
import os
import time
from shapetest import *

_ZcenError = "ZcenError"

def zcen_ (x, i = 0) :
#  zcen_(x, i) does the same thing as in Yorick: x(...,zcen_,...)
#  where zcen_ is the ith subscript. (works for up to 5 dimensions).

   if is_scalar (x) :
      raise _ZcenError, "zcen_ must be called with an array."
   dims = shape (x)
   ndims = len (dims)
   if i < 0 or i > ndims - 1 :
      raise _ZcenError, "i <" + `i+1` + \
         "> is out of the range of x's dimensions<" + `ndims` +"."
   if i == 0 :
      newx = (x [0:dims [0]-1] + x [1:dims [0]]) /2.0
   elif i == 1 :
      newx = (x [:, 0:dims [1]-1] + x[:, 1:dims [1]]) / 2.0
   elif i == 2 :
      newx = (x [:, :, 0:dims [2]-1] + x[:, :, 1:dims [2]]) / 2.0
   elif i == 3 :
      newx = (x [:, :, :, 0:dims [3]-1] + x[:, :, :, 1:dims [3]]) / 2.0
   elif i == 4 :
      newx = (x [:, :, :, :, 0:dims [4]-1] + \
              x [:, :, :, :, 0:dims [4]]) / 2.0

   return newx

_DifError = "DifError"

def dif_ (x, i = 0) :
#  dif_(x, i) does the same thing as in Yorick: x(...,dif_,...)
#  where dif_ is the ith subscript. (works for up to 5 dimensions).

   if is_scalar (x) :
      raise _DifError, "dif_ must be called with an array."
   dims = shape (x)
   ndims = len (dims)
   if i < 0 or i > ndims - 1 :
      raise _DifError, "i <" + `i+1` + \
         "> is out of the range of x's dimensions <" + `ndims` +">."
   if i == 0 :
      newx = x [1:dims [0]] - x [0:dims [0] - 1]
   elif i == 1 :
      newx = x [:, 1:dims [1]] - x[:, 0:dims [1] - 1]
   elif i == 2 :
      newx = x [:, :, 1:dims [2]] - x [:, :, 0:dims [2] - 1]
   elif i == 3 :
      newx = x [:, :, :, 1:dims [3]] - x [:, :, :, 0:dims [3] - 1]
   elif i == 4 :
      newx = x [:, :, :, :, 1:dims [4]] - x [:, :, :, :, 0:dims [4] - 1]
   return newx

def maxelt_ (*x) :
#  maxelt_ accepts a sequence of one or more possible multi-dimensional
#  objects and computes their maximum. In principle these can be of
#  arbitrary complexity, since the routine recurses. 
   if len (x) == 0 :
      return None
   elif len (x) == 1 :
      z = x [0]
      if is_scalar (z) :
         return z
      if len (shape (z)) >= 1 :
         zz = array (z)
         return maximum.reduce (ravel (zz))
   else :
      maxelt = maxelt_ (x [0])
      for i in range (1, len (x)) :
         maxelt = max (maxelt, maxelt_ (x [i]))
      return maxelt

def minelt_ (*x) :
#  minelt_ accepts a sequence of one or more possible multi-dimensional
#  objects and computes their minimum. In principle these can be of
#  arbitrary complexity, since the routine recurses.
   if len (x) == 0 :
      return None
   elif len (x) == 1 :
      z = x [0]
      if is_scalar (z) :
         return z
      if len (shape (z)) >= 1 :
         zz = array (z)
         return minimum.reduce (ravel (zz))
   else :
      minelt = minelt_ (x [0])
      for i in range (1, len (x)) :
         minelt = min (minelt, minelt_ (x [i]))
      return minelt

def rem_0_ (z) :
#   rem_0_ (z) goes through array z and replaces any zero
#  elements with 1.e-35. Assumes z has one or two dimensions.
   if len (shape (z)) == 1 :
      for i in range (len (z)) :
         z [i] = z [i] + (z [i] == 0.0) * 1.e-35
   elif len (shape (z)) == 2 :
      (k, l) = shape (z)
      for i in range (k) :
         for j in range (l) :
            z [i] [j] = z [i] [j] + (z [i] [j] == 0.0) * 1.e-35

def avg_ (z) :
#  avg_ (z) returns the average of all elements of its array
#  argument.
   zz = array (z, copy = 1 )
   return add.reduce (ravel (zz)) / len (ravel (zz))

def sign_ (x) :
   if x >= 0 :
      return (1)
   else :
      return (- 1)

def timer_ (elapsed, *split) :
#   timer (elapsed) returns a triple consisting of the times
#  [cpu, system, wall]. 
#  timer (elapsed, split) returns a sequence whose first element
#  is [cpu, system, wall] and whose second element is the
#  sum of split and the difference between ththe new and old values
#  of 'elapsed.'

   stime = os.times ( )
   wtime = time.time ( )
   retval = array ( [stime [0], stime [1], wtime], Float )
   if len (split) == 0 :
      return retval
   else :
      return [retval, split [0] + retval - elapsed]

def timer_print (label, split, *other_args) :
#  timer_print (label1, split1 [,label2, split2, ...]) prints
#  out a timing summary for splits accumulated by timer.

   print label, split
   i = 0
   while i < len (other_args) :
      print other_args [i], other_args [i + 1]
      i = i + 2
