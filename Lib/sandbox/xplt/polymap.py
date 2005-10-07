# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *

class Polymap :
   """This class is for Gist and possibly other graphics packages;
   it creates a class consisting of filled polygons, which can then
   be given to and plotted by a graph2d object.

   x = Polymap ( <keyword arguments> ) will create a Polymap object
   with coordinates and various other characteristics. The keyword
   arguments are:
      x = <sequence of floating point values>
      y = <sequence of floating point values>
          These are the coordinates of the vertices of the
          polygons. The way this data is sewt up, vertices of
          adjacent polygons will be repeated.
      n = <sequence of integer values>
          The entries in this array tell you how many vertices
          each polygon has. The first n[0] entries in x and y define
          the first polygon, the next n[1] entries define the second,
          etc. The sum of all the entried in n is the length of
          vectors x and y.
      z = <sequence of numerical or character values> (this vector
          is the same length as n) tells how to denote the polygons.
          numbers are interpolated into a palette, characters are used
          as is.
      hide = 0/1 (1 to hide this part of the graph)
      label = <string> label for this part of the graph.
   """

   def type (self) :
       return PolyMapType

   _PolyMapSpecError = "PolyMapSpecError"

   def __init__ ( self , * kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       if not keywords.has_key ("x") or not keywords.has_key ("y") or \
          not keywords.has_key ("z") or not keywords.has_key ("n" ) :
          raise self._PolyMapSpecError , \
             "You need all four keywords x, y, z, and n."
       self.x = keywords ["x"]
       self.y = keywords ["y"]
       self.z = keywords ["z"]
       self.n = keywords ["n"]
       total_size = 0
       if len (self.z) != len (self.n) :
          raise self._PolyMapSpecError , \
             "z and n are not the same size (" + `len (self.z)` + " vs " + \
              `len (self.n)` + ")."
       for i in range ( len(self.n) ) :
          total_size = total_size + self.n [i]
       if len (self.x) != total_size or len (self.y) != total_size :
          raise self._PolyMapSpecError , \
             "Either x or y or both are not the right size (" + \
             `total_size` + ")."
       if keywords.has_key ("hide") :
          self.hide = keywords ["hide"]
       else:
          self.hide = 0
       if keywords.has_key ("label") :
          self.label = keywords ["label"]
       else:
          self.label = " "

   def new ( self, ** keywords ) :
       """ new (...keyword arguments...) allows you to reuse a
       previously existing curve.
       """
       self.__init__ ( keywords )

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       curve characteristics. No error checking is done.
       """
       for k in keywords.keys ():
           if k == "type" :
              self.line_type = keywords ["type"]
           else :
              setattr (self, k, keywords [k])
