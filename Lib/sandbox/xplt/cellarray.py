# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *

class CellArray :
   """This class is for Gist and possibly other graphics packages;
   it creates a class consisting of a Cellarray,  which can then
   be given to and plotted by a graph2d object.

   x = CellArray ( <keyword arguments> ) will create a CellArray
   object. The keyword arguments are:
      z = <2d sequence of numeric or character values>
      specifies the coloring to be given to the CellArray. If
      numeric, the values are interpolated onto a palette. If
      character, the values are used as is, whatever that means.
      THE ARGUMENT z IS REQUIRED.
      x0, y0 -- floating point scalars; if present, the coordinates
      of the lower left corner of the cell array. The default is (0.,0.).
      These coordinates are optional, but if they are present then
      x1, y1 must be also (see below).
      x1, y1 -- floating point scalars; if present, the coordinates
      of the upper right corner of the cell array. If these optional
      keywords are missing, then x0, y0 must be also missing, and the
      call array will be plotted in a 1X1 square with lower left
      corner at the origin.
      hide = 0/1 (1 to hide this part of the graph)
      label = <string> label for this part of the graph.
   """

   _CellArraySpecError = "CellArraySpecError"

   def type (self) :
       return CellArrayType

   def __init__ ( self , * kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       self.line_type = "unspecified"
       if not keywords.has_key ("z") :
          raise _CellArraySpecError , \
             "CellArray requires z keyword argument."
       self.z = keywords ["z"]
       if keywords.has_key ("x1") and ( not keywords.has_key ("x0") or
          not keywords.has_key ("y0") ) :
          raise _CellArraySpecError , \
             "CellArray: if x1, y1 are present then x0, y0 must be also."
       if keywords.has_key ("x1") :
          self.x1 = keywords ["x1"]
          if not keywords.has_key ("y1") :
             raise _CellArraySpecError , \
                "CellArray: if x1 is present then y1 must be also."
          self.y1 = keywords ["y1"]
       else :
          self.x1 = None
          self.y1 = None
       if keywords.has_key ("x0") :
          self.x0 = keywords ["x0"]
          if not keywords.has_key ("y0") :
             raise _CellArraySpecError , \
                "CellArray: if x0 is present then y0 must be also."
          self.y0 = keywords ["y0"]
       else :
          self.x0 = None
          self.y0 = None
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
       previously existing cellarray.
       """
       self.__init__ ( keywords )

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       cellarray characteristics. No error checking is done.
       """
       for k in keywords.keys ():
           if k == "type" :
              self.line_type = keywords ["type"]
           else :
              setattr (self, k, keywords [k])
