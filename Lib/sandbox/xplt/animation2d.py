# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
from shapetest import *
from graftypes import *

class Animation2d :
   """
   anim = Animation2d ( <keyword arguments> ) will create a class which
   supplies the controls for an animation. For the most part, the
   keyword arguments are the names of functions which perform various
   parts of the calculations for the animation, under control of a
   Plotter object. Here the keywords are:
     initialize = <name of an initialization function>. This function
        should have one argument, the name of the instantiation 'anim',
        and when called should initialize any internal variables needed
        before beginning to comput the animation.
     calculations = <name or list of names of calculations of coordinates>
        this routine (or these routines) are called from within a loop
        in the Plotter(s) associated with anim. Each of the calculations
        routines should have 'anim' as the argument. They should compute
        the current values of anim.x and anim.y, the coordinates of the
        curve(s) in this step of the animation. They start with whatever
        was 'initialize'd, then in subsequent calls, use 'update'd values.
        If more than one calculation is specified, then a plot command
        will be issued after each one.
     update = <name of routine to update the variables used in calculations>.
        This updates (increments, decrements) variables used in calculating
        the frames. Like the other functions, has 'anim' as its sole
        argument.
   The following keywords are scalar values:
     animation = 0/1 (If 1, supplies a smoother animation for, e. g.,
        making a movie. default value 1.)
     nsteps = number of animation steps desired. Default: 100.
     color = <value> where <value> is an integer from 0 to 63
        representing an entry in a color chart, or a
        common color name like "red", "blue", "background", etc.
     In the interest of speed, other keywords relating to curve type,
     thickness, etc., are currently not allowed.
   """

   def type (self) :
      return Animation2dType

   _AnimationError = "AnimationError"

   def __init__ ( self, * kwds , ** keywords ) :
       if len (kwds) == 1 :
          keywords = kwds[0]
       self.line_type = "unspecified"
       if keywords.has_key ( "initialize" ) :
          self.initialize = keywords ["initialize"]
       else :
          raise self._AnimationError , \
             "Animation2d: keyword 'initialize' missing."
       if keywords.has_key ( "calculations" ) :
          self.calculations = keywords ["calculations"]
          if is_scalar (self.calculations) :
             self.calculations = [self.calculations]
       else :
          raise self._AnimationError , \
             "Animation2d: keyword 'calculations' missing."
       if keywords.has_key ( "update" ) :
          self.update = keywords ["update"]
       else :
          raise self._AnimationError , \
             "Animation2d: keyword 'update' missing."
       if keywords.has_key ( "animation" ) :
          self.animation = keywords ["animation"]
       else :
          self.animation = 1
       if keywords.has_key ( "nsteps" ) :
          self.nsteps = keywords ["nsteps"]
       else :
          self.nsteps = 100
       if keywords.has_key ( "color" ) :
          self.color = keywords ["color"]
       else :
          self.color = "fg"

   def new ( self, ** keywords ) :
       """ new (...keyword arguments...) allows you to reuse a
       previously existing animation.
       """
       del self.color, self.nsteps, self.animation, self.calculations, \
           self.initialize, self.update
       self.__init__ ( keywords )

   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       animation characteristics. No error checking is done.
       """
       for k in keywords.keys ():
           if k == "type" :
              self.line_type = keywords ["type"]
           else :
              setattr (self, k, keywords [k])
