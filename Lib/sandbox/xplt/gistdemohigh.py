# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
from GistPlotter import *
GraphicsError = "GraphicsError"
import os
try:
   graphics = os.environ["PYGRAPH"]
except KeyError:
   graphics = "Gist"
 
if graphics [0:3] == "Nar" :
   import NarPlotter
elif graphics == "Gist" :
   import GistPlotter
else :
   raise GraphicsError , \
      graphics + " is an unknown graphics package. Check PYGRAPH " + \
         "environment variable."

def span (lb, ub, n) :
   if n < 2: raise ValueError, "Third argument must be at least 2."
   b = lb
   a = (ub - lb) / (n - 1.0)
   return map(lambda x, A=a, B=b: A*x+B,range(n))

def a3 (lb, ub, n) :
   return reshape (array(n*span(lb,ub,n), Float), (n,n))

def mag ( *args ) :
   r = 0
   for i in range (len (args)) :
      r = r + args[i] * args[i]
   return sqrt (r)

def a2 (lb, ub, n) :
   return reshape (array ((n - 1) * spanz (lb, ub, n),  Float), (n-1, n-1))

def spanz (lb, ub, n) :
   if n < 3 : raise ValueError, "3rd argument must be at least 3"
   c = 0.5 * (ub - lb) / (n - 1.0)
   b = lb + c
   a = (ub -c -b) / (n - 2.0)
   return map (lambda x, A = a, B = b: A * x + B, range (n - 1))

def init ( self ) :
    self.t = 2*pi*arange (400, typecode = Float) / 399.0
    self.na1 = 1
    self.nb1 = 5
    self.na2 = 2
    self.nb2 = 7
    self.rc1 = 40.
    self.rc2 = 160.
    self.size = 40.
    self.phase = self.theta = 0.
    self.dtheta = pi / (self.nsteps - 1)
    self.dphase = 2 * pi / (self.nsteps - 1)

def calc1 ( self ) :
    self.cost = cos (self.theta)
    self.sint = sin (self.theta)
    self.x = self.rc1 * self.cost + self.size * cos (self.na1 * self.t)
    self.y = self.rc1 * self.sint + self.size * sin (self.nb1 * self.t +
       self.phase)

def calc2 ( self ) :
    self.x = self.rc2 * self.cost + self.size * cos (self.na2 * self.t)
    self.y = self.rc2 * self.sint + self.size * sin (self.nb2 * self.t +
       self.phase)

def incr  ( self ) :
    self.theta = self.theta + self.dtheta
    self.phase = self.phase + self.dphase

from curve import Curve
from graph2d import *
def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return
 
print "\nComprehensive graphics test (sort of) for Python interface"
print "with Narcisse and (mainly) Gist. This mirrors Lee's low level"
print "routine. Each frame will be described at the terminal; compare"
print "what you see with the description, then hit <RETURN> to see"
print "the next test, or ^C (ctrl-C) to escape.\n"
print "Type 'demo()' to begin.\n"

from quadmesh import *
from cellarray import *
from lines import *
from animation2d import *


def demo ( ) :
   c1 = Curve ( y = [0,1] , marks = 1 , marker = "A")
   pl = Plotter ( dpi = 75 )
   g1 = Graph2d ( c1, plotter = pl , titles = "Curve marked with A" ,
               title_colors = "blue")
   g1.plot ( )
   print "\nA small (75 dpi) window with line marked A from (0,0) to (1,1)."
   paws ( )
#  pl = Plotter ( dpi = 100 )
#  g1 = Graph2d ( c1, plotter = pl , titles = "Curve marked with A" ,
#              title_colors = "fg")
#  g1.plot ( )
#  print "\nA large (100 dpi) window with line marked A from (0,0) to (1,1)."
#  paws ( )
   c2 = Curve ( y = [1,0] , marks = 1 , marker = "B")
   g1.add (c2)
   g1.change (titles = "New curve marked B.")
   g1.plot ( )
   print "\nAdded line marked B from (0,1) to (1,0) in previous plot."
   paws ( )
   c1.set ( x = [1,2] )
   c2.set ( x = [1,2] )
   g1.change (axis_scales = "loglin", titles = "Same, x axis now logarithmic.")
   g1.plot ( )
   print "\nSame, x axis now logarithmic."
   paws ( )
   g1.change (x_axis_scale = "lin", titles = "x axis back to linear.")
   g1.plot ( )
   print "\nSame, x axis now linear again."
   paws ( )
   g1.change(axis_limits=[[1.2,1.8],[0.2,0.8]],
          titles="Limits now 1.2<x<1.8, 0.2<y<0.8.")
   g1.plot ( )
   print "\nLimits now 1.2<x<1.8, 0.2<y<0.8."
   paws ( )
#  g1.change(y_axis_limits=[0.4,0.6],
#         titles="Limits now 1.2<x<1.8, 0.4<y<0.6")
#  g1.plot ( )
#  print "\nLimits now 1.2<x<1.8, 0.4<y<0.6"
#  paws ( )
   g1.change(axis_limits="defaults",
          titles="Limits now back to extreme values.")
   g1.plot ( )
   print "\nLimits now back to extreme values."
   paws ( )
   x=10*pi*arange(200, typecode = Float)/199.0
   c1 = Curve ( x = x , y = sin(x),marks = 1, marker= "A")
   g1.delete (2)
   g1.delete (1)
   g1.add (c1)
   g1.change (titles = "Five cycles of a sine wave, marked A.")
   g1.plot ( )
   print "\nFive cycles of a sine wave."
   paws ( )
   c1.set(marks = 0, width = 6 , type = "dash")
   g1.change (titles = "A mark gone, bold dashed.")
   g1.plot ( )
   print "\nTurn off A marker, plot curve as bold dashed."
   paws ( )
   x=2*pi*arange(200, typecode = Float)/199.0
   crvs = []
   for i in range (1,7) :
      r = 0.5*i -(5-0.5*i)*cos(x)
      s = `i`
      crvs.append(Curve(y=r*sin(x),x=r*cos(x),marks=0, color=-4-i,label=s))
   g1=Graph2d(crvs,plotter = pl,titles="Nested cardioids in colors")
   g1.plot ( )
   print "\nNested cardioids in primary and secondary colors."
   paws()
   g1.change(titles = ["colors","cardioids","Nested","in many"],
          title_colors = ["red","yellow","magenta","blue"])
   g1.plot ( )
   print "\nThe same with four titles in many colors."
   paws()
   mks = ['.','+','*','o','x','A']
   for i in range (1,7) :
      crvs[i-1].set (color="fg", type = 0, marker = mks[i-1])
   g1.change (titles = "Marked curves", title_colors = "fg")
   g1.plot ( )
   print "\nChanges color to foreground, types to no lines."
   print "Markers are point, plus, asterisk, circle, cross, and A."
   paws()
   for i in range (1,6) :
      crvs[i-1].set (marks = 0, type = i-1)
   crvs[i].set (type = 1, width = 4)
   g1.change (titles = "Different curve types.")
   print "\nChanges line types to solid, dash, dot, dashdot, dashdotdot."
   print "Outermost cardioid becomes a thick solid curve."
   g1.plot()
   paws ()
   ValueError = "ValueError"
   x=a3(-1,1,26)
   y=transpose (x)
   zz=x+1j*y
   zz=5.*zz/(5.+zz*zz)
   xx=zz.real
   yy=zz.imaginary
   q1 = QuadMesh ( x=array(xx, copy = 1), y=array(yy, copy = 1) )
   g1.new(q1, plotter = pl, titles = "Round mesh with bites out of sides",
      axis_limits = [[min(ravel(xx)),max(ravel(xx))],
                     [min(ravel(yy)),max(ravel(yy))]])
   g1.plot()
   print "\nQuadrilateral mesh -- round with bites out of its sides."
   paws (  )
   q1.set (boundary = 1, vx = x+.5, vy=y-.5,ktype = "none", ltype = "none")
   g1.change(titles = "Vector field on mesh.")
   g1.plot ( )
   print "\nVelocity vectors. Try zooming and panning with mouse."
   paws ( )
   q1.set (vx = None, vy = None, z = mag (x+.5,y-.5),marks = 1,
      boundary=1,boundary_type="dash",marker="A")
   g1.change(titles = "Lettered contour plot, dashed boundary.")
   g1.plot ( )
   print "\nLettered contours, with boundary dashed."
   paws ( )
   q1.set (filled=1, type = "dash", color = "white")
   g1.plot ( )
   print "\nFilled lettered dashed contours, with boundary dashed."
   paws ()
   q1.set(marks = 0,boundary_type="solid",levels=array([0.5,1.0,1.5]),
       filled=1, contours = 1, type = "dash", color = "white",
       width = 3)
   g1.change(titles = "Filled mesh, dashed contours")
   g1.plot ( )
   print "\nFilled mesh (color only) with three dashed contours overlaid."
   paws()
   print "\nAfter each new palette is installed, hit <RETURN> for the next."
   pal=["heat.gp","stern.gp","rainbow.gp","gray.gp","yarg.gp","earth.gp"]
   for i in range(6):
      g1.quick_plot(color_card = pal[i],titles = "Palette " + pal[i])
   #  g1.plot ( )
      paws ( )
   pal = ["vg.gs", "boxed.gs", "vgbox.gs", "nobox.gs", "work.gs"]
   q1.set (marks = 1,boundary_type="dash",levels = 8, filled = 0, color = "fg")
   for i in range (5) :
      pl.set_style (pal [i])
##    g1 = Graph2d (q1, plotter = pl, titles = "Style Name: " + pal[i])
      g1.change (plotter = pl, titles = "Style Name: " + pal[i])
      g1.plot ( )
      paws ( )
   x=a3(-6,6,200)
   y=transpose (x)
   r=mag(y,x)
   theta=arctan2(y,x)
   funky=cos(r)**2*cos(3*theta)
   c1 = CellArray(z=funky)
   g1.replace (1,c1)
   g1.change (color_card = "earth.gp",
           titles ="Cell array, three cycles in theta,r",
           axis_limits = "defaults")
   g1.plot()
   paws()

   theta = a2 (0, 2*pi, 18)
   x = cos (theta)
   y = sin (theta)

   ln = Lines (x0 = x, y0 = y, x1 = transpose (x), y1 = transpose (y))
   g1 = Graph2d (ln, plotter = pl, xyequal = 1, titles = "Seventeen pointed star.")
   g1.plot ()
   print "\nAll 17 pointed stars."
   paws ( )

   print "\nNext we will see a lissajous-style animation."
   print "First run without animation mode."
   print "Second run with animation mode."


   anim = Animation2d ( initialize = init, calculations = [calc1, calc2],
                     update = incr, animation = 0, nsteps = 200 )

   g1 = Graph2d ( anim , plotter = pl)
   g1.plot ( )
   paws ( )
   anim.set (animation = 1)
   g1.plot ( )

