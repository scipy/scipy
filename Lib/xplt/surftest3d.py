# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *

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

if graphics [0:3] == "Nar" :
   print "This is a test of the Python interface to the Limeil Lab graphics"
   print "package, Narcisse. You  need Narcisse to be running. Fire it up by"
   print "doing setenv PORT_SERVEUR 0, typing  /dist/basis/Narcisse/bin/Narcisse &,"
   print "and then do another senetv PORT_SERVEUR to the port number which"
   print "appears in the Narcisse GUI."
else :
   print "This is a test of the Python interface to the Gist graphics package."

print "Invoke function demo () to run."

from plane import *
from surface import Surface
from graph3d import Graph3d
from mesh3d import *
from MLab import rand

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

def demo () :
#  vsf = 0.
#  c = 1
   s = 1000.
   kmax = 25
   lmax = 35
   # The following computations define an interesting 3d surface.
   xr = multiply.outer (
      arange (1, kmax + 1, typecode = Float), ones (lmax, Float))
   yr = multiply.outer (
      ones (kmax, Float), arange (1, lmax + 1, typecode = Float))
   zt = 5. + xr + .2 * rand (kmax, lmax)   # ranf (xr)
   rt = 100. + yr + .2 * rand (kmax, lmax)   # ranf (yr)
   z = s * (rt + zt)
   z = z + .02 * z * rand (kmax, lmax)   # ranf (z)
   ut = rt / sqrt (rt ** 2 + zt ** 2)
   vt = zt / sqrt (rt ** 2 + zt ** 2)
   ireg =  multiply.outer ( ones (kmax), ones (lmax))
   ireg [0:1, 0:lmax] = 0
   ireg [0:kmax, 0:1] = 0
   ireg [1:15, 7:12] = 2
   ireg [1:15, 12:lmax] = 3
   ireg [3:7, 3:7] = 0
   freg = ireg + .2 * (1. - rand (kmax, lmax))  # ranf (ireg))
   freg = array (freg, Float)
   #rt [4:6, 4:6] = -1.e8
   z [3:10, 3:12] = z [3:10, 3:12] * .9
   z [5, 5] = z [5, 5] * .9
   z [17:22, 15:18] = z [17:22, 15:18] * 1.2
   z [16, 16] = z [16, 16] * 1.1
   
   s1 = Surface (z = z, mask = "max", opt_3d = ["wm", "i3"])
   g1 = Graph3d ( s1 , titles = "Surface with contour lines",
                  xyequal = 1., hardcopy = "talk.ps",
                  theta = 45., phi = 10., roll = 0.)
   g1.plot ( )
   paws ( )
   s1.set (opt_3d = ["wm", "f3"])
   g1.change (titles = "Flat mode")
   g1.plot ()
   paws ()
   s1.set (opt_3d = ["wm", "s3"])
   g1.change (titles = "Smooth mode")
   g1.plot ()
   paws ()
   s1.new (x = xr, y = yr, z = z - z [kmax/2, lmax/2],
           mask = "max", opt_3d = ["wm" , "i3"])
   s11 = sslice (s1, Plane (array ( [0., 0., 1.]), array ( [0., 0., 0.])))
   g1.replace (1, s11)
## g1.change ( titles = "Part of surface above xy plane", phi = 30.,
##          y_factor = 2.0,
##          axis_limits = [[0., 0.],[0., 0.], [0., 100000.]])
   g1.change ( titles = "Part of surface above xy plane",
            y_factor = 2.0,
            axis_limits = [[0., 0.],[0., 0.], [0., 100000.]])
   g1.plot ( )
   paws ( )
   s1.set (opt_3d = ["wm", "f3"])
   s11 = sslice (s1, Plane (array ( [0., 0., 1.]), array ( [0., 0., 0.])))
   g1.replace (1, s11)
   g1.change (titles = "Flat mode")
   g1.plot ()
   paws ()
   s1.set (opt_3d = ["wm", "s3"])
   s11 = sslice (s1, Plane (array ( [0., 0., 1.]), array ( [0., 0., 0.])))
   g1.replace (1, s11)
   g1.change (titles = "Smooth mode")
   g1.plot ()
   paws ()
   s1.set (z = z, c = freg, opt_3d = ["wm", "s4"])
## s1 = Surface (x = xr, y = yr, z = z, c = freg, opt_3d = ["wm", "s4"])
   g1.replace (1, s1)
   g1.change ( titles = "Surface colored by mesh values", phi = 20.,
           xyequal = 1,
           axis_limits = [[0., 0.], [0., 0.], [0., 0.], [0., 0.]])
   g1.plot ( )
   paws ( )
   xr1 = zeros ( (16, 7), Float)
   yr1 = zeros ( (16, 7), Float)
   z1 = zeros ( (16, 7), Float)
   zs1 = zeros ( (16, 7), Float)
   for i in range (16) :
       for j in range (7) :
           xr1 [i, j] = xr [i, j+6]
           yr1 [i, j] = yr [i, j+6]
           z1 [i, j] = z [i, j+6]
           zs1 [i, j] = freg [i, j+6]
   s1.set (x = xr1, y = yr1, z = z1, c = zs1)
   g1.change ( titles = "Region 2 colored by mesh values", phi = 10.)
   g1.plot ( )
   paws ( )
   s1.new (x = xr1, y = yr1, z = z1, opt_3d = ["wm", "s3"], mask = "max")
   g1.change ( titles = "Region 2 with mesh and contours", phi = 10.)
   g1.plot ( )
   paws ()
   s1.set (opt_3d = "s3")
   g1.change ( titles = "Region 2 with contours alone")
   g1.plot ()
   paws ( )
   zs1 = zeros ( (16, lmax - 7), Float)
   for i in range (16) :
       for j in range (lmax - 7) :
           zs1 [i, j] = z [i, j+6]
   s1.new (z = zs1, opt_3d = ["wm", "s3"], mask = "max")
   g1.change ( titles = "Regions 2 and 3, mesh and contours", 
            theta = 70., phi = 10., roll = 0.)
   g1.plot ( )
   paws ( )
   s1.set (opt_3d = "s3")
   g1.change ( titles = "Regions 2 and 3, contours alone")
   g1.plot ( )
