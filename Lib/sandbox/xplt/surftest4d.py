# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
from MLab import rand
from surface import *
from graph3d import *
from mesh3d import *
def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

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

def demo () :
   vsf = 0.
   c = 1
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
    
   # Sombrero function
   x = arange (-20, 21, typecode = Float)
   y = arange (-20, 21, typecode = Float)
   z = zeros ( (41, 41), Float)
   for i in range (0, 41):
       for j in range (0, 41):
           r = sqrt (x [i] ** 2 + y [j] ** 2) + 1e-6
           z [i, j] = sin (r) / r
   s1 = Surface (z = z, opt_3d = "s3", mask = "max")
   g1 = Graph3d ( s1 , titles = "Sombrero function",
                  y_factor = 2.0) #, phi = 30., theta = 45.)
   g1.plot ( )
   paws ( )
   s1.set (opt_3d = "w3")
   g1.plot ()
   paws ( )
   foo = zeros ( (kmax, lmax), Float)
   for k in range (kmax) :
      foo [k, 0:lmax] = log (arange (1, lmax + 1, typecode = Float))
      xxx = exp (k / 5.)
      for l in range (lmax) :
         foo [k, l] = xxx * foo [k, l]
   s1.new (x = xr, y = yr, z = foo, opt_3d = "w3")
   g1.change ( titles = "Exponential surface, z logarithmic", 
            # phi = 30., theta = 45.,
            axis_scales = ["lin", "lin", "log"],
            z_contours_scale = "log", y_factor = 1.0)
   g1.plot ( )
   paws ( )
   xr = array (xr)
   yr = array (yr)
   zz = xr ** 2 + yr ** 2
   g1.change ( titles = "Exponential surface, z linear", phi = 45.,
            axis_scales = ["lin", "lin", "lin"] )
   g1.plot ( )
   paws ( )
   s1.new (z = zz, opt_3d = "w3")
   g1.change ( z_contours_scale = "lin",
            titles = "Graph of xr**2 + yr**2", phi = 35.)
   g1.plot ( )
   paws ( )
   zs = zeros ( (50, 25), Float)
   xs = zeros ( (50, 25), Float)
   ys = zeros ( (50, 25), Float)
   ctheta = zeros (50, Float)
   stheta = zeros (50, Float)
   cphi = zeros (25, Float)
   sphi = zeros (25, Float)
   pi = 3.1415926535
   for i in range (50):
      ctheta [i] = cos (2. * i / 49. * pi)
      stheta [i] = sin (2. * i / 49. * pi)
   for i in range (25):
      cphi [i] = cos (pi * i / 24.)
      sphi [i] = sin (pi * i / 24.)
   for i in range (25):
      zs [0:50, i:i + 1] = cphi [i]
      for j in range (50):
         xs [j, i] = sphi [i] * ctheta [j]
         ys [j, i] = sphi [i] * stheta [j]
   hx = xs  / 1.6
   hy = ys / 1.6
   hz = zs / 1.6
   # Plot a sphere with a random distribution of colors
   s1.new (mask = "sort", opt_3d = "s4", x = hx, y = hy, z = zs,
           c = rand (50, 25)) # ranf (zs))
   g1.replace ( 1, s1 )
   g1.change_plot (titles = "Randomly colored sphere", send = 1)
   paws ()
   s1 = Surface (mask = "sort", opt_3d = "s4", x = hx, y = hy,
                 z = zs, c = rand (50, 25))
   ##m1 = Mesh3d (color_card = "random", opt_3d = "f4", mask = "min",
   ##             x = array ( [0.5, 1.2, 1.4, 1.05], Float),
   ##             y = array ( [1.2, 0.5, 1.4, 1.05], Float),
   ##             z = array ( [0., 0., 0., 1.2], Float),
   ##             c = array ( [0.4, 0.6, 0.6, 0.6], Float),
   ##             avs = 1, tet = [ 1, array ([[3,0,1,2]],Int) ])
   m1 = Slice (array ( [3, 3, 3, 3], Int),
               array ( [
                        [0.5, 1.2, 0.], [1.2, 0.5, 0.], [1.4, 1.4, 0.],
                        [0.5, 1.2, 0.], [1.4, 1.4, 0.], [1.05, 1.05, 1.2],
                        [0.5, 1.2, 0.], [1.05, 1.05, 1.2], [1.2, 0.5, 0.],
                        [1.2, 0.5, 0.], [1.05, 1.05, 1.2], [1.4, 1.4, 0.]
                     ] ), array ( [0.4, 0.5, 0.6, 0.7], Float))
   g1.new ([s1, m1], link = 1,
                   titles = "Randomly colored sphere and solid tetrahedron",
                   axis_limits = [[-0.62, 1.4],
                   [-0.62, 1.4], [-1.0, 1.5]],
                   send = 1, y_factor = 2.0) # phi = 55, theta = 30,
   g1.plot ( )
   #pl.set_color_card (8) #random
   #pl.set_3d_options ("f4") #flat 4d
   #pl.set_grid_type ("none") #don't redraw axes
   # Note to myself: When the mesh plotter is done, use it for the
   # following. Python should never call a plotter directly.
   #pl.plot_surface (array ( [0.5, 1.2, 1.4, 1.05], Float),
   #                 array ( [1.2, 0.5, 1.4, 1.05], Float),
   #                 array ( [0., 0., 0., 1.2], Float),
   #                 array ( [0.4, 0.6, 0.6, 0.6], Float),
   #                 array ( [3, 6, 9, 12, 1, 2, 3, 0, 2,
   #                          3, 0, 1, 3, 0, 1, 2], Int), 4)
   #pl.send_graph ()
   paws ()
   s1 = Surface ( mask = "sort", opt_3d = "s4", x = xs, y = ys,
                  z = zs + 2., c = zs)
   s2 = Surface ( mask = "sort", opt_3d = "wm", x = hx, y = hy, z = hz, c = zs)
   # g1.replace (1, s1)
   # g1.add ( s2 )
   # g1.change_plot ( link = 1 , titles = "Imploding Sphere in R3",
   #                  phi = 70.0, theta = 30.0,
   #                  axis_limits  =  [[-1., 1.], [-1., 1.], [-0.7, 3.0]])
   g1.new ( [s1, s2], link = 1 , titles = "Two Spheres in R3",
                      axis_limits  =  [[-1., 1.], [-1., 1.], [-0.7, 3.0]],
                      x_factor = 1.7)
   ##                 phi = 70.0, theta = 30.0,
   g1.plot ( )
