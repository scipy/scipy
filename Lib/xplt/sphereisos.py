# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
from PR import *
from mesh3d import *
from graph3d import *
from plane import *
from GistPlotter import *

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

print "Test of slices through an imploding sphere. Type demo () to see the test."
paws ()

def demo () :
   f = PR ('./bills_plot')
   n_nodes = f.NumNodes
   n_z = f.NodesOnZones
   x = f.XNodeCoords
   y = f.YNodeCoords
   z = f.ZNodeCoords
   c = f.ZNodeVelocity
   n_zones = f.NumZones
   # Put vertices in right order for Gist
## n_z = transpose (
##    take (transpose (n_z), array ( [0, 1, 3, 2, 4, 5, 7, 6])))
   
   m1 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1, hex = [n_zones, n_z])
   
   s1 = sslice (m1, .9 * max (c), varno = 1)
   s11 = sslice (m1, -1.8, varno = 1)
## pyz = Plane (array ([1., 0., 0.], Float ), zeros (3, Float))
   pyz = Plane (array ([1., 0., 0.], Float ), array ( [0.0001, 0., 0.], Float))
   pxz = Plane (array ([0., 1., 0.], Float ), array ( [0., 0.0001, 0.], Float))
   p2 = Plane (array ([1., 0., 0.], Float ), array ( [0.35, 0., 0.], Float))
   
   s2 = sslice (m1, pyz, varno = 1, opt_3d = ["wm", "s4"])
   s22 = sslice (m1, p2, varno = 1, opt_3d = ["wm", "s4"])
   s23 = sslice (m1, pxz, varno = 1, opt_3d = ["wm", "s4"])
   
## pl = GistPlotter.Plotter (" ", hcp = "talk.ps", dump = 1, dpi = 75)
## g1 = Graph3d( [s2, s22, s23], plotter = pl, color_card = "rainbowhls",
##      opt_3d = ["wm", "s4"], mask = "min", color_bar = 1, split = 0)
   g1 = Graph3d( [s2, s22, s23], color_card = "rainbowhls",
        opt_3d = ["wm", "s4"], mask = "min", color_bar = 1, split = 0,
        hardcopy = "talk.ps")
   g1.plot ()
   paws ()
   g1.change (split = 1)
   g1.plot ()
   paws ()
   g1.set_surface_list ( [s11, s2, s22, s23])
   g1.plot ()
   paws ()
   s2 = sslice (m1, .9 * min (c), varno = 1, opt_3d = "w3")
   s3 = sslice (m1, .5 * min (c), varno = 1, opt_3d = "w3")
   s4 = sslice (s2, -pyz, opt_3d = "w3")
   s5 = sslice (m1, .5 * max (c), varno = 1, opt_3d = "w3")
   s6 = sslice (s5, -pyz, opt_3d = "w3")
   g1.set_surface_list ( [s1, s2, s6])
   g1.plot ()
   paws ()
   # Do multiple slices
   slice_list = []
   for i in range (8) :
      sl = sslice (m1, .9 * min (c) + i * (.9 * max (c) - .9 * min (c)) / 8., 
         varno = 1, opt_3d = "none")
      slice_list.append (sslice (sl, pxz))
   g1.set_surface_list ( slice_list)
   g1.plot ()
   paws ()
