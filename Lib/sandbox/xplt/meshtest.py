# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
import NarPlotter
from MLab import rand
from graph3d import *
from surface import *

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return
 
print "This is a test of the Python interface to the Limeil Lab graphics"
print "package, Narcisse. You  need Narcisse to be running. Fire it up by"
print "doing setenv PORT_SERVEUR 0, typing  /dist/basis/Narcisse/bin/Narcisse &,"
print "and then do another senetv PORT_SERVEUR to the port number which"
print "appears in the Narcisse GUI."
 
paws ()
 
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
ut = rt/sqrt (rt ** 2 + zt ** 2)
vt = zt/sqrt (rt ** 2 + zt ** 2)
ireg =  multiply.outer ( ones (kmax), ones (lmax))
ireg [0:1, 0:lmax]=0
ireg [0:kmax, 0:1]=0
ireg [1:15, 7:12]=2
ireg [1:15, 12:lmax]=3
ireg [3:7, 3:7]=0
freg=ireg + .2 * (1. - rand (kmax, lmax))  # ranf (ireg))
freg=array (freg, Float)
#rt [4:6, 4:6] = -1.e8
z [3:10, 3:12] = z [3:10, 3:12] * .9
z [5, 5] = z [5, 5] * .9
z [17:22, 15:18] = z [17:22, 15:18] * 1.2
z [16, 16] = z [16, 16] * 1.1
 
s1 = Surface (x = zt, y = rt, z = freg, mask = "max")
g1 = Graph3d (s1, titles = "Rotating mesh", theta = 45., phi = 5., roll = 0.0)
g1.plot ( )
for i in range (6) :
   g1.quick_plot (roll = (i + 1) * 22.5)
paws ()
g1.quick_plot (surface = 1, opt_3d = "i3", titles = ["(1) contour lines" ,
        "Different options of viewing mesh"], text = "Region 2",
        text_color = "black", text_size = 20, text_pos = [.20, .35])
paws ( )
g1.quick_plot (surface = 1, opt_3d = "s3", titles = ["(2) smooth coloring",
        "Different options of viewing mesh"])
paws ( )
g1.quick_plot (surface = 1, opt_3d = "w3", titles = ["(3) colored mesh",
        "Different options of viewing mesh"])
paws ( )
s2 = Surface (x = zt + vt, y = rt + vt, z = freg + ut,
              mask = "sort", opt_3d = "w3")
g1 = Graph3d ( [s1, s2], connect = 1, text = " ", roll = 0.0,
               titles = "Mesh follows velocity field one time step")
g1.plot ( )
paws ( )
# The following is good but comment it out if you want to save time.
for i in range (8) :
   g1.quick_plot (titles = ["Mesh follows velocity field one time step",
         "(Animation)"], phi = 5. + 10. * (i + 1))
   paws ()
for i in range (4) :
   g1.quick_plot (theta = 45. - 10. * (i + 1))
   paws ( )
