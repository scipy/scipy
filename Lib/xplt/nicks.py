# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
from surface import *
from graph3d import *
def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return
from PFB import *
f = PFB ('nicks_plot')
x = f.xout
y = f.yout
z = f.rhoout
s1 = Surface (x = y, y = x, z = z, opt_3d = "s3", mask = "max")
g1 = Graph3d (s1, theta = 11, phi = 45, titles = "Density plot: reflecting shock wave")
g1.plot ( )
paws ( )
g1.quick_plot (surface = 1, opt_3d = "w3")
paws ( )
vx = f.vxout
vy = f.vyout
massfrac = f.massfrac1out
 
f.close ()
s1.set (z = vx, mask = "max", opt_3d = "w3")
g1.change_plot (send = 1, theta = 160, phi = 45, titles = "x velocity: reflecting shock wave")
paws ()
s1.set (z = vy, mask = "max", opt_3d = "w3")
g1.change_plot (send = 1, theta = 11, phi = 45, titles = "y velocity: reflecting shock wave")
paws ()
s1.set (z = massfrac)
g1.change_plot (send = 1, theta = 170, phi = 45, titles = "Mass fraction: reflecting shock wave")
paws ()
f = PFB ('mikes_plot')
x8 = f.x8
y8 = f.y8
z8 = f.z8
s2 = Surface (x = x8, y = y8, z = z8, opt_3d = "s3", mask = "max",
     color_card = "redgreen")
g2 = Graph3d (s2, phi = 26, theta = 80, titles = "3-D square mode")
g2.plot ( )
f.close ( )
