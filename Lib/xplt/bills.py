# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from PR import *
from Numeric import *
from mesh3d import Mesh3d
from graph3d import Graph3d

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return
f = PR ('./bills_plot')

n_nodes = f.NumNodes
n_z = f.NodesOnZones
x = f.XNodeCoords
y = f.YNodeCoords
z = f.ZNodeCoords
c = f.ZNodeVelocity
n_zones = f.NumZones

m1 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1, hex = [n_zones, n_z],
           opt_3d = "w4", mask = "none")
g1 = Graph3d (m1,
     titles = ["Vertical component of velocity","Imploding Sphere"])
g1.plot ( )
paws ( )

g1.quick_plot (surface = 1, mask = "max")
paws ( )

# Now we're going to plot it with all cells missing which have an
# x coordinate greater than 0.0.
n_zones_used = 0
zones_not_used = [] # subscripts of zones to ignore
zones_used = []
for i in range (n_zones) :
   nz = n_z [i]
   used = 1
   for j in range (8) :
      if x [nz [j]] > 0.005 :
         zones_not_used.append (i)
         used = 0
         break
   if used == 1 :
      zones_used.append (i)
      n_zones_used = n_zones_used + 1
         
new_n_z = zeros ( (n_zones_used, 8), Int32)
for i in range (n_zones_used) :
   new_n_z [i] = n_z [zones_used [i]]

m1.new (x = x, y = y, z = z, c = c, avs = 1, hex = [n_zones_used, new_n_z],
        mask = "max", opt_3d = "w4")
# Uncomment below when we take the front face away
g2 = Graph3d (m1,
     titles = ["Vertical component of velocity","Imploding Sphere"])
g2.plot ( )
paws ( )

g2.quick_plot (surface = 1, opt_3d = ["wm", "f4"], mask = "max")
paws ( )

g2.quick_plot (surface = 1, opt_3d = ["s4", "wm"])
paws ( )

f.close ( )

# Now the general mesh....

f = PR ("ball.s0001")
ZLss = f.ZLstruct_shapesize
ZLsc = f.ZLstruct_shapecnt
ZLsn = f.ZLstruct_nodelist
x = f.sap_mesh_coord0
y = f.sap_mesh_coord1
z = f.sap_mesh_coord2
c = f.W_vel_data
# Now we need to convert this information to avs-style data
istart = 0 # beginning index into ZLstruct_nodelist
NodeError = "NodeError"
ntet = 0
nhex = 0
npyr = 0
nprism = 0
nz_tet = []
nz_hex = []
nz_pyr = []
nz_prism = []
for i in range (4) :
   if ZLss [i] == 4 : # TETRAHEDRON
      nz_tet = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
               (ZLsc [i], ZLss [i]))
      ntet = ZLsc [i]
      istart = istart + ZLss [i] * ZLsc [i]
   elif ZLss[i] == 5 : # PYRAMID
      nz_pyr = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
               (ZLsc [i], ZLss [i]))
      npyr = ZLsc [i]
      # Now reorder the points (bill has the apex last instead of first)
      for ip in range (npyr) :
         tmp = nz_pyr [ip, 4]
         for jp in range (4) :
            nz_pyr [ip, 4 - jp] = nz_pyr [ip, 3 - jp]
         nz_pyr [ip, 0] = tmp
      istart = istart + ZLss [i] * ZLsc [i]
   elif ZLss[i] == 6 : # PRISM
      nz_prism = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
               (ZLsc [i], ZLss [i]))
      nprism = ZLsc [i]
      # now reorder the points (bill has a square face first)
      for ip in range (nprism) :
         tmp = nz_prism [ip, 1]
         tmpp = nz_prism [ip, 2]
         nz_prism [ip, 1] = nz_prism [ip, 4]
         nz_prism [ip, 2] = nz_prism [ip, 3]
         nz_prism [ip, 3] = tmp
         nz_prism [ip, 4] = nz_prism [ip, 5]
         nz_prism [ip, 5] = tmpp
      istart = istart + ZLss [i] * ZLsc [i]
   elif ZLss[i] == 8 : # HEXAHEDRON
      nz_hex = reshape (ZLsn [istart: istart + ZLss [i] * ZLsc [i]],
               (ZLsc [i], ZLss [i]))
      nhex = ZLsc [i]
      istart = istart + ZLss [i] * ZLsc [i]
   else :
      raise NodeError, `ZLss[i]` + "is an incorrect number of nodes."

m1 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1,
           hex = [nhex, nz_hex] ,
           pyr = [npyr, nz_pyr] ,
           tet = [ntet, nz_tet] ,
           prism = [nprism, nz_prism] , mask = "max", opt_3d = ["s4","wm"])
m2 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1,
           hex = [nhex, nz_hex] , mask = "max", opt_3d = ["s4","wm"])
m3 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1,
           pyr = [npyr, nz_pyr] , mask = "max", opt_3d = ["s4","wm"])
m4 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1,
           tet = [ntet, nz_tet] ,  mask = "max", opt_3d = ["s4","wm"])
m5 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1,
           prism = [nprism, nz_prism] , mask = "max", opt_3d = ["s4","wm"])

g1 = Graph3d (m5)
g1.plot ()
paws ()
g1 = Graph3d (m4)
g1.plot ()
paws ()
g1 = Graph3d (m3)
g1.plot ()
paws ()
g1 = Graph3d (m2)
g1.plot ()
paws ()
g1 = Graph3d (m1)
g1.plot ()
paws ()
# Eliminate all cells with an x-coordinate > 0
# TETRAHEDRA
tetinc = []
for i in range (ntet) :
   if x [nz_tet [i, 0]] <= 0.1 and x [nz_tet [i, 1]] <= 0.1 and \
      x [nz_tet [i, 2]] <= 0.1 and x [nz_tet [i, 3]] <= 0.1 :
      tetinc.append (i)
ntet = len (tetinc)
for i in range (ntet) :
   nz_tet [i] = nz_tet [tetinc [i]]
# HEXAHEDRA
hexinc = []
for i in range (nhex) :
   if x [nz_hex [i, 0]] <= 0.1 and x [nz_hex [i, 1]] <= 0.1 and \
      x [nz_hex [i, 2]] <= 0.1 and x [nz_hex [i, 3]] <= 0.1 and \
      x [nz_hex [i, 4]] <= 0.1 and x [nz_hex [i, 5]] <= 0.1 and \
      x [nz_hex [i, 6]] <= 0.1 and x [nz_hex [i, 7]] <= 0.1 :
      hexinc.append (i)
nhex = len (hexinc)
for i in range (nhex) :
   nz_hex [i] = nz_hex [hexinc [i]]
# PRISMS
prisminc = []
for i in range (nprism) :
   if x [nz_prism [i, 0]] <= 0.1 and x [nz_prism [i, 1]] <= 0.1 and \
      x [nz_prism [i, 2]] <= 0.1 and x [nz_prism [i, 3]] <= 0.1 and \
      x [nz_prism [i, 4]] <= 0.1 and x [nz_prism [i, 5]] <= 0.1 :
      prisminc.append (i)
nprism = len (prisminc)
for i in range (nprism) :
   nz_prism [i] = nz_prism [prisminc [i]]
# PYRAMIDS
pyrinc = []
for i in range (npyr) :
   if x [nz_pyr [i, 0]] <= 0.1 and x [nz_pyr [i, 1]] <= 0.1 and \
      x [nz_pyr [i, 2]] <= 0.1 and x [nz_pyr [i, 3]] <= 0.1 and \
      x [nz_pyr [i, 4]] <= 0.1 :
      pyrinc.append (i)
npyr = len (pyrinc)
for i in range (npyr) :
   nz_pyr [i] = nz_pyr [pyrinc [i]]

m1 = Mesh3d (x = x, y = y, z = z, c = c, avs = 1,
           hex = [nhex, nz_hex] ,
           pyr = [npyr, nz_pyr] ,
           tet = [ntet, nz_tet] ,
           prism = [nprism, nz_prism] , mask = "max", opt_3d = ["w4","wm"])
g1 = Graph3d (m1)
g1.plot ()
g1.quick_plot ( mesh = 1 , opt_3d = ["s4","wm"])

f.close ( )
