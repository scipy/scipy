# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

helpdict = {
   "movie" : """
   movie (draw_frame, time_limit = 120., min_interframe = 0.0,
      bracket_time = array ([2., 2.], Float ), lims = None, timing = 0)
     runs a movie based on the given DRAW_FRAME function.  The movie
     stops after a total elapsed time of TIME_LIMIT seconds, which
     defaults to 60 (one minute), or when the DRAW_FRAME function
     returns zero.

     note: All but the first argument are keyword arguments, with
     defaults as shown.

     def draw_frame(1) :
       # Input argument i is the frame number.
       # draw_frame should return non-zero if there are more
       # frames in this movie.  A zero return will stop the
       # movie.
       # draw_frame must NOT include any fma command if the
       # making_movie variable is set (movie sets this variable
       # before calling draw_frame)

     If MIN_INTERFRAME is specified, a pause will be added as
     necessary to slow down the movie.  MIN_INTERFRAME is a time
     necessary to slow down the movie.  MIN_INTERFRAME is a time
     in seconds (default 0).

     The keyword bracket_time= (again a time in seconds) can be
     used to adjust the duration of the pauses after the first
     and last frames.  It may also be a two element array [beg, end].
     If the pause at the end is greater than five seconds, you will
     be prompted to explain that hitting <RETURN> will abort the final
     pause. Well, the Python version does not have this capability.

     timing = 1 enables a timing printout for your movie.

     If every frame of your movie has the same limits, use the
     limits command to fix the limits before you call movie.

   BUG:  If you hit <RETURN> to start a movie early, it will not
         pause at the end of the movie at all.  You probably should
         not use long initial pauses.
   """,
   "movie_stats" :
   """
   movie_stats ( ) or movie_stats ( timing )
     prints statistics from the last movie command, or from the
     command which produced TIMING.  TIMING is the contents of the
     movie_timing external variable after the movie command completes.
   """,
   "'pl3d'" :
   """
   General overview of module pl3d:

   (1) Viewing transform machinery.  Arguably the simplest model
       is the CAD/CAM notion that the object you see is oriented
       as you see it in the current picture.  You can then move
       it left, right, up, down, or toward or away from you,
       or you can rotate it about any of the three axes (horizontal,
       vertical, or out of the screen).  The xyz coordinates of the
       object remains unchanged throughout all of this, but this
       object coordinate system changes relative to the fixed
       xyz of the viewer, in which x is always to the right, y is
       up, and z is directed out of the screen.  Initially, the
       two coordinate systems coincide.

       rot3 (xangle,yangle,zangle)
         Rotate the object about viewer's x-axis by xangle, then
         about viewer's y-axis by yangle, then about viewer's
         z-axis by zangle
       mov3 (xchange,ychange,zchange)
         Move the object by the specified amounts.

       setz3 (zcamera)
         The "camera" is located at (0,0,zcamera) in the viewer's
         coordinate system, looking in the minus-z direction.
         Initially, zcamera is very large, and the magnification
         factor is correspondingly large, giving an isometric view.
         Decreasing zcamera makes the perspective more extreme.
         If parts of the object are behind the camera, strange things
         may happen.

       undo3 ()
       undo3 (n)
         Undo the last N (default 1) viewpoint commands (rot3, mov3,
         or setz3).  Up to 100 viewpoint changes are remembered.
       viewpoint= save3()
       ...
       restore3 (viewpoint)
         The current viewpoint transformation can be saved and later
         restored.

       gnomon (on_off)
         Toggle the gnomon (a simple display showing the orientation
         of the xyz axes of the object).
   """,
   "set_draw3_" :
   """
   set_draw3_ ( 0 | 1 ) is used to set the global draw3_,
   which controls whether the function draw3 actually shows a drawing.
   """,
   "rot3" :
   """
   rot3 (xa, ya, za)
   rotate the current 3D plot by XA about viewer's x-axis,
   YA about viewer's y-axis, and ZA about viewer's z-axis.
   SEE ALSO: orient3, mov3, aim3, setz3, undo3, save3, restore3, light3
   """,
   "mov3" :
   """
   mov3 ( [xa [, ya [, za]]])
   move the current 3D plot by XA along the viewer's x axis,
   YA along the viewer's y axis, and ZA along the viewer's z axis.
   SEE ALSO: rot3, orient3, setz3, undo3, save3, restore3, light3
   """,
   "aim3" :
   """
   aim3 ( [xa [, ya [, za]]])
   move the current 3D plot to put the point (XA, YA, ZA) in object
   coordinates at the point (0, 0, 0) -- the aim point -- in the
   viewer's coordinates. If any of the XA, YA, or ZA is nil, it defaults
   SEE ALSO: mov3, rot3, orient3, setz3, undo3, save3, restore3, light3
   to zero.
   """,
   "setz3" :
   """
   setz3 ( [zc] )
   Set the camera position to z = ZC (x = y = 0) in the viewer's coordinate
   system. If zc is None, set the camera to infinity (default).
   SEE ALSO: rot3, orient3, undo3, save3, restore3, light3
   """,
   "orient3" :
   """
   orient3 ( [phi = val1, theta = val2] )
   Set the orientation of the object to (PHI, THETA). Orientations
   are a subset of the possible rotation matrices in which the z axis
   of the object appears vertical on the screen (that is, the object
   z axis projects onto the viewer y axis). The THETA angle is the
   angle from the viewer y axis to the object z axis, positive if
   the object z axis is tilted towards you (toward viewer +z). PHI is
   zero when the object x axis coincides with the viewer x axis. If
   neither PHI nor THETA is specified, PHI defaults to - pi / 4 and
   THETA defaults to pi / 6. If only PHI is specified, THETA remains
   unchanged, unless the current THETA is near pi / 2, in which case
   THETA returns to pi / 6, or unless the current orientation does
   not have a vertical z axis, in which case THETA returns to its
   default.
   Unlike rot3, orient3 is not a cumulative operation.
   SEE ALSO: rot3, mov3, aim3, save3, restore3, light3
   """,
   "save3" :
   """
   view = save3 ( )
     Save the current 3D viewing transformation and lighting.
     Actually, this doesn't save anything; it returns a copy
     of the current 3D viewing transformation and lighting, so
     that the user can put it aside somewhere.
   SEE ALSO: restore3, rot3, mov3, aim3, light3
   """,
   "restore3" :
   """
   restore3 ( view )
   Restore a previously saved 3D viewing transformation and lighting.
   If view is missing, rotate object to viewer's coordinate system.
   SEE ALSO: restore3, rot3, mov3, aim3, light3
   """,
   "light3" :
   """
   light3 (ambient=a_level,
                    diffuse=d_level,
                    specular=s_level,
                    spower=n,
                    sdir=xyz)
     Sets lighting properties for 3D shading effects.
     A surface will be shaded according to its to its orientation
     relative to the viewing direction.

     The ambient level A_LEVEL is a light level (arbitrary units)
     that is added to every surface independent of its orientation.

     The diffuse level D_LEVEL is a light level which is proportional
     to cos(theta), where theta is the angle between the surface
     normal and the viewing direction, so that surfaces directly
     facing the viewer are bright, while surfaces viewed edge on are
     unlit (and surfaces facing away, if drawn, are shaded as if they
     faced the viewer).

     The specular level S_LEVEL is a light level proportional to a high
     power spower=N of 1+cos(alpha), where alpha is the angle between
     the specular reflection angle and the viewing direction.  The light
     source for the calculation of alpha lies in the direction XYZ (a
     3 element vector) in the viewer's coordinate system at infinite
     distance.  You can have ns light sources by making S_LEVEL, N, and
     XYZ (or any combination) be vectors of length ns (3-by-ns in the
     case of XYZ).  (See source code for specular_hook function
     definition if powers of 1+cos(alpha) aren't good enough for you.)

     With no arguments, return to the default lighting.

   EXAMPLES:
     light3 ( diffuse=.1, specular=1., sdir=[0,0,-1])
       (dramatic "tail lighting" effect)
     light3 ( diffuse=.5, specular=1., sdir=[1,.5,1])
       (classic "over your right shoulder" lighting)
     light3 ( ambient=.1,diffuse=.1,specular=1.,
             sdir=[[0,0,-1],[1,.5,1]],spower=[4,2])
       (two light sources combining previous effects)
   SEE ALSO: rot3, save3, restore3
   """,
   "get3_light" :
   """
   get3_light(xyz, nxyz)
      or get3_light(xyz)

     return 3D lighting for polygons with vertices XYZ.  If NXYZ is
     specified, XYZ should be sum(nxyz)-by-3, with NXYZ being the
     list of numbers of vertices for each polygon (as for the plfp
     function).  If NXYZ is not specified, XYZ should be a quadrilateral
     mesh, ni-by-nj-by-3 (as for the plf function).  In the first case,
     the return value is len (NXYZ) long; in the second case, the
     return value is (ni-1)-by-(nj-1).

     The parameters of the lighting calculation are set by the
     light3 function.

     SEE ALSO: light3, set3_object, get3_normal, get3_centroid
   """,
   "get3_normal" :
   """
     get3_normal(xyz, nxyz)
         or get3_normal(xyz)

     return 3D normals for polygons with vertices XYZ.  If NXYZ is
     specified, XYZ should be sum(nxyz)-by-3, with NXYZ being the
     list of numbers of vertices for each polygon (as for the plfp
     function).  If NXYZ is not specified, XYZ should be a quadrilateral
     mesh, ni-by-nj-by-3 (as for the plf function).  In the first case,
     the return value is len(NXYZ)-by-3; in the second case, the
     return value is (ni-1)-by-(nj-1)-by-3.

     The normals are constructed from the cross product of the lines
     joining the midpoints of two edges which as nearly quarter the
     polygon as possible (the medians for a quadrilateral).  No check
     is made that these not be parallel; the returned "normal" is
     [0,0,0] in that case.  Also, if the polygon vertices are not
     coplanar, the "normal" has no precisely definable meaning.

     SEE ALSO: get3_centroid, get3_light
   """,
   "get3_centroid" :
   """
     get3_centroid(xyz, *nxyz)
         or get3_centroid(xyz)

     return 3D centroids for polygons with vertices XYZ.  If NXYZ is
     specified, XYZ should be sum(nxyz)-by-3, with NXYZ being the
     list of numbers of vertices for each polygon (as for the plfp
     function).  If NXYZ is not specified, XYZ should be a quadrilateral
     mesh, ni-by-nj-by-3 (as for the plf function).  In the first case,
     the return value is len(NXYZ) in length; in the second case, the
     return value is (ni-1)-by-(nj-1)-by-3.

     The centroids are constructed as the mean value of all vertices
     of each polygon.

     SEE ALSO: get3_normal, get3_light
   """,
   "get3_xy" :
   """
     get3_xy (xyz)
         or get3_xy(xyz, 1)

     Given anything-by-3 coordinates XYZ, return X and Y in viewer's
     coordinate system (set by rot3, mov3, orient3, etc.).  If the
     second argument is present and non-zero, also return Z (for use
     in sort3d or get3_light, for example).  If the camera position
     has been set to a finite distance with setz3, the returned
     coordinates will be tangents of angles for a perspective
     drawing (and Z will be scaled by 1/zc).
     Unlike the Yorick version, this function returns a 3-by-anything
     array of coordinates.
     Actually, what it returns is a 3-by-anything python array, whose
     0th element is the x array, whose 1th element is the y array, and
     whose 2th element is the z array if asked for.
     I believe that x, y, and z can be either 1d or 2d, so this
     routine is written in two cases.
   """,
   "undo3" :
   """
     undo3 ()
         or undo3 (n)
     Undo the effects of the last N (default 1) rot3, orient3, mov3, aim3,
     setz3, or light3 commands.
   """,
   "set3_object" :
   """
     set3_object (drawing_function, [arg1,arg2,...])

     set up to trigger a call to draw3, adding a call to the
     3D display list of the form:

        DRAWING_FUNCTION ( [ARG1, ARG2, ...]))

     When draw3 calls DRAWING_FUNCTION, the external variable draw3_
     will be non-zero, so DRAWING_FUNCTION can be written like this:

     def drawing_function(arg) :
      
       if (draw3_) :
          arg1= arg [0]
          arg1= arg [1]
          ...
          ...<calls to get3_xy, sort3d, get3_light, etc.>...
          ...<calls to graphics functions plfp, plf, etc.>...
          return

       ...<verify args>...
       ...<do orientation and lighting independent calcs>...
       set3_object (drawing_function, [arg1,arg2,...])

   SEE ALSO: get3_xy, get3_light, sort3d
   """,
   "window3" :
   """
   window3 ( ) or window3 (n)
   initialize style="nobox.gs" window for 3D graphics
   """,
   "sort3d" :
   """
   sort3d(z, npolys)
     given Z and NPOLYS, with len(Z)==sum(npolys), return
     a 2-element list [LIST, VLIST] such that Z[VLIST] and NPOLYS[LIST] are
     sorted from smallest average Z to largest average Z, where
     the averages are taken over the clusters of length NPOLYS.
     Within each cluster (polygon), the cyclic order of Z[VLIST]
     remains unchanged, but the absolute order may change.

     This sorting order produces correct or nearly correct order
     for a plfp command to make a plot involving hidden or partially
     hidden surfaces in three dimensions.  It works best when the
     polys form a set of disjoint closed, convex surfaces, and when
     the surface normal changes only very little between neighboring
     polys.  (If the latter condition holds, then even if sort3d
     mis-orders two neighboring polys, their colors will be very
     nearly the same, and the mistake won't be noticeable.)  A truly
     correct 3D sorting routine is impossible, since there may be no
     rendering order which produces correct surface hiding (some polys
     may need to be split into pieces in order to do that).  There
     are more nearly correct algorithms than this, but they are much
     slower.
   SEE ALSO: get3_xy
   """,
   "draw3" :
   """
      draw3 (called_as_idler = 0, lims = None):
   Draw the current 3d display list.
   Ordinarily triggered automatically when the drawing changes.
   """,
   "gnomon" :
   """
   gnomon ()
      or gnomon (onoff)
     Toggle the gnomon display. If on is present and non-zero,
     turn on the gnomon. If zero, turn it off.

     The gnomon shows the X, Y, and Z axis directions in the
     object coordinate system. The directions are labeled.
     The gnomon is always infinitely far behind the object
     (away from the camera).

     There is a mirror-through-the-screen-plane ambiguity in the
     display which is resolved in two ways: (1) the (X, Y, Z)
     coordinate system is right-handed, and (2) If the tip of an
     axis projects into the screen, its label is drawn in opposite
     polarity to the other text in the screen.
   """,
   "spin3" :
   """
   spin3 ( ) or spin3 (nframes) os spin3 (nframes, axis)
     Spin the current 3D display list about AXIS over NFRAMES.  Keywords
     tlimit= the total time allowed for the movie in seconds (default 60),
     dtmin= the minimum allowed interframe time in seconds (default 0.0),
     bracket_time= (as for movie function in movie.i), timing = 1 if
     you want timing measured and printed out, 0 if not.

     The default AXIS is [-1,1,0] and the default NFRAMES is 30.
   SEE ALSO: rot3
   """,
   "plwf" :
   """
   plwf (z)
   or plwf (z, y, x)

     plots a 3-D wire frame of the given Z array, which must have the
     same dimensions as the mesh (X, Y).  If X and Y are not given, they
     default to the first and second indices of Z, respectively.
     The drawing order of the zones is determined by a simple "painter's
     algorithm", which works fairly well if the mesh is reasonably near
     rectilinear, but can fail even then if the viewpoint is chosen to
     produce extreme fisheye perspective effects.  Look at the resulting
     plot carefully to be sure the algorithm has correctly rendered the
     model in each case.

   KEYWORDS: fill   -- optional colors to use (default is to make zones
                       have background color), same dimension options as
                       for z argument to plf function
             shade  -- set non-zero to compute shading from current
                       3D lighting sources
             edges  -- default is 1 (draw edges), but if you provide fill
                       colors, you may set to 0 to supress the edges
             ecolor, ewidth  -- color and width of edges
             cull   -- default is 1 (cull back surfaces), but if you want
                       to see the "underside" of the model, set to 0
             scale  -- by default, Z is scaled to "reasonable" maximum
                       and minimum values related to the scale of (X,Y).
                       This keyword alters the default scaling factor, in
                       the sense that scale=2.0 will produce twice the
                       Z-relief of the default scale=1.0.
             cmax   -- the ambient= keyword in light3 can be used to
                       control how dark the darkest surface is; use this
                       to control how light the lightest surface is
                       the lightwf routine can change this parameter
                       interactively

   SEE ALSO: lightwf, plm, plf, orient3, light3, fma3, window3
   """,
   "lightwf" :
   """
   lightwf (cmax)
     Sets the cmax= parameter interactively, assuming the current
     3D display list contains the result of a previous plwf call.
     This changes the color of the brightest surface in the picture.
     The darkest surface color can be controlled using the ambient=
     keyword to light3.

   SEE ALSO: plwf, light3
   """,
   "xyz_wf" :
   """
   xyz_wf (z, [y, x] [,scale = 1.0])
      returns a 3-by-ni-by-nj array whose 0th entry is x, 1th entry
      is y, and 2th entry is z. z is ni-by-nj. x and y, if present,
      must be the same shape. If not present, integer ranges will
      be used to create an equally spaced coordinate grid in x and y.
      The function which scales the "topography" of z(x,y) is
      potentially useful apart from plwf.
      For example, the xyz array used by plwf can be converted from
      a quadrilateral mesh plotted using plf to a polygon list plotted
      using plfp like this:
        xyz= xyz_wf(z,y,x,scale=scale);
        ni= shape(z)[1];
        nj= shape(z)[2];
        list = ravel (add.outer (
           ravel(add.outer (adders,zeros(nj-1, Int))) +
           arange((ni-1)*(nj-1), typecode = Int),
           array ( [[0, 1], [nj + 1, nj]])))
        xyz=array([take(ravel(xyz[0]),list),
           take(ravel(xyz[1]),list),
           take(ravel(xyz[2]),list)])
        nxyz= ones((ni-1)*(nj-1)) * 4;
      The resulting array xyz is 3-by-(4*(nj-1)*(ni-1)).
      xyz[0:3,4*i:4*(i+1)] are the clockwise coordinates of the
      vertices of cell number i.
   """,
   "plane3" :
   """
   plane3(normal, point)
         or plane3([nx,ny,nz], [px,py,pz])

     returns [nx,ny,nz,pp] for the specified plane.
   """,
   "mesh3" :
   """
    mesh3(x,y,z)
         or mesh3(x,y,z, funcs = [f1,f2,...])
         or mesh3(xyz, funcs = [f1,f2,...])
         or mesh3(nxnynz, dxdydz, x0y0z0, funcs = [f1,f2,...])

     make mesh3 argument for slice3, xyz3, getv3, etc., functions.
     X, Y, and Z are each 3D coordinate arrays.  The optional F1, F2,
     etc. are 3D arrays of function values (e.g. density, temperature)
     which have one less value along each dimension than the coordinate
     arrays.  The "index" of each zone in the returned mesh3 is
     the index in these cell-centered Fi arrays, so every index from
     one through the total number of cells indicates one real cell.
     The Fi arrays can also have the same dimensions as X, Y, or Z
     in order to represent point-centered quantities.

     If X has four dimensions and the length of the first is 3, then
     it is interpreted as XYZ (which is the quantity actually stored
     in the returned cell list).

     If X is a vector of 3 integers, it is interpreted as [nx,ny,nz]
     of a uniform 3D mesh, and the second and third arguments are
     [dx,dy,dz] and [x0,y0,z0] respectively.  (DXDYDZ represent the
     size of the entire mesh, not the size of one cell, and NXNYNZ are
     the number of cells, not the number of points.)

     Added by ZCM 1/13/97: if x, y, and z are one-dimensional of
     the same length and if the keyword verts exists and yields
     an NCELLS by 8 integer array, then we have an unstructured
     rectangular mesh, and the subscripts of cell i's vertices
     are verts[i, 0:8].

     Other sorts of meshes are possible -- a mesh which lives
     in a binary file is an obvious example -- which would need
     different workers for xyz3, getv3, getc3, and iterator3
     iterator3_rect may be more general than the other three;
     as long as the cell dimensions are the car of the list
     which is the 2nd car of m3, it will work.
   """,
   "slice3" :
   """
   slice3 (m3, fslice, nverts, xyzverts)
         or color_values= slice3(m3, fslice, nverts, xyzverts, fcolor)
         or color_values= slice3(m3, fslice, nverts, xyzverts, fcolor, 1)

     slice the 3D mesh M3 using the slicing function FSLICE, returning
     the list [NVERTS, XYZVERTS, color].  Note that it is impossible to
     pass arguments as addresses, as yorick does in this routine.
     NVERTS is the number of vertices in each polygon of the slice, and
     XYZVERTS is the 3-by-sum(NVERTS) list of polygon vertices.  If the
     FCOLOR argument is present, the values of that coloring function on
     the polygons are returned as the value of the slice3 function
     (numberof(color_values) == numberof(NVERTS) == number of polygons).

     If the slice function FSLICE is a function, it should be of the
     form:
        func fslice(m3, chunk)
     returning a list of function values on the specified chunk of the
     mesh m3.  The format of chunk depends on the type of m3 mesh, so
     you should use only the other mesh functions xyz3 and getv3 which
     take m3 and chunk as arguments.  The return value of fslice should
     have the same dimensions as the return value of getv3; the return
     value of xyz3 has an additional first dimension of length 3.

     If FSLICE is a list of 4 numbers, it is taken as a slicing plane
     with the equation FSLICE(+:1:3)*xyz(+)-FSLICE(4), as returned by
     plane3.

     If FSLICE is a single integer, the slice will be an isosurface for
     the FSLICEth variable associated with the mesh M3.  In this case,
     the keyword value= must also be present, representing the value
     of that variable on the isosurface.

     If FCOLOR is nil, slice3 returns nil.  If you want to color the
     polygons in a manner that depends only on their vertex coordinates
     (e.g.- by a 3D shading calculation), use this mode.

     If FCOLOR is a function, it should be of the form:
        func fcolor(m3, cells, l, u, fsl, fsu, ihist)
     returning a list of function values on the specified cells of the
     mesh m3.  The cells argument will be the list of cell indices in
     m3 at which values are to be returned.  l, u, fsl, fsu, and ihist
     are interpolation coefficients which can be used to interpolate
     from vertex centered values to the required cell centered values,
     ignoring the cells argument.  See getc3 source code.
     The return values should always have dimsof(cells).

     If FCOLOR is a single integer, the slice will be an isosurface for
     the FCOLORth variable associated with the mesh M3.

     If the optional argument after FCOLOR is non-nil and non-zero,
     then the FCOLOR function is called with only two arguments:
        func fcolor(m3, cells)
   """,
   "slice3mesh" :
   """
   slice3mesh returns a triple [nverts, xyzverts, color]
    nverts is no_cells long and the ith entry tells how many
       vertices the ith cell has.
    xyzverts is sum (nverts) by 3 and gives the vertex
       coordinates of the cells in order.
    color, if present, is len (nverts) long and contains
       a color value for each cell in the mesh.

   There are a number of ways to call slice3mesh:
 
      slice3mesh (z, color = None)

   z is a two dimensional array of cell function values, assumed
      to be on a uniform mesh nx by ny cells (assuming z is nx by ny)
      nx being the number of cells in the x direction, ny the number
      in the y direction.
   color, if specified, is either an nx by ny array
      of cell-centered values by which the surface is to
      be colored, or an nx +1 by ny + 1 array of vertex-
      centered values, which will be averaged over each
      cell to give cell-centered values.

      slice3mesh (nxny, dxdy, x0y0, z, color = None)

   In this case, slice3mesh accepts the specification for
   a regular 2d mesh: nxny is the number of cells in the
   x direction and the y direction; x0y0 are the initial
   values of x and y; and dxdy are the increments in the
   two directions. z is the height of a surface above
   the xy plane and must be dimensioned nx + 1 by ny + 1.
   color, if specified, is as above.

      slice3mesh (x, y, z, color = None)

   z is as above, an nx by ny array of function values
   on a mesh of the same dimensions. There are two choices
   for x and y: they can both be one-dimensional, dimensioned
   nx and ny respectively, in which case they represent a
   mesh whose edges are parallel to the axes; or else they
   can both be nx by ny, in which case they represent a
   general quadrilateral mesh.
   color, if specified, is as above.
   """,
   "iterator3" :
   """
   iterator3 (m3)
   iterator3 (m3, chunk, clist)
   iterator3_rect (m3)
   iterator3_rect (m3, chunk, clist)
   iterator3_irreg (m3)
   iterator3_irreg (m3, chunk, clist)

   The iterator3 functions combine three distinct operations:
   (1) If only the M3 argument is given, return the initial
       chunk of the mesh.  The chunk will be no more than
       chunk3_limit cells of the mesh.
   (2) If only M3 and CHUNK are given, return the next CHUNK,
       or None if there are no more chunks.
   (3) If M3, CHUNK, and CLIST are all specified, return the
       absolute cell index list corresponding to the index list
       CLIST of the cells in the CHUNK.
       Do not increment the chunk in this case.

   The form of the CHUNK argument and return value for cases (1)
   and (2) is not specified, but it must be recognized by the
   xyz3 and getv3 functions which go along with this iterator3.
   (For case (3), CLIST and the return value are both ordinary
   index lists.)
   In the irregular case, it is guaranteed that the returned chunk
   consists of only one type of cell (tetrahedra, hexahedra,
   pyramids, or prisms).
   """,
   "getv3" :
   """
   getv3(i, m3, chunk)

     return vertex values of the Ith function attached to 3D mesh M3
     for cells in the specified CHUNK.  The CHUNK may be a list of
     cell indices, in which case getv3 returns a 2x2x2x(dimsof(CHUNK))
     list of vertex coordinates.  CHUNK may also be a mesh-specific data
     structure used in the slice3 routine, in which case getv3 may
     return a (ni)x(nj)x(nk) array of vertex values.  For meshes which
     are logically rectangular or consist of several rectangular
     patches, this is up to 8 times less data, with a concomitant
     performance advantage.  Use getv3 when writing slicing functions
     for slice3.
   getv3_rect(i, m3, chunk) does the job for a regular rectangular
     mesh.
   getv3_irreg (i, m3, chunk) :
     for an irregular mesh, returns a 3-list whose elements are:
     (1) the function values for the ith function on the vertices of the
     given chunk. (The function values must have the same dimension
     as the coordinates; there is no attempt to convert zone-centered
     values to vertex-centered values.)
     (2) an array of relative cell numbers within the list of cells
     of this type.
     (3) a number that can be added to these relative numbers to gives
     the absolute cell numbers for correct access to their coordinates
     and function values.
   """,
   "getc3" :
   """
   getc3(i, m3, chunk)
         or getc3(i, m3, clist, l, u, fsl, fsu, cells)

     return cell values of the Ith function attached to 3D mesh M3
     for cells in the specified CHUNK.  The CHUNK may be a list of
     cell indices, in which case getc3 returns a (dimsof(CHUNK))
     list of vertex coordinates.  CHUNK may also be a mesh-specific data
     structure used in the slice3 routine, in which case getc3 may
     return a (ni)x(nj)x(nk) array of vertex values.  There is no
     savings in the amount of data for such a CHUNK, but the gather
     operation is cheaper than a general list of cell indices.
     Use getc3 when writing colorng functions for slice3.

     If CHUNK is a CLIST, the additional arguments L, U, FSL, and FSU
     are vertex index lists which override the CLIST if the Ith attached
     function is defined on mesh vertices.  L and U are index lists into
     the (dimsof(CLIST))x2x2x2 vertex value array, say vva, and FSL
     and FSU are corresponding interpolation coefficients; the zone
     centered value is computed as a weighted average of involving these
     coefficients.  The CELLS argument is required by histogram to do
     the averaging.  See the source code for details.
     By default, this conversion (if necessary) is done by averaging
     the eight vertex-centered values.
   getc3_rect (i, m3, chunk, l, u, fsl, fsu, cells) does the job
     for a regular rectangular mesh.
   getc3_irreg (i, m3, chunk, l, u, fsl, fsu, cells) :
      Same thing as getc3_rect, i. e., returns the same type of
      data structure, but from an irregular mesh.
      m3 [1] is a 2-list; m3[1] [0] is an array whose ith element
         is an array of coordinate indices for the ith cell,
         or a list of up to four such arrays.
         m3 [1] [1] is the 3 by nverts array of coordinates.
      m3 [2] is a list of arrays of vertex-centered or cell-centered
         data.
   chunk may be a list, in which case chunk [0] is a 2-sequence
      representing a range of cell indices; or it may be a one-dimensional
      array, in which case it is a nonconsecutive set of cell indices.
      It is guaranteed that all cells indexed by the chunk are the
      same type.
   """,
   "slice2x" :
   """
   slice2x (plane, nverts, xyzverts, values)

     Slice a polygon list, retaining only those polygons or
     parts of polygons on the positive side of PLANE, that is,
     the side where xyz(+)*PLANE(+:1:3)-PLANE(4) > 0.0.
     The NVERTS, VALUES, and XYZVERTS arrays have the meanings of
     the return values from the slice3 function.
     Python returns a sextuple
     [nverts, xyzverts, values, nvertb, xyzvertb, valueb]
     with None in the place of missing or None input arguments.

   slice2_precision= precision
     Controls how slice2 (or slice2x) handles points very close to
     the slicing plane.  PRECISION should be a positive number or zero.
     Zero PRECISION means to clip exactly to the plane, with points
     exactly on the plane acting as if they were slightly on the side
     the normal points toward.  Positive PRECISION means that edges
     are clipped to parallel planes a distance PRECISION on either
     side of the given plane.  (Polygons lying entirely between these
     planes are completely discarded.)

     Default value is 0.0.
   """,
   "slice2" :
   """
   slice2 (plane, nverts, xyzverts, values)

     Slice a polygon list, retaining only those polygons or
     parts of polygons on the positive side of PLANE, that is,
     the side where xyz(+)*PLANE(+:1:3)-PLANE(4) > 0.0.
     The NVERTS, VALUES, and XYZVERTS arrays have the meanings of
     the return values from the slice3 function.
     Python returns a sextuple
     [nverts, xyzverts, values, nvertb, xyzvertb, valueb]
     with None in the place of missing or None input arguments.
     It is legal for the VALUES argument to be None (e.g.- if there
     is no fcolor function).

     In order to plot two intersecting slices, one could
     slice (for example) the horizontal plane twice (slice2x) -
     first with the plane of the vertical slice, then with minus
     that same plane.  Then, plot first the back part of the
     slice, then the vertical slice, then the front part of the
     horizontal slice.  Of course, the vertical plane could
     be the one to be sliced, and "back" and "front" vary
     depending on the view point, but the general idea always
     works.

   slice2_precision= precision
     Controls how slice2 (or slice2x) handles points very close to
     the slicing plane.  PRECISION should be a positive number or zero.
     Zero PRECISION means to clip exactly to the plane, with points
     exactly on the plane acting as if they were slightly on the side
     the normal points toward.  Positive PRECISION means that edges
     are clipped to parallel planes a distance PRECISION on either
     side of the given plane.  (Polygons lying entirely between these
     planes are completely discarded.)

     Default value is 0.0.
   """,
   "pl3surf" :
   """
   pl3surf (nverts, xyzverts)
         or pl3surf (nverts, xyzverts, values)

     Perform simple 3D rendering of an object created by slice3
     (possibly followed by slice2).  NVERTS and XYZVERTS are polygon
     lists as returned by slice3, so XYZVERTS is sum(NVERTS)-by-3,
     where NVERTS is a list of the number of vertices in each polygon.
     If present, the VALUES should have the same length as NVERTS;
     they are used to color the polygon.  If VALUES is not specified,
     the 3D lighting calculation set up using the light3 function
     will be carried out.  Keywords cmin= and cmax= as for plf, pli,
     or plfp are also accepted.  (If you do not supply VALUES, you
     probably want to use the ambient= keyword to light3 instead of
     cmin= here, but cmax= may still be useful.)
   """,
   "pl3tree" :
   """
   pl3tree (nverts, xyzverts = None, values = None, plane = None,
      cmin = None, cmax = None)
   
     Add the polygon list specified by NVERTS (number of vertices in
     each polygon) and XYZVERTS (3-by-sum(NVERTS) vertex coordinates)
     to the currently displayed b-tree.  If VALUES is specified, it
     must have the same dimension as NVERTS, and represents the color
     of each polygon.  If VALUES is not specified, the polygons
     are assumed to form an isosurface which will be shaded by the
     current 3D lighting model; the isosurfaces are at the leaves of
     the b-tree, sliced by all of the planes.  If PLANE is specified,
     the XYZVERTS must all lie in that plane, and that plane becomes
     a new slicing plane in the b-tree.

     Each leaf of the b-tree consists of a set of sliced isosurfaces.
     A node of the b-tree consists of some polygons in one of the
     planes, a b-tree or leaf entirely on one side of that plane, and
     a b-tree or leaf on the other side.  The first plane you add
     becomes the root node, slicing any existing leaf in half.  When
     you add an isosurface, it propagates down the tree, getting
     sliced at each node, until its pieces reach the existing leaves,
     to which they are added.  When you add a plane, it also propagates
     down the tree, getting sliced at each node, until its pieces
     reach the leaves, which it slices, becoming the nodes closest to
     the leaves.

     tree is a 4-element list like this:
      [plane, back_tree, inplane_leaf, front_tree]
      plane= tree [0]  is None if this is just a leaf
                       in which case, only inplane_leaf is not None
      back_tree= tree [1]    is the part behind plane
      inplane_leaf= tree [2] is the part in the plane itself
      front_tree= tree [3]   is the part in front of plane

     This structure is relatively easy to plot, since from any
     viewpoint, a node can always be plotted in the order from one
     side, then the plane, then the other side.

     This routine assumes a "split palette"; the colors for the
     VALUES will be scaled to fit from color 0 to color 99, while
     the colors from the shading calculation will be scaled to fit
     from color 100 to color 199.  (If VALUES is specified as a char
     array, however, it will be used without scaling.)
     You may specifiy a cmin= or cmax= keyword to affect the
     scaling; cmin is ignored if VALUES is not specified (use the
     ambient= keyword from light3 for that case).
   """,
   "split_palette" :
   """
   split_palette
         or split_palette ("palette_name.gp")
     split the current palette or the specified palette into two
     parts; colors 0 to 99 will be a compressed version of the
     original, while colors 100 to 199 will be a gray scale.
   """,
   "split_bytscl" :
   """
   split_bytscl (x, upper, cmin = None, cmax = None)
     as bytscl function, but scale to the lower half of a split
     palette (0-99, normally the color scale) if the second parameter
     is zero or nil, or the upper half (100-199, normally the gray
     scale) if the second parameter is non-zero.
   """,
   "xyz3" :
   """
   xyz3 (m3, chunk)

     return vertex coordinates for CHUNK of 3D mesh M3.  The CHUNK
     may be a list of cell indices, in which case xyz3 returns a
     (dimsof(CHUNK))x3x2x2x2 list of vertex coordinates.  CHUNK may
     also be a mesh-specific data structure used in the slice3
     routine, in which case xyz3 may return a 3x(ni)x(nj)x(nk)
     array of vertex coordinates.  For meshes which are logically
     rectangular or consist of several rectangular patches, this
     is up to 8 times less data, with a concomitant performance
     advantage.  Use xyz3 when writing slicing functions or coloring
     functions for slice3.
   """,
   "to_corners3" :
   """
   to_corners3(list, nj, nk)
     convert an array of cell indices in an (ni-1)-by-(NJ-1)-by-(NK-1)
     logically rectangular grid of cells into the list of
     len(LIST)-by-2-by-2-by-2 cell corner indices in the
     corresponding ni-by-NJ-by-NK array of vertices.
     Note that this computation in Yorick gives an absolute offset
     for each cell quantity in the grid. In Yorick it is legal to
     index a multidimensional array with an absolute offset. In
     Python it is not. However, an array can be flattened if
     necessary.
     Other changes from Yorick were necessitated by row-major
     order and 0-origin indices, and of course the lack of
     Yorick array facilities.
   """,
   "'yorick'" :
   """
   The yorick module supplies Python versions of some common
   yorick functions: zcen_, dif_, maxelt_, minelt_, rem_0_,
   avg_, timer_, timer_print.
   """,
   "zcen_" :
   """
   zcen_(x, i) does the same thing as in Yorick: x(...,zcen,...)
   where zcen is the ith subscript. (works for up to 5 dimensions).
   Namely, the elements along the ith dimension of x are replaced
   by the averages of adjacent pairs, and the dimension decreases
   by one. Remember that Python sunscripts are counted from 0.
   """,
   "dif_" :
   """
   dif_(x, i) does the same thing as in Yorick: x(...,dif_,...)
   where dif_ is the ith subscript. (works for up to 5 dimensions).
   Namely, the elements along the ith dimension of x are replaced
   by the differences of adjacent pairs, and the dimension decreases
   by one. Remember that Python sunscripts are counted from 0.
   """,
   "maxelt_" :
   """
   maxelt_ accepts a sequence of one or more possible multi-dimensional
   objects and computes their maximum. In principle these can be of
   arbitrary complexity, since the routine recurses.
   """,
   "minelt_" :
   """
   minelt_ accepts a sequence of one or more possible multi-dimensional
   objects and computes their minimum. In principle these can be of
   arbitrary complexity, since the routine recurses.
   """,
   "rem_0_" :
   """
   rem_0_ (z) goes through array z and replaces any zero
   elements with 1.e-35. Assumes z has one or two dimensions.
   """,
   "avg_" :
   """
   avg_ (z) returns the average of all elements of its array
   argument.
   """,
   "timer_" :
   """
   timer (elapsed) returns a triple consisting of the times
   [cpu, system, wall].
   timer (elapsed, split) returns a sequence whose first element
   is [cpu, system, wall] and whose second element is the
   sum of split and the difference between ththe new and old values
   of 'elapsed.'
   """,
   "timer_print" : 
   """
   timer_print (label1, split1 [,label2, split2, ...]) prints
   out a timing summary for splits accumulated by timer_.
   """
}

from string import *

def gist3dhelp ( routine ) :
   routine = `routine` [0:min (len (`routine`), 58)]
   key = None
   lenkey = 0
   for k in helpdict.keys () :
      if find (routine, k) >= 0 :
         if len (k) > lenkey :
            key = k
            lenkey = len (k)
   if key is not None:
      print helpdict [key]
   else :
      print "No help available for " + routine + "."
   return
