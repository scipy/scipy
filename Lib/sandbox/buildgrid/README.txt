    Build Grid,  Version 0.5, Jun 24, 2006
    --------------------------------------
    
Build regular grids from scattered 2D data.
Python extension module.

Summary
========  
Input:
    a set of arbitrary displaced points in 2D plane,
    grid geometry: x/y min/max and cell size, 
    a few (optional) control parameters 

Output:
    to a python list or to file

Method: 
    Put every input data point into nearest grid node (take average if more
    than one), this is main source of errors. Apply empirical method of 
    building: make sequence of averagings with full interpolation. Use 
    distances from data points to remove field values far from data points.
    The method is fair for large data volumes (and grid sizes). 
    The method is irrelevant if (x,y) must be treated as exact values,
    or the surface should have breaks or spikes of any kind.
    The surface has no special predefined properties like derivatives,
    curvatures, etc, but it is as smooth as possible in this approach.


Details
========

Installation:
    The module depends only on Python.h and standard python libraries.
    Compiled (as C, not C++) with MinGW Dev.Studio 2.05, MS VC++ 6.0,
    (Windows XP) and gcc-4 (ubuntu). To install on Windows one can just
    rename .dll to a directory on the path as BuildGrid.pyd. To compile 
    and install on Linux one can use setup.py. 

Tests:
    Python test script is BuildGrid_Test.py. To set/change test parameters
    see the script. There are two surfaces, first is almost a plane,
    second is like cos(a*x*x) to make it impossible to interpolate
    the surface at the far corners of grid. Some differences between
    built surface and exact surface are seen in terminal window.
    I used QuikGrid from http://www.PerspectiveEdge.com to see the
    results in contours and colors.
    On my laptop (1.6 Ghz) the module produces about 10**5 output nodes
    per second (without trimming, there is a weak dependence on input 
    data volume, precision, etc).

General Limitations:
    - minimal grid sizes in nodes: 9 (lesser has no sense in the method)
    - maximal grid sizes in nodes: about 10000 (currently 10915),
        that needs ~300 MB of memory to build
    - minimal cell size: 0.001 (hardcoded to avoid zero dividing)
    - all sizes, coordinates, distances should be in same units
    - at least one input triple (x,y,z) must be given
    - method is irrelevant when exact values of x and y must be used
        or the surface should have breaks or spikes of any kind
    - in grid nodes numeration, (0,0) is node with coordinates (xmin,ymin),
        x values and column numbers increase to East,
        y values and row numbers increase to North
    - output list of grid node values is indexed from bottom to top row,
        from left to right column in the rows
    - output data take at least 8*ncol*nrow bytes (double values),
        maximal memory consumption during the building is about 
        30*ncol*nrow bytes

Input data and parameters:
    - xyz-file of text lines with x,y,z values, space separated, or
    - lists of x,y,z data. Lists are separate for x, y, z data. 
      No special order of (x,y,z) triples is supposed for input data. 
    - grid boundaries are set from (xmin,ymin), (square) cell size (step)
      and grid sizes in nodes.
    - grid building method: 'Good' (default) or 'Best' (not implemented)
    - acceptable error level - absolute or relative. Default is 0.0 (Note:
      it is not an 'interpolation'). If absolute error value > 0.0, it will
      be used, otherwise if relative error value > 0.0, then it will be
      used, otherwise error will be equal to (almost) 0.0.  Relative error
      value is used as fraction of input data standard deviation.
    - trimming distance: throw out (i.e. replace by non-values) all the 
      values with distances from sources larger than trimdistance, starting 
      from grid boundaries (i.e. internal nodes may be untouched). This
      parameter is ignored if its value < step.
    - unvalue, use for output, to mark nodes without values (default: 12345678.9)

Python interface (registered functions):

    filestatistics(xyzfilename)
        The file is a sequence of text lines, each with x,y,z values,
        separated by spaces. 
		Return (number_of_values, 
                x_minimum, x_maximum, y_minimum, y_maximum, 
                s_minimum, s_maximum, s_average, s_standard)
		Raise IOError -- could not open file.
    

    fromfile(
        xyzfile=filename,       - input file of text lines: x y z
		nx=xnodes, ny=ynodes,   - grid sizes in nodes
		step=cellsize,          - cell sizes, cells are square
		xmin=minX, ymin=minY,   - left bottom, x - east, y - north
        # Optional parameters:
		method='Good',          - or 'Best' (not implemented)
		trimming=trimdistance,  - distance to trim, ignored if trimming < step
		unvalue=non_value,      - caller unvalue, ignored if no trimming
        abserror=AbsError,      - acceptable error if > 0
        relerror=RelError,      - the same, ignored if abserror > 0.0
                                  is fraction of standard (if abserror==0)
		)
        Return grid as list, rows from bottom, columns from left
		Raise ValueError/IOError for errors with file/memory
    
   fromxyz(
        xdata=x, ydata=y, zdata=z, - lists of input data (of same length)
		nx=xnodes, ny=ynodes,
		step=cellsize,       
		xmin=minX, ymin=Ymin,
        # Optional parameters are the same as for fromfile()
		method='Good',         
		trimming=trimdistance, 
		unvalue=non_value,    
        abserror=AbsError,    
        relerror=RelError,    
        )
		Return grid as list, rows from bottom, columns from left
		Raise ValueError/IOError for errors with values/memory
        
   tofile(
        outfile=filename,       - output file name
		griddata=xyz,           - list as returned by fromfile() or fromxyz()
		nx=xnodex, ny=ynodes,   - grid sizes in nodes
        filetype=grid_type,     - 'xyz' ('gxf'/'grd'/... not implemented) 
        # Optional parameters:
		step=cellsize,          - cell sizes, default: 1.0
		xmin=minX, ymin=Ymin,   - left bottom, x - east, y - north
                                - default: xmin=ymin=0.0
		unvalue=non_value,      - caller unvalue, default: 12345678.9
                                - grid nodes with unvalue will not be put to 'xyz' file
        )
		Return 0 for Ok
		Raise ValueError/IOError for errors with file/memory

ToDo:
    The module can be extended and improved in some directions:
    - more accurate calculations of errors in grid nodes - try to
        reach limits in each given node
    - more precise computing input data in nodes (add mrrgins, decrease cells)
    - various input/output data formats
    - more information about input and intermediate data
    - the same method can be used to build 3D grids, s(x,y,z,...)

Eugene Druker, 2006 
eugene.druker@gmail.com

