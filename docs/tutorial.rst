Overview
--------

The interpolate package provides tools for interpolating and extrapolating new data points from a set known set of data points.  Intepolate provides both a functional interface that is flexible and easy to use as well as an object oriented interface that can be more efficient and flexible for some cases.  It is able to interpolate and extrapolate in 1D, 2D, and even N dimensions.[fixme: 1D only right now]  

For 1D interpolation, it handles linear and spline(cubic, quadratic, and quintic) for both uniformly and non-uniformly spaced data points "out of the box."  Users can control the behavior of values that fall outside of the range of interpolation either by When new values fall outside of the range of interpolation data, the tools can be   

For 2D interpolation, 

1D Interpolation with the Functional Interface
----------------------------------------------

The following example uses the 'interp1d' function to linearly interpolate a sin curve from a sparse set of values.::

	# start up ipython for our examples.
	$ ipython -pylab
	
	In [1]: from interpolate import interp1d
	
	# Create our "known" set of 5 points with the x values in one array and the y values in another.
	In [2]: x = linspace(0,2*pi,5)
	In [3]: y = sin(x)
	
	# Now interpolate from these x,y values to create a more dense set of new_x, new_y values.
	In [4]: new_x = linspace(0,2*pi, 21)
	In [5]: new_y = interp1d(x,y, new_x)
	
	# Plot the results using matplotlib. [note examples assume you are running in ipython -pylab]
	In [6]: plot(x,y,'ro', new_x, new_y, 'b-')
	
.. image:: interp1d_linear_simple.png

 