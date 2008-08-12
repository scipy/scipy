""" ND Interpolation wrapping using code from NDImage"""

from numpy import array, arange, NaN
import numpy as np
import _nd_image # extension module with interpolation code

def interpNd(data, coordinates, starting_coords=None, spacings=None, kind='linear',out=NaN):
    """ A function for interpolation of 1D, real-valued data.
        
        Parameters
        -----------
            
            data -- NumPy array (N-dimensional) or list of lists
                indicates the known values of the function.
            
            coordinates -- array or list
                To interpolate at a set of L points, array must be NxL, where
                each column denotes a point.  If only one point is desired, its
                coordinates may be entered as either a list or an array.
                
        Optional Arguments
        -------------------
        
            starting_coords -- array or list
                indicates the point in space
                whose value is given by data[0, ..., 0].
                Defaults to being all zeros.
                
            spacings -- array or list
                jth component gives spacing
                of points along the jth axis.  Defaults
                to being all ones.
        
            kind -- A string or integer
                Indicates what interpolation method to perform on
                points within the region of interpolation
                
                0, 'block' -- block interpolation based on interval midpoints
                1, 'linear' -- linear interpolation
                2, 'quadratic' -- spline order 2 interpolation
                3, 'cubic' -- cubic spline interpolation
                4, 'quartic' -- 4th order spline interpolation
                5, 'quintic' -- 5th order spine interpolation
                
            out -- string or NaN
                Indicates how to extrapolate values at points outside 
                the region of interpolation.
            
                NaN -- return NaN for all points out of range
                'nearest' -- return value at nearest valid point
                'constant' -- returns 0
                'wrap' -- points over one boundary wrap around to the other side
                'reflect' -- out-of-bounds points are reflected into the valid region
                
        Example
        --------
        
            >>> import numpy as np
            >>> from interpolate import interpNd
            >>> boring_data = np.ones((5,5,5))
            >>> nd.interpNd(boring_data, np.array([[2.3], [1.0], [3.9]]))
            1.0
    """
    return InterpolateNd(data = data,
                                    starting_coords = starting_coords,
                                    spacings = spacings,
                                    kind = kind,
                                    out = out
                                    )(coordinates)

class InterpolateNd:
    """ A callable class for interpolation of 1D, real-valued data.
    
        Under the hood, interpolation is always done with splines.  Those
        of order 2 or higher have the boundary condition of zero derivative.
        This is a potential shortcoming of the package; natural splines would be better.
        
        Parameters
        -----------
            
            data -- NumPy array (N-dimensional) or list of lists
                indicates the known values of the function.
                
        Optional Arguments
        -------------------
        
            starting_coords -- array or list
                indicates the point in space
                whose value is given by data[0, ..., 0].
                Defaults to being all zeros.
                
            spacings -- array or list
                jth component gives spacing
                of points along the jth axis.  Defaults
                to being all ones.
        
            kind -- A string or integer
                Indicates what interpolation method to perform on
                points within the region of interpolation
                
                0, 'block' -- block interpolation based on interval midpoints
                1, 'linear' -- linear interpolation
                2, 'quadratic' -- spline order 2 interpolation
                3, 'cubic' -- cubic spline interpolation
                4, 'quartic' -- 4th order spline interpolation
                5, 'quintic' -- 5th order spine interpolation
                
            out -- string or NaN
                Indicates how to extrapolate values at points outside 
                the region of interpolation.
            
                NaN -- return NaN for all points out of range
                'nearest' -- return value at nearest valid point
                'constant' -- returns 0
                'wrap' -- points over one boundary wrap around to the other side
                'reflect' -- out-of-bounds points are reflected into the valid region
                
        Example
        --------
        
            >>> import numpy as np
            >>> from interpolate import InterpolateNd
            >>> boring_data = np.ones((5,5,5))
            >>> nd.InterpolateNd(boring_data)( np.array([[2.3], [1.0], [3.9]]) )
            1.0
    """
    # Under the hood this works by interpolating entries in an array, rather than
    # values at points in space.  Known data are stored as an array; info on spacings
    # and starting coords describes the mapping from the array to ND space.  When the
    # function is called, the coordinates are first converted into array values.
    # The array interpolation is handled by the C extension _nd_image.
    def __init__(self, data, starting_coords =None, spacings = None, 
                        kind='linear', out=NaN):
        """ Parameters
            -----------
                data = array or list of lists
                
            Optional Parameters
            -------------------
                starting_coords = None, list, 1D array or 2D (nx1) array
                        spatial coordinate where data[0,..., 0] is measured.
                        Defaults to zeros.
                spacings = None, list, 1D array or 2D (nx1) array
                        spacings[i] is spacing along axis i
                        defaults to all ones
                kind = string or integer
                        indicates type of interpolation to perform
                out = string in 'nearest', 'wrap', 'reflect', 'mirror', 'constant'
                        or just NaN
        """
        
        # checking format of input
        data = array(data)
        
        # for proper processing later, starting_coords and spacings must be of shape (data.ndim, 1)
        if starting_coords == None:
            starting_coords = np.zeros(( data.ndim, 1 ))
        else:
            starting_coords = array(starting_coords)
            assert starting_coords.size == data.ndim, "There must be one element of \
                            starting_coords per data dimension.  Size mismatch."
            starting_coords = np.reshape(starting_coords, (data.ndim, 1))
        if spacings == None:
            spacings = np.ones(( data.ndim, 1 ))
        else:
            spacings = array(spacings)
            assert starting_coords.size == data.ndim, "There must be one element of \
                            starting_coords per data dimension"
            spacings = np.reshape(spacings, (data.ndim, 1))
        
        # determining the order
        order_dict = \
            { 0:0,
                '0':0,
                'block':0,
                1:1,
                '1':1,
                'linear':1,
                2:2,
                '2':2,
                'quadratic':2,
                'quad':2,
                3:3,
                '3':3,
                'spline':3,
                'cubic':3,
                'natural':3,
                4:4,
                '4':4,
                'quartic':4,
                5:5,
                '5':5,
                'quintic':5,
                'quint':5,
                }
        if order_dict.has_key(str(kind).lower()):
            self.order = order_dict[str(kind).lower()]
        elif isinstance(kind, int):
            raise ValueError, "Only spline orders 0, 1, ..., 5 are supported"
        else:
            raise ValueError, "argument kind = %s not recognized" % str(kind)
                
        # Spline pre-filtering
        # This step is done because it is required by the ndimage code that I'm scavenging.
        # I don't fully understand why it must do this, and that's a problem.  But empirically
        # this step is needed in order to get good-looking data.
        if self.order >1:
            self._data_array = self._spline_filter(data, self.order)
        else:
            self._data_array = data
        
        # storing relevant data
        self.ndim = self._data_array.ndim
        self.shape_as_column_vec = np.reshape(np.shape(self._data_array), 
                                                                (self.ndim, 1))
        self._spacings = spacings
        self._min_coords = starting_coords
        self.out = out
        
    def __call__(self, coordinates):
        """ coordinates is an n x L array, where n is the dimensionality of the data
            and L is number of points.  That is, each column of coordinates
            indicates a point at which to interpolate.
        """
        
        # format checking
        coordinates = array(coordinates)
        if coordinates.ndim == 1: # passed in a single point
            coordinates = np.reshape(coordinates, ( self.ndim, 1))
        assert coordinates.ndim == 2, "Coordinates must be 1 or 2 dimensional"
        n, num_points = coordinates.shape
        assert n == self.ndim, "The first dimension of the input \
                must be as long as the dimensionality of the space"
        
        # converting from points in ND space to array indices
        indices = (coordinates - self._min_coords)/self._spacings
        
        if self.out in ['nearest', 'wrap', 'reflect', 'mirror', 'constant']:
            # out of bounds can be performed by _interpolate_array_entry
            result = self._interpolate_array_entry(self._data_array, indices, self.order, out = self.out)
        else:
            # need to return NaN when entry is out of bounds
            in_bounds_mask = self._index_in_bounds(indices)
            in_bounds = indices[:, in_bounds_mask]
            out_bounds = indices[:, ~in_bounds_mask]
            
            result = np.zeros(num_points)
            result[in_bounds_mask] = \
                self._interpolate_array_entry(self._data_array, indices[:,in_bounds_mask], self.order)
            result[~in_bounds_mask] = NaN
            
        return result
        
    
    def _interpolate_array_entry(self, data_array, indices, order, out='nearest'):
        """ indices is nxL matrix, where n is data_array.ndim
            returns array of length L giving interpolated entries.
        """
        
        extrap_code_register = { 'nearest':0,
                                        'wrap': 1,
                                        'reflect':2,
                                        'mirror':3,
                                        'constant':4,
                                        }
        
        n, L = np.shape(indices)
        
        output = np.zeros( L , dtype=np.float64 ) # place to store the data

        # geometric transform takes data_array, interpolates its values at indices, and
        # stores those values in output.  Other parameters give details of interpolation method.
        _nd_image.geometric_transform(data_array, None, indices, None, None, \
               output, order, extrap_code_register[out], 0.0, None, None)
               
        return output
        
    def _index_in_bounds(self, indices):
        """ return an array of bools saying which
            points are in interpolation bounds
        """
        
        # entry is 1 if that coordinate of a point is in its bounds
        index_in_bounds = (0 <= indices) & \
                                    (indices <= self.shape_as_column_vec)
        
        # for each point, number of coordinates that are in bounds
        num_indices_in_bounds = np.sum(index_in_bounds, axis=0)
        
        # True if each coordinate for the point is in bounds
        return num_indices_in_bounds == self.ndim

    def _spline_filter(self, data, order = 3):
        """ This step is required by the extension module.  I (Field Cady)
            do not fully understand why.
        
            Multi-dimensional spline filter.

            Note: The multi-dimensional filter is implemented as a sequence of
            one-dimensional spline filters. The intermediate arrays are stored
            in the same data type as the output. Therefore, for output types
            with a limited precision, the results may be imprecise because
            intermediate results may be stored with insufficient precision.
            
            order must be from 0 to 5
        """
                                                        
        # allocate memory in which to put the output
        output = np.zeros( np.shape(data) , dtype=np.float64 )
        
        if order in [0, 1]:
            output = data
        else:
            for axis in range(data.ndim):
                _nd_image.spline_filter1d(data, order, axis, output)
                data = output
        return output