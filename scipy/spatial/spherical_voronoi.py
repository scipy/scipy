"""
Spherical Voronoi Code

.. versionadded:: 0.17.0

"""
#
# Copyright (C)  Tyler Reddy, Ross Hemsley, Edd Edmondson,
#                Nikolai Nowaczyk, Joe Pitt-Francis, 2015.
#
# Distributed under the same BSD license as Scipy.
#

import numpy as np
import numpy.matlib
import scipy
import math

__all__ = ['SphericalVoronoi']


def calc_circumcenters(tetrahedrons):
    """ Calculates the cirumcenters of the circumspheres of tetrahedrons.

    An implementation based on
    http://mathworld.wolfram.com/Circumsphere.html

    Parameters
    ----------
    tetrahedrons : an array of shape (N, 4, 3)
        consisting of N tetrahedrons defined by 4 points in 3D

    Returns
    ----------
    circumcenters : an array of shape (N, 3)
        consisting of the N circumcenters of the tetrahedrons in 3D

    """

    num = tetrahedrons.shape[0]
    a = np.concatenate((tetrahedrons, np.ones((num, 4, 1))), axis=2)

    sums = np.sum(tetrahedrons ** 2, axis=2)
    d = np.concatenate((sums[:, :, np.newaxis], a), axis=2)

    dx = np.linalg.det(np.delete(d, 1, axis=2))
    dy = -np.linalg.det(np.delete(d, 2, axis=2))
    dz = np.linalg.det(np.delete(d, 3, axis=2))

    nominator = np.vstack((dx, dy, dz))
    denominator = np.matlib.repmat(2 * np.linalg.det(a), 3, 1)
    return (nominator / denominator).transpose()


def convert_cartesian_array_to_spherical_array(coord_array,angle_measure='radians'):
    '''Take shape (N,3) cartesian coord_array and return an array of the same shape in spherical polar form (r, theta, phi). Based on StackOverflow response: http://stackoverflow.com/a/4116899
    use radians for the angles by default, degrees if angle_measure == 'degrees' '''
    spherical_coord_array = np.zeros(coord_array.shape)
    xy = coord_array[...,0]**2 + coord_array[...,1]**2
    spherical_coord_array[...,0] = np.sqrt(xy + coord_array[...,2]**2)
    spherical_coord_array[...,1] = np.arctan2(coord_array[...,1], coord_array[...,0])
    spherical_coord_array[...,2] = np.arccos(coord_array[...,2] / spherical_coord_array[...,0])
    if angle_measure == 'degrees':
        spherical_coord_array[...,1] = np.degrees(spherical_coord_array[...,1])
        spherical_coord_array[...,2] = np.degrees(spherical_coord_array[...,2])
    return spherical_coord_array

def calculate_surface_area_of_a_spherical_Voronoi_polygon(array_ordered_Voronoi_polygon_vertices,sphere_radius):
    '''Calculate the surface area of a polygon on the surface of a sphere. Based on equation provided here: http://mathworld.wolfram.com/LHuiliersTheorem.html
    Decompose into triangles, calculate excess for each'''
    #have to convert to unit sphere before applying the formula
    spherical_coordinates = convert_cartesian_array_to_spherical_array(array_ordered_Voronoi_polygon_vertices)
    spherical_coordinates[...,0] = 1.0
    array_ordered_Voronoi_polygon_vertices = convert_spherical_array_to_cartesian_array(spherical_coordinates)
    #handle nearly-degenerate vertices on the unit sphere by returning an area close to 0 -- may be better options, but this is my current solution to prevent crashes, etc.
    #seems to be relatively rare in my own work, but sufficiently common to cause crashes when iterating over large amounts of messy data
    if scipy.spatial.distance.pdist(array_ordered_Voronoi_polygon_vertices).min() < (10 ** -7):
        return 10 ** -8
    else:
        n = array_ordered_Voronoi_polygon_vertices.shape[0]
        #point we start from
        root_point = array_ordered_Voronoi_polygon_vertices[0]
        totalexcess = 0
        #loop from 1 to n-2, with point 2 to n-1 as other vertex of triangle
        # this could definitely be written more nicely
        b_point = array_ordered_Voronoi_polygon_vertices[1]
        root_b_dist = calculate_haversine_distance_between_spherical_points(root_point, b_point, 1.0)
        for i in 1 + np.arange(n - 2):
            a_point = b_point
            b_point = array_ordered_Voronoi_polygon_vertices[i+1]
            root_a_dist = root_b_dist
            root_b_dist = calculate_haversine_distance_between_spherical_points(root_point, b_point, 1.0)
            a_b_dist = calculate_haversine_distance_between_spherical_points(a_point, b_point, 1.0)
            s = (root_a_dist + root_b_dist + a_b_dist) / 2
            totalexcess += 4 * math.atan(math.sqrt( math.tan(0.5 * s) * math.tan(0.5 * (s-root_a_dist)) * math.tan(0.5 * (s-root_b_dist)) * math.tan(0.5 * (s-a_b_dist))))
        return totalexcess * (sphere_radius ** 2)

def calculate_surface_area_of_planar_polygon_in_3D_space(array_ordered_Voronoi_polygon_vertices):
    '''Based largely on: http://stackoverflow.com/a/12653810
    Use this function when spherical polygon surface area calculation fails (i.e., lots of nearly-coplanar vertices and negative surface area).'''
    #unit normal vector of plane defined by points a, b, and c
    def unit_normal(a, b, c):
        x = np.linalg.det([[1,a[1],a[2]],
             [1,b[1],b[2]],
             [1,c[1],c[2]]])
        y = np.linalg.det([[a[0],1,a[2]],
             [b[0],1,b[2]],
             [c[0],1,c[2]]])
        z = np.linalg.det([[a[0],a[1],1],
             [b[0],b[1],1],
             [c[0],c[1],1]])
        magnitude = (x**2 + y**2 + z**2)**.5
        return (x/magnitude, y/magnitude, z/magnitude)

    #area of polygon poly
    def poly_area(poly):
        '''Accepts a list of xyz tuples.'''
        assert len(poly) >= 3, "Not a polygon (< 3 vertices)."
        total = np.zeros(3)
        N = len(poly)
        for i in range(N):
            vi1 = poly[i]
            vi2 = poly[(i+1) % N]
            prod = np.cross(vi1, vi2)
            total += prod
        result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
        return abs(result/2)

    list_vertices = list(map(tuple, array_ordered_Voronoi_polygon_vertices))
    planar_polygon_surface_area = poly_area(list_vertices)
    return planar_polygon_surface_area

def calculate_and_sum_up_inner_sphere_surface_angles_Voronoi_polygon(array_ordered_Voronoi_polygon_vertices,sphere_radius):
    '''Takes an array of ordered Voronoi polygon vertices (for a single generator) and calculates the sum of the inner angles on the sphere surface. The resulting value is theta in the equation provided here: http://mathworld.wolfram.com/SphericalPolygon.html '''
    num_vertices_in_Voronoi_polygon = array_ordered_Voronoi_polygon_vertices.shape[0] #the number of rows == number of vertices in polygon
    spherical_polar_polygon_vertices = convert_cartesian_array_to_spherical_array(array_ordered_Voronoi_polygon_vertices)
    min_vertex_radius = spherical_polar_polygon_vertices[...,0].min()
    assert sphere_radius - min_vertex_radius < 0.1, "The minimum projected Voronoi vertex r value should match the sphere_radius of {sphere_radius}, but got {r_min}.".format(sphere_radius=sphere_radius,r_min=min_vertex_radius)
    #two edges (great circle arcs actually) per vertex are needed to calculate tangent vectors / inner angle at that vertex
    current_vertex_index = 0
    list_Voronoi_poygon_angles_radians = []
    while current_vertex_index < num_vertices_in_Voronoi_polygon:
        current_vertex_coordinate = array_ordered_Voronoi_polygon_vertices[current_vertex_index]
        if current_vertex_index == 0:
            previous_vertex_index = num_vertices_in_Voronoi_polygon - 1
        else:
            previous_vertex_index = current_vertex_index - 1
        if current_vertex_index == num_vertices_in_Voronoi_polygon - 1:
            next_vertex_index = 0
        else:
            next_vertex_index = current_vertex_index + 1
        #try using the law of cosines to produce the angle at the current vertex (basically using a subtriangle, which is a common strategy anyway)
        current_vertex = array_ordered_Voronoi_polygon_vertices[current_vertex_index] 
        previous_vertex = array_ordered_Voronoi_polygon_vertices[previous_vertex_index]
        next_vertex = array_ordered_Voronoi_polygon_vertices[next_vertex_index] 
        a = calculate_Vincenty_distance_between_spherical_points(current_vertex,next_vertex,sphere_radius)
        b = calculate_Vincenty_distance_between_spherical_points(next_vertex,previous_vertex,sphere_radius)
        c = calculate_Vincenty_distance_between_spherical_points(previous_vertex,current_vertex,sphere_radius)
        pre_acos_term = (math.cos(b) - math.cos(a)*math.cos(c)) / (math.sin(a)*math.sin(c))
        current_vertex_inner_angle_on_sphere_surface = math.acos(pre_acos_term)
        list_Voronoi_poygon_angles_radians.append(current_vertex_inner_angle_on_sphere_surface)
        current_vertex_index += 1

    if abs(pre_acos_term) > 1.0:
        theta = 0
    else:
        theta = np.sum(list_Voronoi_poygon_angles_radians)
    return theta 

def calculate_haversine_distance_between_spherical_points(cartesian_array_1,cartesian_array_2,sphere_radius):
    '''Calculate the haversine-based distance between two points on the surface of a sphere. Should be more accurate than the arc cosine strategy. See, for example: http://en.wikipedia.org/wiki/Haversine_formula'''
    spherical_array_1 = convert_cartesian_array_to_spherical_array(cartesian_array_1)
    spherical_array_2 = convert_cartesian_array_to_spherical_array(cartesian_array_2)
    lambda_1 = spherical_array_1[1]
    lambda_2 = spherical_array_2[1]
    phi_1 = spherical_array_1[2]
    phi_2 = spherical_array_2[2]
    #we rewrite the standard Haversine slightly as long/lat is not the same as spherical coordinates - phi differs by pi/4
    spherical_distance = 2.0 * sphere_radius * math.asin(math.sqrt( ((1 - math.cos(phi_2-phi_1))/2.) + math.sin(phi_1) * math.sin(phi_2) * ( (1 - math.cos(lambda_2-lambda_1))/2.)  ))
    return spherical_distance

def calculate_Vincenty_distance_between_spherical_points(cartesian_array_1,cartesian_array_2,sphere_radius):
    '''Apparently, the special case of the Vincenty formula (http://en.wikipedia.org/wiki/Great-circle_distance) may be the most accurate method for calculating great-circle distances.'''
    spherical_array_1 = convert_cartesian_array_to_spherical_array(cartesian_array_1)
    spherical_array_2 = convert_cartesian_array_to_spherical_array(cartesian_array_2)
    lambda_1 = spherical_array_1[1]
    lambda_2 = spherical_array_2[1]
    phi_1 = spherical_array_1[2]
    phi_2 = spherical_array_2[2]
    delta_lambda = abs(lambda_2 - lambda_1)
    delta_phi = abs(phi_2 - phi_1)
    radian_angle = math.atan2( math.sqrt( (math.sin(phi_2)*math.sin(delta_lambda))**2 + (math.sin(phi_1)*math.cos(phi_2) - math.cos(phi_1)*math.sin(phi_2)*math.cos(delta_lambda)  )**2 ),  (math.cos(phi_1) * math.cos(phi_2) + math.sin(phi_1) * math.sin(phi_2) * math.cos(delta_lambda) ) )
    spherical_distance = sphere_radius * radian_angle
    return spherical_distance

def convert_spherical_array_to_cartesian_array(spherical_coord_array,angle_measure='radians'):
    '''Take shape (N,3) spherical_coord_array (r,theta,phi) and return an array of the same shape in cartesian coordinate form (x,y,z). Based on the equations provided at: http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#From_spherical_coordinates
    use radians for the angles by default, degrees if angle_measure == 'degrees' '''
    cartesian_coord_array = np.zeros(spherical_coord_array.shape)
    #convert to radians if degrees are used in input (prior to Cartesian conversion process)
    if angle_measure == 'degrees':
        spherical_coord_array[...,1] = np.deg2rad(spherical_coord_array[...,1])
        spherical_coord_array[...,2] = np.deg2rad(spherical_coord_array[...,2])
    #now the conversion to Cartesian coords
    cartesian_coord_array[...,0] = spherical_coord_array[...,0] * np.cos(spherical_coord_array[...,1]) * np.sin(spherical_coord_array[...,2])
    cartesian_coord_array[...,1] = spherical_coord_array[...,0] * np.sin(spherical_coord_array[...,1]) * np.sin(spherical_coord_array[...,2])
    cartesian_coord_array[...,2] = spherical_coord_array[...,0] * np.cos(spherical_coord_array[...,2])
    return cartesian_coord_array

def generate_random_array_spherical_generators(num_generators,sphere_radius,prng_object):
    '''Recoded using standard uniform selector over theta and acos phi, http://mathworld.wolfram.com/SpherePointPicking.html
    Same as in iPython notebook version'''
    u = prng_object.uniform(low=0,high=1,size=num_generators)
    v = prng_object.uniform(low=0,high=1,size=num_generators)
    theta_array = 2 * math.pi * u
    phi_array = np.arccos((2*v - 1.0))
    r_array = sphere_radius * np.ones((num_generators,))
    spherical_polar_data = np.column_stack((r_array,theta_array,phi_array))
    cartesian_random_points = convert_spherical_array_to_cartesian_array(spherical_polar_data)
    return cartesian_random_points


class SphericalVoronoi:
    """ Voronoi diagrams on the surface of a sphere.

    .. versionadded:: 0.17.0

    Attributes
    ----------
        points : double array of shape (npoints, 3)
            the points in 3D to generate the Voronoi diagram from
        radius : double
            radius of the sphere
            Default: None (forces estimation, which is less precise)
        center : double array of shape (3,)
            center of the sphere
            Default: None (assumes sphere is centered at origin)
        vertices : double array of shape (nvertices, 3)
            Voronoi vertices corresponding to points
        regions : list of list of integers of shape (npoints, _ )
            the n-th entry is a list consisting of the indices
            of the vertices belonging to the n-th point in points

    Notes
    ----------
    The spherical Voronoi diagram algorithm proceeds as follows. The Convex
    Hull of the input points (generators) is calculated, and is equivalent to
    their Delaunay triangulation on the surface of the sphere [Caroli]_.
    A 3D Delaunay tetrahedralization is obtained by including the origin of
    the coordinate system as the fourth vertex of each simplex of the Convex
    Hull. The circumcenters of all tetrahedra in the system are calculated and
    projected to the surface of the sphere, producing the Voronoi vertices.
    The Delaunay tetrahedralization neighbour information is then used to
    order the Voronoi region vertices around each generator. The latter
    approach is substantially less sensitive to floating point issues than
    angle-based methods of Voronoi region vertex sorting.

    The surface area of spherical polygons is calculated by decomposing them
    into triangles and using L'Huilier's Theorem to calculate the spherical
    excess of each triangle [Weisstein]_. The sum of the spherical excesses is
    multiplied by the square of the sphere radius to obtain the surface area
    of the spherical polygon. For nearly-degenerate spherical polygons an area
    of approximately 0 is returned by default, rather than attempting the
    unstable calculation.

    Empirical assessment of spherical Voronoi algorithm performance suggests
    quadratic time complexity (loglinear is optimal, but algorithms are more
    challenging to implement). The reconstitution of the surface area of the
    sphere, measured as the sum of the surface areas of all Voronoi regions,
    is closest to 100 % for larger (>> 10) numbers of generators.

    References
    ----------

    .. [Caroli] Caroli et al. Robust and Efficient Delaunay triangulations of
                points on or close to a sphere. Research Report RR-7004, 2009.
    .. [Weisstein] "L'Huilier's Theorem." From MathWorld -- A Wolfram Web
                Resource. http://mathworld.wolfram.com/LHuiliersTheorem.html

    See Also
    --------
    Voronoi : Conventional Voronoi diagrams in N dimensions.

    Examples
    --------

    Produce a Voronoi diagram for a pseudo-random set of points on the
    unit sphere:

    >>> import matplotlib
    >>> import matplotlib.pyplot as plt
    >>> import matplotlib.colors as colors
    >>> from mpl_toolkits.mplot3d import Axes3D
    >>> from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    >>> import numpy as np
    >>> import scipy as sp
    >>> from scipy.spatial import SphericalVoronoi
    >>> import scipy.spatial.spherical_voronoi
    >>> #pin down the pseudo random number generator (prng) object to avoid certain pathological generator sets
    >>> prng = np.random.RandomState(117) #otherwise, would need to filter the random data to ensure Voronoi diagram is possible
    >>> #produce 1000 random points on the unit sphere using the above seed
    >>> random_coordinate_array = scipy.spatial.spherical_voronoi.generate_random_array_spherical_generators(1000,1.0,prng)
    >>> #produce the Voronoi diagram data
    >>> voronoi_instance = SphericalVoronoi(random_coordinate_array,1.0)
    >>> dictionary_voronoi_polygon_vertices = voronoi_instance.calc_vertices()
    >>> #plot the Voronoi diagram
    >>> fig = plt.figure()
    >>> fig.set_size_inches(4,4)
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> ax.set_title('Voronoi Diagram on Sphere')
    >>> for generator_index, voronoi_region in dictionary_voronoi_polygon_vertices.items():
    ...    random_color = colors.rgb2hex(sp.rand(3))
    ...    #fill in the Voronoi region (polygon) that contains the generator:
    ...    polygon = Poly3DCollection([voronoi_region],alpha=1.0)
    ...    polygon.set_color(random_color)
    ...    ax.add_collection3d(polygon)
    >>> ax.set_xlim(-1,1);ax.set_ylim(-1,1);ax.set_zlim(-1,1);
    (-1, 1)
    (-1, 1)
    (-1, 1)
    >>> ax.set_xticks([-1,1]);ax.set_yticks([-1,1]);ax.set_zticks([-1,1]); #doctest: +ELLIPSIS
    [<matplotlib.axis.XTick object at 0x...>, <matplotlib.axis.XTick object at 0x...>]
    [<matplotlib.axis.XTick object at 0x...>, <matplotlib.axis.XTick object at 0x...>]
    [<matplotlib.axis.XTick object at 0x...>, <matplotlib.axis.XTick object at 0x...>]
    >>> plt.tick_params(axis='both', which='major', labelsize=6)
    >>> plt.tight_layout()
    >>> plt.show()

    Now, calculate the surface areas of the Voronoi region polygons and verify that the reconstituted surface area is sensible:

    >>> import math
    >>> dictionary_voronoi_polygon_surface_areas = voronoi_instance.voronoi_region_surface_areas_spherical_surface()
    >>> theoretical_surface_area_unit_sphere = 4 * math.pi
    >>> reconstituted_surface_area_Voronoi_regions = sum(dictionary_voronoi_polygon_surface_areas.itervalues())
    >>> percent_area_recovery = round((reconstituted_surface_area_Voronoi_regions / theoretical_surface_area_unit_sphere) * 100., 5)
    >>> print percent_area_recovery
    99.91979

    For completeness, produce the Delaunay triangulation on the surface of the unit sphere for the same data set:

    >>> Delaunay_triangles = voronoi_instance.delaunay_triangulation_spherical_surface()
    >>> Delaunay_triangles = random_coordinate_array[Delaunay_triangles]
    >>> fig2 = plt.figure()
    >>> fig2.set_size_inches(4,4)
    >>> ax = fig2.add_subplot(111, projection='3d')
    >>> ax.set_title('Delaunay Triangulation on Sphere')
    >>> for triangle_coordinate_array in Delaunay_triangles:
    ...     m = ax.plot(triangle_coordinate_array[...,0],triangle_coordinate_array[...,1],triangle_coordinate_array[...,2],c='r',alpha=0.1)
    ...     connecting_array = np.delete(triangle_coordinate_array,1,0)
    ...     n = ax.plot(connecting_array[...,0],connecting_array[...,1],connecting_array[...,2],c='r',alpha=0.1)
    >>> o = ax.scatter(random_coordinate_array[...,0],random_coordinate_array[...,1],random_coordinate_array[...,2],c='k',lw=0,s=0.9)
    >>> ax.set_xlim(-1,1);ax.set_ylim(-1,1);ax.set_zlim(-1,1);
    (-1, 1)
    (-1, 1)
    (-1, 1)
    >>> ax.set_xticks([-1,1]);ax.set_yticks([-1,1]);ax.set_zticks([-1,1]); #doctest: +ELLIPSIS
    [<matplotlib.axis.XTick object at 0x...>, <matplotlib.axis.XTick object at 0x...>]
    [<matplotlib.axis.XTick object at 0x...>, <matplotlib.axis.XTick object at 0x...>]
    [<matplotlib.axis.XTick object at 0x...>, <matplotlib.axis.XTick object at 0x...>]
    >>> plt.tick_params(axis='both', which='major', labelsize=6)
    >>> plt.tight_layout()
    >>> plt.show()
    """

    def __init__(self, points, radius=None, center=None):
        # translate generators such that sphere center is at origin
        if np.any(center):
            self.points = points - center
        else:
            self.points = points

        # estimate radius if not known
        if radius:
            self.radius = radius
        else:
            self.radius = np.average(
                scipy.spatial.distance.cdist(self.points,
                                             np.zeros((3,))[np.newaxis, :]))
        self.vertices = None
        self.regions = None

    def delaunay_triangulation_spherical_surface(self):
        '''Delaunay tessellation of the points on the surface of the sphere. This is simply the 3D convex hull of the points. Returns a shape (N,3,3) array of points representing the vertices of the Delaunay triangulation on the sphere (i.e., N three-dimensional triangle vertex arrays).'''
        hull = scipy.spatial.ConvexHull(self.points)
        array_points_vertices_Delaunay_triangulation = hull.simplices
        return array_points_vertices_Delaunay_triangulation

    def calc_vertices(self):
        """
        Returns a dictionary with the sorted (non-intersecting) polygon
        vertices for the Voronoi regions associated with each generator index.
        A dictionary entry would be structured as follows:
            {generator_index : array_polygon_vertices, ...}

        This algorithm was discussed at PyData London 2015 by
        Tyler Reddy, Ross Hemsley and Nikolai Nowaczyk
        """

        #step 1: perform 3D Delaunay triangulation on data set
        #        (here ConvexHull can also be used, and is faster)
        tri = scipy.spatial.ConvexHull(self.points)

        #step 2: add the origin to each of the simplices in tri to get the
        #        same tetrahedrons we'd have gotten from Delaunay
        #        tetrahedralization
        tetrahedrons = tri.points[tri.simplices]
        tetrahedrons = np.insert(tetrahedrons, 3, np.zeros((1, 3)), axis=1)

        #step 3: produce circumcenters of tetrahedrons from 3D Delaunay
        circumcenters = calc_circumcenters(tetrahedrons)

        #step 4: project tetrahedron circumcenters to the surface of the sphere
        #        to produce the Voronoi vertices
        lengths = scipy.spatial.distance.cdist(circumcenters, np.zeros((1, 3)))
        self.vertices = (self.radius / np.abs(lengths)) * circumcenters

        self.regions = [[k for k in range(0, len(tri.simplices))
                         if n in tri.simplices[k]]
                        for n in range(0, len(self.points))]

        #step 5: use the Delaunay tetrahedralization neighbour information to
        #        connect the Voronoi vertices around the generators, to
        #        produce the Voronoi regions
        return self._sort_vertices(tetrahedrons, tri)

    def _sort_vertices(self, tetrahedrons, tri):
        dictionary_sorted_Voronoi_point_coordinates_for_each_generator = {}
        array_tetrahedra = tetrahedrons
        generator_index_array = np.arange(self.points.shape[0])
        filter_tuple = np.where((np.expand_dims(tri.simplices, -1) == generator_index_array).any(axis=1))
        dictionary_generators_and_triangle_indices_containing_those_generators = {}
        for generator_index, triangle_index in zip(filter_tuple[1], filter_tuple[0]):
            dictionary_generators_and_triangle_indices_containing_those_generators.setdefault(generator_index, []).append(triangle_index)
        generator_index = 0
        for generator in tri.points[:-1]:
            indices_of_triangles_surrounding_generator = dictionary_generators_and_triangle_indices_containing_those_generators[generator_index]
            #pick any one of the triangles surrounding the generator and pick a non-generator vertex
            first_tetrahedron_index = indices_of_triangles_surrounding_generator[0]
            first_tetrahedron = array_tetrahedra[first_tetrahedron_index]
            first_triangle = first_tetrahedron[:-1,...]
            #pick one of the two non-generator vertices in the first triangle
            indices_non_generator_vertices_first_triangle = np.unique(np.where(first_triangle != generator)[0])
            ordered_list_tetrahedron_indices_surrounding_current_generator = [first_tetrahedron_index]
            #determine the appropriate ordering of Voronoi vertices to close the Voronoi region (polygon) by traversing the Delaunay neighbour data structure from scipy
            vertices_remaining = len(indices_of_triangles_surrounding_generator) - 1
            #choose the neighbour opposite the first non-generator vertex of the first triangle
            neighbour_tetrahedral_index = tri.neighbors[first_tetrahedron_index][indices_non_generator_vertices_first_triangle[0]]
            ordered_list_tetrahedron_indices_surrounding_current_generator.append(neighbour_tetrahedral_index)
            vertices_remaining -= 1

            #for all subsequent triangles it is the common non-generator vertex with the previous neighbour that should be used to propagate the connection chain to the following neighbour
            #the common vertex with the previous neighbour is the the vertex of the previous neighbour that was NOT used to locate the current neighbour
            #since there are only two candidate vertices on the previous neighbour and I've chosen to use the vertex with index 0, the remaining vertex on the previous neighbour is the non-generator vertex with index 1
            common_vertex_coordinate = first_triangle[indices_non_generator_vertices_first_triangle[1]]
            while vertices_remaining > 0:
                current_tetrahedron_index = ordered_list_tetrahedron_indices_surrounding_current_generator[-1]
                current_tetrahedron_coord_array = array_tetrahedra[current_tetrahedron_index]
                current_triangle_coord_array = current_tetrahedron_coord_array[:-1,...]
                indices_candidate_vertices_current_triangle_excluding_generator = np.unique(np.where(current_triangle_coord_array != generator)[0])
                array_candidate_vertices = current_triangle_coord_array[indices_candidate_vertices_current_triangle_excluding_generator]
                current_tetrahedron_index_for_neighbour_propagation = np.unique(np.where(current_tetrahedron_coord_array == common_vertex_coordinate)[0])
                next_tetrahedron_index_surrounding_generator = tri.neighbors[current_tetrahedron_index][current_tetrahedron_index_for_neighbour_propagation][0]
                common_vertex_coordinate = array_candidate_vertices[array_candidate_vertices != common_vertex_coordinate] #for the next iteration
                ordered_list_tetrahedron_indices_surrounding_current_generator.append(next_tetrahedron_index_surrounding_generator)
                vertices_remaining -= 1
            dictionary_sorted_Voronoi_point_coordinates_for_each_generator[generator_index] = self.vertices[ordered_list_tetrahedron_indices_surrounding_current_generator]
            generator_index += 1
        return dictionary_sorted_Voronoi_point_coordinates_for_each_generator

    def voronoi_region_surface_areas_spherical_surface(self):
        '''Returns a dictionary with the estimated surface areas of the Voronoi region polygons corresponding to each generator (original data point) index. An example dictionary entry: `{generator_index : surface_area, ...}`.'''
        dictionary_sorted_Voronoi_point_coordinates_for_each_generator = self.calc_vertices()
        dictionary_Voronoi_region_surface_areas_for_each_generator = {}
        for generator_index, Voronoi_polygon_sorted_vertex_array in dictionary_sorted_Voronoi_point_coordinates_for_each_generator.items():
            current_Voronoi_polygon_surface_area_on_sphere = calculate_surface_area_of_a_spherical_Voronoi_polygon(Voronoi_polygon_sorted_vertex_array,self.radius)
            assert current_Voronoi_polygon_surface_area_on_sphere > 0, "Obtained a surface area of zero for a Voronoi region."
            dictionary_Voronoi_region_surface_areas_for_each_generator[generator_index] = current_Voronoi_polygon_surface_area_on_sphere
        return dictionary_Voronoi_region_surface_areas_for_each_generator
