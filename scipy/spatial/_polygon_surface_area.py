"""

Polygon Surface Area Code

.. versionadded:: 1.0.0

"""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.spatial.distance import pdist, cdist
from . import _surface_area
from six.moves import xrange
import math

#
# Copyright (C)  James Nichols and Tyler Reddy
# Contributions from Tyler Reddy are owned by LANL
# and respect the scipy distribution license
#
# Distributed under the same BSD license as Scipy.

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

def _vertex_index_strider(index, num_vertices):
    # handle the wrapping / iteration over
    # polygon vertices in either CW or CCW
    # sort order
    forward_index = index + 1
    backward_index = index - 1
    if forward_index >= num_vertices:
        forward_index = 0
    return forward_index, backward_index

def _planar_polygon_area(vertices):
    num_vertices = vertices.shape[0]
    area_sum = 0
    for i in xrange(num_vertices):
        forward_index, backward_index = _vertex_index_strider(i, num_vertices)
        delta_x = (vertices[forward_index][0] -
                   vertices[backward_index][0])
        area_sum += delta_x * vertices[i][1]
    area = -0.5 * area_sum
    return area

def _slerp(start_coord,
           end_coord,
           n_pts):
    # spherical linear interpolation between points
    # on great circle arc
    # see: https://en.wikipedia.org/wiki/Slerp#Geometric_Slerp
    # NOTE: could we use scipy.interpolate.RectSphereBivariateSpline instead?
    omega = np.arccos(np.dot(start_coord, end_coord))
    t_values = np.linspace(0, 1, n_pts)
    new_pts = []
    for t in t_values:
        new_pt = (((np.sin((1 - t) * omega) / np.sin(omega)) * start_coord) +
                  ((np.sin(t * omega) / np.sin(omega)) * end_coord))
        new_pts.append(new_pt)
    return np.array(new_pts)

def _spherical_polygon_area(vertices, radius, discretizations):
    num_vertices = vertices.shape[0]
    area_sum = 0

    for i in xrange(num_vertices):
        new_pts = _slerp(vertices[i], vertices[i-1], discretizations)

        lambda_range = np.arctan2(new_pts[...,1], new_pts[...,0])
        phi_range = np.arcsin((new_pts[...,2]))
        # NOTE: early stage handling of antipodes
        if phi_range[0] == phi_range[-1]:
            phi_range[:] = phi_range[0]
        area_element = 0
        for j in xrange(discretizations - 1):
            delta_lambda = (lambda_range[j+1] -
                            lambda_range[j])

            # at the + / - pi transition point
            # of the unit circle
            # add or subtract 2 * pi to the delta
            # based on original sign
            if delta_lambda > np.pi:
                delta_lambda -= 2 * np.pi
            elif delta_lambda < (-np.pi):
                delta_lambda += 2 * np.pi

            second_term = 2 + np.sin(phi_range[j]) + np.sin(phi_range[j+1])
            area_element += (delta_lambda * second_term * (radius ** 2) * 0.5)
        area_sum += area_element
    area = area_sum
    return area

def pole_in_polygon(vertices):
    # determine if the North or South Pole is contained
    # within the spherical polygon defined by vertices
    # the implementation is based on discussion
    # here: https://blog.element84.com/determining-if-a-spherical-polygon-contains-a-pole.html
    # and the course calculation component of the aviation
    # formulary here: http://www.edwilliams.org/avform.htm#Crs

    #TODO: handle case where initial point is within
    # i.e., expected machine precision of a pole

    #TODO: handle case where both the N & S poles
    # are contained within a given polygon

    course_angle_sum = 0
    lambda_range = np.arctan2(vertices[...,1], vertices[...,0])
    phi_range = np.arcsin((vertices[...,2]))
    for i in range(vertices.shape[0] - 1):
        if i == (vertices.shape[0] - 1):
            next_index = 0
        else:
            next_index = i + 1

        lambda_1 = lambda_range[i]
        lambda_2 = lambda_range[next_index]
        phi_1 = phi_range[i]
        phi_2 = phi_range[next_index]
        delta_lambda = lambda_1 - lambda_2

        term_1 = math.sin(delta_lambda) * math.cos(phi_2)
        term_2 = math.cos(phi_1) * math.cos(phi_2)
        term_3 = math.sin(phi_1) * math.cos(phi_2) * math.cos(delta_lambda)
        result = math.atan2(term_1, term_2 - term_3)
        course_angle = result % (2 * math.pi)
        course_angle_sum += course_angle

    if abs(course_angle_sum) == 360:
        # there's a pole in the polygon
        return 1
    else:
        # no single pole in the polygon
        return 0

def poly_area(vertices, radius=None, threshold=1e-21,
              cython=None, discretizations=500):
    # calculate the surface area of a planar or spherical polygon
    # crude pure Python implementation for handling a single
    # polygon at a time
    # based on JPL Publication 07-3 by Chamberlain and Duquette (2007)
    # for planar polygons we currently still require x,y,z coords
    # can just set i.e., z = 0 for all vertices
    num_vertices = vertices.shape[0]
    # require that a planar or spherical triangle is
    # the smallest possible input polygon
    if num_vertices < 3:
        err_str = "An input polygon must have at least 3 vertices."
        raise ValueError(err_str)

    min_vertex_dist = pdist(vertices).min()
    if min_vertex_dist < threshold:
        err_str = '''Duplicate vertices detected based on minimum
                     distance {min_vertex_dist} and threshold value
                     {threshold}.'''.format(min_vertex_dist=min_vertex_dist,
                                            threshold=threshold)
        raise ValueError(err_str)


    if radius is not None: # spherical polygons
        if radius <= threshold:
            err_str = 'radius must be > {threshold}'.format(threshold=threshold)
            raise ValueError(err_str)

        # if any two *consecutive* vertices in the spherical polygon
        # are antipodes, there's an infinite set of
        # geodesics connecting them, and we cannot
        # possibly guess the appropriate surface area
        dist_matrix = cdist(vertices, vertices)
        # TODO: use a threshold for the floating
        # point dist comparison here
        matches = np.where(dist_matrix == (2. * radius))
        # can't handle consecutive matches that are antipodes
        if np.any(np.abs(matches[0] - matches[1]) == 1):
            raise ValueError("Consecutive antipodal vertices are ambiguous.")

        # a great circle can only occur if the origin lies in a plane
        # with all vertices of the input spherical polygon
        # the input spherical polygon must have at least three vertices
        # to begin with (i.e., spherical triangle as smallest possible
        # spherical polygon) so we can safely assume (and check above)
        # three input vertices plus the origin for our test here

        # points lie in a common plane if the determinant below
        # is zero (informally, see:
        # https://math.stackexchange.com/a/684580/480070 
        # and http://mathworld.wolfram.com/Plane.html eq. 18)

        # have to iterate through the vertices in chunks
        # to preserve the three-point form of the plane-check
        # determinant
        current_vert = 0
        plane_counts = 0
        plane_failures = 0
        while current_vert <= (num_vertices - 3):
            candidate_plane = np.concatenate((np.zeros((1,3)),
                                              vertices[current_vert:current_vert + 3]))
            one_pads = np.ones((1, candidate_plane.shape[0]))

            candidate_plane = np.concatenate((candidate_plane.T,
                                              one_pads))

            if np.linalg.det(candidate_plane) == 0:
                plane_counts += 1
            else:
                # if any of the vertices aren't on the plane with
                # the origin we can safely break out
                plane_failures += 1
                break

            current_vert += 1

        if plane_failures == 0:
            # we have a great circle so
            # return hemisphere area
            return 2. * np.pi * (radius ** 2) 

        # normalize vertices to unit sphere
        vertices = convert_cartesian_array_to_spherical_array(vertices)
        vertices[...,0] = 1.
        vertices = convert_spherical_array_to_cartesian_array(vertices)

        if cython is None:
            area = _spherical_polygon_area(vertices, radius, discretizations)
        else: # cython code for spherical polygon SA
            area = _spherical_polygon_area(vertices, radius, discretizations)
            #area = _surface_area.spherical_polygon_area(vertices, radius,
                                                        #discretizations)
    else: # planar polygon
        if cython is None:
            area = _planar_polygon_area(vertices)
        else: # cython code for planar polygon SA
            area = _surface_area.planar_polygon_area(vertices)
    
    return abs(area)






