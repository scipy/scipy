# implements algorithm found at 
# http://vcg.isti.cnr.it/publications/papers/dewall.pdf

# WARNING
# This code is grotesquely inefficient and messy.
# The goal is to have it be 1) correct, and 2) understandable
# The next phase is to rewrite it using arrays and efficient algorithms.
# After that, the slow parts should be translated into C.

# In particular, calculation of the circumcircle is a purely
# mathematical operation that really should be made into C.

import numpy as np
from numpy.linalg import norm

eps =1e-5
compare_first_elem = lambda t1, t2 : cmp(t1[0], t2[0])
point_in_list = lambda p, L : np.any([ (p==elem).all() for elem in L ])
face_in_list = lambda f, L :  np.any([ np.alltrue([point_in_list(p,face) for p in f]) for face in L])
compare_pointlists = lambda L1, L2 : np.alltrue([ point_in_list(l1, L2) for l1 in L1]) and \
                                                    np.alltrue([ point_in_list(l2, L1) for l2 in L2])

def dewall (P, #set of points 
                AFL = [], # list of faces: (d-1)face list
                ):
    
    # checking input
    assert isinstance(P, list)
    if len(P)>0:
        assert isinstance(P[0], np.ndarray)
    assert isinstance(AFL, list)
    if len(AFL)>0: 
        assert isinstance(AFL[0],tuple)
        assert isinstance(AFL[0][0], list)
        assert isinstance(AFL[0][0][0], np.ndarray)
        assert isinstance(AFL[0][1], np.ndarray)
    
    # base case
    if len(P) == 0: return []
    if len(P) <= len(P[0]): return [] # better to return [] whenever points are co-hyper-planar
    
    # lists of active faces
    # elem = ( [list of d points] , outvec )
    AFL_alpha = []
    AFL1 = []
    AFL2 = []

    # list of simplices
    Sigma= []
    
    alpha = select_alpha(P, AFL)
                            
    # divide points into two sets separated by alpha
    P1, P2 = pointset_partition(P, alpha) # both lists of points
        
    # Simplex Wall Construction
    just_starting = False #source of problem?
    if len(AFL) == 0: # This is only executed once, at the start of the algorithm
        just_starting = True
        first_simplex = make_first_simplex(P, alpha)
        AFL = [ (face, get_out_vec(face,first_simplex))\
                for face in faces(first_simplex)] # d+1 of them
        Sigma.append(first_simplex)
    for face, outvec in AFL:
        if is_intersected(face, alpha): # not counting as intersected
            AFL_alpha.append((face, \
                            get_out_vec(face,first_simplex) if just_starting \
                            else outvec))
        if is_subset(face, P1):
            AFL1.append((face,outvec))
        if is_subset(face, P2):
            AFL2.append((face,outvec))
    while len(AFL_alpha) != 0:
        
        face, outvec = AFL_alpha.pop()
        
        if outvec is not None:
            outward_points = filter( lambda p: (np.dot(p,outvec)>np.dot(face[0],outvec)),\
                                            P)
        else:
            outward_points = []#filter( lambda p: np.all([not point_in_list(p,vertex) for vertex in face]), P)
                
        t = make_simplex(face, outward_points) # make only over outer half space
        
        if t is not None:
            Sigma.append(t)
            # for f0 != f in faces(t) , ie for new outward faces
            for f0 in filter(lambda f: not face_in_list(f,[face]), faces(t)):
                new_pair = (f0, get_out_vec(f0, t))
                if is_intersected(f0, alpha):
                    # continue building wall out in that direction
                    AFL_alpha = update(new_pair, AFL_alpha)
                    #np.random.shuffle(AFL_alpha)
                if is_subset(f0, P1):
                    AFL1 = update(new_pair, AFL1)
                if is_subset(f0, P2):
                    AFL2 = update(new_pair, AFL2)
                    
    # now Sigma contains all simplices that intersect alpha
    
    # Recursive Triangulation
    if len(AFL2)!=0:
        Sigma = Sigma + dewall(P2,AFL2)
    if len(AFL1)!=0:
        Sigma = Sigma + dewall(P1,AFL1) # not doing this branch of recursion
    return Sigma

def select_alpha(P, AFL):
    # dividing plane.  This must divide points into 2 non-empty parts
    # must intersect active face if there is one
    # must not intersect any points
    
    d = len(P[0])
    
    # plane through avg of cluster.  Guarantees separation
    # must also make sure intersects face in AFL
    if len(AFL) != 0:
        mid = sum(AFL[0][0])/(d)
    else:
        mid = sum(P)/len(P)
    if norm(mid)==0:
        direction = np.random.random_sample((d))
        alpha = ( direction/norm(direction), 0)
    else: 
        alpha =(mid/norm(mid), norm(mid))
        
    return alpha

def get_out_vec(face, simplex): #SEEMS GOOD
    # face is a face of simplex
    # returns vector normal to face pointing out of simplex
    assert len(face)==len(face[0]), "face needs d points"
    assert len(simplex)-len(face)==1, "simplex needs one more point"
    
    d = len(face[0]) # dimension
    
    # get point in simplex that's not in face
    # vector must face away from it
    other_point = filter(lambda p: not point_in_list(p, face),
                                simplex)[0]
    
    # determinants of submatrices give components of a vector
    # that is orthogonal to each row
    vecs_along_face = np.array([point-face[0] for point in face[1:]])
    normal_vector = np.zeros(d)
    
    # normal vector must be orthogonal to each element of matrix_of_diffs
    # components given by determinants of submatrices
    for i in range(d):
        normal_vector[i] = ((-1)**i)*\
                np.linalg.det(vecs_along_face[:,np.arange(d)!=i])
    unit_normal = normal_vector/norm(normal_vector)
    
    # return unit_normal point in correct direction
    if np.dot(other_point, unit_normal) < np.dot(face[0], unit_normal):
        return unit_normal
    else:
        return -1*unit_normal

def is_subset(S1, S2):
    # both are lists of arrays
    return np.alltrue([ point_in_list(s1, S2) for s1 in S1])

def update (face_pair, face_pair_list):
    # returns face_list with face_pair added if it wasn't there, else deleted
    face, outvec = face_pair
    face_list = [face for face, outvec in face_pair_list]
    if face_in_list(face, face_list):
        f_not_equal_face = lambda face_pair :  not np.alltrue([ point_in_list(p, face_pair[0]) for p in face ])
        face_pair_list = filter(f_not_equal_face, face_pair_list)
    else:
        face_pair_list.append(face_pair)
    return face_pair_list
        
def pointset_partition(P, alpha): #WORKS
    P1 = [p for p in P if np.dot(p,alpha[0])<  alpha[1]]
    P2 = [p for p in P if np.dot(p,alpha[0])>=alpha[1]]
    return P1, P2
    
def is_intersected(f, alpha): #WORKS
    assert isinstance(f, list), "the face must be a list: "+str(f)
    assert isinstance(alpha, tuple), "alpha must be tuple: "+str(alpha)
    assert isinstance(alpha[0],np.ndarray), "normal vector: "+str(alpha[0])
    assert isinstance(f[0],np.ndarray), "point in f: "+str(f[0])
    list_of_dot_prods = [np.dot(alpha[0], point) for point in f]
    all_dots_nonneg = [ prod >= alpha[1] for prod in list_of_dot_prods]
    all_dots_nonpos = [ prod <= alpha[1] for prod in list_of_dot_prods]
    return not (np.all(all_dots_nonneg) or np.all(all_dots_nonpos))
    
def faces(simplex): #WORKS
    # given simplex (as list of points) returns list of faces (each face a list of points)
    faces = []
    for point in simplex:
        faces.append( [p for p in simplex if (p!=point).any()] )
    return faces

def make_first_simplex(P, alpha): #WORKS
    # alpha = unit_normal_vec, distance
    # assume no points on plane

    d = len(P[0])
    unit_normal, distance = alpha
    points_and_distances = [(np.dot(unit_normal, point), point) for point in P]

    points_and_dist_from_alpha = [ (norm(distance-dist), point) for (dist, point) in points_and_distances]

    first_point = sorted(points_and_dist_from_alpha,
                                cmp = compare_first_elem \
                                )[0][1] # closest to alpha
                                                    
    possible_second_pts = [ (circumcircle([first_point, point])[1], point) \
                                        for point in P if (point != first_point).any() and \
                                        (np.dot(unit_normal, first_point)-distance)*(np.dot(unit_normal, point)-distance)<=0 \
                                    ]
    second_point = sorted(possible_second_pts, \
                                    cmp = compare_first_elem \
                                    )[0][1]
    simplex = [first_point, second_point]
    for i in range(d-1):
        radii_of_circumcircles = [(circumcircle(  copy_list(simplex)+[point.copy()]  )[1], point) \
                                        for point in P if not point_in_list(point, simplex) ]
        new_point = sorted(radii_of_circumcircles, \
                                    cmp = compare_first_elem \
                                    )[0][1]
        simplex.append(new_point)

    return simplex
    
def copy_list(list_of_arrays): #WORKS
    # returns a list of copies of the arrays
    # use if want to modify an array list without 
    result = [array.copy() for array in list_of_arrays]
    return result
    

        
def make_simplex(f,P): #WORKS
    # must be only in the outer halfspace
    
    # returns the simlex
    
    # clean up by making sure only points not in hyperplane(f) are in P
        
    valid_points = [p for p in P if not point_in_list(p,f)]
    if len(valid_points) == 0: 
        return None
    delaunay_distances = [(delaunay_dist(f,p) ,p) for p in valid_points]
    new_point = sorted(delaunay_distances, cmp = compare_first_elem)[0][1]
    simplex = f+[new_point]
    return simplex

def delaunay_dist(f, p): #WORKS
    # |distance| = radius of circumcircles
    # sign depends on whether center is on same side as p
    
    # FIXME : don't let pathological stuff come in
    
    normal, distance = spanning_hyperplane(f)
    center, radius = circumcircle( [elem.copy() for elem in f]+[p] ) #need copy of f and add to it
    if (np.dot(normal,center) >= distance and np.dot(normal,p) >= distance) \
        or (np.dot(normal,center) <= distance and np.dot(normal,p) <= distance) \
        or radius == np.Inf:
        return radius
    else:
        return -1*radius
        
def spanning_hyperplane(f): #WORKS
    # f=list of points
    # return outward unit normal to hyperplane they define
    # and distance to hyperplane from the origin
    d = len(f[0])
    
    # determinants of submatrices give components of a vector
    # that is orthogonal to each row
    matrix_of_diffs = np.array([point-f[0] for point in f[1:]])
    normal_vector = np.zeros(d)
    
    # normal vector must be orthogonal to each element of matrix_of_diffs
    # components given by determinants of submatrices
    for i in range(d):
        normal_vector[i] = ((-1)**i)*np.linalg.det(\
                                matrix_of_diffs[:,np.arange(d)!=i])
    unit_normal = normal_vector/norm(normal_vector)
    distance = np.dot(unit_normal, f[0])
    
    # want a positive distance
    if distance < 0:
        return -1*unit_normal, -1*distance
    return unit_normal, distance
    
def circumcircle(Points): #WORKS
    # P = list of distinct points
    # n distinct points, with n<= d+1
    # returns (center=array, radius)
    
    # FIXME : account for collinear/same and error check
    
    # calculates center by forming system of linear constraints
    # on its components, and calculates radius after that
    
    d = len(Points[0]) # dimensionality of space
    n = len(Points)
    
    # need "infinitely big" circle in this case
    if not linearly_independent(Points):
        return np.Inf*np.ones(d), np.Inf
    
    matrix = np.zeros((d,d), dtype=np.float64) # d constraints to find the center
    vector = np.zeros((d,1), dtype=np.float64)
    
    difference_vecs = [(Points[i]-Points[0])/norm(Points[i]-Points[0]) \
                                for i in range(1,n) if norm(Points[i]-Points[0])>eps] # n-1 difference vectors.  Point along hyperplane
    n = len(difference_vecs) + 1 # correct for possible degeneracies in the input
        
    # for each P[0] -- P[i] line segment, center is along the middle of it
    for i in range(n-1):
        matrix[i,:] = difference_vecs[i]
        vector[i,0] = np.dot( difference_vecs[i] , (Points[0]+Points[i+1])/2)
    
    # center is in the hyperplane defined by the points
    # ie its componeny orthogonal to difference_vecs is same as that of all points
    orthog_comp = orthogonal_complement(difference_vecs)
    m = len(orthog_comp) # m+n = d
    if m+(n-1) != d:
        assert False, "messed up.  m=%i, and n=%i " % (m, n)
    for i in range(m):
        matrix[n-1+i,:] = orthog_comp[i]
        vector[n-1+i,0] = np.dot(Points[0], orthog_comp[i])
    
    # calculating and returning
    center = np.linalg.solve(matrix, vector)
    radius = norm(center-np.reshape(Points[0],(d,1)))
    return np.reshape(center,(d)), radius
    
def orthogonal_complement(P): #WORKS
    # return orthogonal complement of P of unit length
    # in = list of 1D arrays
    # out = list of 1D arrays
    
    # nothing here can be pathological
    
    P = copy_list(P) # this way won't harm original P
    d = len(P[0]) # dimensionality
    n = len(P) # number of points.  <= d
    
    # Graham-Schmidt orthonormalization of input vectors
    for j in range(n):
        for k in range(j):
            P[j] = P[j] - np.dot(P[k],P[j])*P[k]
        if norm(P[j]) > eps:
            P[j]=P[j]/norm(P[j])
    
    # remove linear dependencies and normalize
    P = np.array([point/norm(point) for point in P if norm(point)>eps])
    
    orthog = np.eye(d,d)
    for j in range(d):
        # remove component of Oj in P
        for k in range(len(P)):
            orthog[j,:] = orthog[j,:] - np.dot(P[k,:],orthog[j,:])*P[k,:]
            
        # remove component of Oj in Oi with i<j
        for m in range(j):
            orthog[j,:] = orthog[j,:] - np.dot(orthog[m,:],orthog[j,:])*orthog[m,:]
            
        # normalize if vector not linearly dependent on previous
        if norm(orthog[j,:]) > eps:
            orthog[j,:] = orthog[j,:]/norm(orthog[j,:])
    
    # return each non-zero vector
    return [np.reshape(orthog[j,:],(d)) for j in range(d) if norm(orthog[j,:])>eps]
    
def linearly_independent(P):
    # return true if the vectors are linearly independent
    # needs n points
    d=len(P[0])
    if len(P)>d: return False
    matrix = np.array(P).reshape((d,len(P)))
    return np.rank(matrix)==len(P)
    
    
    
    
    