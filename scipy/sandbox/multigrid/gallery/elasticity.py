"""Linear Elasticity Examples"""

__all__ = [ 'linear_elasticity' ]

from scipy import array, matrix, ones, zeros, arange, empty, \
        hstack, vstack, tile, ravel, mgrid, concatenate, \
        cumsum
from scipy.linalg import inv, det
from scipy.sparse import coo_matrix, bsr_matrix
   

def linear_elasticity( grid, spacing=None, E=1e5, nu=0.3, format=None):
    if len(grid) == 2:
        return q12d( grid, spacing=spacing, E=E, nu=nu, format=format )
    else:
        raise NotImplemented,'no support for grid=%s' % str(grid)

def q12d( grid, spacing=None, E = 1e5, nu = 0.3, dirichlet_boundary = True, format=None):
    """ Q1 elements in 2 dimensions """
    X,Y = grid
    
    if X < 1 or Y < 1:
        raise ValueError,'invalid grid shape'

    if dirichlet_boundary:
        X += 1
        Y += 1

    pts = mgrid[0:X+1,0:Y+1]
    pts = hstack((pts[0].T.reshape(-1,1) - X/2.0, pts[1].T.reshape(-1,1) - Y/2.0)) 
    
    if spacing is None:
        DX,DY = 1,1
    else:
        DX,DY = spacing
        pts = [DX,DY]

    #compute local stiffness matrix
    lame = E * nu / ((1 + nu) * (1 - 2*nu)) # Lame's first parameter
    mu   = E / (2 + 2*nu)                   # shear modulus
    vertices = array([[0,0],[DX,0],[DX,DY],[0,DY]])
    K = stima4(vertices, lame, mu ) 
    
    nodes = arange( (X+1)*(Y+1) ).reshape(X+1,Y+1)
    LL = nodes[:-1,:-1]
    I  = (2*LL).repeat( K.size ).reshape(-1,8,8)
    J  = I.copy()
    I += tile( [0,1,2,3, 2*X + 4, 2*X + 5, 2*X + 2, 2*X + 3], (8,1) )
    J += tile( [0,1,2,3, 2*X + 4, 2*X + 5, 2*X + 2, 2*X + 3], (8,1) ).T
    V  = tile( K, (X*Y,1) )
     
    I = ravel(I)
    J = ravel(J)
    V = ravel(V)

    A = coo_matrix( (V,(I,J)), shape=(pts.size,pts.size) ).tocsr() #sum duplicates
    A = A.tobsr(blocksize=(2,2))

    del I,J,V,LL,nodes

    B = zeros(( 2 * (X+1)*(Y+1), 3))
    B[0::2,0] = 1
    B[1::2,1] = 1
    B[0::2,2] = -pts[:,1]
    B[1::2,2] =  pts[:,0]

    if dirichlet_boundary:
        mask = zeros((X+1, Y+1),dtype='bool')
        mask[1:-1,1:-1] = True
        mask = ravel(mask)
        data = zeros( ((X-1)*(Y-1),2,2) )
        data[:,0,0] = 1
        data[:,1,1] = 1
        indices = arange( (X-1)*(Y-1) )
        indptr = concatenate((array([0]),cumsum(mask)))
        P = bsr_matrix((data,indices,indptr),shape=(2*(X+1)*(Y+1),2*(X-1)*(Y-1)))
        Pt = P.T
        A = P.T * A * P

        B = Pt * B

    return A.asformat(format),B 

def stima4(vertices, lame, mu):
    """local stiffness matrix for two dimensional elasticity on a square element

    Material Parameters:
        - lame  : Lame's first parameter
        - mu : shear modulus

    Note:
        Vertices should be listed in counter-clockwise order:
            [3]----[2]
             |      |
             |      |
            [0]----[1]
        
        Degrees of freedom are enumerated as follows:
            [x=6,y=7]----[x=4,y=5]
                |            |
                |            |
            [x=0,y=1]----[x=2,y=3]
    """

    M    = lame + 2*mu # P-wave modulus

    R_11 = matrix([[  2, -2, -1,  1],
                   [ -2,  2,  1, -1],
                   [ -1,  1,  2, -2],
                   [  1, -1, -2,  2]]) / 6.0
    
    R_12 = matrix([[  1,  1, -1, -1],
                   [ -1, -1,  1,  1],
                   [ -1, -1,  1,  1],
                   [  1,  1, -1, -1]]) / 4.0
    
    R_22 = matrix([[  2,  1, -1, -2],
                   [  1,  2, -2, -1],
                   [ -1, -2,  2,  1],
                   [ -2, -1,  1,  2]]) / 6.0
     
    F = inv( vstack( (vertices[1] - vertices[0], vertices[3] - vertices[0]) ) )
    
    K = zeros((8,8)) # stiffness matrix
    
    E = F.T * matrix([[M, 0],[0, mu]]) * F
    K[0::2,0::2] = E[0,0] * R_11 + E[0,1] * R_12 + E[1,0] * R_12.T + E[1,1] * R_22
    
    E = F.T * matrix([[mu, 0],[0, M]]) * F
    K[1::2,1::2] = E[0,0] * R_11 + E[0,1] * R_12 + E[1,0] * R_12.T + E[1,1] * R_22
    
    E = F.T * matrix([[0, mu],[lame, 0]]) * F;
    K[1::2,0::2] = E[0,0] * R_11 + E[0,1] * R_12 + E[1,0] * R_12.T + E[1,1] * R_22
    
    K[0::2,1::2] = K[1::2,0::2].T
    
    K /= det(F)

    return K


   
   
   
