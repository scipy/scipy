
from scipy import *
from pydec import *
from pydec.multigrid import *
from pydec.multigrid.discrete_laplacian import boundary_hierarchy, discrete_laplacian_solver, hodge_solver

from scipy.sandbox.multigrid import smoothed_aggregation_solver
from scipy.sandbox.multigrid.utils import expand_into_blocks


## Load mesh from file
mesh_path = '../../../../../hirani_group/wnbell/meshes/'
#mesh = read_mesh(mesh_path + 'rocket/rocket.xml')
#mesh = read_mesh(mesh_path + 'genus3/genus3_168k.xml')
#mesh = read_mesh(mesh_path + 'genus3/genus3_455k.xml')
mesh = read_mesh(mesh_path + '/torus/torus.xml')
for i in range(3):
    mesh['vertices'],mesh['elements'] = loop_subdivision(mesh['vertices'],mesh['elements'])
cmplx = simplicial_complex(mesh['vertices'],mesh['elements'])

## Construct mesh manually
#bitmap = ones((60,60),dtype='bool')
#bitmap[1::10,1::10] = False
#bitmap[100:150,100:400] = False
#cmplx = regular_cube_complex(regular_cube_mesh(bitmap))



def whitney_innerproduct_cache(cmplx,k):
    h = hash(cmplx.vertices.tostring()) ^ hash(cmplx.simplices.tostring()) ^ hash(k)

    filename = "/home/nathan/.pydec/cache/whitney_" + str(h) + ".mtx"

    try:
        import pydec
        M = pydec.io.read_array(filename)
    except:
        import pydec
        M = whitney_innerproduct(cmplx,k)
        pydec.io.write_array(filename,M)

    return M



def cube_innerproduct_cache(cmplx,k):
    h = hash(cmplx.mesh.bitmap.tostring()) ^ hash(cmplx.mesh.bitmap.shape) ^ hash(k)

    filename = "/home/nathan/.pydec/cache/cube_" + str(h) + ".mtx"

    try:
        import pydec
        M = pydec.io.read_array(filename)
    except:
        import pydec
        M = regular_cube_innerproduct(cmplx,k)
        pydec.io.write_array(filename,M)

    return M



#solve d_k d_k problem for all reasonable k
#from pylab import semilogy,show,xlabel,ylabel,legend,ylim,xlim
#from matplotlib.font_manager import fontManager, FontProperties

cochain_complex = cmplx.cochain_complex()

for i in [1]: #range(len(cochain_complex)-1):
    print "computing mass matrix"

    if isinstance(cmplx,simplicial_complex):
        Mi = whitney_innerproduct_cache(cmplx,i+1)
    else:
        Mi = regular_cube_innerproduct(cmplx,i+1)

    ##print "constructing solver"
    ##ss = discrete_laplacian_solver(cochain_complex,len(cochain_complex)-i-1,innerproduct=Mi)
    ##print ss
    ##
    ##print "solving"
    ##x,res = ss.solve(b=zeros(ss.A.shape[0]),x0=rand(ss.A.shape[0]),return_residuals=True)

    bh = boundary_hierarchy(cochain_complex)
    while len(bh) < 3:
        bh.coarsen()
    print repr(bh)

    N = len(cochain_complex) - 1

    B =  bh[0][N - i].B

    A = B.T.tocsr() * B
    #A = B.T.tocsr() * Mi * B

    constant_prolongators = [lvl[N - i].I for lvl in bh[:-1]]

    if i == 0:
        candidates = None
    else:
        #candidates = [ones(A.shape[0])]

        #TODO test
        candidates = []
        for coord in range(mesh['vertices'].shape[1]):
            candidates.append( bh[0][N-i+1].B * mesh['vertices'][:,coord] )

        K = len(candidates)

        constant_prolongators = [constant_prolongators[0]] + \
                [expand_into_blocks(T,K,1).tocsr() for T in constant_prolongators[1:] ]


    ml = smoothed_aggregation_solver(A,candidates,aggregation=constant_prolongators)
    #ml = smoothed_aggregation_solver(A,candidates)

    x = rand(A.shape[0])
    b = zeros_like(x)
    #b = A*rand(A.shape[0])

    if True:
        x_sol,residuals = ml.solve(b,x0=x,maxiter=50,tol=1e-12,return_residuals=True)
    else:
        residuals = []
        def add_resid(x):
            residuals.append(linalg.norm(b - A*x))
        A.psolve = ml.psolve
        x_sol = linalg.cg(A,b,x0=x,maxiter=30,tol=1e-12,callback=add_resid)[0]


    residuals = array(residuals)/residuals[0]
    avg_convergence_ratio = residuals[-1]**(1.0/len(residuals))
    print "average convergence ratio",avg_convergence_ratio
    print "last convergence ratio",residuals[-1]/residuals[-2]

    print residuals




##candidates = None
##blocks = None
##
##
##
##A = io.mmread('tests/sample_data/elas30_A.mtx').tocsr()
##candidates = io.mmread('tests/sample_data/elas30_nullspace.mtx')
##candidates = [ array(candidates[:,x]) for x in range(candidates.shape[1]) ]
##blocks = arange(A.shape[0]/2).repeat(2)
##
##ml = smoothed_aggregation_solver(A,candidates,blocks=blocks,epsilon=0,max_coarse=10,max_levels=10)
###ml = ruge_stuben_solver(A)
##
##x = rand(A.shape[0])
###b = zeros_like(x)
##b = A*rand(A.shape[0])
##
##if True:
##    x_sol,residuals = ml.solve(b,x0=x,maxiter=30,tol=1e-12,return_residuals=True)
##else:
##    residuals = []
##    def add_resid(x):
##        residuals.append(linalg.norm(b - A*x))
##    A.psolve = ml.psolve
##    x_sol = linalg.cg(A,b,x0=x,maxiter=25,tol=1e-12,callback=add_resid)[0]
##
##
##residuals = array(residuals)/residuals[0]
##avg_convergence_ratio = residuals[-1]**(1.0/len(residuals))
##print "average convergence ratio",avg_convergence_ratio
##print "last convergence ratio",residuals[-1]/residuals[-2]
##
##print residuals
##
