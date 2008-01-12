
from scipy import *
from scipy.sparse import *
from pydec import *
from pydec.multigrid.discrete_laplacian import boundary_hierarchy, discrete_laplacian_solver, hodge_solver

from scipy.sandbox.multigrid import smoothed_aggregation_solver,multigridtools,multilevel_solver
from scipy.sandbox.multigrid.adaptive import adaptive_sa_solver
from scipy.sandbox.multigrid.sa import sa_smoothed_prolongator
from scipy.sandbox.multigrid.utils import expand_into_blocks


## Load mesh from file
mesh_path = '../../../../../hirani_group/wnbell/meshes/'
#mesh = read_mesh(mesh_path + 'rocket/rocket.xml')
#mesh = read_mesh(mesh_path + 'genus3/genus3_168k.xml')
#mesh = read_mesh(mesh_path + 'genus3/genus3_455k.xml')
#mesh = read_mesh(mesh_path + '/torus/torus.xml')
mesh = read_mesh(mesh_path + '/sq14tri/sq14tri.xml')
for i in range(5):
    mesh['vertices'],mesh['elements'] = loop_subdivision(mesh['vertices'],mesh['elements'])
cmplx = simplicial_complex(mesh['vertices'],mesh['elements'])

## Construct mesh manually
#bitmap = ones((60,60),dtype='bool')
#bitmap[1::10,1::10] = False
#bitmap[100:150,100:400] = False
#cmplx = regular_cube_complex(regular_cube_mesh(bitmap))

def curl_curl_prolongator(D_nodal,vertices):
    if not isspmatrix_csr(D_nodal):
        raise TypeError('expected csr_matrix')

    A = D_nodal.T.tocsr() * D_nodal
    aggs = multigridtools.sa_get_aggregates(A.shape[0],A.indptr,A.indices)

    num_edges = D_nodal.shape[0]
    num_basis = vertices.shape[1]
    num_aggs  = aggs.max() + 1

    # replace with CSR + eliminate duplicates
    #indptr  = (2*num_basis) * arange(num_edges+1)
    ## same same
    #csr_matrix((data,indices,indptr),shape=(num_edges,num_aggs))

    row  = arange(num_edges).repeat(2*num_basis)
    col  = (num_basis*aggs[D_nodal.indices]).repeat(num_basis)
    col = col.reshape(-1,num_basis) + arange(num_basis)
    col = col.reshape(-1)
    data = tile(0.5 * (D_nodal*vertices),(1,2)).reshape(-1)

    return coo_matrix((data,(row,col)),shape=(num_edges,num_basis*num_aggs)).tocsr()





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


    dimension = mesh['vertices'].shape[1]

    if True:

        d0 = cmplx[0].d
        d1 = cmplx[1].d

        #A = (d1.T.tocsr() * d1 + d0 * d0.T.tocsr()).astype('d')
        A = (d1.T.tocsr()  * d1).astype('d')

        P = curl_curl_prolongator(d0,mesh['vertices'])

        num_blocks = P.shape[1]/dimension
        blocks = arange(num_blocks).repeat(dimension)

        P = sa_smoothed_prolongator(A,P,epsilon=0,omega=4.0/3.0)

        PAP = P.T.tocsr() * A * P

        candidates = None
        candidates = zeros((num_blocks,dimension,dimension))
        for n in range(dimension):
            candidates[:,n,n] = 1.0
        candidates = candidates.reshape(-1,dimension)

        ml = smoothed_aggregation_solver(PAP,epsilon=0.0,candidates=candidates,blocks=blocks)
        #A = PAP
        ml = multilevel_solver([A] + ml.As, [P] + ml.Ps)
    else:

        bh = boundary_hierarchy(cochain_complex)
        while len(bh) < 3:
            bh.coarsen()
        print repr(bh)

        N = len(cochain_complex) - 1

        B =  bh[0][N - i].B

        A = (B.T.tocsr() * B).astype('d')
        #A = B.T.tocsr() * Mi * B

        constant_prolongators = [lvl[N - i].I for lvl in bh[:-1]]

        method = 'aSA'

        if method == 'RS':
            As = [A]
            Ps = []
            for T in constant_prolongators:
                Ps.append( sa_smoothed_prolongator(As[-1],T,epsilon=0.0,omega=4.0/3.0) )
                As.append(Ps[-1].T.tocsr() * As[-1] * Ps[-1])
            ml = multilevel_solver(As,Ps)

        else:
            if method == 'BSA':
                if i == 0:
                    candidates = None
                else:
                    candidates = cmplx[0].d * mesh['vertices']
                    K = candidates.shape[1]

                    constant_prolongators = [constant_prolongators[0]] + \
                            [expand_into_blocks(T,K,1).tocsr() for T in constant_prolongators[1:] ]

                    ml = smoothed_aggregation_solver(A,candidates,aggregation=constant_prolongators)
            elif method == 'aSA':
                asa = adaptive_sa_solver(A,aggregation=constant_prolongators,max_candidates=dimension,epsilon=0.0)
                ml = asa.solver
            else:
                raise ValuerError,'unknown method'

        #ml = smoothed_aggregation_solver(A,candidates)

    #x = d0 * mesh['vertices'][:,0]
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

        from pydec import cg
        x_sol = cg(A,b,x0=x,maxiter=40,tol=1e-8,callback=add_resid)[0]


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
