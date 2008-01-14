import numpy,scipy,scipy.sparse
from numpy import sqrt, ravel, diff, zeros, zeros_like, inner, concatenate, \
                  asarray, hstack, ascontiguousarray, isinf, dot
from numpy.random import randn
from scipy.sparse import csr_matrix,coo_matrix

from relaxation import gauss_seidel
from multilevel import multilevel_solver
from sa import sa_constant_interpolation,sa_fit_candidates
from utils import approximate_spectral_radius,hstack_csr,vstack_csr,diag_sparse



def augment_candidates(AggOp, old_Q, old_R, new_candidate):
    #TODO update P and A also

    K = old_R.shape[1]

    #determine blocksizes
    if new_candidate.shape[0] == old_Q.shape[0]:
        #then this is the first prolongator
        old_bs = (1,K)
        new_bs = (1,K+1)
    else:
        old_bs = (K,K)
        new_bs = (K+1,K+1)

    AggOp = expand_into_blocks(AggOp,new_bs[0],1).tocsr() #TODO switch to block matrix


    # tentative prolongator
    #TODO USE BSR
    Q_indptr  = (K+1)*AggOp.indptr
    Q_indices = ((K+1)*AggOp.indices).repeat(K+1)
    for i in range(K+1):
        Q_indices[i::K+1] += i
    Q_data = zeros((AggOp.indptr[-1]/new_bs[0],) + new_bs)
    Q_data[:,:old_bs[0],:old_bs[1]] = old_Q.data.reshape((-1,) + old_bs)   #TODO BSR change

    # coarse candidates
    R = zeros((AggOp.shape[1],K+1,K+1))
    R[:,:K,:K] = old_R.reshape(-1,K,K)

    c = new_candidate.reshape(-1)[diff(AggOp.indptr) == 1]  #eliminate DOFs that aggregation misses
    threshold = 1e-10 * abs(c).max()   # cutoff for small basis functions

    X = csr_matrix((c,AggOp.indices,AggOp.indptr),shape=AggOp.shape)

    #orthogonalize X against previous
    for i in range(K):
        old_c = ascontiguousarray(Q_data[:,:,i].reshape(-1))
        D_AtX = csr_matrix((old_c*X.data,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X
        R[:,i,K] = D_AtX
        X.data -= D_AtX[X.indices] * old_c

    #normalize X
    D_XtX = csr_matrix((X.data**2,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of X.T * X
    col_norms = sqrt(D_XtX)
    mask = col_norms < threshold  # find small basis functions
    col_norms[mask] = 0           # and set them to zero

    R[:,K,K] = col_norms      # store diagonal entry into R

    col_norms = 1.0/col_norms
    col_norms[mask] = 0
    X.data *= col_norms[X.indices]
    Q_data[:,:,-1] = X.data.reshape(-1,new_bs[0])

    Q_data = Q_data.reshape(-1)  #TODO BSR change
    R = R.reshape(-1,K+1)

    Q = csr_matrix((Q_data,Q_indices,Q_indptr),shape=(AggOp.shape[0],(K+1)*AggOp.shape[1]))

    return Q,R






def smoothed_prolongator(P,A):
    #just use Richardson for now
    #omega = 4.0/(3.0*approximate_spectral_radius(A))
    #return P - omega*(A*P)
    #return P  #TEST

    D = diag_sparse(A)
    D_inv_A = diag_sparse(1.0/D)*A
    omega = 4.0/(3.0*approximate_spectral_radius(D_inv_A))
    print "spectral radius",approximate_spectral_radius(D_inv_A) #TODO remove this
    D_inv_A *= omega

    return P - D_inv_A*P



def sa_hierarchy(A,B,AggOps):
    """
    Construct multilevel hierarchy using Smoothed Aggregation
        Inputs:
          A  - matrix
          Ps - list of constant prolongators
          B  - "candidate" basis function to be approximated
        Ouputs:
          (As,Ps,Ts) - tuple of lists
                  - As - [A, Ts[0].T*A*Ts[0], Ts[1].T*A*Ts[1], ... ]
                  - Ps - smoothed prolongators
                  - Ts - tentative prolongators
    """
    As = [A]
    Ps = []
    Ts = []
    Bs = [B]

    for AggOp in AggOps:
        P,B = sa_fit_candidates(AggOp,B)
        I   = smoothed_prolongator(P,A)
        A   = I.T.tocsr() * A * I
        As.append(A)
        Ts.append(P)
        Ps.append(I)
        Bs.append(B)
    return As,Ps,Ts,Bs


class adaptive_sa_solver:
    def __init__(self, A, blocks=None, aggregation=None, max_levels=10, max_coarse=100,\
                 max_candidates=1, mu=5, epsilon=0.1):

        self.A = A

        self.Rs = []

        if self.A.shape[0] <= max_coarse:
            raise ValueError,'small matrices not handled yet'

        #first candidate
        x,AggOps = self.__initialization_stage(A, blocks = blocks, \
                                               max_levels = max_levels, \
                                               max_coarse = max_coarse, \
                                               mu = mu, epsilon = epsilon, \
                                               aggregation = aggregation )

        #create SA using x here
        As,Ps,Ts,Bs = sa_hierarchy(A,x,AggOps)

        for i in range(max_candidates - 1):
            x = self.__develop_new_candidate(As,Ps,Ts,Bs,AggOps,mu=mu)

            #TODO which is faster?
            As,Ps,Ts,Bs = self.__augment_cycle(As,Ps,Ts,Bs,AggOps,x)

            #B = hstack((Bs[0],x))
            #As,Ps,Ts,Bs = sa_hierarchy(A,B,AggOps)

        #improve candidates?
        if True:
            print "improving candidates"
            B = Bs[0]
            for i in range(max_candidates):
                B = B[:,1:]
                As,Ps,Ts,Bs = sa_hierarchy(A,B,AggOps)
                x = self.__develop_new_candidate(As,Ps,Ts,Bs,AggOps,mu=mu)
                B = hstack((B,x))
            As,Ps,Ts,Bs = sa_hierarchy(A,B,AggOps)

        self.Ts = Ts
        self.solver = multilevel_solver(As,Ps)
        self.AggOps = AggOps
        self.Bs = Bs

    def __initialization_stage(self,A,blocks,max_levels,max_coarse,mu,epsilon,aggregation):
        if aggregation is not None:
            max_coarse = 0
            max_levels = len(aggregation) + 1

        # aSA parameters
        # mu      - number of test relaxation iterations
        # epsilon - minimum acceptable relaxation convergence factor

        #step 1
        A_l = A
        x   = randn(A_l.shape[0],1)
        skip_f_to_i = False

        #step 2
        gauss_seidel(A_l,x,zeros_like(x),iterations=mu,sweep='symmetric')

        #step 3
        #TODO test convergence rate here

        As = [A]
        AggOps = []
        Ps = []

        while len(AggOps) + 1 < max_levels and A_l.shape[0] > max_coarse:
            if aggregation is None:
                W_l   = sa_constant_interpolation(A_l,epsilon=0,blocks=blocks) #step 4b
            else:
                W_l   = aggregation[len(AggOps)]
            P_l,x = sa_fit_candidates(W_l,x)                             #step 4c
            I_l   = smoothed_prolongator(P_l,A_l)                          #step 4d
            A_l   = I_l.T.tocsr() * A_l * I_l                              #step 4e

            blocks = None #not needed on subsequent levels

            print 'A_l.shape',A_l.shape
            AggOps.append(W_l)
            Ps.append(I_l)
            As.append(A_l)

            if A_l.shape <= max_coarse:  break

            if not skip_f_to_i:
                print "."
                x_hat = x.copy()                                                   #step 4g
                gauss_seidel(A_l,x,zeros_like(x),iterations=mu,sweep='symmetric')  #step 4h
                x_A_x = dot(x.T,A_l*x)
                if (x_A_x/dot(x_hat.T,A_l*x_hat))**(1.0/mu) < epsilon:             #step 4i
                    print "sufficient convergence, skipping"
                    skip_f_to_i = True
                    if x_A_x == 0:
                        x = x_hat  #need to restore x

        #update fine-level candidate
        for A_l,I in reversed(zip(As[1:],Ps)):
            gauss_seidel(A_l,x,zeros_like(x),iterations=mu,sweep='symmetric')         #TEST
            x = I * x
        gauss_seidel(A,x,zeros_like(x),iterations=mu,sweep='symmetric')         #TEST

        return x,AggOps  #first candidate,aggregation


    def __develop_new_candidate(self,As,Ps,Ts,Bs,AggOps,mu):
        A = As[0]

        x = randn(A.shape[0],1)
        b = zeros_like(x)

        x = multilevel_solver(As,Ps).solve(b, x0=x, tol=1e-10, maxiter=mu)

        #TEST FOR CONVERGENCE HERE

        temp_Ps = []
        temp_As = [A]

        def make_bridge(P,K):
            indptr = P.indptr[:-1].reshape(-1,K-1)
            indptr = hstack((indptr,indptr[:,-1].reshape(-1,1)))
            indptr = indptr.reshape(-1)
            indptr = hstack((indptr,indptr[-1:])) #duplicate last element
            return csr_matrix((P.data,P.indices,indptr),shape=(K*P.shape[0]/(K-1),P.shape[1]))

        for i in range(len(As) - 2):
            T,R = augment_candidates(AggOps[i], Ts[i], Bs[i+1], x)

            P = smoothed_prolongator(T,A)
            A = P.T.tocsr() * A * P

            temp_Ps.append(P)
            temp_As.append(A)

            #TODO USE BSR (K,K) -> (K,K-1)
            bridge = make_bridge(Ps[i+1],R.shape[1])

            solver = multilevel_solver( [A] + As[i+2:], [bridge] + Ps[i+2:] )

            x = R[:,-1].reshape(-1,1)
            x = solver.solve(zeros_like(x), x0=x, tol=1e-8, maxiter=mu)

        for A,P in reversed(zip(temp_As,temp_Ps)):
            x = P * x
            gauss_seidel(A,x,zeros_like(x),iterations=mu,sweep='symmetric')

        return x

    def __augment_cycle(self,As,Ps,Ts,Bs,AggOps,x):
        A = As[0]

        new_As = [A]
        new_Ps = []
        new_Ts = []
        new_Bs = [ hstack((Bs[0],x)) ]

        for i in range(len(As) - 1):
            T,R = augment_candidates(AggOps[i], Ts[i], Bs[i+1], x)

            P = smoothed_prolongator(T,A)
            A = P.T.tocsr() * A * P

            new_As.append(A)
            new_Ps.append(P)
            new_Ts.append(T)
            new_Bs.append(R)

            x = R[:,-1].reshape(-1,1)

        return new_As,new_Ps,new_Ts,new_Bs


if __name__ == '__main__':
    from scipy import *
    from utils import diag_sparse
    from multilevel import poisson_problem1D,poisson_problem2D

    blocks = None
    aggregation = None

    #A = poisson_problem2D(200,1e-2)
    #aggregation = [ sa_constant_interpolation(A*A*A,epsilon=0.0) ]

    #A = io.mmread("tests/sample_data/laplacian_41_3dcube.mtx").tocsr()
    #A = io.mmread("laplacian_40_3dcube.mtx").tocsr()
    #A = io.mmread("/home/nathan/Desktop/9pt/9pt-100x100.mtx").tocsr()
    #A = io.mmread("/home/nathan/Desktop/BasisShift_W_EnergyMin_Luke/9pt-5x5.mtx").tocsr()


    #D = diag_sparse(1.0/sqrt(10**(12*rand(A.shape[0])-6))).tocsr()
    #A = D * A * D

    A = io.mmread("tests/sample_data/elas30_A.mtx").tocsr()
    blocks = arange(A.shape[0]/2).repeat(2)

    from time import clock; start = clock()
    asa = adaptive_sa_solver(A,max_candidates=3,mu=5,blocks=blocks,aggregation=aggregation)
    print "Adaptive Solver Construction: %s seconds" % (clock() - start); del start

    scipy.random.seed(0)  #make tests repeatable
    x = randn(A.shape[0])
    b = A*randn(A.shape[0])
    #b = zeros(A.shape[0])


    print "solving"
    if False:
        x_sol,residuals = asa.solver.solve(b,x0=x,maxiter=20,tol=1e-12,return_residuals=True)
    else:
        residuals = []
        def add_resid(x):
            residuals.append(linalg.norm(b - A*x))
        A.psolve = asa.solver.psolve
        x_sol = linalg.cg(A,b,x0=x,maxiter=30,tol=1e-12,callback=add_resid)[0]

    residuals = array(residuals)/residuals[0]

    print "residuals ",residuals
    print "mean convergence factor",(residuals[-1]/residuals[0])**(1.0/len(residuals))
    print "last convergence factor",residuals[-1]/residuals[-2]

    print
    print asa.solver

    print "constant Rayleigh quotient",dot(ones(A.shape[0]),A*ones(A.shape[0]))/float(A.shape[0])

    def plot2d_arrows(x):
        from pylab import figure,quiver,show
        x = x.reshape(-1)
        N = (len(x)/2)**0.5
        assert(2 * N * N == len(x))
        X = linspace(-1,1,N).reshape(1,N).repeat(N,0).reshape(-1)
        Y = linspace(-1,1,N).reshape(1,N).repeat(N,0).T.reshape(-1)

        dX = x[0::2]
        dY = x[1::2]

        figure()
        quiver(X,Y,dX,dY)
        show()

    def plot2d(x):
        from pylab import pcolor,figure,show
        figure()
        pcolor(x.reshape(sqrt(len(x)),sqrt(len(x))))
        show()


    for c in asa.Bs[0].T:
        #plot2d(c)
        plot2d_arrows(c)
        print "candidate Rayleigh quotient",dot(c,A*c)/dot(c,c)


    ##W = asa.AggOps[0]*asa.AggOps[1]
    ##pcolor((W * rand(W.shape[1])).reshape((200,200)))
