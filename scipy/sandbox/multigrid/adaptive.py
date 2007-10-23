import numpy,scipy,scipy.sparse
from numpy import sqrt, ravel, diff, zeros, zeros_like, inner, concatenate, \
                  asarray, hstack
from scipy.sparse import csr_matrix,coo_matrix

from relaxation import gauss_seidel
from multilevel import multilevel_solver
from sa import sa_constant_interpolation,sa_fit_candidates
from utils import approximate_spectral_radius,hstack_csr,vstack_csr,expand_into_blocks



def augment_candidates(AggOp, old_Q, old_R, new_candidate, level):
    K = old_R.shape[1]
    
    #determine blocksizes
    if level == 0:
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
    Q_data = zeros((AggOp.shape[0]/new_bs[0],) + new_bs)
    Q_data[:,:new_bs[0],:new_bs[1]] = old_Q.data.reshape((-1,) + old_bs)   #TODO BSR change 

    # coarse candidates
    R = zeros(((K+1)*AggOp.shape[1],K+1)) 
    for i in range(K):
        R[i::K+1,:K] = old_R[i::K,:] 

    c = new_candidate.reshape(-1)[diff(AggOp.indptr) == 1]  #eliminate DOFs that aggregation misses
    X = csr_matrix((c,AggOp.indices,AggOp.indptr),dims=AggOp.shape)

    #orthogonalize X against previous
    for i in range(K):
        old_c = ascontiguousarray(Q_data[:,:,i].reshape(-1))
        D_AtX = csr_matrix((old_c*X.data,X.indices,X.indptr),dims=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X            
        R[i::K+1,-1] = D_AtX
        X.data -= D_AtX[X.indices] * old_c
     
    #normalize X
    D_XtX = csr_matrix((X.data**2,X.indices,X.indptr),dims=X.shape).sum(axis=0).A.flatten() #same as diagonal of X.T * X            
    col_norms = sqrt(D_XtX)
    R[K::K+1,-1] = col_norms

    col_norms = 1.0/col_norms
    col_norms[isinf(col_norms)] = 0
    X.data *= col_norms[X.indices]
    last_column = Q_data[:,:,-1].reshape(-1)
    last_column = X.data.reshape(-1)
    Q_data = Q_data.reshape(-1)  #TODO BSR change

    Q = csr_matrix((Q_data,Q_indices,Q_indptr),dims=(AggOp.shape[0],(K+1)*AggOp.shape[1]))

    return Q,R
    



def orthonormalize_prolongator(P_l,x_l,W_l,W_m):
    """
     
    """

    #candidate prolongator (assumes every value from x is used)  #TODO permit gaps
    X = csr_matrix((x_l,W_l.indices,W_l.indptr),dims=W_l.shape,check=False)  
    
    R = (P_l.T.tocsr() * X)  # R has at most 1 nz per row
    X = X - P_l*R            # othogonalize X against P_l
    
    #TODO DROP REDUNDANT COLUMNS FROM P (AND R?) HERE (NULL OUT R ACCORDINGLY?)    
    #TODO REMOVE CORRESPONDING COLUMNS FROM W_l AND ROWS FROM A_m ALSO
    W_l_new = W_l
    W_m_new = W_m

    #normalize surviving columns of X
    col_norms = ravel(sqrt(csr_matrix((X.data*X.data,X.indices,X.indptr),dims=X.shape,check=False).sum(axis=0)))
    print "zero cols",sum(col_norms == 0)
    print "small cols",sum(col_norms < 1e-8)
    Xcopy = X.copy()
    X.data *= (1.0/col_norms)[X.indices]

    P_l_new = hstack_csr(P_l,X)


    #check orthonormality
    print "norm(P.T*P - I) ",scipy.linalg.norm((P_l_new.T * P_l_new - scipy.sparse.spidentity(P_l_new.shape[1])).data)
    #assert(scipy.linalg.norm((P_l_new.T * P_l_new - scipy.sparse.spidentity(P_l_new.shape[1])).data)<1e-8)
    
    x_m = zeros(P_l_new.shape[1],dtype=x_l.dtype)
    x_m[:P_l.shape[1]][diff(R.indptr).astype('bool')] = R.data
    x_m[P_l.shape[1]:] = col_norms

    print "||x_l - P_l*x_m||",scipy.linalg.norm(P_l_new* x_m - x_l) #see if x_l is represented exactly

    return P_l_new,x_m,W_l,W_m


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
 


def sa_hierarchy(A,Ws,B):
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

    for W in Ws:
        P,B = sa_fit_candidates(W,B)
        I   = smoothed_prolongator(P,A)  
        A   = I.T.tocsr() * A * I
        As.append(A)
        Ts.append(P)
        Ps.append(I)
        Bs.append(B)
    return As,Ps,Ts,Bs
         
def make_bridge(I,N):
    tail = I.indptr[-1].repeat(N - I.shape[0])
    ptr = concatenate((I.indptr,tail))
    return csr_matrix((I.data,I.indices,ptr),dims=(N,I.shape[1]),check=False)

class adaptive_sa_solver:
    def __init__(self,A,blocks=None,options=None,max_levels=10,max_coarse=100,\
                        max_candidates=1,mu=5,epsilon=0.1):

        self.A = A

        self.Rs = [] 
        
        if self.A.shape[0] <= max_coarse:
            raise ValueError,'small matrices not handled yet'
       
        #first candidate 
        x,AggOps = self.__initialization_stage(A, blocks = blocks, \
                                               max_levels = max_levels, \
                                               max_coarse = max_coarse, \
                                               mu = mu, epsilon = epsilon) 
        
        Ws = AggOps

        fine_candidates = x

        #create SA using x here
        As,Ps,Ts,self.candidates = sa_hierarchy(A,AggOps,fine_candidates)
        #self.candidates = [x]

        for i in range(max_candidates - 1):
            #x = self.__develop_candidate_OLD(A,As,Ps,Ts,Ws,AggOps,mu=mu)    
            #As,Ps,Ts,Ws = self.__augment_cycle(A,As,Ts,Ws,AggOps,x)  
            #self.candidates.append(x)
            
            #TODO which is faster?
            x = self.__develop_candidate(As,Ps,Ts,AggOps,self.candidates,mu=mu)    
            fine_candidates = hstack((fine_candidates,x))
            As,Ps,Ts,self.candidates = sa_hierarchy(A,AggOps,fine_candidates)

        self.Ts = Ts 
        self.solver = multilevel_solver(As,Ps)
        self.AggOps = AggOps
                

    def __initialization_stage(self,A,blocks,max_levels,max_coarse,mu,epsilon):
        AggOps = []
        Ps     = []

        # aSA parameters 
        # mu      - number of test relaxation iterations
        # epsilon - minimum acceptable relaxation convergence factor 
        
        #step 1
        A_l = A
        x   = scipy.rand(A_l.shape[0],1)
        skip_f_to_i = False
        
        #step 2
        gauss_seidel(A_l,x,zeros_like(x),iterations=mu,sweep='symmetric')

        #step 3
        #TODO test convergence rate here 
      
        As = [A]

        while len(AggOps) + 1 < max_levels  and A_l.shape[0] > max_coarse:
            W_l   = sa_constant_interpolation(A_l,epsilon=0,blocks=blocks) #step 4b
            P_l,x = sa_fit_candidates(W_l,x)                             #step 4c  
            I_l   = smoothed_prolongator(P_l,A_l)                          #step 4d
            A_l   = I_l.T.tocsr() * A_l * I_l                              #step 4e
           
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


    def __develop_candidate(self,As,Ps,Ts,AggOps,candidates,mu):
        A = As[0]
        
        x = A*scipy.rand(A.shape[0],1)
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
            return csr_matrix((P.data,P.indices,indptr),dims=(K*P.shape[0]/(K-1),P.shape[1]))

        for i in range(len(As) - 2):            
            #TODO test augment_candidates against fit candidates
            if i == 0:
                temp = hstack((candidates[0],x))
            else:
                K = candidates[i].shape[1] 
                temp = zeros((x.shape[0]/(K+1),K+1,K + 1))
                temp[:,:-1,:-1] = candidates[i].reshape(-1,K,K)
                temp[:,:,-1] = x.reshape(-1,K+1,1)
                temp = temp.reshape(-1,K+1)
            T_,R_ = sa_fit_candidates(AggOps[i],temp)
            #print "T - T_",(T - T_).data.max()
            #assert((T - T_).data.max() < 1e-10)
            #assert((R - R_).data.max() < 1e-10)
            T,R = T_,R_
            #TODO end test

            #T,R = augment_candidates(AggOps[i], Ts[i], candidates[i+1], x, i)
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
        
        #for P in reversed(temp_Ps):
        #    x = P*x
        
        return x


##    def __develop_candidate_OLD(self,A,As,Ps,Ts,Ws,AggOps,mu):
##        #scipy.random.seed(0)  #TEST
##        x = scipy.rand(A.shape[0])
##        b = zeros_like(x)
##
##        solver = multilevel_solver(As,Ps)
##
##        x = solver.solve(b, x0=x, tol=1e-10, maxiter=mu)
##    
##        #TEST FOR CONVERGENCE HERE
## 
##        A_l,P_l,W_l,x_l = As[0],Ts[0],Ws[0],x
##        
##        temp_Ps = []
##        for i in range(len(As) - 2):
##            P_l_new, x_m, W_l_new, W_m_new = orthonormalize_prolongator(P_l, x_l, W_l, AggOps[i+1])    
## 
##            I_l_new = smoothed_prolongator(P_l_new,A_l)
##            A_m_new = I_l_new.T.tocsr() * A_l * I_l_new
##            bridge = make_bridge(Ps[i+1],A_m_new.shape[0]) 
## 
##            temp_solver = multilevel_solver( [A_m_new] + As[i+2:], [bridge] + Ps[i+2:] )
## 
##            for n in range(mu):
##                x_m = temp_solver.solve(zeros_like(x_m), x0=x_m, tol=1e-8, maxiter=1)
## 
##            temp_Ps.append(I_l_new)
## 
##            W_l = vstack_csr(Ws[i+1],W_m_new)  #prepare for next iteration
##            A_l = A_m_new
##            x_l = x_m
##            P_l = make_bridge(Ts[i+1],A_m_new.shape[0])
## 
##        x = x_l
##        for I in reversed(temp_Ps):
##            x = I*x
##        
##        return x
           

    def __augment_cycle(self,A,As,Ts,Ws,AggOps,x):
        #make a new cycle using the new candidate
        A_l,P_l,W_l,x_l = As[0],Ts[0],AggOps[0],x            
        
        new_As,new_Ps,new_Ts,new_Ws = [A],[],[],[AggOps[0]]

        for i in range(len(As) - 2):
            P_l_new, x_m, W_l_new, W_m_new = orthonormalize_prolongator(P_l, x_l, W_l, AggOps[i+1])

            I_l_new = smoothed_prolongator(P_l_new,A_l)
            A_m_new = I_l_new.T.tocsr() * A_l * I_l_new
            W_m_new = vstack_csr(Ws[i+1],W_m_new)

            new_As.append(A_m_new)
            new_Ws.append(W_m_new)
            new_Ps.append(I_l_new)
            new_Ts.append(P_l_new)

            #prepare for next iteration
            W_l = W_m_new 
            A_l = A_m_new
            x_l = x_m
            P_l = make_bridge(Ts[i+1],A_m_new.shape[0]) 
        
        P_l_new, x_m, W_l_new, W_m_new = orthonormalize_prolongator(P_l, x_l, W_l, csr_matrix((P_l.shape[1],1)))
        I_l_new = smoothed_prolongator(P_l_new,A_l)
        A_m_new = I_l_new.T.tocsr() * A_l * I_l_new

        new_As.append(A_m_new)
        new_Ps.append(I_l_new)
        new_Ts.append(P_l_new)

        return new_As,new_Ps,new_Ts,new_Ws



from scipy import *
from utils import diag_sparse
from multilevel import poisson_problem1D,poisson_problem2D

blocks = None

A = poisson_problem2D(100)
#A = io.mmread("tests/sample_data/laplacian_41_3dcube.mtx").tocsr()
#A = io.mmread("laplacian_40_3dcube.mtx").tocsr()
#A = io.mmread("/home/nathan/Desktop/9pt/9pt-100x100.mtx").tocsr()
#A = io.mmread("/home/nathan/Desktop/BasisShift_W_EnergyMin_Luke/9pt-5x5.mtx").tocsr()


#D = diag_sparse(1.0/sqrt(10**(12*rand(A.shape[0])-6))).tocsr()
#A = D * A * D

A = io.mmread("tests/sample_data/elas30_A.mtx").tocsr()
blocks = arange(A.shape[0]/2).repeat(2)

asa = adaptive_sa_solver(A,max_candidates=4,mu=5)
scipy.random.seed(0)  #make tests repeatable
x = rand(A.shape[0])
#b = A*rand(A.shape[0])
b = zeros(A.shape[0])


print "solving"
if True:
    x_sol,residuals = asa.solver.solve(b,x0=x,maxiter=25,tol=1e-7,return_residuals=True)
else:
    residuals = []
    def add_resid(x):
        residuals.append(linalg.norm(b - A*x))
    A.psolve = asa.solver.psolve
    x_sol = linalg.cg(A,b,x0=x,maxiter=20,tol=1e-12,callback=add_resid)[0]

residuals = array(residuals)/residuals[0]

print "residuals ",residuals
print "mean convergence factor",(residuals[-1]/residuals[0])**(1.0/len(residuals))
print "last convergence factor",residuals[-1]/residuals[-2]

print
print asa.solver

print "constant Rayleigh quotient",dot(ones(A.shape[0]),A*ones(A.shape[0]))/float(A.shape[0])

for c in asa.candidates[0].T:
    print "candidate Rayleigh quotient",dot(c,A*c)/dot(c,c)



##W = asa.AggOps[0]*asa.AggOps[1]
##pcolor((W * rand(W.shape[1])).reshape((200,200)))

def plot2d_arrows(x):
    from pylab import quiver
    x = x.reshape(-1)
    N = (len(x)/2)**0.5
    assert(2 * N * N == len(x))
    X = linspace(-1,1,N).reshape(1,N).repeat(N,0).reshape(-1)
    Y = linspace(-1,1,N).reshape(1,N).repeat(N,0).T.reshape(-1)

    dX = x[0::2]
    dY = x[1::2]

    quiver(X,Y,dX,dY)

def plot2d(x):
    from pylab import pcolor
    pcolor(x.reshape(sqrt(len(x)),sqrt(len(x))))
    
