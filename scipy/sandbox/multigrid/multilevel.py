from scipy import ones,zeros,rand,array,array_split,hstack,transpose,sum,ones_like,sqrt,concatenate
from scipy.sparse import spidentity,csr_matrix,coo_matrix
from numpy.linalg import norm
from numpy import zeros_like,arange,inner,diff,ravel

import pydec
from pydec import diag_sparse,inf_norm, mls_polynomial_coeffs,polynomial_smoother
     
from multigrid import sa_interpolation,rs_interpolation,sa_constant_interpolation,sa_no_threshold
import multigrid
import multigridtools
from relaxation import gauss_seidel,jacobi 

import scipy
import numpy


##    import scipy.sandbox.arpack as arpack
##    eigs,vecs = arpack.eigen(A,maxiter=10)
##    raise ValueError
##    return eigs.max()

    
def avg_work_per_digit(ml_solver,residuals):
    digits = numpy.log(residuals[0]/residuals[-1])/numpy.log(10)
    return (ml_solver.cycle_complexity() * len(residuals)) / digits


def avg_convergence_rate(residuals):
    return (residuals[-1]/residuals[0]) ** (1.0/len(residuals))


def asym_work_per_digit(ml_solver,residuals):
    digits = numpy.log(residuals[-2]/residuals[-1])/numpy.log(10)
    return (ml_solver.cycle_complexity()) / digits




class coarse_grid_solver:
    def __init__(self,A,options):
        self.opts = options
        
        self.A = A

        solver = self.opts['coarse: type'] 
        
        if solver == 'pinv':
            self.pinv = scipy.linalg.pinv(self.A.todense())
            self.nnz  = self.pinv.size
            self.__solve = lambda b : numpy.dot(self.pinv,b)
        elif solver == 'pinv2':
            self.pinv = scipy.linalg.pinv2(self.A.todense())
            self.nnz  = self.pinv.size
            self.__solve = lambda b : numpy.dot(self.pinv,b)
        elif solver == 'splu':
            import scipy.linsolve.umfpack as um
            self.umfpack = um.UmfpackContext()
            self.umfpack.numeric( self.A )
            self.nnz  = self.umfpack.info[um.umfDefines['UMFPACK_LU_ENTRIES']]
            self.__solve = lambda b : self.umfpack.solve( um.UMFPACK_A, self.A, b, autoTranspose = True )
        elif solver in ['bicg','bicgstab','cg','cgs','gmres','qmr']:
            #self.__solve = lambda b : scipy.linalg.cg(self.A,b,tol=1e-12,maxiter=100)[0]
            #it_solver = getattr(scipy.linalg.iterative,solver)
            
            it_solver = pydec.numerical.iterative.cg
            self.__solve = lambda b : it_solver(self.A,b,tol=1e-12)[0]
        else:
            raise ValueError,('unknown solver: %s' % solver)
            
    def solve(self,b):
        #M = self.A.todense()
        #val,vec = scipy.linalg.eig(M)
        #pet = vec[:,val < 1e-8][:,0]
        #print pet
        #return self.__solve(b) + pet
        return self.__solve(b)

    def nnz(self):
        return self.nnz


class multilevel_solver(list):
    class grid_data:
        pass

    class options(dict):        
        def __repr__(self):
            keys  = sorted([k for k in self.keys() if ':' not in k])
            keys += sorted([k for k in self.keys() if ':'     in k])
            
            output = "solver options:\n"
            for k in keys:
               output += "   %-25s %-30s\n" % (k,self[k])
            return output
        def sub_options(self,sub_opt):
            """
            Filter options with a given prefix

            Example:
               opts.sub_options('smoother:')
               
            """
            return dict([ (k,v) for (k,v) in self.iteritems() if k.startswith(sub_opt)])
            
    def __init__(self,A,options=None):
        assert(False) #should not instantiated

            
    def __repr__(self):
        output = '%s\n'% type(self).__name__
        output += 'Number of Levels:     %d     (max: %d)\n' % (len(self),self.opts['max levels'])
        output += 'Operator Complexity: %6.3f\n' % self.operator_complexity()
        output += 'Grid Complexity:     %6.3f\n' % self.grid_complexity()
        output += 'Cycle Complexity:    %6.3f\n' % self.cycle_complexity()

        total_nnz =  sum([lvl.A.nnz for lvl in self])
        
        for lvl,data in enumerate(self):
            output += '   [level %2d]  unknowns: %10d  nnz: %5.2f%%\n' % (lvl,data.A.shape[1],(100*float(data.A.nnz)/float(total_nnz)))

        #output += '\n' + repr(self.opts)
        return output
    


    def operator_complexity(self):
        """number of nonzeros on all levels / number of nonzeros on the finest level"""
        return sum([lvl.A.nnz for lvl in self])/float(self[0].A.nnz)
    def grid_complexity(self):
        """number of unknowns on all levels / number of unknowns on the finest level"""
        return sum([lvl.A.shape[0] for lvl in self])/float(self[0].A.shape[0])
    def cycle_complexity(self):
        """total FLOPs in one MG cycle / FLOPs in single smoother sweep on the finest level"""
        return self.cycle_flops()/float(self[0].A.nnz)
    def cycle_flops(self):
        """total FLOPs in one MG cycle"""
        total_flops = 0

        gamma  = self.opts['cycle: gamma']
        passes = self.opts['smoother: passes']

        if self.opts['smoother: type'] in ['jacobi','symmetric gauss-seidel','richardson']:
            passes *= 2
            passes += 1 #residual computation
            
        if self.opts['smoother: type'] in ['polynomial']:
            print "poly degree:",len(self.opts['smoother: omega'][-1])
            passes *= 2*len(self.opts['smoother: omega'][-1])
            #residual computation already factored in
        

        for n,lvl in enumerate(self):
            total_flops += (gamma**n)*lvl.A.nnz*passes

        #account for iterative solver using this as a preconditioner
        if self.opts['solver: type'] != 'standalone':
            total_flops += self.A.nnz

        return total_flops
    
    def solve(self,b, x0=None, tol=1e-5, maxiter=100, callback=None, return_residuals=False, precond=False):

        if x0 is None:
            x = zeros(b.shape,max(self.A.dtype,b.dtype))
        else:
            x = x0.copy()


        #was invoked as a preconditioner
        if precond:
	    #return b #no precond
            self.__solve(0,x,b)
            return x


        if self.opts['solver: type'] == 'standalone':
            residuals = [norm(b-self[0].A*x,2)]

            while len(residuals) <= maxiter and residuals[-1]/residuals[0] > tol:
                self.__solve(0,x,b)

                residuals.append(scipy.linalg.norm(b-self[0].A*x,2))

                if callback is not None:
                    callback(x)

        else:
            #using acceleration
            
            #residuals = [scipy.linalg.norm(b-self[0].A*x,2)]            
            #callback = lambda x_k : residuals.append(scipy.linalg.norm(b-self[0].A*x_k,2))
            #solver = getattr(scipy.linalg.iterative,self.opts['solver: type'])

            assert(self.opts['solver: type'] == 'cg') #only CG supported now
            solver = pydec.iterative.cg

            mtx = self[0].A
            mtx.psolve = lambda b : self.solve(b,precond=True)
            
            x,residuals = solver(mtx,b,x0=x,tol=tol,maxiter=maxiter,callback=callback)

        if return_residuals:
            return x,residuals
        else:
            return x




    def __smooth(self,lvl,x,b,which):
        smoother_type   = self.opts['smoother: type']
        smoother_passes = self.opts['smoother: passes']

        A = self[lvl].A
        
        if smoother_type == 'jacobi':
            omega = self.opts['smoother: omega'][lvl]
            jacobi(A,x,b,iterations=smoother_passes,omega=omega)
        elif smoother_type == 'richardson':
            omega = self.opts['smoother: omega'][lvl]            
            x += omega*(b - A*x)
        elif smoother_type == 'polynomial':
            coeffs = self.opts['smoother: omega'][lvl]
            polynomial_smoother(A,x,b,coeffs)
        elif smoother_type == 'symmetric gauss-seidel':
            if which == 'pre':
                gauss_seidel(A,x,b,iterations=smoother_passes,sweep="forward")
            else:
                gauss_seidel(A,x,b,iterations=smoother_passes,sweep="backward")
        else:
            raise ValueError,'unknown smoother'
        
    def __solve(self,lvl,x,b):

        if len(self) == 1:
            x[:] = self[0].coarse_solver.solve(b)
            return
        
        A = self[lvl].A

        self.__smooth(lvl,x,b,which='pre')
            
        residual = b - A*x

        coarse_x = zeros((self[lvl+1].A.shape[0]))
        coarse_b = self[lvl].P.T * residual
        
        if lvl == len(self) - 2:
            coarse_x[:] = self[-1].coarse_solver.solve(coarse_b)
        else:                
            for i in range(self.opts['cycle: gamma']):
                self.__solve(lvl+1,coarse_x,coarse_b)
                
        x += self[lvl].P * coarse_x

        self.__smooth(lvl,x,b,which='post')
        
        




class scalar_solver(multilevel_solver):
    def __init__(self,A,options=None):
        self.A = A
        
        if options is None:
            self.opts = scalar_solver.default_options()
        else:
            self.opts = options

        self.__construct_hierarchy()

    def default_options():
        opts = multilevel_solver.options()
        opts['max levels']            = 8
        opts['cycle: gamma']          = 1
        opts['coarse: type']          = 'splu'
        opts['coarse: max size']      = 2000
        opts['aggregation: type']     = 'SA'
        opts['aggregation: epsilon']  = 0.05
        opts['smoother: passes']      = 1
        opts['smoother: type']        = 'symmetric gauss-seidel'
#        opts['smoother: type']        = 'jacobi'
        opts['solver: type']          = 'cg'        
        return opts
    default_options = staticmethod(default_options)

    def __construct_hierarchy(self):
        A = self.A
        
        agg_type   = self.opts['aggregation: type']
        max_levels = self.opts['max levels']
        max_coarse = self.opts['coarse: max size']

        while len(self) < max_levels  and A.shape[0] > max_coarse:
            self.append(self.grid_data())
            
            if agg_type == 'SA':
                P,I = sa_interpolation(A)
            elif agg_type == 'RS':
                P = rs_interpolation(A)
            else:
                raise ValueError,'unknown aggregation type: %s' % agg_type

            self[-1].A = A
            self[-1].P = P

            A = (P.T.tocsr() * A) * P

        self.append(self.grid_data())
        
        self[-1].coarse_solver = coarse_grid_solver(A,self.opts.sub_options('coarse:'))
        self[-1].A = A
        if self.opts['smoother: type'] == 'jacobi':
            omegas = []
            for lvl in self:
                A = lvl.A
                D_inv = diag_sparse(1.0/diag_sparse(A))
    
                D_inv_A  = D_inv * A
                omegas.append((4.0/3.0)/inf_norm(D_inv_A))
            self.opts['smoother: omega'] = omegas
        

class multilevel_solver2:
    def __init__(self,As,Ps,options=None):
        self.As = As
        self.Ps = Ps
        self.ops = options
            
    def solve(self,b, x0=None, tol=1e-5, maxiter=100, callback=None, return_residuals=False):

        if x0 is None:
            x = zeros(b.shape,max(self.A.dtype,b.dtype))
        else:
            x = array(x0)

        self.__solve(0,x,b)

        return x
        
    def __solve(self,lvl,x,b):

        A = self.As[lvl]
        
        if len(self.As) == 1:
            x[:] = scipy.linalg.solve(A.todense(),b)
            return x


        self.__smooth(lvl,x,b,which='pre')
            
        residual = b - A*x

        coarse_x = zeros((self.As[lvl+1].shape[0]))
        coarse_b = self.Ps[lvl].T * residual
        
        if lvl == len(self.As) - 2:
            pass
            coarse_x[:] = scipy.linalg.solve(self.As[-1].todense(),coarse_b)
            #coarse_x[:] = self[-1].coarse_solver.solve(coarse_b)   #next level is coarsest 
        else:   
            self.__solve(lvl+1,coarse_x,coarse_b)
                
        x += self.Ps[lvl] * coarse_x

        self.__smooth(lvl,x,b,which='post')


    def __smooth(self,lvl,x,b,which):
        A = self.As[lvl]
        if which == 'pre':
            gauss_seidel(A,x,b,iterations=1,sweep="forward")
        else:
            gauss_seidel(A,x,b,iterations=1,sweep="backward")


def inf_norm(A):
    return abs(A).sum(axis=1).max()  #max abs row sum

def fit_candidate(I,x):
    """
    For each aggregate in I (i.e. each column of I) compute vector R and 
    sparse matrix Q (having the sparsity of I) such that the following holds:

        Q*R = x     and   Q^T*Q = I

    In otherwords, find a prolongator Q with orthonormal columns so that
    x is represented exactly on the coarser level by R.
    """
    Q = csr_matrix((x.copy(),I.indices,I.indptr),dims=I.shape,check=False)
    R = sqrt(numpy.ravel(csr_matrix((x*x,I.indices,I.indptr),dims=I.shape,check=False).sum(axis=0)))  #column 2-norms  
    Q.data *= (1.0/R)[Q.indices]
    print "norm(Q*R - x)",linalg.norm(Q*R - x)
    return Q,R


def scaled_columns_csr(A,scales):
    scales = numpy.ravel(scales)
    A = A.copy()
    A.data *= scales[A.indices]
    return A

def orthonormalize_candidate(I,x,basis):
    Px = csr_matrix((x,I.indices,I.indptr),dims=I.shape,check=False) 
    Rs = []
    #othogonalize columns of Px against other candidates 
    for b in basis:
        Pb = csr_matrix((b,I.indices,I.indptr),dims=I.shape,check=False)
        R  = ravel(csr_matrix((Pb.data*Px.data,I.indices,I.indptr),dims=I.shape,check=False).sum(axis=0)) # columnwise projection of Px on Pb
        Px.data -= R[I.indices] * Pb.data   #subtract component in b direction
        Rs.append(R)

    #filter columns here, set unused cols to 0, add to mask
    
    #normalize columns of Px
    R = ravel(csr_matrix((x**x,I.indices,I.indptr),dims=I.shape,check=False).sum(axis=0))
    Px.data *= (1.0/R)[I.indices]
    Rs.append(R.reshape(-1,1))
    return Rs

def hstack_csr(A,B):
    #OPTIMIZE THIS
    assert(A.shape[0] == B.shape[0])
    A = A.tocoo()
    B = B.tocoo()
    I = concatenate((A.row,B.row))
    J = concatenate((A.col,B.col+A.shape[1]))
    V = concatenate((A.data,B.data))
    return coo_matrix((V,(I,J)),dims=(A.shape[0],A.shape[1]+B.shape[1])).tocsr()


def vstack_csr(A,B):
    #OPTIMIZE THIS
    assert(A.shape[1] == B.shape[1])
    A = A.tocoo()
    B = B.tocoo()
    I = concatenate((A.row,B.row+A.shape[0]))
    J = concatenate((A.col,B.col))
    V = concatenate((A.data,B.data))
    return coo_matrix((V,(I,J)),dims=(A.shape[0]+B.shape[0],A.shape[1])).tocsr()
    


def orthonormalize_prolongator(P_l,x_l,W_l,W_m):
    """
     
    """
    X = csr_matrix((x_l,W_l.indices,W_l.indptr),dims=W_l.shape,check=False)  #candidate prolongator (assumes every value from x is used)
    
    R = (P_l.T.tocsr() * X)  # R has at most 1 nz per row
    X = X - P_l*R            # othogonalize X against P_l
    
    #DROP REDUNDANT COLUMNS FROM P (AND R?) HERE (NULL OUT R ACCORDINGLY?)    
    #REMOVE CORRESPONDING COLUMNS FROM W_l AND ROWS FROM A_m ALSO
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



def prolongation_smoother(A):
    omega = (4.0/3.0)/inf_norm(A)
    S = (spidentity(A.shape[0]).T  - omega*A)
    return S


def smoothed_prolongator(P,A):
    #just use Richardson for now
    omega = 4.0/(3.0*inf_norm(A))
    return P - omega*(A*P)
 



def sa_hierarchy(A,Ws,x):
    """
    Construct multilevel hierarchy using Smoothed Aggregation
        Inputs:
          A  - matrix
          Is - list of constant prolongators
          x  - "candidate" basis function to be approximated
        Ouputs:
          (As,Is,Ps) - tuple of lists
                  - As - [A, Ps[0].T*A*Ps[0], Ps[1].T*A*Ps[1], ... ]
                  - Is - smoothed prolongators 
                  - Ps - tentative prolongators
    """
    Ps = []
    Is = []
    As = [A]

    for W in Ws:
        P,x = fit_candidate(W,x)
        I   = smoothed_prolongator(P,A)  
        A   = I.T.tocsr() * A * I
        As.append(A)
        Ps.append(P)
        Is.append(I)
    return As,Is,Ps
         
def make_bridge(I,N):
    tail = I.indptr[-1].repeat(N - I.shape[0])
    ptr = concatenate((I.indptr,tail))
    return csr_matrix((I.data,I.indices,ptr),dims=(N,I.shape[1]),check=False)

class adaptive_sa_solver:
    def __init__(self,A,options=None):
        self.A = A

        self.Rs = [] 
        self.__construct_hierarchy(A)

    def __construct_hierarchy(self,A):
        #if self.A.shape[0] <= self.opts['coarse: max size']:
        #    raise ValueError,'small matrices not handled yet'
        
        x,AggOps = self.__initialization_stage(A) #first candidate
        Ws = AggOps
        
        #x[:] = 1 #TEST 

        self.candidates = [x]
        
        #create SA using x here
        As,Is,Ps = sa_hierarchy(A,Ws,x)

        for i in range(0):
            x           = self.__develop_candidate(A,As,Is,Ps,Ws,AggOps)
            #x[:]  = arange(x.shape[0])
            #x[x.shape[0]/2:] = 2*x[x.shape[0]/2] - x[x.shape[0]/2:]
            As,Is,Ps,Ws = self.__augment_cycle(A,As,Ps,Ws,AggOps,x)
             
            self.candidates.append(x)
        
        #As,Is,Ps = sa_hierarchy(A,AggOps,x)  #TESTING
        self.Ps = Ps 
        self.solver = multilevel_solver2(As,Is)
                


    def __develop_candidate(self,A,As,Is,Ps,Ws,AggOps):
        x = rand(A.shape[0])
        b = zeros_like(x)

        #x[:] = 1 #TEST

        mu = 5
 
        solver = multilevel_solver2(As,Is)

        for n in range(mu):
            x = solver.solve(b, x0=x, tol=1e-8, maxiter=1)
        #TEST FOR CONVERGENCE HERE
 
        A_l,P_l,W_l,x_l = As[0],Ps[0],Ws[0],x
        
        temp_Is = []
        for i in range(len(As) - 2):
            P_l_new, x_m, W_l_new, W_m_new = orthonormalize_prolongator(P_l, x_l, W_l, AggOps[i+1])    
 
            I_l_new = smoothed_prolongator(P_l_new,A_l)
            A_m_new = I_l_new.T.tocsr() * A_l * I_l_new
            bridge = make_bridge(Is[i+1],A_m_new.shape[0]) 
 
            temp_solver = multilevel_solver2( [A_m_new] + As[i+2:], [bridge] + Is[i+2:] )
 
            for n in range(mu):
                x_m = temp_solver.solve(zeros_like(x_m), x0=x_m, tol=1e-8, maxiter=1)
 
            temp_Is.append(I_l_new)
 
            W_l = vstack_csr(Ws[i+1],W_m_new)  #prepare for next iteration
            A_l = A_m_new
            x_l = x_m
            P_l = make_bridge(Ps[i+1],A_m_new.shape[0])
 
        x = x_l
        for I in reversed(temp_Is):
            x = I*x
        
        return x
           

    def __augment_cycle(self,A,As,Ps,Ws,AggOps,x):
        #As,Is,Ps,Ws = self.__augment_cycle(A,Ps,Ws,AggOps,x)

        #make a new cycle using the new candidate
        A_l,P_l,W_l,x_l = As[0],Ps[0],AggOps[0],x            
        
        new_As,new_Is,new_Ps,new_Ws = [A],[],[],[AggOps[0]]

        for i in range(len(As) - 2):
            P_l_new, x_m, W_l_new, W_m_new = orthonormalize_prolongator(P_l, x_l, W_l, AggOps[i+1])

            I_l_new = smoothed_prolongator(P_l_new,A_l)
            A_m_new = I_l_new.T.tocsr() * A_l * I_l_new
            W_m_new = vstack_csr(Ws[i+1],W_m_new)

            new_As.append(A_m_new)
            new_Ws.append(W_m_new)
            new_Is.append(I_l_new)
            new_Ps.append(P_l_new)

            #prepare for next iteration
            W_l = W_m_new 
            A_l = A_m_new
            x_l = x_m
            P_l = make_bridge(Ps[i+1],A_m_new.shape[0]) 
        
        P_l_new, x_m, W_l_new, W_m_new = orthonormalize_prolongator(P_l, x_l, W_l, csr_matrix((P_l.shape[1],1)))
        I_l_new = smoothed_prolongator(P_l_new,A_l)
        A_m_new = I_l_new.T.tocsr() * A_l * I_l_new

        new_As.append(A_m_new)
        new_Is.append(I_l_new)
        new_Ps.append(P_l_new)

        return new_As,new_Is,new_Ps,new_Ws


    def __initialization_stage(self,A):
        max_levels = 10
        max_coarse = 50

        AggOps = []
        Is     = []

        # aSA parameters 
        mu      = 5    # number of test relaxation iterations
        epsilon = 0.1  # minimum acceptable relaxation convergence factor 

        #step 1
        A_l = A
        x   = scipy.rand(A_l.shape[0])
        skip_f_to_i = False
        
        #step 2
        b = zeros_like(x)
        gauss_seidel(A_l,x,b,iterations=mu)
        #step 3
        #test convergence rate here 
       
        while len(AggOps) < max_levels  and A_l.shape[0] > max_coarse:
            W_l   = sa_constant_interpolation(A_l)                        #step 4b 
            #W_l   = sa_no_threshold(A_l)                        #step 4b  TEST
            P_l,x = fit_candidate(W_l,x)                                  #step 4c  
            I_l   = smoothed_prolongator(P_l,A_l)                         #step 4d
            A_l   = I_l.T.tocsr() * A_l * I_l                             #step 4e
            
            AggOps.append(W_l)
            Is.append(I_l)

            if A_l.shape <= max_coarse:  break

            if not skip_f_to_i:
                print "."
                x_hat = x.copy()                                        #step 4g
                gauss_seidel(A_l,x,zeros_like(x),iterations=mu)         #step 4h
                x_A_x = inner(x,A_l*x) 
                if (x_A_x/inner(x_hat,A_l*x_hat))**(1.0/mu) < epsilon:  #step 4i
                    print "sufficient convergence, skipping"
                    skip_f_to_i = True
                    if x_A_x == 0:       
                        x = x_hat  #need to restore x
        
        #update fine-level candidate
        for I in reversed(Is):
            x = I * x
                
        #gauss_seidel(A,x,zeros_like(x),iterations=mu)    #TEST
    
        #x[:] = 1 #TEST
    
        return x,AggOps  #first candidate,aggregation



from scipy import *
from pydec import diag_sparse
from multigrid import poisson_problem,poisson_problem1D
#A = poisson_problem(100).T
A = poisson_problem1D(100).T
D = diag_sparse(1.0/sqrt(10**(12*rand(A.shape[0])-6))).tocsr()
A = D * A * D
#A = A*A
#A = io.mmread("nos2.mtx").tocsr()
asa = adaptive_sa_solver(A)
x = rand(A.shape[0])
b = zeros_like(x)

resid = []

for n in range(50):
    x = asa.solver.solve(b,x)
    resid.append(linalg.norm(A*x))




