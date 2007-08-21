from scipy import ones,zeros,rand,array,array_split,hstack,transpose,sum,ones_like,sqrt
from scipy.sparse import spidentity
from numpy.linalg import norm

import pydec
from pydec import diag_sparse,inf_norm, mls_polynomial_coeffs,polynomial_smoother
     
from multigrid import sa_interpolation,rs_interpolation
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
        


