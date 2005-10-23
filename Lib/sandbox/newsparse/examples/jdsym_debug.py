import math
import Numeric
import spmatrix, itsolvers, jdsym, precon, superlu

test = 1

if test == 1:
    # Test: compare K=eye with K=None
    #
    #   results are not the same, ritz-value in 2nd iterations differ

    path = '/local/home/geus/matrices'
    A = spmatrix.ll_mat_from_mtx(path + 'edge6x3x5_A.mtx')
    M = spmatrix.ll_mat_from_mtx(path + 'edge6x3x5_B.mtx')
    n = A.shape[0]
    sigma = 25.0

    I = spmatrix.ll_mat(n, n)
    for i in xrange(n):
        I[i,i] = 1.0
    Keye = precon.jacobi(I)
    
    k_conv, lmbd, Q, it, it_inner = jdsym.jdsym(A.to_sss(), M.to_sss(), None, 5, sigma, 1e-10, 15, itsolvers.qmrs,
                                                jmin=5, jmax=10, eps_tr=1e-4, toldecay=2.0, linitmax=200, clvl=1, strategy=1)

    k_conv, lmbd, Q, it, it_inner  = jdsym.jdsym(A.to_sss(), M.to_sss(), Keye, 5, sigma, 1e-10, 15, itsolvers.qmrs,
                                                 jmin=5, jmax=10, eps_tr=1e-4, toldecay=2.0, linitmax=200, clvl=1, strategy=1)

elif test == 2:
    
    # Test 2: K = diag(A - sigma*M), test diagonal prec using Matlab

    Asigma = A.copy()
    Asigma.shift(-sigma, M)
    K = precon.jacobi(Asigma.to_sss())
    
    b = Numeric.ones(n, 'd')
    x = Numeric.zeros(n, 'd')
    K.precon(b, x)
    print 'norm(idiag) = %.16g' % (math.sqrt(Numeric.dot(x, x)), )

    k_conv, lmbd, Q, it, it_inner  = jdsym.jdsym(A.to_sss(), M.to_sss(), K, 5, sigma, 1e-10, 150, itsolvers.qmrs,
                                       jmin=5, jmax=10, eps_tr=1e-4, toldecay=2.0, linitmax=200, clvl=1, strategy=1)

elif test == 3:
    
    Asigma = A.copy()
    Asigma.shift(-sigma, M)
    K = precon.ssor(Asigma.to_sss())
    k_conv, lmbd, Q, it, it_inner  = jdsym.jdsym(A.to_sss(), M.to_sss(), K, 5, sigma, 1e-10, 150, itsolvers.qmrs,
                                       jmin=5, jmax=10, eps_tr=1e-4, toldecay=2.0, linitmax=200, clvl=1, strategy=1)
    
elif test == 4:
    
    Asigma = A.copy()
    Asigma.shift(-sigma, M)
    K = precon.jacobi(Asigma.to_sss())
    k_conv, lmbd, Q, it, it_inner  = jdsym.jdsym(A.to_sss(), M.to_sss(), K, 5, sigma, 1e-10, 150, itsolvers.cgs,
                                       jmin=5, jmax=10, eps_tr=1e-4, toldecay=2.0, linitmax=200, clvl=1,
                                       strategy=1, optype=1)

elif test == 5:

    # time jdtest -e1 -k1 -o linsolver=QMRS,optype=SYM,jmin=5,jmax=10,linitmax=1000,jdtol=1e-6,strategy=1 1.4 1 ~/matrices/cop18_el5_A.mtx ~/matrices/cop18_el5_M.mtx

    path = '/home/geus/matrices/'
    A = spmatrix.ll_mat_from_mtx(path + 'cop18_el5_A.mtx')
    M = spmatrix.ll_mat_from_mtx(path + 'cop18_el5_M.mtx')
    n = A.shape[0]
    sigma = 1.4
    
    Asigma = A.copy()
    Asigma.shift(-sigma, M)
    K = precon.jacobi(Asigma.to_sss())
    
    k_conv, lmbd, Q, it, it_inner  = jdsym.jdsym(A.to_sss(), M.to_sss(), K, 1, sigma, 1e-6, 150, itsolvers.qmrs,
                                       jmin=5, jmax=10, eps_tr=1e-4, toldecay=2.0, linitmax=1000, clvl=1,
                                       strategy=1)
    
    print k_conv, lmbd, it, it_inner

elif test == 6:

    class prec2lev:
        def __init__(self, A, n11):
            n = A.shape[0]
            self.n11 = n11
            self.lu11 = superlu.factorize(A[:n11,:n11].to_csr(), diag_pivot_thresh=0)
            self.K22 = precon.ssor(A[n11:,n11:].to_sss())
            self.shape = (n, n)
        def precon(self, x, y):
            self.lu11.solve(x[:n11], y[:n11])
            self.K22.precon(x[n11:], y[n11:])

    print 'Loading matrices...'
    path = '/home/geus/matrices/'
    A = spmatrix.ll_mat_from_mtx(path + 'cop18_el5_A.mtx')
    M = spmatrix.ll_mat_from_mtx(path + 'cop18_el5_M.mtx')
    n = A.shape[0]
    sigma = 1.4
    n11 = 4688

    print 'Constructing preconditioner...'
    Asigma = A.copy()
    Asigma.shift(-sigma, M)
    K = prec2lev(Asigma, n11)

    k_conv, lmbd, Q, it, it_inner  = jdsym.jdsym(A.to_sss(), M.to_sss(), K, 2, sigma, 1e-8, 300, itsolvers.qmrs,
                                       jmin=10, jmax=25, eps_tr=1e-3, toldecay=1.5, linitmax=3000, clvl=1,
                                       strategy=1)

    print k_conv, lmbd, it, it_inner

