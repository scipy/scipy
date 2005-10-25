import os
import spmatrix, itsolvers, jdsym, precon

path = os.path.join(os.environ['HOME'], 'matrices')
A = spmatrix.ll_mat_from_mtx(os.path.join(path, 'edge6x3x5_A.mtx'))
M = spmatrix.ll_mat_from_mtx(os.path.join(path, 'edge6x3x5_B.mtx'))

sigma = 25.0
Asigma = A.copy()
Asigma.shift(-sigma, M)
K = precon.jacobi(Asigma.to_sss())
del Asigma

##n = A.shape[0]
##I = spmatrix.ll_mat(n, n)
##for i in xrange(n):
##    I[i,i] = 1.0
##K = precon.jacobi(I)

k_conv, lmbd, Q, it, it_inner  = \
        jdsym.jdsym(A.to_sss(), M.to_sss(), K, 5, sigma, 1e-10, 150, itsolvers.qmrs,
                    jmin=5, jmax=10, eps_tr=1e-4, toldecay=2.0, linitmax=200, clvl=1, strategy=1)
print k_conv, lmbd, it, it_inner


##path = '/homes/geus/jdbsym/test/'
##A = spmatrix.ll_mat_from_mtx(path + 'sameh_10000.mtx').to_sss()

##K = precon.ssor(A)
##k_conv, lmbd, Q, it, it_inner = jdsym.jdsym(A, None, K, 4, 0.0, 1e-10, 50, itsolvers.qmrs,
##                                   jmax=20, eps_tr=1e-4, toldecay=2.0, linitmax=60, clvl=1)


##path = '../spmatrix/'
##A_ll = spmatrix.ll_mat_from_mtx(path + 'cop18_el5_A.mtx')
##M_ll = spmatrix.ll_mat_from_mtx(path + 'cop18_el5_M.mtx')

##A = A_ll.to_sss()
##M = M_ll.to_sss()

##k_conv, lmbd, Q, it, it_inner = jdsym.jdsym(A, M, None, 5, 5.0, 1e-4, 50, itsolvers.minres,
##                                   linitmax=3000, strategy=1, clvl=1)


##path = '/homes/geus/jdbsym/test/'
##A = spmatrix.ll_mat_from_mtx(path + 'node4x3x1_A.mtx').to_sss()
##M = spmatrix.ll_mat_from_mtx(path + 'node4x3x1_M.mtx').to_sss()

##K = precon.jacobi(A)
##k_conv, lmbd, Q, it, it_inner = jdsym.jdsym(A, M, K, 15, 0.0, 1e-10, 150, itsolvers.qmrs,
##                                   jmin=15, jmax=30, eps_tr=1e-4, toldecay=2.0, linitmax=60, clvl=1)

##print lmbd
