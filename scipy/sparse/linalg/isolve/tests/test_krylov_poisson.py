#########################################################################
# File Name: test_krylov_poisson.py
# Author: Gang Zhao
# Mail: zhaog6@lsec.cc.ac.cn
# Created Time: 02/04/2021 Thursday 10:35:15
#########################################################################

###
 # Add new CG, CGS, improved-CGS, BiCG, BiCGSTAB, CGNR, CGNE, GMRES
 # Test for the benchmark problem: 
 #      1. Poisson equations with homogeneous Dirichlet boundary condition
 #         - \Delta u = f
 #                  u = 0
 ##

import numpy as np
import time
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg as scipy_cg
from scipy.sparse.linalg import cgs as scipy_cgs
from scipy.sparse.linalg import bicg as scipy_bicg
from scipy.sparse.linalg import bicgstab as scipy_bicgstab
from scipy.sparse.linalg import gmres as scipy_gmres
from scipy.sparse.linalg.isolve.utils import make_system
from scipy.sparse.linalg.isolve.krylov import *

#------------------
#   Loading
#------------------
nn = 128
print("Loading Poisson...")
cpu = time.time()
mat = np.loadtxt('data/poisson_mat_'+str(nn)+'x'+str(nn)+'.dat')
b = np.loadtxt('data/poisson_rhs_'+str(nn)+'x'+str(nn)+'.dat')
cpu = time.time() - cpu
print(" --- global data loaded (in {})".format(cpu))

#------------------------
#   Building CSR
#------------------------
print("Generating CSR...")
cpu = time.time()
I = mat[:, 0] - 1
J = mat[:, 1] - 1
V = mat[:, 2]
A = csr_matrix((V, (I, J)), shape=(len(b), len(b)))
cpu = time.time() - cpu
print(" --- global CSR created (in {})".format(cpu))
print("--------------------------------------------------")

#------------------------
#   SciPy/CG
#------------------------
print("#------------------------")
print("#   SciPy/CG")
print("#------------------------")
print("SciPy/Original CG: beginning solve...")
cpu = time.time()
(x, flag) = scipy_cg(A, b, tol=1e-8)
cpu = time.time() - cpu
print(" --- system solved with SciPy/Original CG (in {} s)".format(cpu))
print("SciPy/CG: beginning solve...")
(y, flag) = cg(A, b, tol=1e-8)
e = x - y
print("||e|| = {}".format(np.linalg.norm(e)))
r = b - A * y
print("||r||/||r_0|| = {}".format(np.linalg.norm(r)/np.linalg.norm(b)))
print("")

#------------------------
#   SciPy/CGS
#------------------------
print("#------------------------")
print("#   SciPy/CGS")
print("#------------------------")
print("SciPy/Original CGS: beginning solve...")
cpu = time.time()
(x, flag) = scipy_cgs(A, b, tol=1e-8)
cpu = time.time() - cpu
print(" --- system solved with SciPy/Original CGS (in {} s)".format(cpu))
print("SciPy/CGS: beginning solve...")
(y, flag) = cgs(A, b, tol=1e-8)
e = x - y
print("||e|| = {}".format(np.linalg.norm(e)))
r = b - A * y
print("||r||/||r_0|| = {}".format(np.linalg.norm(r)/np.linalg.norm(b)))
print("")

#-------------------------
#   SciPy/Improved CGS
#-------------------------
print("#------------------------")
print("#   SciPy/Improved CGS")
print("#------------------------")
print("SciPy/Original CGS: beginning solve...")
cpu = time.time()
(x, flag) = scipy_cgs(A, b, tol=1e-8)
cpu = time.time() - cpu
print(" --- system solved with SciPy/Original CGS (in {} s)".format(cpu))
print("SciPy/Improved CGS: beginning solve...")
(y, flag) = icgs(A, b, tol=1e-8)
e = x - y
print("||e|| = {}".format(np.linalg.norm(e)))
r = b - A * y
print("||r||/||r_0|| = {}".format(np.linalg.norm(r)/np.linalg.norm(b)))
print("")

#-------------------------
#   SciPy/BiCG
#-------------------------
print("#------------------------")
print("#   SciPy/BiCG")
print("#------------------------")
print("SciPy/Original BiCG: beginning solve...")
cpu = time.time()
(x, flag) = scipy_bicg(A, b, tol=1e-8)
cpu = time.time() - cpu
print(" --- system solved with SciPy/Original BiCG (in {} s)".format(cpu))
print("SciPy/BiCG: beginning solve...")
(y, flag) = bicg(A, b, tol=1e-8)
e = x - y
print("||e|| = {}".format(np.linalg.norm(e)))
r = b - A * y
print("||r||/||r_0|| = {}".format(np.linalg.norm(r)/np.linalg.norm(b)))
print("")

#-------------------------
#   SciPy/BiCGSTAB
#-------------------------
print("#------------------------")
print("#   SciPy/BiCGSTAB")
print("#------------------------")
print("SciPy/Original BiCGSTAB: beginning solve...")
cpu = time.time()
(x, flag) = scipy_bicgstab(A, b, tol=1e-8)
cpu = time.time() - cpu
print(" --- system solved with SciPy/Original BiCGSTAB (in {} s)".format(cpu))
print("SciPy/BiCGSTAB: beginning solve...")
(y, flag) = bicgstab(A, b, tol=1e-8)
e = x - y
print("||e|| = {}".format(np.linalg.norm(e)))
r = b - A * y
print("||r||/||r_0|| = {}".format(np.linalg.norm(r)/np.linalg.norm(b)))
print("")

#-------------------------
#   SciPy/CGNR
#-------------------------
print("#------------------------")
print("#   SciPy/CGNR")
print("#------------------------")
print("SciPy/CGNR: beginning solve...")
(x, flag) = cgnr(A, b, tol=1e-8)
print("")

#-------------------------
#   SciPy/CGNE
#-------------------------
print("#------------------------")
print("#   SciPy/CGNE")
print("#------------------------")
print("SciPy/CGNE: beginning solve...")
(x, flag) = cgne(A, b, tol=1e-8)
print("")

#-------------------------
#   SciPy/GMRES
#-------------------------
print("#------------------------")
print("#   SciPy/GMRES")
print("#------------------------")
print("SciPy/Original GMRES: beginning solve...")
cpu = time.time()
(x, flag) = scipy_gmres(A, b, tol=1e-8)  # restart = 20 by default
cpu = time.time() - cpu
print(" --- system solved with SciPy/Original GMRES (in {} s)".format(cpu))
print("SciPy/GMRES: beginning solve...")
(y, flag) = gmres(A, b, restart=20, tol=1e-8)
e = x - y
print("||e|| = {}".format(np.linalg.norm(e)))
r = b - A * y
print("||r||/||r_0|| = {}".format(np.linalg.norm(r)/np.linalg.norm(b)))
