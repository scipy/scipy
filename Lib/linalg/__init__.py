""" Linear algebra routines.

!!   There is some fortran/Python confusion here in places.  After
!!   clapack is fixed, get rid of confusion.
"""

from Numeric import *
#import cblas
#import clapack
#import fblas
import flapack

#-------------------------------------------------------------------
# Borrowed From LinearAlgebra.py
#-------------------------------------------------------------------

LinAlgError = 'LinAlgError'

# Helper routines
#lapack_type = {'f': 0, 'd': 1, 'F': 2, 'D': 3}
#lapack_letter = ['s', 'd', 'c', 'z']
_array_kind = {'i':0, 'l': 0, 'f': 0, 'd': 0, 'F': 1, 'D': 1}
_array_precision = {'i': 1, 'l': 1, 'f': 0, 'd': 1, 'F': 0, 'D': 1}
_array_type = [['f', 'd'], ['F', 'D']]

def _common_type(*arrays):
    kind = 0
#    precision = 0
#   force higher precision in lite version
    precision = 1
    for a in arrays:
        t = a.typecode()
        kind = max(kind, _array_kind[t])
        precision = max(precision, _array_precision[t])
    return _array_type[kind][precision]

def _assert_squareness(*arrays):
    for a in arrays:
        if max(a.shape) != min(a.shape):
            raise LinAlgError, 'Array must be square'
    
#-------------------------------------------------------------------
# End borrowed From LinearAlgebra.py
#-------------------------------------------------------------------

fortran_transpose = {}
fortran_transpose[0] = 'N'
fortran_transpose[1] = 'T'
fortran_transpose[2] = 'C'

c_transpose = {}
c_transpose[0] = 111
c_transpose[1] = 112
c_transpose[2] = 113

SingularError = 'SingularError'

def solve(a,b,overwrite=0,fortran=0):
    a_lu,pivots = lu_factor(a,overwrite,fortran=fortran) # fortran for now
    x =  lu_solve(a_lu,pivots,b,fortran=fortran) # fortran for now
    return x
    
def lu_factor(a, overwrite = 0,fortran=0):
    """ Factor the general dense matrix, a, into lower and upper traingular 
        matrices such that

                    a = PLU
                    
        where P is a permutation matrix, L is the lower triangular
        with unit diagonal elements and U is the the upper traingular.      
              
        lu_factor returns the factored array in matrix form with L in the 
        lower half and U in the upper half of the matrix.  A vector of pivots
        indices used during factorization is also returned (see XXX) 
        
        If overwrite = 0 (default), the contents of a are unaffected.  If 
        overwrite = 1, the array, a, is overwritten with the factored array.  
        This choice requires half the memory which can be helfpul for large 
        arrays.
        
        If fortran = 1, a is assumed stored in Fortran's column major ordering.
        If fortran = 0 (default), a is assumed stored in Python and C's 
        row major storage format. 
                
        lu_factor uses LAPACK's ?getrf() routine.
        
        Note: Currently only general dense matrices are supported, but their is
              no reason that, if a were a matrix class that specified its storage
              (sparse, symmetric, packed, etc.) that this routine couldn't handle
              these cases also.
    """
    # should look into using an LAPACK (BLAS) copy routine?
    if not overwrite:  aa = array(a,copy=1)
    else:              aa = a                       
    
    if fortran:        lapack = flapack
    else:              lapack = clapack
    
    atype = aa.typecode()
    
    # choose function based on type.  Unknown types will result in
    # using double percision.    
    if atype == Float32:     pivots,err = lapack.sgetrf(aa)
    elif atype == Float64:   pivots,err = lapack.dgetrf(aa)
    elif atype == Complex32: pivots,err = lapack.cgetrf(aa)
    elif atype == Complex64: pivots,err = lapack.zgetrf(aa)
    else:                    pivots,err = lapack.dgetrf(aa)
    
    if err < 0:
        raise TypeError, ("Argument %d to lapack's ?getrf() has an"  \
                          "illegal value." % abs(err))
    if err > 0:
        raise LinAlgError, ("Diagonal entry %d is zero.  Solving" \
                            " will cause a divide by zero error" % err)        
    return aa,pivots

def lu_solve(a_lu,pivots,b,a_transpose = 1,fortran=0):    
    """ Solve Ax = b. A must first be processed by lu_factor().
    
        a_lu    -- an NxN matrix. See a_transpose for more information.
        pivots  -- the partial pivot information about a_lu.  It, like a_lu,
                   is output from lu_factor(). 
        b       -- length N array. It is the right-hand-side of the equation.
        a_transpose -- 0 if a_lu is the lu factorization of A
                       1 (default) a_lu is the lu factorization of transpose(A) 
                       2 a_lu is the lu factorization of conjugate-transpose(A)
                       !! tranpose shouldn't be an option.
        fortran  --  0 (default), a is assumed stored in Python and C's row 
                       major storage format. 
                     1, a is assumed stored in Fortran's column major ordering.
                   
        
                       
        Note: We default to transpose(A) because this, combined with the default
              call to lu_factor() implicitly handles the C/Fortran array storage
              issue.  Fix this when CLAPACK is working!
    """    
    atype = _common_type(a_lu,b)
    _assert_squareness(a_lu)
    assert(len(a_lu) == len(b))     
    #how to handle casting of arrays?
       
    if fortran:
        lapack = flapack
        trans = fortran_transpose[a_transpose]
    else:
        lapack = clapack
        trans = c_transpose[a_transpose]
    
    if atype == Float32:   x,err = lapack.sgetrs(a_lu,pivots,b,trans)
    if atype == Float64:   x,err = lapack.dgetrs(a_lu,pivots,b,trans)
    if atype == Complex32: x,err = lapack.cgetrs(a_lu,pivots,b,trans)
    if atype == Complex64: x,err = lapack.zgetrs(a_lu,pivots,b,trans)
    if err < 0:
        raise TypeError, ("Argument %d to lapack's ?getrs() has an " \
                          "illegal value." % err)
    if err > 0:
        raise TypeError, ("Unknown error occured int ?getrs(): " \
                          "error code = %d" % err)
    return x
