#from base import Struct, pause
import scipy.base as nm
import scipy.sparse as sp
import _umfpack as _um
import re

__doc__ = """
Interface to the UMFPACK library.
Routines for symbolic and numeric LU factorization of sparse
matrices and for solving systems of linear equations with sparse matrices.

Tested with UMFPACK V4.4 (Jan. 28, 2005),
Copyright (c) 2005 by Timothy A. Davis.  All Rights Reserved.
UMFPACK homepage: http://www.cise.ufl.edu/research/sparse/umfpack

Contains: UmfpackContext class

Use 'print UmfpackContext().funs' to see all UMFPACK library functions the
module exposes, if you need something not covered by the examples below.

Examples:
=========

Assuming this module imported as um

Sparse matrix in CSR or CSC format: mtx
Righthand-side: rhs
Solution: sol

# Contruct the solver.
umfpack = um.UmfpackContext() # Use default 'di' family of UMFPACK routines.

# One-shot solution.
sol = umfpack( um.UMFPACK_A, mtx, rhs, autoTranspose = True )
# same as:
sol = umfpack.linsolve( um.UMFPACK_A, mtx, rhs, autoTranspose = True )

-or-

# Make LU decomposition.
umfpack.numeric( mtx )
...
# Use already LU-decomposed matrix.
sol1 = umfpack( um.UMFPACK_A, mtx, rhs1, autoTranspose = True )
sol2 = umfpack( um.UMFPACK_A, mtx, rhs2, autoTranspose = True )
# same as:
sol1 = umfpack.solve( um.UMFPACK_A, mtx, rhs1, autoTranspose = True )
sol2 = umfpack.solve( um.UMFPACK_A, mtx, rhs2, autoTranspose = True )

-or-

# Make symbolic decomposition.
umfpack.symbolic( mtx0 )
# Print statistics.
umfpack.report_symbolic()

...

# Make LU decomposition of mtx1 which has same structure as mtx0.
umfpack.numeric( mtx1 )
# Print statistics.
umfpack.report_numeric()

# Use already LU-decomposed matrix.
sol1 = umfpack( um.UMFPACK_A, mtx1, rhs1, autoTranspose = True )

...

# Make LU decomposition of mtx2 which has same structure as mtx0.
umfpack.numeric( mtx2 )
sol2 = umfpack.solve( um.UMFPACK_A, mtx2, rhs2, autoTranspose = True )

# Print all statistics.
umfpack.report_info()

Setting control parameters:
===========================
Assuming this module imported as um:

List of control parameter names is accessible as 'um.umfControls' - their
meaning and possible values are described in the UMFPACK documentation.
To each name corresponds an attribute of the 'um' module, such as,
for example 'um.UMFPACK_PRL' (controlling the verbosity of umfpack report
functions).  These attributes are in fact indices into the control array
- to set the corresponding control array value, just do the following:

umfpack = um.UmfpackContext()
umfpack.control[um.UMFPACK_PRL] = 4 # Let's be more verbose.

--
Author: Robert Cimrman
"""

##
# 30.11.2005, c
def updateDictWithVars( adict, module, pattern, group = None ):
    match = re.compile( pattern ).match

    for name in [ii for ii in vars( module ).keys()
                 if match( ii )]:
        if group is not None:
            outName = match( name ).group( group )
        else:
            outName = name

        adict[outName] = module.__dict__[name]

    return adict

##
# Export UMFPACK constants from _um.
umfDefines = updateDictWithVars( {}, _um, 'UMFPACK_.*' )
locals().update( umfDefines )

##
# How to list these automagically?
umfControls = [
    'UMFPACK_PRL',
    'UMFPACK_DENSE_ROW',
    'UMFPACK_DENSE_COL',
    'UMFPACK_BLOCK_SIZE',
    'UMFPACK_STRATEGY',
    'UMFPACK_2BY2_TOLERANCE',
    'UMFPACK_FIXQ',
    'UMFPACK_AMD_DENSE',
    'UMFPACK_AGGRESSIVE',
    'UMFPACK_PIVOT_TOLERANCE',
    'UMFPACK_ALLOC_INIT',
    'UMFPACK_SYM_PIVOT_TOLERANCE',
    'UMFPACK_SCALE',
    'UMFPACK_FRONT_ALLOC_INIT',
    'UMFPACK_DROPTOL',
    'UMFPACK_IRSTEP',
    'UMFPACK_COMPILED_WITH_BLAS',
    'UMFPACK_COMPILED_FOR_MATLAB',
    'UMFPACK_COMPILED_WITH_GETRUSAGE',
    'UMFPACK_COMPILED_IN_DEBUG_MODE',
    'UMFPACK_STRATEGY_AUTO',
    'UMFPACK_STRATEGY_UNSYMMETRIC',
    'UMFPACK_STRATEGY_2BY2',
    'UMFPACK_STRATEGY_SYMMETRIC',
    'UMFPACK_SCALE_NONE',
    'UMFPACK_SCALE_SUM',
    'UMFPACK_SCALE_MAX',
]

umfInfo = [
    'UMFPACK_STATUS',
    'UMFPACK_NROW',
    'UMFPACK_NCOL',
    'UMFPACK_NZ',
    'UMFPACK_SIZE_OF_UNIT',
    'UMFPACK_SIZE_OF_INT',
    'UMFPACK_SIZE_OF_LONG',
    'UMFPACK_SIZE_OF_POINTER',
    'UMFPACK_SIZE_OF_ENTRY',
    'UMFPACK_NDENSE_ROW',
    'UMFPACK_NEMPTY_ROW',
    'UMFPACK_NDENSE_COL',
    'UMFPACK_NEMPTY_COL',
    'UMFPACK_SYMBOLIC_DEFRAG',
    'UMFPACK_SYMBOLIC_PEAK_MEMORY',
    'UMFPACK_SYMBOLIC_SIZE',
    'UMFPACK_SYMBOLIC_TIME',
    'UMFPACK_SYMBOLIC_WALLTIME',
    'UMFPACK_STRATEGY_USED',
    'UMFPACK_ORDERING_USED',
    'UMFPACK_QFIXED',
    'UMFPACK_DIAG_PREFERRED',
    'UMFPACK_PATTERN_SYMMETRY',
    'UMFPACK_NZ_A_PLUS_AT',
    'UMFPACK_NZDIAG',
    'UMFPACK_SYMMETRIC_LUNZ',
    'UMFPACK_SYMMETRIC_FLOPS',
    'UMFPACK_SYMMETRIC_NDENSE',
    'UMFPACK_SYMMETRIC_DMAX',
    'UMFPACK_2BY2_NWEAK',
    'UMFPACK_2BY2_UNMATCHED',
    'UMFPACK_2BY2_PATTERN_SYMMETRY',
    'UMFPACK_2BY2_NZ_PA_PLUS_PAT',
    'UMFPACK_2BY2_NZDIAG',
    'UMFPACK_COL_SINGLETONS',
    'UMFPACK_ROW_SINGLETONS',
    'UMFPACK_N2',
    'UMFPACK_S_SYMMETRIC',
    'UMFPACK_NUMERIC_SIZE_ESTIMATE',
    'UMFPACK_PEAK_MEMORY_ESTIMATE',
    'UMFPACK_FLOPS_ESTIMATE',
    'UMFPACK_LNZ_ESTIMATE',
    'UMFPACK_UNZ_ESTIMATE',
    'UMFPACK_VARIABLE_INIT_ESTIMATE',
    'UMFPACK_VARIABLE_PEAK_ESTIMATE',
    'UMFPACK_VARIABLE_FINAL_ESTIMATE',
    'UMFPACK_MAX_FRONT_SIZE_ESTIMATE',
    'UMFPACK_MAX_FRONT_NROWS_ESTIMATE',
    'UMFPACK_MAX_FRONT_NCOLS_ESTIMATE',
    'UMFPACK_NUMERIC_SIZE',
    'UMFPACK_PEAK_MEMORY',
    'UMFPACK_FLOPS',
    'UMFPACK_LNZ',
    'UMFPACK_UNZ',
    'UMFPACK_VARIABLE_INIT',
    'UMFPACK_VARIABLE_PEAK',
    'UMFPACK_VARIABLE_FINAL',
    'UMFPACK_MAX_FRONT_SIZE',
    'UMFPACK_MAX_FRONT_NROWS',
    'UMFPACK_MAX_FRONT_NCOLS',
    'UMFPACK_NUMERIC_DEFRAG',
    'UMFPACK_NUMERIC_REALLOC',
    'UMFPACK_NUMERIC_COSTLY_REALLOC',
    'UMFPACK_COMPRESSED_PATTERN',
    'UMFPACK_LU_ENTRIES',
    'UMFPACK_NUMERIC_TIME',
    'UMFPACK_UDIAG_NZ',
    'UMFPACK_RCOND',
    'UMFPACK_WAS_SCALED',
    'UMFPACK_RSMIN',
    'UMFPACK_RSMAX',
    'UMFPACK_UMIN',
    'UMFPACK_UMAX',
    'UMFPACK_ALLOC_INIT_USED',
    'UMFPACK_FORCED_UPDATES',
    'UMFPACK_NUMERIC_WALLTIME',
    'UMFPACK_NOFF_DIAG',
    'UMFPACK_ALL_LNZ',
    'UMFPACK_ALL_UNZ',
    'UMFPACK_NZDROPPED',
    'UMFPACK_IR_TAKEN',
    'UMFPACK_IR_ATTEMPTED',
    'UMFPACK_OMEGA1',
    'UMFPACK_OMEGA2',
    'UMFPACK_SOLVE_FLOPS',
    'UMFPACK_SOLVE_TIME',
    'UMFPACK_SOLVE_WALLTIME',
    'UMFPACK_ORDERING_COLAMD',
    'UMFPACK_ORDERING_AMD',
    'UMFPACK_ORDERING_GIVEN',
]

umfStatus = {
    UMFPACK_OK : 'UMFPACK_OK',
    UMFPACK_WARNING_singular_matrix : 'UMFPACK_WARNING_singular_matrix',
    UMFPACK_WARNING_determinant_underflow : 'UMFPACK_WARNING_determinant_underflow',
    UMFPACK_WARNING_determinant_overflow : 'UMFPACK_WARNING_determinant_overflow',
    UMFPACK_ERROR_out_of_memory : 'UMFPACK_ERROR_out_of_memory',
    UMFPACK_ERROR_invalid_Numeric_object : 'UMFPACK_ERROR_invalid_Numeric_object',
    UMFPACK_ERROR_invalid_Symbolic_object : 'UMFPACK_ERROR_invalid_Symbolic_object',
    UMFPACK_ERROR_argument_missing : 'UMFPACK_ERROR_argument_missing',
    UMFPACK_ERROR_n_nonpositive : 'UMFPACK_ERROR_n_nonpositive',
    UMFPACK_ERROR_invalid_matrix : 'UMFPACK_ERROR_invalid_matrix',
    UMFPACK_ERROR_different_pattern : 'UMFPACK_ERROR_different_pattern',
    UMFPACK_ERROR_invalid_system : 'UMFPACK_ERROR_invalid_system',
    UMFPACK_ERROR_invalid_permutation : 'UMFPACK_ERROR_invalid_permutation',
    UMFPACK_ERROR_internal_error : 'UMFPACK_ERROR_internal_error',
    UMFPACK_ERROR_file_IO : 'UMFPACK_ERROR_file_IO',
}

umfSys = [
    UMFPACK_A,
    UMFPACK_At,
    UMFPACK_Aat,
    UMFPACK_Pt_L,
    UMFPACK_L,
    UMFPACK_Lt_P,
    UMFPACK_Lat_P,
    UMFPACK_Lt,
    UMFPACK_U_Qt,
    UMFPACK_U,
    UMFPACK_Q_Ut,
    UMFPACK_Q_Uat,
    UMFPACK_Ut,
    UMFPACK_Uat,
]    

# Real, complex.
umfSys_transposeMap = [
    {UMFPACK_A : UMFPACK_At,
     UMFPACK_At : UMFPACK_A,
     UMFPACK_Aat : UMFPACK_A},
    {UMFPACK_A : UMFPACK_Aat,
     UMFPACK_Aat : UMFPACK_A}
]

umfFamilyTypes = {'di' : int, 'dl' : long, 'zi' : int, 'zl' : long}
umfRealTypes = ('di', 'dl')
umfComplexTypes = ('zi', 'zl')


##
# 02.01.2005
class Struct( object ):
    # 03.10.2005, c
    # 26.10.2005
    def __init__( self, **kwargs ):
        if kwargs:
            self.__dict__.update( kwargs )
        
    # 08.03.2005
    def __str__( self ):
        ss = "%s\n" % self.__class__
        for key, val in self.__dict__.iteritems():
            if (issubclass( self.__dict__[key].__class__, Struct )):
                ss += "  %s:\n    %s\n" % (key, self.__dict__[key].__class__)
            else:
                aux = "\n" + str( val )
                aux = aux.replace( "\n", "\n    " );
                ss += "  %s:\n%s\n" % (key, aux[1:])
        return( ss.rstrip() )

##
# 30.11.2005, c
class UmfpackContext( Struct ):

    ##
    # 30.11.2005, c
    # 01.12.2005
    def __init__( self, family = 'di', **kwargs ):
        """family .. family of UMFPACK functions ('di', 'dl', 'zi', 'zl')"""
        Struct.__init__( self, **kwargs )

        if family not in umfFamilyTypes.keys():
            raise TypeError, 'wrong family: %s' % family

        self.family = family
        self.control = nm.zeros( (UMFPACK_CONTROL, ), dtype = nm.double )
        self.info = nm.zeros( (UMFPACK_INFO, ), dtype = nm.double )
        self._symbolic = None
        self._numeric = None
        self.mtx = None
         
        ##
        # Functions corresponding to <family> are stored in self.funs.
        pattern = 'umfpack_' + family + '_(.*)'
        fn = updateDictWithVars( {}, _um, pattern, group = 1 )
        self.funs = Struct( **fn )

        self.funs.defaults( self.control )
        self.control[UMFPACK_PRL] = 3

    ##
    # 30.11.2005, c
    def strControl( self ):
        maxLen = max( [len( name ) for name in umfControls] )
        format = '%%-%ds : %%d' % maxLen
        aux = [format % (name, self.control[umfDefines[name]])
               for name in umfControls if umfDefines.has_key( name )]
        return '\n'.join( aux )

    ##
    # 01.12.2005, c
    def strInfo( self ):
        maxLen = max( [len( name ) for name in umfInfo] )
        format = '%%-%ds : %%d' % maxLen
        aux = [format % (name, self.info[umfDefines[name]])
               for name in umfInfo if umfDefines.has_key( name )]
        return '\n'.join( aux )

    ##
    # 30.11.2005, c
    # 01.12.2005
    def _getIndx( self, mtx ):

        ##
        # Should check types of indices to correspond to familyTypes.
        if sp.isspmatrix_csc( mtx ):
            indx = mtx.rowind
            self.isCSR = 1
        elif sp.isspmatrix_csr( mtx ):
            indx = mtx.colind
            self.isCSR = 0
        else:
            raise TypeError, 'must be a CSC/CSR matrix'

        return indx

    ##
    # 30.11.2005, c
    # 01.12.2005
    def symbolic( self, mtx ):
        """Symbolic object (symbolic LU decomposition) computation for a given
        sparsity pattern."""
        self.free_symbolic()
        
        indx = self._getIndx( mtx )
        status, self._symbolic\
                = self.funs.symbolic( mtx.shape[0], mtx.shape[1],
                                      mtx.indptr, indx, mtx.data,
                                      self.control, self.info )
##         print status, self._symbolic

        if status != UMFPACK_OK:
            raise RuntimeError, '%s failed with %s' % (self.funs.symbolic,
                                                       umfStatus[status])
            
        self.mtx = mtx

    ##
    # 30.11.2005, c
    # 01.12.2005
    # 02.12.2005
    def numeric( self, mtx ):
        """Numeric object (LU decomposition) computation using the
        symbolic decomposition. The symbolic decomposition is (re)computed
        if necessary."""

        self.free_numeric()

        if self._symbolic is None:
            self.symbolic( mtx )

        indx = self._getIndx( mtx )

        failCount = 0
        while 1:
            status, self._numeric\
                    = self.funs.numeric( mtx.indptr, indx, mtx.data,
                                         self._symbolic,
                                         self.control, self.info )
##             print status, self._numeric

            if status != UMFPACK_OK:
                if status == UMFPACK_WARNING_singular_matrix:
                    print 'warning: singular matrix'
                    break
                elif status in (UMFPACK_ERROR_different_pattern,
                                UMFPACK_ERROR_invalid_Symbolic_object):
                    # Try again.
                    print 'warning: recomputing symbolic'
                    self.symbolic( mtx )
                    failCount += 1
                else:
                    failCount += 100
            else:
                break
            if failCount >= 2:
                raise RuntimeError, '%s failed with %s' % (self.funs.numeric,
                                                           umfStatus[status])

    ##
    # 14.12.2005, c
    def report_symbolic( self ):
        """Print information about the symbolic object. Output depends on
        self.control[UMFPACK_PRL]."""
        self.funs.report_symbolic( self._symbolic, self.control )

    ##
    # 14.12.2005, c
    def report_numeric( self ):
        """Print information about the numeric object. Output depends on
        self.control[UMFPACK_PRL]."""
        self.funs.report_numeric( self._numeric, self.control )

    ##
    # 14.12.2005, c
    def report_info( self ):
        """Print all status information."""
        self.funs.report_info( self.control, self.info )

    ##
    # 30.11.2005, c
    # 01.12.2005
    def free_symbolic( self ):
        if self._symbolic is not None:
            self.funs.free_symbolic( self._symbolic )
            self._symbolic = None
            self.mtx = None

    ##
    # 30.11.2005, c
    # 01.12.2005
    def free_numeric( self ):
        if self._numeric is not None:
            self.funs.free_numeric( self._numeric )
            self._numeric = None
            self.free_symbolic()

    ##
    # 30.11.2005, c
    def free( self ):
        self.free_symbolic()
        self.free_numeric()

    ##
    # 30.11.2005, c
    # 01.12.2005
    # 02.12.2005
    def solve( self, sys, mtx, rhs, autoTranspose = False ):
        """Solution of system of linear equation using the Numeric object."""
        if sys not in umfSys:
            raise ValueError, 'sys must be in' % umfSys

        if autoTranspose and self.isCSR:
            ##
            # UMFPACK uses CSC internally...
            if self.family in umfRealTypes: ii = 0
            else: ii = 1
            if umfSys_transposeMap[ii].has_key( sys ):
                sys = umfSys_transposeMap[ii][sys]
            else:
                raise RuntimeError, 'autoTranspose ambiguous, switch it off'

        if self._numeric is not None:
            if self.mtx is not mtx:
                raise ValueError, 'must be called with same matrix as numeric()'
        else:
            raise RuntimeError, 'numeric() not called'

        sol = nm.zeros( (mtx.shape[1],), dtype = nm.double )
        indx = self._getIndx( mtx )
        status = self.funs.solve( sys, mtx.indptr, indx, mtx.data, sol, rhs,
                                  self._numeric, self.control, self.info )
        #self.funs.report_info( self.control, self.info )
        #pause()
        if status != UMFPACK_OK:
            if status == UMFPACK_WARNING_singular_matrix:
                ## Change inf, nan to zeros.
                print 'zeroing nan and inf entries...'
                sol[~nm.isfinite( sol )] = 0.0
            else:
                raise RuntimeError, '%s failed with %s' % (self.funs.solve,
                                                           umfStatus[status])
        return sol

    ##
    # 30.11.2005, c
    # 01.12.2005
    def linsolve( self, sys, mtx, rhs, autoTranspose = False ):
        """One-shot solution of system of linear equation. Reuses Numeric
        object if possible."""

        if sys not in umfSys:
            raise ValueError, 'sys must be in' % umfSys

        if self._numeric is None:
            self.numeric( mtx )
        else:
            if self.mtx is not mtx:
                self.numeric( mtx )

        sol = self.solve( sys, mtx, rhs, autoTranspose )
        self.free_numeric()

        return sol

    ##
    # 30.11.2005, c
    # 01.12.2005
    def __call__( self, sys, mtx, rhs, autoTranspose = False ):
        """Uses solve() or linsolve() depending on the presence of
        the Numeric object."""

        if self._numeric is not None:
            return self.solve( sys, mtx, rhs, autoTranspose )
        else:
            return self.linsolve( sys, mtx, rhs, autoTranspose )
