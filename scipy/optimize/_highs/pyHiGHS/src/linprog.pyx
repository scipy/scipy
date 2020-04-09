# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE, tmpfile

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, make_unique

cimport numpy as np
import numpy as np
#from scipy.sparse import csc_matrix
#from scipy.optimize import OptimizeResult
from warnings import warn

from HConst cimport (
    # Constants
    HIGHS_CONST_I_INF as _HIGHS_CONST_I_INF,
    HIGHS_CONST_INF as _HIGHS_CONST_INF,
    HIGHS_CONST_TINY as _HIGHS_CONST_TINY,
    HIGHS_CONST_ZERO as _HIGHS_CONST_ZERO,
    HIGHS_THREAD_LIMIT as _HIGHS_THREAD_LIMIT,

    # Verbosity levels
    ML_DETAILED as _ML_DETAILED,
    ML_NONE as _ML_NONE,
    ML_VERBOSE as _ML_VERBOSE,
    ML_MINIMAL as _ML_MINIMAL,

    # HighsBasisStatus
    HighsBasisStatus,
    LOWER,
    BASIC,
    UPPER,
    ZERO,
    NONBASIC,
    SUPER,

    # Solvers
    SOLVER_OPTION_SIMPLEX as _SOLVER_OPTION_SIMPLEX,
    SOLVER_OPTION_CHOOSE as _SOLVER_OPTION_CHOOSE,
    SOLVER_OPTION_IPM as  _SOLVER_OPTION_IPM)
from Highs cimport Highs
from HighsLp cimport (
    HighsSolution,
    HighsBasis,

    # Model Statuses
    HighsModelStatus,
    HighsModelStatusNOTSET as _HighsModelStatusNOTSET,
    HighsModelStatusLOAD_ERROR as _HighsModelStatusLOAD_ERROR,
    HighsModelStatusMODEL_ERROR as _HighsModelStatusMODEL_ERROR,
    HighsModelStatusMODEL_EMPTY as _HighsModelStatusMODEL_EMPTY,
    HighsModelStatusPRESOLVE_ERROR as _HighsModelStatusPRESOLVE_ERROR,
    HighsModelStatusSOLVE_ERROR as _HighsModelStatusSOLVE_ERROR,
    HighsModelStatusPOSTSOLVE_ERROR as _HighsModelStatusPOSTSOLVE_ERROR,
    HighsModelStatusPRIMAL_INFEASIBLE as _HighsModelStatusPRIMAL_INFEASIBLE,
    HighsModelStatusPRIMAL_UNBOUNDED as _HighsModelStatusPRIMAL_UNBOUNDED,
    HighsModelStatusOPTIMAL as _HighsModelStatusOPTIMAL,
    HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND as _HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
    HighsModelStatusREACHED_TIME_LIMIT as _HighsModelStatusREACHED_TIME_LIMIT,
    HighsModelStatusREACHED_ITERATION_LIMIT as _HighsModelStatusREACHED_ITERATION_LIMIT)
from HighsInfo cimport HighsInfo
from highs_c_api cimport Highs_passLp

# Export some constants to Python
HIGHS_CONST_I_INF = _HIGHS_CONST_I_INF
HIGHS_CONST_INF = _HIGHS_CONST_INF
HIGHS_CONST_TINY = _HIGHS_CONST_TINY
HIGHS_CONST_ZERO = _HIGHS_CONST_ZERO
HIGHS_THREAD_LIMIT = _HIGHS_THREAD_LIMIT
ML_NONE = _ML_NONE
ML_VERBOSE = _ML_VERBOSE
ML_DETAILED = _ML_DETAILED
ML_MINIMAL = _ML_MINIMAL
SOLVER_OPTION_CHOOSE = _SOLVER_OPTION_CHOOSE
SOLVER_OPTION_SIMPLEX = _SOLVER_OPTION_SIMPLEX
SOLVER_OPTION_IPM = _SOLVER_OPTION_IPM

HighsModelStatusNOTSET = <int>_HighsModelStatusNOTSET
HighsModelStatusLOAD_ERROR = <int>_HighsModelStatusLOAD_ERROR
HighsModelStatusMODEL_ERROR = <int>_HighsModelStatusMODEL_ERROR
HighsModelStatusMODEL_EMPTY = <int>_HighsModelStatusMODEL_EMPTY
HighsModelStatusPRESOLVE_ERROR = <int>_HighsModelStatusPRESOLVE_ERROR
HighsModelStatusSOLVE_ERROR = <int>_HighsModelStatusSOLVE_ERROR
HighsModelStatusPOSTSOLVE_ERROR = <int>_HighsModelStatusPOSTSOLVE_ERROR
HighsModelStatusPRIMAL_INFEASIBLE = <int>_HighsModelStatusPRIMAL_INFEASIBLE
HighsModelStatusPRIMAL_UNBOUNDED = <int>_HighsModelStatusPRIMAL_UNBOUNDED
HighsModelStatusOPTIMAL = <int>_HighsModelStatusOPTIMAL
HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND = <int>_HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND
HighsModelStatusREACHED_TIME_LIMIT = <int>_HighsModelStatusREACHED_TIME_LIMIT
HighsModelStatusREACHED_ITERATION_LIMIT = <int>_HighsModelStatusREACHED_ITERATION_LIMIT


cdef int Highs_call(int numcol, int numrow, int numnz, double* colcost,
                    double* collower, double* colupper, double* rowlower,
                    double* rowupper, int* astart, int* aindex, double* avalue,
                    double* colvalue, double* coldual, double* rowvalue,
                    double* rowdual, int* colbasisstatus, int* rowbasisstatus,
                    int* modelstatus, Highs & highs):
    # cdef Highs highs
    cdef int status = Highs_passLp(&highs, numcol, numrow, numnz, colcost, collower, colupper,
                                   rowlower, rowupper, astart, aindex, avalue)

    # Customize sense : MIN or MAX
    # This API is not currently working, do it manually in caller
    # highs.changeObjectiveSense(sense)

    if (status != 0):
        return status
    status = <int>highs.run()

    # See how we did
    cdef int model_status = <int>highs.getModelStatus()
    cdef int scaled_model_status = <int>highs.getModelStatus(True);
    if model_status != scaled_model_status:
        if scaled_model_status == HighsModelStatusOPTIMAL:
            # The scaled model has been solved to optimality, but not the
            # unscaled model, flag this up, but report the scaled model
            # status
            model_status = scaled_model_status

    cdef unique_ptr[HighsSolution] solution
    cdef HighsBasis basis
    if (status == 0 and model_status in [
            HighsModelStatusOPTIMAL,
            HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
            HighsModelStatusREACHED_TIME_LIMIT,
            HighsModelStatusREACHED_ITERATION_LIMIT]):
        solution = make_unique[HighsSolution](highs.getSolution())
        basis = highs.getBasis()

        # Set the modelstatus for return
        modelstatus[0] = <int>highs.getModelStatus()

        for ii in range(numcol):
            colvalue[ii] = solution.get().col_value[ii]
            coldual[ii] = solution.get().col_dual[ii]
            colbasisstatus[ii] = <int>basis.col_status[ii]

        for ii in range(numrow):
            rowvalue[ii] = solution.get().row_value[ii]
            rowdual[ii] = solution.get().row_dual[ii]
            rowbasisstatus[ii] = <int>basis.row_status[ii]
    else:
        # Let them know something was wrong with model
        modelstatus[0] = model_status

        # Set all statuses to a custom unset value
        for ii in range(numcol):
            colbasisstatus[ii] = -1
        for ii in range(numrow):
            rowbasisstatus[ii] = -1

    return status

cdef apply_options(dict options, Highs & highs):
    '''Take options from dictionary and apply to HiGHS object.'''

    # Send logging to dummy file to get rid of output from stdout
    cdef FILE * f
    if options.get('message_level', None) == ML_NONE:
        f = tmpfile()
        highs.setHighsLogfile(f)

    # Do all the ints
    for opt in [
            'max_threads',
            'message_level',
            'min_threads',
            'simplex_crash_strategy',
            'simplex_dual_edge_weight_strategy',
            'simplex_iteration_limit',
            'simplex_primal_edge_weight_strategy',
            'simplex_strategy',
            'simplex_update_limit',
            'small_matrix_value']:
        val = options.get(opt, None)
        if val is not None:
            highs.setHighsOptionValueInt(opt.encode(), val)

    # Do all the doubles
    for opt in [
            'dual_feasibility_tolerance',
            'dual_objective_value_upper_bound',
            'infinite_bound',
            'infinite_cost',
            'large_matrix_value',
            'primal_feasibility_tolerance',
            'small_matrix_value',
            'time_limit']:
        val = options.get(opt, None)
        if val is not None:
            highs.setHighsOptionValueDbl(opt.encode(), val)

    # Do all the strings
    for opt in ['solver']:
        val = options.get(opt, None)
        if val is not None:
            highs.setHighsOptionValueStr(opt.encode(), val.encode())

    # Do all the bool to strings
    for opt in ['parallel', 'presolve']:
        val = options.get(opt, None)
        if val is not None:
            if val:
                val0 = b'on'
            else:
                val0 = b'off'
            highs.setHighsOptionValueStr(opt.encode(), val0)

def highs_wrapper(
        double[::1] c,
        A,
        double[::1] rhs,
        double[::1] lhs=None,
        double[::1] lb=None,
        double[::1] ub=None,
        dict options=None):
    '''Solve linear programs using HiGHS [1]_.

    Assume problems of the form:

        MIN/MAX c.T @ x
        s.t. lhs <= A @ x <= rhs
             lb <= x <= ub

    Default is MIN (for MAX set `sense=-1`).

    Parameters
    ----------
    c : 1-D array, (n,)
        Array of objective value coefficients.
    A : 2-D array, (m, n)
        Sparse (or dense) matrix of constraint coefficients.
    rhs : 1-D array, (m,)
        Array of right hand side values of the inequality constraints.
    lhs : 1-D array (or None), (m,)
        Array of left hand side values of the inequality constraints.
        If `lhs=None`, then an array of `-inf` is assumed.
    lb : 1-D array (or None), (n,)
        Lower bounds on solution variables x.  If `lb=None`, then an
        array of all `0` is assumed.
    ub : 1-D array (or None), (n,)
        Upper bounds on solution variables x.  If `ub=None`, then an
        array of `inf` is assumed.
    options : dict
        A dictionary of solver options with the following fields:

            - dual_feasibility_tolerance : double
                Dual feasibility tolerance
            - dual_objective_value_upper_bound : double
                Upper bound on objective value for dual simplex:
                algorithm terminates if reached
            - infinite_bound : double
                Limit on abs(constraint bound): values larger than
                this will be treated as infinite
            - infinite_cost : double
                Limit on cost coefficient: values larger than this
                will be treated as infinite.
            - large_matrix_value : double
                Upper limit on abs(matrix entries): values larger than
                this will be treated as infinite
            - max_threads : int
                Maximum number of threads in parallel execution.
            - message_level : int {0, 1, 2, 4}
                Verbosity level, corresponds to:

                    - `0`: ML_NONE
                    - `1`: ML_VERBOSE
                    - `2`: ML_DETAILED
                    - `4`: ML_MINIMAL

            - min_threads : int
                Minimum number of threads in parallel execution.
            - parallel : bool
                Run the solver in serial (False) or parallel (True).
            - presolve : bool
                Run the presolve or not (or if `None`, then choose).
            - primal_feasibility_tolerance : double
                Primal feasibility tolerance.
            - sense : int {1, -1}
                `sense=1` corresponds to the MIN problem, `sense=-1`
                corresponds to the MAX problem.
            - simplex_crash_strategy : int {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
                Strategy for simplex crash: off / LTSSF / Bixby (0/1/2).
                Default is `0`.  Corresponds to the following:

                    - `0`: `SIMPLEX_CRASH_STRATEGY_OFF`
                    - `1`: `SIMPLEX_CRASH_STRATEGY_LTSSF_K`
                    - `2`: `SIMPLEX_CRASH_STRATEGY_BIXBY`
                    - `3`: `SIMPLEX_CRASH_STRATEGY_LTSSF_PRI`
                    - `4`: `SIMPLEX_CRASH_STRATEGY_LTSF_K`
                    - `5`: `SIMPLEX_CRASH_STRATEGY_LTSF_PRI`
                    - `6`: `SIMPLEX_CRASH_STRATEGY_LTSF`
                    - `7`: `SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS`
                    - `8`: `SIMPLEX_CRASH_STRATEGY_BASIC`
                    - `9`: `SIMPLE_CRASH_STRATEGY_TEST_SING`

            - simplex_dual_edge_weight_strategy : int {0, 1, 2, 3, 4}
                Strategy for simplex dual edge weights:
                Dantzig / Devex / Steepest Edge. Corresponds
                to the following:

                    - `0`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG`
                    - `1`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX`
                    - `2`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH`
                    - `3`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE`
                    - `4`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL`

            - simplex_iteration_limit : int
                Iteration limit for simplex solver.

            - simplex_primal_edge_weight_strategy : int {0, 1}
                Strategy for simplex primal edge weights:
                Dantzig / Devex.  Corresponds to the following:

                    - `0`: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG`
                    - `1`: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX`

            - simplex_strategy : int {0, 1, 2, 3, 4}
                Strategy for simplex solver. Default: 1. Corresponds
                to the following:

                    - `0`: `SIMPLEX_STRATEGY_MIN`
                    - `1`: `SIMPLEX_STRATEGY_DUAL`
                    - `2`: `SIMPLEX_STRATEGY_DUAL_TASKS`
                    - `3`: `SIMPLEX_STRATEGY_DUAL_MULTI`
                    - `4`: `SIMPLEX_STRATEGY_PRIMAL`

            - simplex_update_limit : int
                Limit on the number of simplex UPDATE operations.
            - small_matrix_value : double
                Lower limit on abs(matrix entries): values smaller
                than this will be treated as zero.
            - solver : str {'simplex', 'ipm'}
                Choose which solver to use.  If `solver='simplex'`
                and `parallel=True` then PAMI will be used.
            - time_limit : double
                Max number of seconds to run the solver for.
            - solution_file : str
                Solution file
            - write_solution_to_file : bool
                Write the primal and dual solution to a file
            - write_solution_pretty : bool
                Write the primal and dual solution in a pretty
                (human-readable) format

        See [2]_ for a list of all options.

    Returns
    -------
    res : dict

        - col_basis_status : dict
            Key: `'statuses'` contains `n` status codes corresponding
                 to the `n` columns.
            Key: `'messages'` contains the `n` messages corresponding
                 to each status.
        - col_dual : 1-D array, (n,)
            The dual solution.
        - col_value : 1-D array, (n,)
            Solution variables.
        - crossover_nit : int
            Number of iterations taken to transform the interior
            solution produced by barrier into a basic solution
        - dual_status : dict
            Key: `'status'` contains the dual solution status code.
            Key: `'message'` contains the corresponding message.
        - fun : double
            The final objective value.
        - ipm_nit : int
            Number of iterations taken by IPM (interior-point solver).
        - max_dual_infeasibility : double
        - max_primal_infeasibility : double
        - model_status : dict
            Key: `'status'` contains the status code of the LP model.
            Key: `'message'` contains the corresponding message.
        - num_dual_infeasibilities : int
        - num_primal_infeasibilities : int
        - primal_status : dict
            Key: `'status'` contains the primal solution status code.
            Key: `'message'` contains the corresponding message.
        - row_basis_status : dict
            Key: `'statuses'` contains `m` status codes corresponding
                 to the `m` rows.
            Key: `'messages'` contains the `m` messages corresponding
                 to each status.
        - row_dual : 1-D array, (m,)
        - simplex_nit : int
            Number of iterations taken by the simplex solver.
        - sum_dual_infeasibilities : double
        - sum_primal_infeasibilities : double

    Notes
    -----
    If `options['write_solution_to_file']` is `True` but
    `options['solution_file']` is unset or `''`, then the solution
    will be printed to `stdout`.

    References
    ----------
    .. [1] https://www.maths.ed.ac.uk/hall/HiGHS
    .. [2] https://www.maths.ed.ac.uk/hall/HiGHS/HighsOptions.html
    '''

    # Assume we are working with a sparse matrix for now;
    # Try to decouple scikit-highs and scipy for easier development
    # with both locally installed
    # if not isinstance(A, csc_matrix):
    #     A = csc_matrix(A)

    # Get dimensions of problem
    cdef int numrow = A.shape[0]
    cdef int numcol = A.shape[1]
    cdef int numnz = A.nnz

    # Objective function coefficients
    # Do MIN/MAX conversion here because API not working for HiGHS
    cdef double[::1] cc = c.copy() # doing a copy here -- prefer to have HiGHS do this
    if options.get('sense', 1) == -1:
        for ii in range(numcol):
            cc[ii] *= -1
    cdef double * colcost = &cc[0]

    # Bounds on variables
    cdef double * collower
    cdef double * colupper
    if lb is None:
        # Default is lower bound of 0
        lb = np.zeros(numcol, dtype='double')
    if ub is None:
        # Default is upper bound of inf
        ub = HIGHS_CONST_INF*np.ones(numcol, dtype='double')
    collower = &lb[0]
    colupper = &ub[0]

    # LHS/RHS constraints
    cdef double * rowlower = NULL
    cdef double * rowupper = NULL
    if lhs is None:
        # Default to no LHS (all -Inf)
        lhs = HIGHS_CONST_TINY*np.ones(numrow, dtype='double')

    # Contents of constraint matrices as memoryviews
    cdef int[::1] Aindptr = A.indptr
    cdef int[::1] Aindices = A.indices
    cdef double[::1] Adata = A.data
    cdef int * astart = NULL
    cdef int * aindex = NULL
    cdef double * avalue = NULL

    # Allocate memoryviews to hold results
    cdef double[::1] colvalue = np.empty(numcol, dtype='double')
    cdef double[::1] coldual = np.empty(numcol, dtype='double')
    cdef double[::1] rowvalue = np.empty(numrow, dtype='double')
    cdef double[::1] rowdual = np.empty(numrow, dtype='double')
    cdef double * rowvalue_ptr = NULL
    cdef double * rowdual_ptr = NULL

    # Result status flags
    cdef int[::1] colbasisstatus = np.empty(numcol, dtype=np.int32)
    cdef int[::1] rowbasisstatus = np.empty(numrow, dtype=np.int32)
    cdef int * rowbasisstatus_ptr = NULL
    cdef int modelstatus = 0

    # If we have no rows, then we can't index into the memoryviews to get pointers
    if numrow > 0:
        rowlower = &lhs[0]
        rowupper = &rhs[0]
        rowvalue_ptr = &rowvalue[0]
        rowdual_ptr = &rowdual[0]
        rowbasisstatus_ptr = &rowbasisstatus[0]
    if Aindptr.size > 0:
        astart = &Aindptr[0]
    if Aindices.size > 0:
        aindex = &Aindices[0]
    if Adata.size > 0:
        avalue = &Adata[0]


    # Instantiate the HiGHS object that will hold the model
    cdef Highs highs

    # Apply any options
    apply_options(options, highs)

    # Call the solver
    cdef int ret = Highs_call(
        numcol, numrow, numnz,
        colcost, collower, colupper,
        rowlower, rowupper,
        astart, aindex, avalue,
        &colvalue[0], &coldual[0], rowvalue_ptr, rowdual_ptr,
        &colbasisstatus[0], rowbasisstatus_ptr, &modelstatus,
        highs)

    # Pull info out of out of highs
    cdef HighsInfo info = highs.getHighsInfo()

    # If the model is unset, it means we've encountered an error during optimization
    if modelstatus == 0:
        # It could also mean that the problem is unbounded
        raise RuntimeError("Model failed during optimization! Could be unbounded! Try `presolve=False` and/or `method='simplex'`")

    # Maybe write to file
    if options.get('write_solution_to_file', None):
        outfile = options.get('solution_file', '')
        outpretty = options.get('write_solution_pretty', False)
        highs.writeSolution(outfile.encode(), outpretty)

    # Decode HighsBasisStatus:
    HighsBasisStatusToStr = {
        -1 : 'UNSET', # custom
        <int>LOWER: 'LOWER: (slack) variable is at its lower bound [including fixed variables]',
        <int>BASIC: 'BASIC: (slack) variable is basic',
        <int>UPPER: 'UPPER: (slack) variable is at its upper bound',
        <int>ZERO: 'ZERO: free variable is non-basic and set to zero',
        <int>NONBASIC: 'NONBASIC: nonbasic with no specific bound information - useful for users and postsolve',
        <int>SUPER: 'SUPER: Super-basic variable: non-basic and either free and nonzero or not at a bound. No SCIP equivalent',
    }

    return {
        # From HighsInfo
        'fun': info.objective_function_value,
        'simplex_nit': info.simplex_iteration_count,
        'ipm_nit': info.ipm_iteration_count,
        'crossover_nit': info.crossover_iteration_count,
        'primal_status': {
            'status': info.primal_status,
            'message': highs.highsPrimalDualStatusToString(info.primal_status).decode(),
        },
        'dual_status': {
            'status': info.dual_status,
            'message': highs.highsPrimalDualStatusToString(info.dual_status).decode(),
        },
        'num_primal_infeasibilities': info.num_primal_infeasibilities,
        'max_primal_infeasibility': info.max_primal_infeasibility,
        'sum_primal_infeasibilities': info.sum_primal_infeasibilities,
        'num_dual_infeasibilities': info.num_dual_infeasibilities,
        'max_dual_infeasibility': info.max_dual_infeasibility,
        'sum_dual_infeasibilities': info.sum_dual_infeasibilities,

        # From C API
        'col_value': np.array(colvalue),
        'col_dual': np.array(coldual),
        'row_value': np.array(rowvalue),
        'row_dual': np.array(rowdual),
        'col_basis_status': {
            'statuses': [colbasisstatus[ii] for ii in range(numcol)],
            'messages': [HighsBasisStatusToStr[colbasisstatus[ii]] for ii in range(numcol)],
        },
        'row_basis_status': {
            'statuses': [rowbasisstatus[ii] for ii in range(numrow)],
            'messages': [HighsBasisStatusToStr[rowbasisstatus[ii]] for ii in range(numrow)],
        },
        'model_status': {
            'status': modelstatus,
            'message': highs.highsModelStatusToString(<HighsModelStatus>modelstatus).decode(),
        },
    }



################################ MPS SOLVER ####################################

from libc.stdio cimport FILE

from libcpp cimport bool
from libcpp.memory cimport unique_ptr, allocator, make_unique
from libcpp.string cimport string

from HConst cimport ML_ALWAYS
from HighsOptions cimport HighsOptions
from HighsRuntimeOptions cimport loadOptions
from HighsIO cimport HighsPrintMessage
from HighsLp cimport (
    HighsLp,
    HighsModelStatus)
from HighsStatus cimport (
    HighsStatus,
    HighsStatusToString,
    HighsStatusOK,
    HighsStatusWarning,
    HighsStatusError)
from LoadProblem cimport loadLpFromFile
from HighsInfo cimport HighsInfo
from Highs cimport Highs
from HighsMipSolver cimport (
    HighsMipStatus,
    HighsMipStatuskOptimal,
    HighsMipSolver)

cdef void reportLpStatsOrError(FILE* output, int message_level, const HighsStatus read_status, const HighsLp& lp):
    if read_status == HighsStatusError:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Error loading file\n")
    else:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "LP       : %s\n", lp.model_name_.c_str())
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Rows     : %d\n", lp.numRow_)
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Cols     : %d\n", lp.numCol_)
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Nonzeros : %d\n", lp.Avalue_.size())
        if lp.numInt_:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "Integer  : %d\n", lp.numInt_)

cdef void reportSolvedLpStats(FILE* output, int message_level, const HighsStatus run_status, const Highs& highs):
    cdef string statusname
    cdef HighsModelStatus model_status
    cdef HighsModelStatus scaled_model_status
    cdef HighsInfo highs_info
    cdef double objective_function_value = 0 # initialized but written over for cython to stop complaining
    cdef const HighsOptions * options

    if run_status == HighsStatusError:
        statusname = HighsStatusToString(run_status)
        HighsPrintMessage(output, message_level, ML_ALWAYS, "HiGHS status: %s\n", statusname.c_str())
    else:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "\n")
        model_status = highs.getModelStatus()
        scaled_model_status = highs.getModelStatus(True)
        highs_info = highs.getHighsInfo()
        if model_status != scaled_model_status:
            if scaled_model_status == _HighsModelStatusOPTIMAL:
                HighsPrintMessage(output, message_level, ML_ALWAYS,
                                  "Primal infeasibility: %10.3e (%d)\n",
                                  highs_info.max_primal_infeasibility,
                                  highs_info.num_primal_infeasibilities);
                HighsPrintMessage(output, message_level, ML_ALWAYS,
                                  "Dual   infeasibility: %10.3e (%d)\n",
                                  highs_info.max_dual_infeasibility,
                                  highs_info.num_dual_infeasibilities);
                model_status = scaled_model_status;

        HighsPrintMessage(output, message_level, ML_ALWAYS, "Model   status      : %s\n", highs.highsModelStatusToString(model_status).c_str())
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Simplex   iterations: %d\n", highs_info.simplex_iteration_count);
        if highs_info.ipm_iteration_count:
            HighsPrintMessage(output, message_level, ML_ALWAYS,
                              "IPM       iterations: %d\n",
                              highs_info.ipm_iteration_count)
        if highs_info.crossover_iteration_count:
            HighsPrintMessage(output, message_level, ML_ALWAYS,
                              "Crossover iterations: %d\n",
                              highs_info.crossover_iteration_count)
        if model_status == _HighsModelStatusOPTIMAL:
            highs.getHighsInfoValue("objective_function_value".encode(), objective_function_value)
            HighsPrintMessage(output, message_level, ML_ALWAYS,
                              "Objective value     : %13.6e\n",
                              objective_function_value)

        # Possibly write the solution to a file
        options = &highs.getHighsOptions()
        if options.write_solution_to_file:
            highs.writeSolution(options.solution_file, options.write_solution_pretty)

cdef HighsStatus callLpSolver(const HighsOptions& options, const HighsLp& lp, FILE* output, int message_level, bool run_quiet):
    # Solve LP case.
    cdef Highs highs
    cdef HighsStatus return_status = highs.passHighsOptions(options)
    if return_status != HighsStatusOK:
        if return_status == HighsStatusWarning:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "HighsStatus::Warning return from passHighsOptions\n")
        else:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "In main: fail return from passHighsOptions\n")
        return return_status

    if run_quiet:
        highs.setHighsLogfile(NULL)
        highs.setHighsOutput(NULL)

    cdef HighsStatus init_status = highs.passModel(lp)
    if init_status != HighsStatusOK:
        if init_status == HighsStatusWarning:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "HighsStatus::Warning return setting HighsLp\n")
        else:
            HighsPrintMessage(output, message_level, ML_ALWAYS, "Error setting HighsLp\n")
        return HighsStatusError

    highs.writeHighsOptions("".encode())

    if run_quiet:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "Before calling highs.run()\n")

    # Run HiGHS.
    cdef HighsStatus run_status = highs.run()

    if run_quiet:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "After calling highs.run()\n")

    reportSolvedLpStats(output, message_level, run_status, highs)
    return run_status

cdef HighsStatus callMipSolver(const HighsOptions& options, const HighsLp& lp, FILE* output, int message_level, bool run_quiet):
    #cdef HighsMipSolver solver(options, lp)
    cdef unique_ptr[HighsMipSolver] solver = make_unique[HighsMipSolver](options, lp)
    cdef HighsMipStatus status = solver.get().runMipSolver()
    if status == HighsMipStatuskOptimal:
        return HighsStatusOK
    return HighsStatusError

def linprog_mps(model_file, presolve=None, solver=None, bool run_quiet=True):
    '''Solve linear program described in an MPS model file.

    Parameters
    ----------
    model_file : str
        Filename of uncompressed .MPS file.
    presolve : bool or None, optional
        Whether to run presolve or not. Values correspond to the HiGHS
        options:

            - `True`: `'on'`
            - `False`: `'off'`
            - `None` : `'choose'`

    solver : str or None {'simplex', 'ipm', None}, optional
        Method used to solve the LP. `solver=None` corresponds to the
        HiGHS option of `'choose'`.
    run_quiet : bool, optional
        Diplay lots of info or just some info.
    '''

    # Map some of the inputs to the correct HiGHS options; everything
    # should be a string after this.
    presolve = {
        True: 'on',
        False: 'off',
        None: 'choose',
    }[presolve]
    if solver is None:
        solver = 'choose'

    # Parse the inputs and put into char**
    args = {
        b'--model_file': model_file.encode(),
        b'--solver': solver.encode(),
        b'--presolve': presolve.encode(),
    }
    cdef allocator[char *] ptr_al
    cdef unique_ptr[char *] argv
    argv.reset(ptr_al.allocate(len(args)*2+1))
    argv.get()[0] = 'highs' # name of program in argv[0]
    for ii, (k, v) in enumerate(args.items()):
        argv.get()[2*ii+1] = k
        argv.get()[2*ii+2] = v

    # Load user options.
    cdef HighsOptions options
    cdef bool options_ok = loadOptions(len(args)*2+1, argv.get(), options)
    if not options_ok:
        return 0

    # Set message level.
    cdef FILE* output = options.output
    cdef int message_level = options.message_level

    #cdef bool run_quiet = True #False
    if run_quiet:
        HighsPrintMessage(output, message_level, ML_ALWAYS, "In main: running highs.run() quietly\n")
    output = options.output
    message_level = options.message_level

    # Load problem.
    cdef HighsLp lp
    cdef HighsStatus read_status = loadLpFromFile(options, lp)
    reportLpStatsOrError(output, message_level, read_status, lp)
    if read_status == HighsStatusError:
        return <int>HighsStatusError

    # Run LP or MIP solver.
    cdef HighsStatus run_status = HighsStatusError
    if not options.mip:
        run_status = callLpSolver(options, lp, output, message_level, run_quiet)
    else:
        run_status = callMipSolver(options, lp, output, message_level, run_quiet)

    return <int>run_status
