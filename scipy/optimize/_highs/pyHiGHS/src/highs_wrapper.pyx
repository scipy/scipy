# distutils: language=c++
# cython: language_level=3

import logging

from libc.stdio cimport FILE, tmpfile
from libcpp.memory cimport unique_ptr

from HConst cimport (
    ML_NONE,
    PrimalDualStatusSTATUS_FEASIBLE_POINT,
)
from Highs cimport Highs
from HighsStatus cimport (
    HighsStatus,
    HighsStatusError,
    HighsStatusWarning,
    HighsStatusOK,
)
from HighsLp cimport (
    HighsLp,
    HighsSolution,
    HighsModelStatus,
    HighsModelStatusNOTSET,
    HighsModelStatusLOAD_ERROR,
    HighsModelStatusMODEL_ERROR,
    HighsModelStatusMODEL_EMPTY,
    HighsModelStatusPRESOLVE_ERROR,
    HighsModelStatusSOLVE_ERROR,
    HighsModelStatusPOSTSOLVE_ERROR,
    HighsModelStatusPRIMAL_INFEASIBLE,
    HighsModelStatusPRIMAL_UNBOUNDED,
    HighsModelStatusOPTIMAL,
    HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
    HighsModelStatusREACHED_TIME_LIMIT,
    HighsModelStatusREACHED_ITERATION_LIMIT,
)
from HighsInfo cimport HighsInfo

cdef apply_options(dict options, Highs & highs):
    '''Take options from dictionary and apply to HiGHS object.'''

    # Send logging to dummy file to get rid of output from stdout
    cdef FILE * f
    if options.get('message_level', None) == ML_NONE:
        f = tmpfile()
        highs.setHighsLogfile(f)
        highs.setHighsOutput(f)

    # Do all the ints
    for opt in set([
            'allowed_simplex_cost_scale_factor',
            'allowed_simplex_matrix_scale_factor',
            'dual_simplex_cleanup_strategy',
            'ipm_iteration_limit',
            'keep_n_rows',
            'max_threads',
            'message_level',
            'min_threads',
            'simplex_crash_strategy',
            'simplex_dual_edge_weight_strategy',
            'simplex_dualise_strategy',
            'simplex_iteration_limit',
            'simplex_permute_strategy',
            'simplex_price_strategy',
            'simplex_primal_edge_weight_strategy',
            'simplex_scale_strategy',
            'simplex_strategy',
            'simplex_update_limit',
            'small_matrix_value',
    ]):
        val = options.get(opt, None)
        if val is not None:
            highs.setHighsOptionValueInt(opt.encode(), val)

    # Do all the doubles
    for opt in set([
            'dual_feasibility_tolerance',
            'dual_objective_value_upper_bound',
            'dual_simplex_cost_perturbation_multiplier',
            'dual_steepest_edge_weight_log_error_threshhold',
            'infinite_bound',
            'infinite_cost',
            'large_matrix_value',
            'primal_feasibility_tolerance',
            'simplex_initial_condition_tolerance',
            'small_matrix_value',
            'time_limit'
    ]):
        val = options.get(opt, None)
        if val is not None:
            highs.setHighsOptionValueDbl(opt.encode(), val)

    # Do all the strings
    for opt in set(['solver']):
        val = options.get(opt, None)
        if val is not None:
            highs.setHighsOptionValueStr(opt.encode(), val.encode())

    # Do all the bool to strings
    for opt in set([
            'parallel',
            'presolve',
    ]):
        val = options.get(opt, None)
        if val is not None:
            if val:
                val0 = b'on'
            else:
                val0 = b'off'
            highs.setHighsOptionValueStr(opt.encode(), val0)

    # Do the actual bools
    for opt in set([
            'less_infeasible_DSE_check',
            'less_infeasible_DSE_choose_row',
            'mps_parser_type_free',
            'run_as_hsol',
            'simplex_initial_condition_check',
            'use_original_HFactor_logic',
    ]):
        val = options.get(opt, None)
        if val is not None:
            highs.setHighsOptionValueBool(opt.encode(), val)

def highs_wrapper(
        double[::1] c,
        int[::1] astart,
        int[::1] aindex,
        double[::1] avalue,
        double[::1] lhs,
        double[::1] rhs,
        double[::1] lb,
        double[::1] ub,
        dict options):
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
    astart : 1-D array
    aindex : 1-D array
    avalue : 1-D array
    lhs : 1-D array (or None), (m,)
        Array of left hand side values of the inequality constraints.
        If `lhs=None`, then an array of `-inf` is assumed.
    rhs : 1-D array, (m,)
        Array of right hand side values of the inequality constraints.
    lb : 1-D array (or None), (n,)
        Lower bounds on solution variables x.  If `lb=None`, then an
        array of all `0` is assumed.
    ub : 1-D array (or None), (n,)
        Upper bounds on solution variables x.  If `ub=None`, then an
        array of `inf` is assumed.
    options : dict
        A dictionary of solver options with the following fields:

            - allowed_simplex_cost_scale_factor : int
                Undocumented advanced option.
            - allowed_simplex_matrix_scale_factor : int
                Undocumented advanced option.
            - dual_feasibility_tolerance : double
                Dual feasibility tolerance
            - dual_objective_value_upper_bound : double
                Upper bound on objective value for dual simplex:
                algorithm terminates if reached
            - dual_simplex_cleanup_strategy : int
                Undocumented advanced option.
            - dual_simplex_cost_perturbation_multiplier : double
                Undocumented advanced option.
            - dual_steepest_edge_weight_log_error_threshhold : double
                Undocumented advanced option.
            - infinite_bound : double
                Limit on abs(constraint bound): values larger than
                this will be treated as infinite
            - infinite_cost : double
                Limit on cost coefficient: values larger than this
                will be treated as infinite.
            - ipm_iteration_limit : int
                Iteration limit for interior-point solver.
            - keep_n_rows : int {-1, 0, 1}
                Undocumented advanced option.

                    - `-1`: KEEP_N_ROWS_DELETE_ROWS
                    - `0`: KEEP_N_ROWS_DELETE_ENTRIES
                    - `1`: KEEP_N_ROWS_KEEP_ROWS

            - large_matrix_value : double
                Upper limit on abs(matrix entries): values larger than
                this will be treated as infinite
            - less_infeasible_DSE_check : bool
                Undocumented advanced option.
            - less_infeasible_DSE_choose_row : bool
                Undocumented advanced option.
            - max_threads : int
                Maximum number of threads in parallel execution.
            - message_level : int {0, 1, 2, 4, 7}
                Verbosity level, corresponds to:

                    - `0`: ML_NONE
                    - `1`: ML_VERBOSE
                    - `2`: ML_DETAILED
                    - `4`: ML_MINIMAL
                    - `7`: ML_ALWAYS

            - min_threads : int
                Minimum number of threads in parallel execution.
            - mps_parser_type_free : bool
                Use free format MPS parsing.
            - parallel : bool
                Run the solver in serial (False) or parallel (True).
            - presolve : bool
                Run the presolve or not (or if `None`, then choose).
            - primal_feasibility_tolerance : double
                Primal feasibility tolerance.
            - run_as_hsol : bool
                Undocumented advanced option.
            - sense : int {1, -1}
                `sense=1` corresponds to the MIN problem, `sense=-1`
                corresponds to the MAX problem. TODO: NOT IMPLEMENTED
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

            - simplex_dualise_strategy : int
                Undocumented advanced option.
            - simplex_dual_edge_weight_strategy : int {0, 1, 2, 3, 4}
                Strategy for simplex dual edge weights:
                Dantzig / Devex / Steepest Edge. Corresponds
                to the following:

                    - `0`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG`
                    - `1`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX`
                    - `2`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH`
                    - `3`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE`
                    - `4`: `SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL`

            - simplex_initial_condition_check : bool
                Undocumented advanced option.
            - simplex_initial_condition_tolerance : double
                Undocumented advanced option.
            - simplex_iteration_limit : int
                Iteration limit for simplex solver.
            - simplex_permute_strategy : int
                Undocumented advanced option.
            - simplex_price_strategy : int
                Undocumented advanced option.
            - simplex_primal_edge_weight_strategy : int {0, 1}
                Strategy for simplex primal edge weights:
                Dantzig / Devex.  Corresponds to the following:

                    - `0`: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG`
                    - `1`: `SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX`

            - simplex_scale_strategy : int {0, 1, 2, 3, 4, 5}
                Strategy for scaling before simplex solver:
                off / on (0/1)

                    - `0`:  `SIMPLEX_SCALE_STRATEGY_OFF`
                    - `1`: `SIMPLEX_SCALE_STRATEGY_HIGHS`
                    - `2`: `SIMPLEX_SCALE_STRATEGY_HIGHS_FORCED`
                    - `3`: `SIMPLEX_SCALE_STRATEGY_HIGHS_015`
                    - `4`: `SIMPLEX_SCALE_STRATEGY_HIGHS_0157`
                    - `5`: `SIMPLEX_SCALE_STRATEGY_HSOL`

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
            - solution_file : str
                Solution file
            - solver : str {'simplex', 'ipm'}
                Choose which solver to use.  If `solver='simplex'`
                and `parallel=True` then PAMI will be used.
            - time_limit : double
                Max number of seconds to run the solver for.
            - use_original_HFactor_logic : bool
                Undocumented advanced option.
            - write_solution_to_file : bool
                Write the primal and dual solution to a file
            - write_solution_pretty : bool
                Write the primal and dual solution in a pretty
                (human-readable) format

        See [2]_ for a list of all options.

    Returns
    -------
    res : dict

        If model_status is one of OPTIMAL,
        REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND, REACHED_TIME_LIMIT,
        REACHED_ITERATION_LIMIT:

            - `status` : int
                Model status code.
            - `message` : str
                Message corresponding to model status code.
            - `x` : list
                Solution variables.
            - `slack` : list
                Slack variables.
            - `lambda` : list
                Lagrange multipliers assoicated with the constraints
                Ax = b.
            - `s` : list
                Lagrange multipliers associated with the constraints
                x >= 0.
            - `fun`
                Final objective value.
            - `simplex_nit` : int
                Number of iterations accomplished by the simplex
                solver.
            - `ipm_nit` : int
                Number of iterations accomplished by the interior-
                point solver.

        If model_status is not one of the above:

            - `status` : int
                Model status code.
            - `message` : str
                Message corresponding to model status code.

    Notes
    -----
    If `options['write_solution_to_file']` is `True` but
    `options['solution_file']` is unset or `''`, then the solution
    will be printed to `stdout`.

    If any iteration limit is reached, no solution will be
    available.

    References
    ----------
    .. [1] https://www.maths.ed.ac.uk/hall/HiGHS
    .. [2] https://www.maths.ed.ac.uk/hall/HiGHS/HighsOptions.html
    '''


    cdef int numcol = c.size
    cdef int numrow = rhs.size
    cdef int numnz = avalue.size

    # Fill up a HighsLp object
    cdef HighsLp lp
    lp.numCol_ = numcol
    lp.numRow_ = numrow

    lp.colCost_.resize(numcol)
    lp.colLower_.resize(numcol)
    lp.colUpper_.resize(numcol)

    lp.rowLower_.resize(numrow)
    lp.rowUpper_.resize(numrow)
    lp.Astart_.resize(numcol + 1)
    lp.Aindex_.resize(numnz)
    lp.Avalue_.resize(numnz)

    # Be careful not index into nothing
    cdef double * colcost_ptr = NULL
    cdef double * collower_ptr = NULL
    cdef double * colupper_ptr = NULL
    cdef double * rowlower_ptr = NULL
    cdef double * rowupper_ptr = NULL
    cdef int * astart_ptr = NULL
    cdef int * aindex_ptr = NULL
    cdef double * avalue_ptr = NULL
    if numrow > 0:
        rowlower_ptr = &lhs[0]
        rowupper_ptr = &rhs[0]
        lp.rowLower_.assign(rowlower_ptr, rowlower_ptr + numrow)
        lp.rowUpper_.assign(rowupper_ptr, rowupper_ptr + numrow)
    else:
        lp.rowLower_.empty()
        lp.rowUpper_.empty()
    if numcol > 0:
        colcost_ptr = &c[0]
        collower_ptr = &lb[0]
        colupper_ptr = &ub[0]
        lp.colCost_.assign(colcost_ptr, colcost_ptr + numcol)
        lp.colLower_.assign(collower_ptr, collower_ptr + numcol)
        lp.colUpper_.assign(colupper_ptr, colupper_ptr + numcol)
    else:
        lp.colCost_.empty()
        lp.colLower_.empty()
        lp.colUpper_.empty()
    if numnz > 0:
        astart_ptr = &astart[0]
        aindex_ptr = &aindex[0]
        avalue_ptr = &avalue[0]
        lp.Astart_.assign(astart_ptr, astart_ptr + numcol + 1)
        lp.Aindex_.assign(aindex_ptr, aindex_ptr + numnz)
        lp.Avalue_.assign(avalue_ptr, avalue_ptr + numnz)
    else:
        lp.Astart_.empty()
        lp.Aindex_.empty()
        lp.Avalue_.empty()

    # Create the options
    cdef Highs highs
    apply_options(options, highs)

    # Make a Highs object and pass it everything
    cdef HighsStatus init_status = highs.passModel(lp)
    if init_status != HighsStatusOK:
        if init_status != HighsStatusWarning:
            logging.warning("Error setting HighsLp");
            return <int> HighsStatusError

    # Solve the fool thing
    cdef HighsStatus run_status = highs.run()

    # Extract what we need from the solution
    #     HighsModelStatus
    #     solution
    #     dual solution
    #     objective value
    #     Number of iterations (simplex or ipm)
    #     sum of primal infeasibilities

    cdef HighsModelStatus model_status = highs.getModelStatus()
    cdef HighsModelStatus scaled_model_status = highs.getModelStatus(True)
    if model_status != scaled_model_status:
        if scaled_model_status == HighsModelStatusOPTIMAL:
            # The scaled model has been solved to optimality, but not the
            # unscaled model, flag this up, but report the scaled model
            # status
            logging.warning('model_status is not optimal, using scaled_model_status instead.')
            model_status = scaled_model_status

    # We might need an info object if we can look up the solution and a place to put solution
    cdef HighsInfo info = highs.getHighsInfo() # it should always be safe to get the info object
    cdef HighsSolution solution

    # If the status is bad, don't look up the solution
    if info.primal_status != PrimalDualStatusSTATUS_FEASIBLE_POINT:
#    if model_status in [
#            HighsModelStatusNOTSET,
#            HighsModelStatusLOAD_ERROR,
#            HighsModelStatusMODEL_ERROR,
#            HighsModelStatusMODEL_EMPTY,
#            HighsModelStatusPRESOLVE_ERROR,
#            HighsModelStatusSOLVE_ERROR,
#            HighsModelStatusPOSTSOLVE_ERROR,
#            HighsModelStatusPRIMAL_INFEASIBLE,
#            HighsModelStatusPRIMAL_UNBOUNDED,
#            HighsModelStatusREACHED_ITERATION_LIMIT,
#    ]:
        return {
            'status': <int> model_status,
            'message': highs.highsModelStatusToString(model_status).decode(),
            'simplex_nit': info.simplex_iteration_count,
            'ipm_nit': info.ipm_iteration_count,
            'fun': info.objective_function_value,
            'crossover_nit': info.crossover_iteration_count,
        }
    # If the model status is such that the solution can be read
    else:
        # Should be safe to read the solution:
        solution = highs.getSolution()
        return {
            'status': <int> model_status,
            'message': highs.highsModelStatusToString(model_status).decode(),

            # Primal solution
            'x': [solution.col_value[ii] for ii in range(numcol)],

            # Ax + s = b => Ax = b - s
            # Note: this is for all constraints (A_ub and A_eq)
            'slack': [rhs[ii] - solution.row_value[ii] for ii in range(numrow)],

            # slacks in HiGHS appear as Ax - s, not Ax + s, so lambda is negated;
            # lambda are the lagrange multipliers associated with Ax=b
            'lambda': [-1*solution.row_dual[ii] for ii in range(numrow)],

            # s are the lagrange multipliers associated with bound conditions
            's': [solution.col_dual[ii] for ii in range(numcol) if solution.col_dual[ii]],

            'fun': info.objective_function_value,
            'simplex_nit': info.simplex_iteration_count,
            'ipm_nit': info.ipm_iteration_count,
            'crossover_nit': info.crossover_iteration_count,
        }

# Export some things
from HConst cimport (
    HIGHS_CONST_I_INF,
    HIGHS_CONST_INF,
)
CONST_I_INF = HIGHS_CONST_I_INF
CONST_INF = HIGHS_CONST_INF

MODEL_STATUS_NOTSET = <int> HighsModelStatusNOTSET
MODEL_STATUS_LOAD_ERROR = <int> HighsModelStatusLOAD_ERROR
MODEL_STATUS_MODEL_ERROR = <int> HighsModelStatusMODEL_ERROR
MODEL_STATUS_MODEL_EMPTY = <int> HighsModelStatusMODEL_EMPTY
MODEL_STATUS_PRESOLVE_ERROR = <int> HighsModelStatusPRESOLVE_ERROR
MODEL_STATUS_SOLVE_ERROR = <int> HighsModelStatusSOLVE_ERROR
MODEL_STATUS_POSTSOLVE_ERROR = <int> HighsModelStatusPOSTSOLVE_ERROR
MODEL_STATUS_PRIMAL_INFEASIBLE = <int> HighsModelStatusPRIMAL_INFEASIBLE
MODEL_STATUS_PRIMAL_UNBOUNDED = <int> HighsModelStatusPRIMAL_UNBOUNDED
MODEL_STATUS_OPTIMAL = <int> HighsModelStatusOPTIMAL
MODEL_STATUS_REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND = <int> HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND
MODEL_STATUS_REACHED_TIME_LIMIT = <int> HighsModelStatusREACHED_TIME_LIMIT
MODEL_STATUS_REACHED_ITERATION_LIMIT = <int> HighsModelStatusREACHED_ITERATION_LIMIT
