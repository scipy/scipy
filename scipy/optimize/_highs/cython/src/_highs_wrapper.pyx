# distutils: language=c++
# cython: language_level=3

import numpy as np
cimport numpy as np
from scipy.optimize import OptimizeWarning
from warnings import warn
import numbers

from libcpp.string cimport string
from libcpp.map cimport map as cppmap
from libcpp.cast cimport reinterpret_cast

from .HConst cimport (
    HIGHS_CONST_INF,

    HighsModelStatus,
    HighsModelStatusNOTSET,
    HighsModelStatusMODEL_ERROR,
    HighsModelStatusOPTIMAL,
    HighsModelStatusREACHED_TIME_LIMIT,
    HighsModelStatusREACHED_ITERATION_LIMIT,

    HighsOptionTypeBOOL,
    HighsOptionTypeINT,
    HighsOptionTypeDOUBLE,
    HighsOptionTypeSTRING,

    HighsBasisStatus,
    HighsBasisStatusLOWER,
    HighsBasisStatusUPPER,

    MatrixFormatkColwise,
    HighsVarType,
)
from .Highs cimport Highs
from .HighsStatus cimport (
    HighsStatus,
    highsStatusToString,
    HighsStatusError,
    HighsStatusWarning,
    HighsStatusOK,
)
from .HighsLp cimport (
    HighsLp,
    HighsSolution,
    HighsBasis,
)
from .HighsInfo cimport HighsInfo
from .HighsOptions cimport (
    HighsOptions,
    OptionRecord,
    OptionRecordBool,
    OptionRecordInt,
    OptionRecordDouble,
    OptionRecordString,
)
from .HighsModelUtils cimport utilBasisStatusToString

np.import_array()

# options to reference for default values and bounds;
# make a map to quickly lookup
cdef HighsOptions _ref_opts
cdef cppmap[string, OptionRecord*] _ref_opt_lookup
cdef OptionRecord * _r = NULL
for _r in _ref_opts.records:
    _ref_opt_lookup[_r.name] = _r


cdef str _opt_warning(string name, val, valid_set=None) noexcept:
    cdef OptionRecord * r = _ref_opt_lookup[name]

    # BOOL
    if r.type == HighsOptionTypeBOOL:
        default_value = (<OptionRecordBool*> r).default_value
        return ('Option "%s" is "%s", but only True or False is allowed. '
                'Using default: %s.' % (name.decode(), str(val), default_value))

    # INT
    if r.type == HighsOptionTypeINT:
        lower_bound = int((<OptionRecordInt*> r).lower_bound)
        upper_bound = int((<OptionRecordInt*> r).upper_bound)
        default_value = int((<OptionRecordInt*> r).default_value)
        if upper_bound - lower_bound < 10:
            int_range = str(set(range(lower_bound, upper_bound + 1)))
        else:
            int_range = '[%d, %d]' % (lower_bound, upper_bound)
        return ('Option "%s" is "%s", but only values in %s are allowed. '
                'Using default: %d.' % (name.decode(), str(val), int_range, default_value))

    # DOUBLE
    if r.type == HighsOptionTypeDOUBLE:
        lower_bound = (<OptionRecordDouble*> r).lower_bound
        upper_bound = (<OptionRecordDouble*> r).upper_bound
        default_value = (<OptionRecordDouble*> r).default_value
        return ('Option "%s" is "%s", but only values in (%g, %g) are allowed. '
                'Using default: %g.' % (name.decode(), str(val), lower_bound, upper_bound, default_value))

    # STRING
    if r.type == HighsOptionTypeSTRING:
        if valid_set is not None:
            descr = 'but only values in %s are allowed. ' % str(set(valid_set))
        else:
            descr = 'but this is an invalid value. %s. ' % r.description.decode()
        default_value = (<OptionRecordString*> r).default_value.decode()
        return ('Option "%s" is "%s", '
                '%s'
                'Using default: %s.' % (name.decode(), str(val), descr, default_value))

    # We don't know what type (should be unreachable)?
    return('Option "%s" is "%s", but this is not a valid value. '
           'See documentation for valid options. '
           'Using default.' % (name.decode(), str(val)))

cdef apply_options(dict options, Highs & highs) noexcept:
    '''Take options from dictionary and apply to HiGHS object.'''

    # Initialize for error checking
    cdef HighsStatus opt_status = HighsStatusOK

    # Do all the ints
    for opt in set([
            'allowed_simplex_cost_scale_factor',
            'allowed_simplex_matrix_scale_factor',
            'dual_simplex_cleanup_strategy',
            'ipm_iteration_limit',
            'keep_n_rows',
            'threads',
            'mip_max_nodes',
            'highs_debug_level',
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
            if not isinstance(val, int):
                warn(_opt_warning(opt.encode(), val), OptimizeWarning)
            else:
                opt_status = highs.setHighsOptionValueInt(opt.encode(), val)
                if opt_status != HighsStatusOK:
                    warn(_opt_warning(opt.encode(), val), OptimizeWarning)
                else:
                    if opt == "threads":
                        highs.resetGlobalScheduler(blocking=True)

    # Do all the doubles
    for opt in set([
            'dual_feasibility_tolerance',
            'dual_objective_value_upper_bound',
            'dual_simplex_cost_perturbation_multiplier',
            'dual_steepest_edge_weight_log_error_threshhold',
            'infinite_bound',
            'infinite_cost',
            'ipm_optimality_tolerance',
            'large_matrix_value',
            'primal_feasibility_tolerance',
            'simplex_initial_condition_tolerance',
            'small_matrix_value',
            'start_crossover_tolerance',
            'time_limit',
            'mip_rel_gap'
    ]):
        val = options.get(opt, None)
        if val is not None:
            if not isinstance(val, numbers.Number):
                warn(_opt_warning(opt.encode(), val), OptimizeWarning)
            else:
                opt_status = highs.setHighsOptionValueDbl(opt.encode(), val)
                if opt_status != HighsStatusOK:
                    warn(_opt_warning(opt.encode(), val), OptimizeWarning)


    # Do all the strings
    for opt in set(['solver']):
        val = options.get(opt, None)
        if val is not None:
            if not isinstance(val, str):
                warn(_opt_warning(opt.encode(), val), OptimizeWarning)
            else:
                opt_status = highs.setHighsOptionValueStr(opt.encode(), val.encode())
                if opt_status != HighsStatusOK:
                    warn(_opt_warning(opt.encode(), val), OptimizeWarning)


    # Do all the bool to strings
    for opt in set([
            'parallel',
            'presolve',
    ]):
        val = options.get(opt, None)
        if val is not None:
            if isinstance(val, bool):
                if val:
                    val0 = b'on'
                else:
                    val0 = b'off'
                opt_status = highs.setHighsOptionValueStr(opt.encode(), val0)
                if opt_status != HighsStatusOK:
                    warn(_opt_warning(opt.encode(), val, valid_set=[True, False]), OptimizeWarning)
            else:
                warn(_opt_warning(opt.encode(), val, valid_set=[True, False]), OptimizeWarning)


    # Do the actual bools
    for opt in set([
            'less_infeasible_DSE_check',
            'less_infeasible_DSE_choose_row',
            'log_to_console',
            'mps_parser_type_free',
            'output_flag',
            'run_as_hsol',
            'run_crossover',
            'simplex_initial_condition_check',
            'use_original_HFactor_logic',
    ]):
        val = options.get(opt, None)
        if val is not None:
            if val in [True, False]:
                opt_status = highs.setHighsOptionValueBool(opt.encode(), val)
                if opt_status != HighsStatusOK:
                    warn(_opt_warning(opt.encode(), val), OptimizeWarning)
            else:
                warn(_opt_warning(opt.encode(), val), OptimizeWarning)


ctypedef HighsVarType* HighsVarType_ptr


def _highs_wrapper(
        double[::1] c,
        int[::1] astart,
        int[::1] aindex,
        double[::1] avalue,
        double[::1] lhs,
        double[::1] rhs,
        double[::1] lb,
        double[::1] ub,
        np.uint8_t[::1] integrality,
        dict options):
    '''Solve linear programs using HiGHS [1]_.

    Assume problems of the form:

        MIN c.T @ x
        s.t. lhs <= A @ x <= rhs
             lb <= x <= ub

    Parameters
    ----------
    c : 1-D array, (n,)
        Array of objective value coefficients.
    astart : 1-D array
        CSC format index array.
    aindex : 1-D array
        CSC format index array.
    avalue : 1-D array
        Data array of the matrix.
    lhs : 1-D array (or None), (m,)
        Array of left hand side values of the inequality constraints.
        If ``lhs=None``, then an array of ``-inf`` is assumed.
    rhs : 1-D array, (m,)
        Array of right hand side values of the inequality constraints.
    lb : 1-D array (or None), (n,)
        Lower bounds on solution variables x.  If ``lb=None``, then an
        array of all `0` is assumed.
    ub : 1-D array (or None), (n,)
        Upper bounds on solution variables x.  If ``ub=None``, then an
        array of ``inf`` is assumed.
    options : dict
        A dictionary of solver options with the following fields:

            - allowed_simplex_cost_scale_factor : int
                Undocumented advanced option.

            - allowed_simplex_matrix_scale_factor : int
                Undocumented advanced option.

            - dual_feasibility_tolerance : double
                Dual feasibility tolerance for simplex.
                ``min(dual_feasibility_tolerance,
                primal_feasibility_tolerance)`` will be used for
                ipm feasibility tolerance.

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

            - ipm_optimality_tolerance : double
                Optimality tolerance for IPM.

            - keep_n_rows : int {-1, 0, 1}
                Undocumented advanced option.

                    - ``-1``: ``KEEP_N_ROWS_DELETE_ROWS``
                    - ``0``: ``KEEP_N_ROWS_DELETE_ENTRIES``
                    - ``1``: ``KEEP_N_ROWS_KEEP_ROWS``

            - large_matrix_value : double
                Upper limit on abs(matrix entries): values larger than
                this will be treated as infinite

            - less_infeasible_DSE_check : bool
                Undocumented advanced option.

            - less_infeasible_DSE_choose_row : bool
                Undocumented advanced option.

            - threads : int
                Maximum number of threads in parallel execution.

            - message_level : int {0, 1, 2, 4, 7}
                Verbosity level, corresponds to:

                    - ``0``: ``ML_NONE``
                        All messaging to stdout is supressed.

                    - ``1``: ``ML_VERBOSE``
                        Includes a once-per-iteration report on simplex/ipm
                        progress and information about each nonzero row and
                        column.

                    - ``2``: ``ML_DETAILED``
                        Includes technical information about progress and
                        events in applying the simplex method.

                    - ``4``: ``ML_MINIMAL``
                        Once-per-solve information about progress as well as a
                        once-per-basis-matrix-reinversion report on progress in
                        simplex or a once-per-iteration report on progress in IPX.

                ``message_level`` behaves like a bitmask, i.e., any
                combination of levels is possible using the bit-or
                operator.

            - mps_parser_type_free : bool
                Use free format MPS parsing.

            - parallel : bool
                Run the solver in serial (False) or parallel (True).

            - presolve : bool
                Run the presolve or not (or if ``None``, then choose).

            - primal_feasibility_tolerance : double
                Primal feasibility tolerance.
                ``min(dual_feasibility_tolerance,
                primal_feasibility_tolerance)`` will be used for
                ipm feasibility tolerance.

            - run_as_hsol : bool
                Undocumented advanced option.

            - run_crossover : bool
                Advanced option. Toggles running the crossover routine
                for IPX.

            - sense : int {1, -1}
                ``sense=1`` corresponds to the MIN problem, ``sense=-1``
                corresponds to the MAX problem. TODO: NOT IMPLEMENTED

            - simplex_crash_strategy : int {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
                Strategy for simplex crash: off / LTSSF / Bixby (0/1/2).
                Default is ``0``.  Corresponds to the following:

                    - ``0``: ``SIMPLEX_CRASH_STRATEGY_OFF``
                    - ``1``: ``SIMPLEX_CRASH_STRATEGY_LTSSF_K``
                    - ``2``: ``SIMPLEX_CRASH_STRATEGY_BIXBY``
                    - ``3``: ``SIMPLEX_CRASH_STRATEGY_LTSSF_PRI``
                    - ``4``: ``SIMPLEX_CRASH_STRATEGY_LTSF_K``
                    - ``5``: ``SIMPLEX_CRASH_STRATEGY_LTSF_PRI``
                    - ``6``: ``SIMPLEX_CRASH_STRATEGY_LTSF``
                    - ``7``: ``SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS``
                    - ``8``: ``SIMPLEX_CRASH_STRATEGY_BASIC``
                    - ``9``: ``SIMPLE_CRASH_STRATEGY_TEST_SING``

            - simplex_dualise_strategy : int
                Undocumented advanced option.

            - simplex_dual_edge_weight_strategy : int {0, 1, 2, 3, 4}
                Strategy for simplex dual edge weights:
                Dantzig / Devex / Steepest Edge. Corresponds
                to the following:

                    - ``0``: ``SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG``
                    - ``1``: ``SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX``
                    - ``2``: ``SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH``
                    - ``3``: ``SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE``
                    - ``4``: ``SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL``

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

                    - ``0``: ``SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG``
                    - ``1``: ``SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX``

            - simplex_scale_strategy : int {0, 1, 2, 3, 4, 5}
                Strategy for scaling before simplex solver:
                off / on (0/1)

                    - ``0``:  ``SIMPLEX_SCALE_STRATEGY_OFF``
                    - ``1``: ``SIMPLEX_SCALE_STRATEGY_HIGHS``
                    - ``2``: ``SIMPLEX_SCALE_STRATEGY_HIGHS_FORCED``
                    - ``3``: ``SIMPLEX_SCALE_STRATEGY_HIGHS_015``
                    - ``4``: ``SIMPLEX_SCALE_STRATEGY_HIGHS_0157``
                    - ``5``: ``SIMPLEX_SCALE_STRATEGY_HSOL``

            - simplex_strategy : int {0, 1, 2, 3, 4}
                Strategy for simplex solver. Default: 1. Corresponds
                to the following:

                    - ``0``: ``SIMPLEX_STRATEGY_MIN``
                    - ``1``: ``SIMPLEX_STRATEGY_DUAL``
                    - ``2``: ``SIMPLEX_STRATEGY_DUAL_TASKS``
                    - ``3``: ``SIMPLEX_STRATEGY_DUAL_MULTI``
                    - ``4``: ``SIMPLEX_STRATEGY_PRIMAL``

            - simplex_update_limit : int
                Limit on the number of simplex UPDATE operations.

            - small_matrix_value : double
                Lower limit on abs(matrix entries): values smaller
                than this will be treated as zero.

            - solution_file : str
                Solution file

            - solver : str {'simplex', 'ipm'}
                Choose which solver to use.  If ``solver='simplex'``
                and ``parallel=True`` then PAMI will be used.

            - start_crossover_tolerance : double
                Tolerance to be satisfied before IPM crossover will
                start.

            - time_limit : double
                Max number of seconds to run the solver for.

            - use_original_HFactor_logic : bool
                Undocumented advanced option.

            - write_solution_to_file : bool
                Write the primal and dual solution to a file

            - write_solution_pretty : bool
                Write the primal and dual solution in a pretty
                (human-readable) format

        See [2]_ for a list of all non-advanced options.

    Returns
    -------
    res : dict

        If model_status is one of OPTIMAL,
        REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND, REACHED_TIME_LIMIT,
        REACHED_ITERATION_LIMIT:

            - ``status`` : int
                Model status code.

            - ``message`` : str
                Message corresponding to model status code.

            - ``x`` : list
                Solution variables.

            - ``slack`` : list
                Slack variables.

            - ``lambda`` : list
                Lagrange multipliers assoicated with the constraints
                Ax = b.

            - ``s`` : list
                Lagrange multipliers associated with the constraints
                x >= 0.

            - ``fun``
                Final objective value.

            - ``simplex_nit`` : int
                Number of iterations accomplished by the simplex
                solver.

            - ``ipm_nit`` : int
                Number of iterations accomplished by the interior-
                point solver.

        If model_status is not one of the above:

            - ``status`` : int
                Model status code.

            - ``message`` : str
                Message corresponding to model status code.

    Notes
    -----
    If ``options['write_solution_to_file']`` is ``True`` but
    ``options['solution_file']`` is unset or ``''``, then the solution
    will be printed to ``stdout``.

    If any iteration limit is reached, no solution will be
    available.

    ``OptimizeWarning`` will be raised if any option value set by
    the user is found to be incorrect.

    References
    ----------
    .. [1] https://highs.dev/
    .. [2] https://www.maths.ed.ac.uk/hall/HiGHS/HighsOptions.html
    '''

    cdef int numcol = c.size
    cdef int numrow = rhs.size
    cdef int numnz = avalue.size
    cdef int numintegrality = integrality.size

    # Fill up a HighsLp object
    cdef HighsLp lp
    lp.num_col_ = numcol
    lp.num_row_ = numrow
    lp.a_matrix_.num_col_ = numcol
    lp.a_matrix_.num_row_ = numrow
    lp.a_matrix_.format_ = MatrixFormatkColwise

    lp.col_cost_.resize(numcol)
    lp.col_lower_.resize(numcol)
    lp.col_upper_.resize(numcol)

    lp.row_lower_.resize(numrow)
    lp.row_upper_.resize(numrow)
    lp.a_matrix_.start_.resize(numcol + 1)
    lp.a_matrix_.index_.resize(numnz)
    lp.a_matrix_.value_.resize(numnz)

    # only need to set integrality if it's not's empty
    cdef HighsVarType * integrality_ptr = NULL
    if numintegrality > 0:
        lp.integrality_.resize(numintegrality)
        integrality_ptr = reinterpret_cast[HighsVarType_ptr](&integrality[0])
        lp.integrality_.assign(integrality_ptr, integrality_ptr + numcol)

    # Explicitly create pointers to pass to HiGHS C++ API;
    # do checking to make sure null memory-views are not
    # accessed (e.g., &lhs[0] raises exception when lhs is
    # empty!)
    cdef:
        double * colcost_ptr = NULL
        double * collower_ptr = NULL
        double * colupper_ptr = NULL
        double * rowlower_ptr = NULL
        double * rowupper_ptr = NULL
        int * astart_ptr = NULL
        int * aindex_ptr = NULL
        double * avalue_ptr = NULL
    if numrow > 0:
        rowlower_ptr = &lhs[0]
        rowupper_ptr = &rhs[0]
        lp.row_lower_.assign(rowlower_ptr, rowlower_ptr + numrow)
        lp.row_upper_.assign(rowupper_ptr, rowupper_ptr + numrow)
    else:
        lp.row_lower_.empty()
        lp.row_upper_.empty()
    if numcol > 0:
        colcost_ptr = &c[0]
        collower_ptr = &lb[0]
        colupper_ptr = &ub[0]
        lp.col_cost_.assign(colcost_ptr, colcost_ptr + numcol)
        lp.col_lower_.assign(collower_ptr, collower_ptr + numcol)
        lp.col_upper_.assign(colupper_ptr, colupper_ptr + numcol)
    else:
        lp.col_cost_.empty()
        lp.col_lower_.empty()
        lp.col_upper_.empty()
        lp.integrality_.empty()
    if numnz > 0:
        astart_ptr = &astart[0]
        aindex_ptr = &aindex[0]
        avalue_ptr = &avalue[0]
        lp.a_matrix_.start_.assign(astart_ptr, astart_ptr + numcol + 1)
        lp.a_matrix_.index_.assign(aindex_ptr, aindex_ptr + numnz)
        lp.a_matrix_.value_.assign(avalue_ptr, avalue_ptr + numnz)
    else:
        lp.a_matrix_.start_.empty()
        lp.a_matrix_.index_.empty()
        lp.a_matrix_.value_.empty()

    # Create the options
    cdef Highs highs
    apply_options(options, highs)

    # Make a Highs object and pass it everything
    cdef HighsModelStatus err_model_status = HighsModelStatusNOTSET
    cdef HighsStatus init_status = highs.passModel(lp)
    if init_status != HighsStatusOK:
        if init_status != HighsStatusWarning:
            err_model_status = HighsModelStatusMODEL_ERROR
            return {
                'status': <int> err_model_status,
                'message': highs.modelStatusToString(err_model_status).decode(),
            }

    # Solve the LP
    cdef HighsStatus run_status = highs.run()
    if run_status == HighsStatusError:
        return {
            'status': <int> highs.getModelStatus(),
            'message': highsStatusToString(run_status).decode(),
        }

    # Extract what we need from the solution
    cdef HighsModelStatus model_status = highs.getModelStatus()

    # We might need an info object if we can look up the solution and a place to put solution
    cdef HighsInfo info = highs.getHighsInfo() # it should always be safe to get the info object
    cdef HighsSolution solution
    cdef HighsBasis basis
    cdef double[:, ::1] marg_bnds = np.zeros((2, numcol))  # marg_bnds[0, :]: lower

    # Failure modes:
    #     LP: if we have anything other than an Optimal status, it
    #         is unsafe (and unhelpful) to read any results
    #    MIP: has a non-Optimal status or has timed out/reached max iterations
    #             1) If not Optimal/TimedOut/MaxIter status, there is no solution
    #             2) If TimedOut/MaxIter status, there may be a feasible solution.
    #                if the objective function value is not Infinity, then the
    #                current solution is feasible and can be returned.  Else, there
    #                is no solution.
    mipFailCondition = model_status not in {
        HighsModelStatusOPTIMAL,
        HighsModelStatusREACHED_TIME_LIMIT,
        HighsModelStatusREACHED_ITERATION_LIMIT,
    } or (model_status in {
        HighsModelStatusREACHED_TIME_LIMIT,
        HighsModelStatusREACHED_ITERATION_LIMIT,
    } and (info.objective_function_value == HIGHS_CONST_INF))
    lpFailCondition = model_status != HighsModelStatusOPTIMAL
    if (highs.getLp().isMip() and mipFailCondition) or (not highs.getLp().isMip() and lpFailCondition):
        return {
            'status': <int> model_status,
            'message': f'model_status is {highs.modelStatusToString(model_status).decode()}; '
                       f'primal_status is {utilBasisStatusToString(<HighsBasisStatus> info.primal_solution_status).decode()}',
            'simplex_nit': info.simplex_iteration_count,
            'ipm_nit': info.ipm_iteration_count,
            'fun': None,
            'crossover_nit': info.crossover_iteration_count,
        }
    # If the model status is such that the solution can be read
    else:
        # Should be safe to read the solution:
        solution = highs.getSolution()
        basis = highs.getBasis()

        # lagrangians for bounds based on column statuses
        for ii in range(numcol):
            if HighsBasisStatusLOWER == basis.col_status[ii]:
                marg_bnds[0, ii] = solution.col_dual[ii]
            elif HighsBasisStatusUPPER == basis.col_status[ii]:
                marg_bnds[1, ii] = solution.col_dual[ii]

        res = {
            'status': <int> model_status,
            'message': highs.modelStatusToString(model_status).decode(),

            # Primal solution
            'x': [solution.col_value[ii] for ii in range(numcol)],

            # Ax + s = b => Ax = b - s
            # Note: this is for all constraints (A_ub and A_eq)
            'slack': [rhs[ii] - solution.row_value[ii] for ii in range(numrow)],

            # lambda are the lagrange multipliers associated with Ax=b
            'lambda': [solution.row_dual[ii] for ii in range(numrow)],
            'marg_bnds': marg_bnds,

            'fun': info.objective_function_value,
            'simplex_nit': info.simplex_iteration_count,
            'ipm_nit': info.ipm_iteration_count,
            'crossover_nit': info.crossover_iteration_count,
        }

        if highs.getLp().isMip():
            res.update({
                'mip_node_count': info.mip_node_count,
                'mip_dual_bound': info.mip_dual_bound,
                'mip_gap': info.mip_gap,
            })

        return res
