# distutils: language=c++
# cython: language_level=3

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
    HighsModelStatus,
    HighsModelStatusOPTIMAL)
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
            if scaled_model_status == HighsModelStatusOPTIMAL:
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
        if model_status == HighsModelStatusOPTIMAL:
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

def linprog_mps(model_file, presolve=None, solver=None, bool run_quiet=False):
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
