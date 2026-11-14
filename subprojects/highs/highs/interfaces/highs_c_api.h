/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_C_API
#define HIGHS_C_API
//
// Welcome to the HiGHS C API!
//
// The simplest way to use HiGHS to solve an LP, MIP or QP from C is
// to pass the problem data to the appropriate method Highs_lpCall,
// Highs_mipCall or Highs_qpCall, and these methods return the
// appropriate solution information
//
// For sophisticated applications, where esoteric solution
// information is needed, or if a sequence of modified models need to
// be solved, use the Highs_create method to generate a pointer to an
// instance of the C++ Highs class, and then use any of a large number
// of models for which this pointer is the first parameter.
//
#include "lp_data/HighsCallbackStruct.h"

// The maximum string length. This corresponds to the `kHighsMaxStringLength`
// const in HConst.h. Note the difference of Max and Maximum.
static const HighsInt kHighsMaximumStringLength = 512;

// Error codes for the return value of functions in the C API. These do not
// appear in the C++ API.
static const HighsInt kHighsStatusError = -1;
static const HighsInt kHighsStatusOk = 0;
static const HighsInt kHighsStatusWarning = 1;

// The values of the `HighsVarType` enum in HConst.h.
static const HighsInt kHighsVarTypeContinuous = 0;
static const HighsInt kHighsVarTypeInteger = 1;
static const HighsInt kHighsVarTypeSemiContinuous = 2;
static const HighsInt kHighsVarTypeSemiInteger = 3;
static const HighsInt kHighsVarTypeImplicitInteger = 4;

// The values of the `HighsOptionType` enum in HConst.h.
static const HighsInt kHighsOptionTypeBool = 0;
static const HighsInt kHighsOptionTypeInt = 1;
static const HighsInt kHighsOptionTypeDouble = 2;
static const HighsInt kHighsOptionTypeString = 3;

// The values of the `HighsInfoType` enum in HConst.h.
static const HighsInt kHighsInfoTypeInt64 = -1;
static const HighsInt kHighsInfoTypeInt = 1;
static const HighsInt kHighsInfoTypeDouble = 2;

// The values of the `ObjSense` enum in HConst.h.
static const HighsInt kHighsObjSenseMinimize = 1;
static const HighsInt kHighsObjSenseMaximize = -1;

// The values of the `MatrixFormat` enum in HConst.h.
static const HighsInt kHighsMatrixFormatColwise = 1;
static const HighsInt kHighsMatrixFormatRowwise = 2;

// The values of the `HessianFormat` enum in HConst.h.
static const HighsInt kHighsHessianFormatTriangular = 1;
static const HighsInt kHighsHessianFormatSquare = 2;

// The values of the `SolutionStatus` enum in HConst.h.
static const HighsInt kHighsSolutionStatusNone = 0;
static const HighsInt kHighsSolutionStatusInfeasible = 1;
static const HighsInt kHighsSolutionStatusFeasible = 2;

// The values of the `BasisValidity` enum in HConst.h
static const HighsInt kHighsBasisValidityInvalid = 0;
static const HighsInt kHighsBasisValidityValid = 1;

// The values of the `HighsPresolveStatus` enum in HConst.h.
static const HighsInt kHighsPresolveStatusNotPresolved = -1;
static const HighsInt kHighsPresolveStatusNotReduced = 0;
static const HighsInt kHighsPresolveStatusInfeasible = 1;
static const HighsInt kHighsPresolveStatusUnboundedOrInfeasible = 2;
static const HighsInt kHighsPresolveStatusReduced = 3;
static const HighsInt kHighsPresolveStatusReducedToEmpty = 4;
static const HighsInt kHighsPresolveStatusTimeout = 5;
static const HighsInt kHighsPresolveStatusNullError = 6;
static const HighsInt kHighsPresolveStatusOptionsError = 7;
static const HighsInt kHighsPresolveStatusOutOfMemory = 8;
static const HighsInt kHighsPresolveStatusException = 9;

// The values of the `HighsModelStatus` enum in HConst.h.
static const HighsInt kHighsModelStatusNotset = 0;
static const HighsInt kHighsModelStatusLoadError = 1;
static const HighsInt kHighsModelStatusModelError = 2;
static const HighsInt kHighsModelStatusPresolveError = 3;
static const HighsInt kHighsModelStatusSolveError = 4;
static const HighsInt kHighsModelStatusPostsolveError = 5;
static const HighsInt kHighsModelStatusModelEmpty = 6;
static const HighsInt kHighsModelStatusOptimal = 7;
static const HighsInt kHighsModelStatusInfeasible = 8;
static const HighsInt kHighsModelStatusUnboundedOrInfeasible = 9;
static const HighsInt kHighsModelStatusUnbounded = 10;
static const HighsInt kHighsModelStatusObjectiveBound = 11;
static const HighsInt kHighsModelStatusObjectiveTarget = 12;
static const HighsInt kHighsModelStatusTimeLimit = 13;
static const HighsInt kHighsModelStatusIterationLimit = 14;
static const HighsInt kHighsModelStatusUnknown = 15;
static const HighsInt kHighsModelStatusSolutionLimit = 16;
static const HighsInt kHighsModelStatusInterrupt = 17;
static const HighsInt kHighsModelStatusMemoryLimit = 18;
static const HighsInt kHighsModelStatusHighsInterrupt = 19;

// The values of the `HighsBasisStatus` enum in HConst.h.
static const HighsInt kHighsBasisStatusLower = 0;
static const HighsInt kHighsBasisStatusBasic = 1;
static const HighsInt kHighsBasisStatusUpper = 2;
static const HighsInt kHighsBasisStatusZero = 3;
static const HighsInt kHighsBasisStatusNonbasic = 4;

// The values of the `HighsCallbackType` enum in HConst.h.
static const HighsInt kHighsCallbackLogging = 0;
static const HighsInt kHighsCallbackSimplexInterrupt = 1;
static const HighsInt kHighsCallbackIpmInterrupt = 2;
static const HighsInt kHighsCallbackMipSolution = 3;
static const HighsInt kHighsCallbackMipImprovingSolution = 4;
static const HighsInt kHighsCallbackMipLogging = 5;
static const HighsInt kHighsCallbackMipInterrupt = 6;
static const HighsInt kHighsCallbackMipGetCutPool = 7;
static const HighsInt kHighsCallbackMipDefineLazyConstraints = 8;
static const HighsInt kHighsCallbackCallbackMipUserSolution = 9;
static const HighsInt kHighsCallbackQpFirstFeasiblePoint = 10;
static const HighsInt kHighsCallbackQpInterrupt = 11;

// Constants for Highs_getCallbackDataOutItem. These do not appear in the C++
// API.
static const char* const kHighsCallbackDataOutLogTypeName = "log_type";
static const char* const kHighsCallbackDataOutRunningTimeName = "running_time";
static const char* const kHighsCallbackDataOutSimplexIterationCountName =
    "simplex_iteration_count";
static const char* const kHighsCallbackDataOutIpmIterationCountName =
    "ipm_iteration_count";
static const char* const kHighsCallbackDataOutPdlpIterationCountName =
    "pdlp_iteration_count";
static const char* const kHighsCallbackDataOutQpasmIterationCountName =
    "qpasm_iteration_count";
static const char* const kHighsCallbackDataOutObjectiveFunctionValueName =
    "objective_function_value";
static const char* const kHighsCallbackDataOutMipNodeCountName =
    "mip_node_count";
static const char* const kHighsCallbackDataOutMipTotalLpIterationsName =
    "mip_total_lp_iterations";
static const char* const kHighsCallbackDataOutMipPrimalBoundName =
    "mip_primal_bound";
static const char* const kHighsCallbackDataOutMipDualBoundName =
    "mip_dual_bound";
static const char* const kHighsCallbackDataOutMipGapName = "mip_gap";
static const char* const kHighsCallbackDataOutMipSolutionName = "mip_solution";
static const char* const kHighsCallbackDataOutCutpoolNumColName =
    "cutpool_num_col";
static const char* const kHighsCallbackDataOutCutpoolNumCutName =
    "cutpool_num_cut";
static const char* const kHighsCallbackDataOutCutpoolNumNzName =
    "cutpool_num_nz";
static const char* const kHighsCallbackDataOutCutpoolStartName =
    "cutpool_start";
static const char* const kHighsCallbackDataOutCutpoolIndexName =
    "cutpool_index";
static const char* const kHighsCallbackDataOutCutpoolValueName =
    "cutpool_value";
static const char* const kHighsCallbackDataOutCutpoolLowerName =
    "cutpool_lower";
static const char* const kHighsCallbackDataOutCutpoolUpperName =
    "cutpool_upper";
static const char* const kHighsCallbackDataOutQpSolutionName = "qp_solution";

// Modified values of the `IisStrategy` enum in HConst.h.
static const HighsInt kHighsIisStrategyLight = 0;
// Forces full IIS calculation as before - ie with the
// kIisStrategyIrreducible = 4 bit set, as well as the
// kIisStrategyFromLp = 2 bit set, and possibly the
// kIisStrategyColPriority = 8 bit set
static const HighsInt kHighsIisStrategyFromLpRowPriority = 6;
static const HighsInt kHighsIisStrategyFromLpColPriority = 14;

// Variants of the `IisBoundStatus` enum in HighsIis.h. The `Status` part has
// been droped from the names.
static const HighsInt kHighsIisBoundFree = 1;
static const HighsInt kHighsIisBoundLower = 2;
static const HighsInt kHighsIisBoundUpper = 3;
static const HighsInt kHighsIisBoundBoxed = 4;

// The values of the `IisStatus` enum in HConst.h.
static const HighsInt kHighsIisStatusNotInConflict = -1;
static const HighsInt kHighsIisStatusMaybeInConflict = 0;
static const HighsInt kHighsIisStatusInConflict = 1;

// The values of the `HessianOracleCallType` enum in HConst.H
static const HighsInt kHighsHessianOracleCallTypeMin = 0;
static const HighsInt kHighsHessianOracleCallTypeEntry = 0;
static const HighsInt kHighsHessianOracleCallTypeColumn = 1;
static const HighsInt kHighsHessianOracleCallTypeProduct = 2;
static const HighsInt kHighsHessianOracleCallTypeMax = 2;

// The values of the `PresolveRuleType` enum in HConst.h.
static const HighsInt kHighsPresolveRuleEmptyRow = 0;
static const HighsInt kHighsPresolveRuleSingletonRow = 1;
static const HighsInt kHighsPresolveRuleRedundantRow = 2;
static const HighsInt kHighsPresolveRuleEmptyCol = 3;
static const HighsInt kHighsPresolveRuleFixedCol = 4;
static const HighsInt kHighsPresolveRuleDominatedCol = 5;
static const HighsInt kHighsPresolveRuleForcingRow = 6;
static const HighsInt kHighsPresolveRuleForcingCol = 7;
static const HighsInt kHighsPresolveRuleFreeColSubstitution = 8;
static const HighsInt kHighsPresolveRuleDoubletonEquation = 9;
static const HighsInt kHighsPresolveRuleDependentEquations = 10;
static const HighsInt kHighsPresolveRuleDependentFreeCols = 11;
static const HighsInt kHighsPresolveRuleAggregator = 12;
static const HighsInt kHighsPresolveRuleParallelRowsAndCols = 13;
static const HighsInt kHighsPresolveRuleSparsify = 14;
static const HighsInt kHighsPresolveRuleProbing = 15;
static const HighsInt kHighsPresolveRuleEnumeration = 16;
static const HighsInt kHighsPresolveRuleDualFixing = 17;
static const HighsInt kHighsPresolveRuleZeroCostSingleton = 18;
static const HighsInt kHighsPresolveRuleColStuffing = 19;
static const HighsInt kHighsPresolveRuleDualFixProbing = 20;
static const HighsInt kHighsPresolveRuleInitialSweep = 21;
static const HighsInt kHighsPresolveRuleFourierMotzkin = 22;
static const HighsInt kHighsPresolveRuleWeaklyDominatedCol = 23;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Formulate and solve a linear program using HiGHS.
 *
 * @param num_col   The number of columns.
 * @param num_row   The number of rows.
 * @param num_nz    The number of nonzeros in the constraint matrix.
 * @param a_format  The format of the constraint matrix as a
 *                  `kHighsMatrixFormat` constant.
 * @param sense     The optimization sense as a `kHighsObjSense` constant.
 * @param offset    The objective constant.
 * @param col_cost  An array of length [num_col] with the column costs.
 * @param col_lower An array of length [num_col] with the column lower bounds.
 * @param col_upper An array of length [num_col] with the column upper bounds.
 * @param row_lower An array of length [num_row] with the row lower bounds.
 * @param row_upper An array of length [num_row] with the row upper bounds.
 * @param a_start   The constraint matrix is provided to HiGHS in compressed
 *                  sparse column form (if `a_format` is
 *                  `kHighsMatrixFormatColwise`, otherwise compressed sparse row
 *                  form). The sparse matrix consists of three arrays,
 *                  `a_start`, `a_index`, and `a_value`. `a_start` is an array
 *                  of length [num_col] containing the starting index of each
 *                  column in `a_index`. If `a_format` is
 *                  `kHighsMatrixFormatRowwise` the array is of length [num_row]
 *                  corresponding to each row.
 * @param a_index   An array of length [num_nz] with indices of matrix entries.
 * @param a_value   An array of length [num_nz] with values of matrix entries.
 *
 * @param col_value      An array of length [num_col], to be filled with the
 *                       primal column solution.
 * @param col_dual       An array of length [num_col], to be filled with the
 *                       dual column solution.
 * @param row_value      An array of length [num_row], to be filled with the
 *                       primal row solution.
 * @param row_dual       An array of length [num_row], to be filled with the
 *                       dual row solution.
 * @param col_basis_status  An array of length [num_col], to be filled with the
 *                          basis status of the columns in the form of a
 *                          `kHighsBasisStatus` constant.
 * @param row_basis_status  An array of length [num_row], to be filled with the
 *                          basis status of the rows in the form of a
 *                          `kHighsBasisStatus` constant.
 * @param model_status      The location in which to place the termination
 *                          status of the model after the solve in the form of a
 *                          `kHighsModelStatus` constant.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_lpCall(const HighsInt num_col, const HighsInt num_row,
                      const HighsInt num_nz, const HighsInt a_format,
                      const HighsInt sense, const double offset,
                      const double* col_cost, const double* col_lower,
                      const double* col_upper, const double* row_lower,
                      const double* row_upper, const HighsInt* a_start,
                      const HighsInt* a_index, const double* a_value,
                      double* col_value, double* col_dual, double* row_value,
                      double* row_dual, HighsInt* col_basis_status,
                      HighsInt* row_basis_status, HighsInt* model_status);

/**
 * Formulate and solve a mixed-integer linear program using HiGHS.
 *
 * The signature of this method is identical to `Highs_lpCall`, except that it
 * has an additional `integrality` argument, and that it is missing the
 * `col_dual`, `row_dual`, `col_basis_status` and `row_basis_status` arguments.
 *
 * @param integrality   An array of length [num_col], containing a
 *                      `kHighsVarType` constant for each column.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_mipCall(const HighsInt num_col, const HighsInt num_row,
                       const HighsInt num_nz, const HighsInt a_format,
                       const HighsInt sense, const double offset,
                       const double* col_cost, const double* col_lower,
                       const double* col_upper, const double* row_lower,
                       const double* row_upper, const HighsInt* a_start,
                       const HighsInt* a_index, const double* a_value,
                       const HighsInt* integrality, double* col_value,
                       double* row_value, HighsInt* model_status);

/**
 * Formulate and solve a quadratic program using HiGHS.
 *
 * The signature of this method is identical to `Highs_lpCall`, except that it
 * has additional arguments for specifying the Hessian matrix.
 *
 * @param q_num_nz  The number of nonzeros in the Hessian matrix.
 * @param q_format  The format of the Hessian matrix in the form of a
 *                  `kHighsHessianStatus` constant. If q_num_nz > 0, this must
 *                  be `kHighsHessianFormatTriangular`.
 * @param q_start   The Hessian matrix is provided to HiGHS as the lower
 *                  triangular component in compressed sparse column form
 *                  (or, equivalently, as the upper triangular component
 *                  in compressed sparse row form). The sparse matrix consists
 *                  of three arrays, `q_start`, `q_index`, and `q_value`.
 *                  `q_start` is an array of length [num_col].
 * @param q_index   An array of length [q_num_nz] with indices of matrix
 *                  entries.
 * @param q_value   An array of length [q_num_nz] with values of matrix entries.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_qpCall(
    const HighsInt num_col, const HighsInt num_row, const HighsInt num_nz,
    const HighsInt q_num_nz, const HighsInt a_format, const HighsInt q_format,
    const HighsInt sense, const double offset, const double* col_cost,
    const double* col_lower, const double* col_upper, const double* row_lower,
    const double* row_upper, const HighsInt* a_start, const HighsInt* a_index,
    const double* a_value, const HighsInt* q_start, const HighsInt* q_index,
    const double* q_value, double* col_value, double* col_dual,
    double* row_value, double* row_dual, HighsInt* col_basis_status,
    HighsInt* row_basis_status, HighsInt* model_status);

/**
 * Create a Highs instance and return the reference.
 *
 * Call `Highs_destroy` on the returned reference to clean up memory
 * allocated for this instance.
 *
 * @returns A pointer to the Highs instance.
 */
void* Highs_create(void);

/**
 * Destroy the model `highs` created by `Highs_create` and free all
 * memory allocated for this instance. Future calls using `highs` are
 * not allowed.
 *
 * Since the scheduler is owned by the thread and not by the Highs
 * instances, its memory cannot be freed by `Highs_destroy`. Hence, to free
 * all memory used by HiGHS on this thread, `Highs_resetThreadScheduler` must
 * also be called.
 *
 * To empty a model without invalidating `highs`, see `Highs_clearModel`.
 *
 * @param highs     A pointer to the Highs instance.
 */
void Highs_destroy(void* highs);

/**
 * Return the HiGHS version number as a string of the form "vX.Y.Z".
 *
 * @returns The HiGHS version as a `char*`.
 */
const char* Highs_version(void);

/**
 * Return the HiGHS major version number.
 *
 * @returns The HiGHS major version number.
 */
HighsInt Highs_versionMajor(void);

/**
 * Return the HiGHS minor version number.
 *
 * @returns The HiGHS minor version number.
 */
HighsInt Highs_versionMinor(void);

/**
 * Return the HiGHS patch version number.
 *
 * @returns The HiGHS patch version number.
 */
HighsInt Highs_versionPatch(void);

/**
 * Return the HiGHS githash.
 *
 * @returns The HiGHS githash.
 */
const char* Highs_githash(void);

/**
 * Read a model from `filename` into `highs`.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The filename to read.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_readModel(void* highs, const char* filename);

/**
 * Write the model in `highs` to `filename`.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The filename to write.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_writeModel(void* highs, const char* filename);

/**
 * Write the presolved model in `highs` to `filename`.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The filename to write.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_writePresolvedModel(void* highs, const char* filename);

/**
 * Reset the options and then call `clearModel`.
 *
 * See `Highs_destroy` to free all associated memory.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_clear(void* highs);

/**
 * Remove all variables and constraints from the model `highs`, but do not
 * invalidate the pointer `highs`. Future calls (for example, adding new
 * variables and constraints) are allowed.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_clearModel(void* highs);

/**
 * Clear all solution data associated with the model.
 *
 * See `Highs_destroy` to clear the model and free all associated memory.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_clearSolver(void* highs);

/**
 * Release all retained memory back to the allocator.
 *
 * Clears all solver state and shrinks internal vectors to free unused
 * capacity. Useful in long-running services that reuse a Highs instance
 * across multiple solves to prevent unbounded RSS growth from heap
 * fragmentation.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_releaseMemory(void* highs);

/**
 * Presolve a model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_presolve(void* highs);

/**
 * Optimize a model. The algorithm used by HiGHS depends on the options that
 * have been set.
 *
 * If this function returns `kHighsStatusError`, the solve failed. If this
 * function returns `kHighsStatusWarning`, something unusual happened, but the
 * model may have been solved correctly. In each case, look at the log for an
 * explanation, and use information such as the model and solution status, and
 * the number of primal and dual infeasibilities to understand the quality of
 * the solution (if one is available).
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_run(void* highs);

/**
 * Postsolve a model using a primal (and possibly dual) solution. The
 * postsolved solution can be retrieved later by calling
 * `Highs_getSolution`.
 *
 * @param highs       A pointer to the Highs instance.
 * @param col_value   An array of length [num_col] with the column solution
 *                    values.
 * @param col_dual    An array of length [num_col] with the column dual
 *                    values, or a null pointer if not known.
 * @param row_dual    An array of length [num_row] with the row dual values,
 *                    or a null pointer if not known.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_postsolve(void* highs, const double* col_value,
                         const double* col_dual, const double* row_dual);

/**
 * Write the solution information (including dual and basis status, if
 * available) to a file.
 *
 * See also: `Highs_writeSolutionPretty`.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The name of the file to write the results to.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_writeSolution(const void* highs, const char* filename);

/**
 * Write the solution information (including dual and basis status, if
 * available) to a file in a human-readable format.
 *
 * The method identical to `Highs_writeSolution`, except that the
 * printout is in a human-readable format.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The name of the file to write the results to.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_writeSolutionPretty(const void* highs, const char* filename);

/**
 * Pass a linear program (LP) to HiGHS in a single function call.
 *
 * The signature of this function is identical to `Highs_passModel`, without the
 * arguments for passing the Hessian matrix of a quadratic program and the
 * integrality vector.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passLp(void* highs, const HighsInt num_col,
                      const HighsInt num_row, const HighsInt num_nz,
                      const HighsInt a_format, const HighsInt sense,
                      const double offset, const double* col_cost,
                      const double* col_lower, const double* col_upper,
                      const double* row_lower, const double* row_upper,
                      const HighsInt* a_start, const HighsInt* a_index,
                      const double* a_value);

/**
 * Pass a mixed-integer linear program (MILP) to HiGHS in a single function
 * call.
 *
 * The signature of function is identical to `Highs_passModel`, without the
 * arguments for passing the Hessian matrix of a quadratic program.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passMip(void* highs, const HighsInt num_col,
                       const HighsInt num_row, const HighsInt num_nz,
                       const HighsInt a_format, const HighsInt sense,
                       const double offset, const double* col_cost,
                       const double* col_lower, const double* col_upper,
                       const double* row_lower, const double* row_upper,
                       const HighsInt* a_start, const HighsInt* a_index,
                       const double* a_value, const HighsInt* integrality);

/**
 * Pass a model to HiGHS in a single function call. This is faster than
 * constructing the model using `Highs_addRow` and `Highs_addCol`.
 *
 * @param highs       A pointer to the Highs instance.
 * @param num_col     The number of columns.
 * @param num_row     The number of rows.
 * @param num_nz      The number of elements in the constraint matrix.
 * @param q_num_nz    The number of elements in the Hessian matrix.
 * @param a_format    The format of the constraint matrix to use in the form of
 *                    a `kHighsMatrixFormat` constant.
 * @param q_format    The format of the Hessian matrix to use in the form of a
 *                    `kHighsHessianFormat` constant.
 * @param sense       The optimization sense in the form of a `kHighsObjSense`
 *                    constant.
 * @param offset      The constant term in the objective function.
 * @param col_cost    An array of length [num_col] with the objective
 *                    coefficients.
 * @param col_lower   An array of length [num_col] with the lower column bounds.
 * @param col_upper   An array of length [num_col] with the upper column bounds.
 * @param row_lower   An array of length [num_row] with the upper row bounds.
 * @param row_upper   An array of length [num_row] with the upper row bounds.
 * @param a_start     The constraint matrix is provided to HiGHS in compressed
 *                    sparse column form (if `a_format` is
 *                    `kHighsMatrixFormatColwise`, otherwise compressed sparse
 *                    row form). The sparse matrix consists of three arrays,
 *                    `a_start`, `a_index`, and `a_value`. `a_start` is an array
 *                    of length [num_col] containing the starting index of each
 *                    column in `a_index`. If `a_format` is
 *                    `kHighsMatrixFormatRowwise` the array is of length
 *                    [num_row] corresponding to each row.
 * @param a_index     An array of length [num_nz] with indices of matrix
 *                    entries.
 * @param a_value     An array of length [num_nz] with values of matrix entries.
 * @param q_start     The Hessian matrix is provided to HiGHS as the lower
 *                    triangular component in compressed sparse column form
 *                    (or, equivalently, as the upper triangular component
 *                    in compressed sparse row form). The sparse matrix consists
 *                    of three arrays, `q_start`, `q_index`, and `q_value`.
 *                    `q_start` is an array of length [num_col]. If the model
 *                    is linear, pass NULL.
 * @param q_index     An array of length [q_num_nz] with indices of matrix
 *                    entries. If the model is linear, pass NULL.
 * @param q_value     An array of length [q_num_nz] with values of matrix
 *                     entries. If the model is linear, pass NULL.
 * @param integrality An array of length [num_col] containing a `kHighsVarType`
 *                    constant for each column.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passModel(void* highs, const HighsInt num_col,
                         const HighsInt num_row, const HighsInt num_nz,
                         const HighsInt q_num_nz, const HighsInt a_format,
                         const HighsInt q_format, const HighsInt sense,
                         const double offset, const double* col_cost,
                         const double* col_lower, const double* col_upper,
                         const double* row_lower, const double* row_upper,
                         const HighsInt* a_start, const HighsInt* a_index,
                         const double* a_value, const HighsInt* q_start,
                         const HighsInt* q_index, const double* q_value,
                         const HighsInt* integrality);

/**
 * Set the Hessian matrix for a quadratic objective.
 *
 * @param highs     A pointer to the Highs instance.
 * @param dim       The dimension of the Hessian matrix. Should be [num_col].
 * @param num_nz    The number of non-zero elements in the Hessian matrix.
 * @param format    The format of the Hessian matrix as a `kHighsHessianFormat`
 *                  constant. This must be `kHighsHessianFormatTriangular`.
 * @param start     The Hessian matrix is provided to HiGHS as the lower
 *                  triangular component in compressed sparse column form
 *                  (or, equivalently, as the upper triangular component
 *                  in compressed sparse row form), using `q_start`, `q_index`,
 *                  and `q_value`.The Hessian matrix is provided to HiGHS as the
 *                  lower triangular component in compressed sparse column form.
 *                  The sparse matrix consists of three arrays, `start`,
 *                  `index`, and `value`. `start` is an array of length
 *                  [num_col] containing the starting index of each column in
 *                  `index`.
 * @param index     An array of length [num_nz] with indices of matrix entries.
 * @param value     An array of length [num_nz] with values of matrix entries.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passHessian(void* highs, const HighsInt dim,
                           const HighsInt num_nz, const HighsInt format,
                           const HighsInt* start, const HighsInt* index,
                           const double* value);

/**
 * Set the Hessian oracle callback method
 *
 * @param highs        A pointer to the Highs instance.
 * @param dim          Dimension of the Hessian
 * @param oracleCall   A pointer to the Hessian oracle callback
 * @param oracle_data  A pointer to the Hessian oracle callback data
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passHessianOracle(void* highs, const HighsInt dim,
                                 HighsCHessianFunctionType oracleCall,
                                 void* oracle_data);

/**
 * Check the Hessian oracle callback method
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_checkHessianOracle(void* highs);

/**
 * Passes multiple linear objective data to HiGHS, clearing any such
 * data already in HiGHS
 *
 * @param highs         A pointer to the Highs instance.
 * @param weight        A pointer to the weights of the linear objective, with
 *                      its positive/negative sign determining whether it is
 *                      minimized or maximized during lexicographic optimization
 * @param offset        A pointer to the objective offsets
 * @param coefficients  A pointer to the objective coefficients
 * @param abs_tolerance A pointer to the absolute tolerances used when
 *                      constructing objective constraints during lexicographic
 *                      optimization
 * @param rel_tolerance A pointer to the relative tolerances used when
 *                      constructing objective constraints during lexicographic
 *                      optimization
 * @param priority      A pointer to the priorities of the objectives during
 *                      lexicographic optimization
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */

HighsInt Highs_passLinearObjectives(const void* highs,
                                    const HighsInt num_linear_objective,
                                    const double* weight, const double* offset,
                                    const double* coefficients,
                                    const double* abs_tolerance,
                                    const double* rel_tolerance,
                                    const HighsInt* priority);

/**
 * Adds linear objective data to HiGHS
 *
 * @param highs         A pointer to the Highs instance.
 * @param weight        The weight of the linear objective, with its
 *                      positive/negative sign determining whether it is
 *                      minimized or maximized during lexicographic
 *                      optimization
 * @param offset        The objective offset
 * @param coefficients  A pointer to the objective coefficients
 * @param abs_tolerance The absolute tolerance used when constructing an
 *                      objective constraint during lexicographic optimization
 * @param rel_tolerance The relative tolerance used when constructing an
 *                      objective constraint during lexicographic optimization
 * @param priority      The priority of this objective during lexicographic
 *                      optimization
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */

HighsInt Highs_addLinearObjective(const void* highs, const double weight,
                                  const double offset,
                                  const double* coefficients,
                                  const double abs_tolerance,
                                  const double rel_tolerance,
                                  const HighsInt priority);

/**
 * Clears any multiple linear objective data in HiGHS
 *
 * @param highs A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */

HighsInt Highs_clearLinearObjectives(const void* highs);
/**
 * Pass the name of a row.
 *
 * @param highs A pointer to the Highs instance.
 * @param row   The row for which the name is supplied.
 * @param name  The name of the row.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passRowName(const void* highs, const HighsInt row,
                           const char* name);

/**
 * Pass the name of a column.
 *
 * @param highs A pointer to the Highs instance.
 * @param col   The column for which the name is supplied.
 * @param name  The name of the column.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passColName(const void* highs, const HighsInt col,
                           const char* name);

/**
 * Pass the name of the model.
 *
 * @param highs A pointer to the Highs instance.
 * @param name  The name of the model.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_passModelName(const void* highs, const char* name);

/**
 * Read the option values from file.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The filename from which to read the option values.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_readOptions(const void* highs, const char* filename);

/**
 * Set a boolean-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     The new value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setBoolOptionValue(void* highs, const char* option,
                                  const HighsInt value);

/**
 * Set an int-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     The new value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setIntOptionValue(void* highs, const char* option,
                                 const HighsInt value);

/**
 * Set a double-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     The new value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setDoubleOptionValue(void* highs, const char* option,
                                    const double value);

/**
 * Set a string-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     The new value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setStringOptionValue(void* highs, const char* option,
                                    const char* value);

/**
 * Get a boolean-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     The location in which the current value of the option should
 *                  be placed.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getBoolOptionValue(const void* highs, const char* option,
                                  HighsInt* value);

/**
 * Get an int-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     The location in which the current value of the option should
 *                  be placed.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getIntOptionValue(const void* highs, const char* option,
                                 HighsInt* value);

/**
 * Get a double-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     The location in which the current value of the option should
 *                  be placed.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getDoubleOptionValue(const void* highs, const char* option,
                                    double* value);

/**
 * Get a string-valued option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param value     A pointer to allocated memory (of at least
 *                  `kMaximumStringLength`) to store the current value of the
 *                  option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getStringOptionValue(const void* highs, const char* option,
                                    char* value);

/**
 * Get the type expected by an option.
 *
 * @param highs     A pointer to the Highs instance.
 * @param option    The name of the option.
 * @param type      A HighsInt in which the corresponding `kHighsOptionType`
 *                  constant should be placed.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getOptionType(const void* highs, const char* option,
                             HighsInt* type);

/**
 * Reset all options to their default value.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_resetOptions(void* highs);

/**
 * Write the current options to file.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The filename to write the options to.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_writeOptions(const void* highs, const char* filename);

/**
 * Write the value of non-default options to file.
 *
 * This is similar to `Highs_writeOptions`, except only options with
 * non-default value are written to `filename`.
 *
 * @param highs     A pointer to the Highs instance.
 * @param filename  The filename to write the options to.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_writeOptionsDeviations(const void* highs, const char* filename);

/**
 * Return the number of options
 *
 * @param highs     A pointer to the Highs instance.
 */
HighsInt Highs_getNumOptions(const void* highs);

/**
 * Get the name of an option identified by index
 *
 * @param highs     A pointer to the Highs instance.
 * @param index     The index of the option.
 * @param name      The name of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getOptionName(const void* highs, const HighsInt index,
                             char** name);

/**
 * Get the current and default values of a bool option
 *
 * @param highs         A pointer to the Highs instance.
 * @param current_value A pointer to the current value of the option.
 * @param default_value A pointer to the default value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getBoolOptionValues(const void* highs, const char* option,
                                   HighsInt* current_value,
                                   HighsInt* default_value);
/**
 * Get the current and default values of a HighsInt option
 *
 * @param highs         A pointer to the Highs instance.
 * @param current_value A pointer to the current value of the option.
 * @param min_value     A pointer to the minimum value of the option.
 * @param max_value     A pointer to the maximum value of the option.
 * @param default_value A pointer to the default value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getIntOptionValues(const void* highs, const char* option,
                                  HighsInt* current_value, HighsInt* min_value,
                                  HighsInt* max_value, HighsInt* default_value);

/**
 * Get the current and default values of a double option
 *
 * @param highs         A pointer to the Highs instance.
 * @param current_value A pointer to the current value of the option.
 * @param min_value     A pointer to the minimum value of the option.
 * @param max_value     A pointer to the maximum value of the option.
 * @param default_value A pointer to the default value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getDoubleOptionValues(const void* highs, const char* option,
                                     double* current_value, double* min_value,
                                     double* max_value, double* default_value);

/**
 * Get the current and default values of a string option
 *
 * @param highs         A pointer to the Highs instance.
 * @param current_value A pointer to the current value of the option.
 * @param default_value A pointer to the default value of the option.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getStringOptionValues(const void* highs, const char* option,
                                     char* current_value, char* default_value);

/**
 * Get an int-valued info value.
 *
 * @param highs     A pointer to the Highs instance.
 * @param info      The name of the info item.
 * @param value     A reference to an integer that the result will be stored in.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getIntInfoValue(const void* highs, const char* info,
                               HighsInt* value);

/**
 * Get a double-valued info value.
 *
 * @param highs     A pointer to the Highs instance.
 * @param info      The name of the info item.
 * @param value     A reference to a double that the result will be stored in.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getDoubleInfoValue(const void* highs, const char* info,
                                  double* value);

/**
 * Get an int64-valued info value.
 *
 * @param highs     A pointer to the Highs instance.
 * @param info      The name of the info item.
 * @param value     A reference to an int64 that the result will be stored in.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getInt64InfoValue(const void* highs, const char* info,
                                 int64_t* value);

/**
 * Get the type expected by an info item.
 *
 * @param highs     A pointer to the Highs instance.
 * @param info      The name of the info item.
 * @param type      A HighsInt in which the corresponding `kHighsOptionType`
 *                  constant is stored.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getInfoType(const void* highs, const char* info, HighsInt* type);

/**
 * Get the primal and dual solution from an optimized model.
 *
 * @param highs      A pointer to the Highs instance.
 * @param col_value  An array of length [num_col], to be filled with primal
 *                   column values.
 * @param col_dual   An array of length [num_col], to be filled with dual column
 *                   values.
 * @param row_value  An array of length [num_row], to be filled with primal row
 *                   values.
 * @param row_dual   An array of length [num_row], to be filled with dual row
 *                   values.
 *
 * @returns A `kHighsStatus` constant indicating whether a primal solution
 * exists.
 */
HighsInt Highs_getSolution(const void* highs, double* col_value,
                           double* col_dual, double* row_value,
                           double* row_dual);

/**
 * Given a linear program with a basic feasible solution, get the column and row
 * basis statuses.
 *
 * @param highs       A pointer to the Highs instance.
 * @param col_status  An array of length [num_col], to be filled with the column
 *                    basis statuses in the form of a `kHighsBasisStatus`
 *                    constant.
 * @param row_status  An array of length [num_row], to be filled with the row
 *                    basis statuses in the form of a `kHighsBasisStatus`
 *                    constant.
 *
 * @returns A `kHighsStatus` constant indicating whether the basis is valid.
 */
HighsInt Highs_getBasis(const void* highs, HighsInt* col_status,
                        HighsInt* row_status);

/**
 * Return the optimization status of the model in the form of a
 * `kHighsModelStatus` constant.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns An integer corresponding to the `kHighsModelStatus` constant
 */
HighsInt Highs_getModelStatus(const void* highs);

/**
 * Indicates whether a dual ray that is a certificate of primal
 * infeasibility currently exists, and (at the expense of solving an
 * LP) gets it if it does not and dual_ray_value is not nullptr.
 *
 * @param highs             A pointer to the Highs instance.
 * @param has_dual_ray      A pointer to a HighsInt to store 1 if a dual ray
 *                          currently exists.
 * @param dual_ray_value    An array of length [num_row] filled with the
 *                          unbounded ray.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getDualRay(const void* highs, HighsInt* has_dual_ray,
                          double* dual_ray_value);

/**
 * Indicates whether a dual unboundedness direction (corresponding to a
 * certificate of primal infeasibility) exists, and (at the expense of
 * solving an LP) gets it if it does not and
 * dual_unboundedness_direction is not nullptr
 *
 * @param highs                                   A pointer to the Highs
 *                                                instance.
 * @param has_dual_unboundedness_direction        A pointer to a HighsInt to
 *                                                store 1 if the dual
 *                                                unboundedness direction
 *                                                exists.
 * @param dual_unboundedness_direction_value      An array of length [num_col]
 *                                                filled with the unboundedness
 *                                                direction.
 */
HighsInt Highs_getDualUnboundednessDirection(
    const void* highs, HighsInt* has_dual_unboundedness_direction,
    double* dual_unboundedness_direction_value);

/**
 * Indicates whether a primal ray that is a certificate of primal
 * unboundedness currently exists, and (at the expense of solving an
 * LP) gets it if it does not and primal_ray_value is not nullptr.
 *
 * @param highs             A pointer to the Highs instance.
 * @param has_primal_ray    A pointer to a HighsInt to store 1 if the primal ray
 *                          exists.
 * @param primal_ray_value  An array of length [num_col] filled with the
 *                          unbounded ray.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getPrimalRay(const void* highs, HighsInt* has_primal_ray,
                            double* primal_ray_value);

/**
 * Get the primal objective function value.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The primal objective function value
 */
double Highs_getObjectiveValue(const void* highs);

/**
 * Get the indices of the rows and columns that make up the basis matrix ``B``
 * of a basic feasible solution.
 *
 * Non-negative entries are indices of columns, and negative entries are
 * `-row_index - 1`. For example, `{1, -1}` would be the second column and first
 * row.
 *
 * The order of these rows and columns is important for calls to the functions:
 *
 *  - `Highs_getBasisInverseRow`
 *  - `Highs_getBasisInverseCol`
 *  - `Highs_getBasisSolve`
 *  - `Highs_getBasisTransposeSolve`
 *  - `Highs_getReducedRow`
 *  - `Highs_getReducedColumn`
 *
 * @param highs             A pointer to the Highs instance.
 * @param basic_variables   An array of size [num_rows], filled with the indices
 *                          of the basic variables.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getBasicVariables(const void* highs, HighsInt* basic_variables);

/**
 * Get a row of the inverse basis matrix ``B^{-1}``.
 *
 * See `Highs_getBasicVariables` for a description of the ``B`` matrix.
 *
 * The arrays `row_vector` and `row_index` must have an allocated length of
 * [num_row]. However, check `row_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         A pointer to the Highs instance.
 * @param row           The index of the row to compute.
 * @param row_vector    An array of length [num_row] in which to store the
 *                      values of the non-zero elements.
 * @param row_num_nz    The number of non-zeros in the row.
 * @param row_index     An array of length [num_row] in which to store the
 *                      indices of the non-zero elements.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getBasisInverseRow(const void* highs, const HighsInt row,
                                  double* row_vector, HighsInt* row_num_nz,
                                  HighsInt* row_index);

/**
 * Get a column of the inverse basis matrix ``B^{-1}``.
 *
 * See `Highs_getBasicVariables` for a description of the ``B`` matrix.
 *
 * The arrays `col_vector` and `col_index` must have an allocated length of
 * [num_row]. However, check `col_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         A pointer to the Highs instance.
 * @param col           The index of the column to compute.
 * @param col_vector    An array of length [num_row] in which to store the
 *                      values of the non-zero elements.
 * @param col_num_nz    The number of non-zeros in the column.
 * @param col_index     An array of length [num_row] in which to store the
 *                      indices of the non-zero elements.

 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getBasisInverseCol(const void* highs, const HighsInt col,
                                  double* col_vector, HighsInt* col_num_nz,
                                  HighsInt* col_index);

/**
 * Compute ``\mathbf{x}=B^{-1}\mathbf{b}`` for a given vector
 * ``\mathbf{b}``.
 *
 * See `Highs_getBasicVariables` for a description of the ``B`` matrix.
 *
 * The arrays `solution_vector` and `solution_index` must have an allocated
 * length of [num_row]. However, check `solution_num_nz` to see how many
 * non-zero elements are actually stored.
 *
 * @param highs             A pointer to the Highs instance.
 * @param rhs               The right-hand side vector ``b``.
 * @param solution_vector   An array of length [num_row] in which to store the
 *                          values of the non-zero elements.
 * @param solution_num_nz   The number of non-zeros in the solution.
 * @param solution_index    An array of length [num_row] in which to store the
 *                          indices of the non-zero elements.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getBasisSolve(const void* highs, const double* rhs,
                             double* solution_vector, HighsInt* solution_num_nz,
                             HighsInt* solution_index);

/**
 * Compute ``\mathbf{x}=B^{-T}\mathbf{b}`` for a given vector
 * ``\mathbf{b}``.
 *
 * See `Highs_getBasicVariables` for a description of the ``B`` matrix.
 *
 * The arrays `solution_vector` and `solution_index` must have an allocated
 * length of [num_row]. However, check `solution_num_nz` to see how many
 * non-zero elements are actually stored.
 *
 * @param highs             A pointer to the Highs instance.
 * @param rhs               The right-hand side vector ``b``
 * @param solution_vector   An array of length [num_row] in which to store the
 *                          values of the non-zero elements.
 * @param solution_num_nz   The number of non-zeros in the solution.
 * @param solution_index    An array of length [num_row] in which to store the
 *                          indices of the non-zero elements.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getBasisTransposeSolve(const void* highs, const double* rhs,
                                      double* solution_vector,
                                      HighsInt* solution_nz,
                                      HighsInt* solution_index);

/**
 * Compute a row of ``B^{-1}A``.
 *
 * See `Highs_getBasicVariables` for a description of the ``B`` matrix.
 *
 * The arrays `row_vector` and `row_index` must have an allocated length of
 * [num_col]. However, check `row_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         A pointer to the Highs instance.
 * @param row           The index of the row to compute.
 * @param row_vector    An array of length [num_col] in which to store the
 *                      values of the non-zero elements.
 * @param row_num_nz    The number of non-zeros in the row.
 * @param row_index     An array of length [num_col] in which to store the
 *                      indices of the non-zero elements.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getReducedRow(const void* highs, const HighsInt row,
                             double* row_vector, HighsInt* row_num_nz,
                             HighsInt* row_index);

/**
 * Compute a column of ``B^{-1}A``.
 *
 * See `Highs_getBasicVariables` for a description of the ``B`` matrix.
 *
 * The arrays `col_vector` and `col_index` must have an allocated length of
 * [num_row]. However, check `col_num_nz` to see how many non-zero elements are
 * actually stored.
 *
 * @param highs         A pointer to the Highs instance.
 * @param col           The index of the column to compute.
 * @param col_vector    An array of length [num_row] in which to store the
*                       values of the non-zero elements.
 * @param col_num_nz    The number of non-zeros in the column.
 * @param col_index     An array of length [num_row] in which to store the
*                       indices of the non-zero elements.

 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getReducedColumn(const void* highs, const HighsInt col,
                                double* col_vector, HighsInt* col_num_nz,
                                HighsInt* col_index);

/**
 * Set a basic feasible solution by passing the column and row basis statuses to
 * the model.
 *
 * @param highs       A pointer to the Highs instance.
 * @param col_status  an array of length [num_col] with the column basis status
 *                    in the form of `kHighsBasisStatus` constants
 * @param row_status  an array of length [num_row] with the row basis status
 *                    in the form of `kHighsBasisStatus` constants
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setBasis(void* highs, const HighsInt* col_status,
                        const HighsInt* row_status);

/**
 * Set a logical basis in the model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setLogicalBasis(void* highs);

/**
 * Set a solution by passing the column and row primal and dual solution values.
 *
 * For any values that are unavailable, pass NULL.
 *
 * @param highs       A pointer to the Highs instance.
 * @param col_value   An array of length [num_col] with the column solution
 *                    values.
 * @param row_value   An array of length [num_row] with the row solution
 *                    values.
 * @param col_dual    An array of length [num_col] with the column dual values.
 * @param row_dual    An array of length [num_row] with the row dual values.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setSolution(void* highs, const double* col_value,
                           const double* row_value, const double* col_dual,
                           const double* row_dual);

/**
 * Set a partial primal solution by passing values for a set of variables
 *
 * @param highs       A pointer to the Highs instance.
 * @param num_entries Number of variables in the set
 * @param index       Indices of variables in the set
 * @param value       Values of variables in the set
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setSparseSolution(void* highs, const HighsInt num_entries,
                                 const HighsInt* index, const double* value);

/**
 * Set the callback method to use for HiGHS
 *
 * @param highs              A pointer to the Highs instance.
 * @param user_callback      A pointer to the user callback
 * @param user_callback_data A pointer to the user callback data
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setCallback(void* highs, HighsCCallbackType user_callback,
                           void* user_callback_data);

/**
 * Start callback of given type
 *
 * @param highs         A pointer to the Highs instance.
 * @param callback_type The type of callback to be started
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_startCallback(void* highs, const HighsInt callback_type);

/**
 * Stop callback of given type
 *
 * @param highs         A pointer to the Highs instance.
 * @param callback_type The type of callback to be stopped
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_stopCallback(void* highs, const HighsInt callback_type);

/**
 * Return the cumulative wall-clock time spent in `Highs_run`.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The cumulative wall-clock time spent in `Highs_run`
 */
double Highs_getRunTime(const void* highs);

/**
 * Reset the clocks in a `highs` model.
 *
 * Each `highs` model contains a single instance of clock that records how much
 * time is spent in various parts of the algorithm. This clock is not reset on
 * entry to `Highs_run`, so repeated calls to `Highs_run` report the cumulative
 * time spent in the algorithm. A side-effect is that this will trigger a time
 * limit termination once the cumulative run time exceeds the time limit, rather
 * than the run time of each individual call to `Highs_run`.
 *
 * As a work-around, call `Highs_zeroAllClocks` before each call to `Highs_run`.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_zeroAllClocks(const void* highs);

/**
 * Add a new column (variable) to the model.
 *
 * @param highs         A pointer to the Highs instance.
 * @param cost          The objective coefficient of the column.
 * @param lower         The lower bound of the column.
 * @param upper         The upper bound of the column.
 * @param num_new_nz    The number of non-zeros in the column.
 * @param index         An array of size [num_new_nz] with the row indices.
 * @param value         An array of size [num_new_nz] with row values.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_addCol(void* highs, const double cost, const double lower,
                      const double upper, const HighsInt num_new_nz,
                      const HighsInt* index, const double* value);

/**
 * Add multiple columns (variables) to the model.
 *
 * @param highs         A pointer to the Highs instance.
 * @param num_new_col   The number of new columns to add.
 * @param costs         An array of size [num_new_col] with objective
 *                      coefficients.
 * @param lower         An array of size [num_new_col] with lower bounds.
 * @param upper         An array of size [num_new_col] with upper bounds.
 * @param num_new_nz    The number of new nonzeros in the constraint matrix.
 * @param starts        The constraint coefficients are given as a matrix in
 *                      compressed sparse column form by the arrays `starts`,
 *                      `index`, and `value`. `starts` is an array of size
 *                      [num_new_cols] with the start index of each row in
 *                      indices and values.
 * @param index         An array of size [num_new_nz] with row indices.
 * @param value         An array of size [num_new_nz] with row values.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_addCols(void* highs, const HighsInt num_new_col,
                       const double* costs, const double* lower,
                       const double* upper, const HighsInt num_new_nz,
                       const HighsInt* starts, const HighsInt* index,
                       const double* value);

/**
 * Add a new variable to the model.
 *
 * @param highs         A pointer to the Highs instance.
 * @param lower         The lower bound of the column.
 * @param upper         The upper bound of the column.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_addVar(void* highs, const double lower, const double upper);

/**
 * Add multiple variables to the model.
 *
 * @param highs         A pointer to the Highs instance.
 * @param num_new_var   The number of new variables to add.
 * @param lower         An array of size [num_new_var] with lower bounds.
 * @param upper         An array of size [num_new_var] with upper bounds.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_addVars(void* highs, const HighsInt num_new_var,
                       const double* lower, const double* upper);

/**
 * Add a new row (a linear constraint) to the model.
 *
 * @param highs         A pointer to the Highs instance.
 * @param lower         The lower bound of the row.
 * @param upper         The upper bound of the row.
 * @param num_new_nz    The number of non-zeros in the row
 * @param index         An array of size [num_new_nz] with column indices.
 * @param value         An array of size [num_new_nz] with column values.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_addRow(void* highs, const double lower, const double upper,
                      const HighsInt num_new_nz, const HighsInt* index,
                      const double* value);

/**
 * Add multiple rows (linear constraints) to the model.
 *
 * @param highs         A pointer to the Highs instance.
 * @param num_new_row   The number of new rows to add
 * @param lower         An array of size [num_new_row] with the lower bounds of
 *                      the rows.
 * @param upper         An array of size [num_new_row] with the upper bounds of
 *                      the rows.
 * @param num_new_nz    The number of non-zeros in the rows.
 * @param starts        The constraint coefficients are given as a matrix in
 *                      compressed sparse row form by the arrays `starts`,
 *                      `index`, and `value`. `starts` is an array of size
 *                      [num_new_rows] with the start index of each row in
 *                      indices and values.
 * @param index         An array of size [num_new_nz] with column indices.
 * @param value         An array of size [num_new_nz] with column values.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_addRows(void* highs, const HighsInt num_new_row,
                       const double* lower, const double* upper,
                       const HighsInt num_new_nz, const HighsInt* starts,
                       const HighsInt* index, const double* value);

/**
 * Ensure that the constraint matrix of the incumbent model is stored
 * column-wise.
 *
 * @param highs         A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_ensureColwise(void* highs);

/**
 * Ensure that the constraint matrix of the incumbent model is stored row-wise.
 *
 * @param highs         A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_ensureRowwise(void* highs);

/**
 * Change the objective sense of the model.
 *
 * @param highs     A pointer to the Highs instance.
 * @param sense     The new optimization sense in the form of a `kHighsObjSense`
 *                  constant.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeObjectiveSense(void* highs, const HighsInt sense);

/**
 * Change the objective offset of the model.
 *
 * @param highs     A pointer to the Highs instance.
 * @param offset    The new objective offset.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeObjectiveOffset(void* highs, const double offset);

/**
 * Change the integrality of a column.
 *
 * @param highs         A pointer to the Highs instance.
 * @param col           The column index to change.
 * @param integrality   The new integrality of the column in the form of a
 *                      `kHighsVarType` constant.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColIntegrality(void* highs, const HighsInt col,
                                    const HighsInt integrality);

/**
 * Change the integrality of multiple adjacent columns.
 *
 * @param highs         A pointer to the Highs instance.
 * @param from_col      The index of the first column whose integrality changes.
 * @param to_col        The index of the last column whose integrality
 *                      changes.
 * @param integrality   An array of length [to_col - from_col + 1] with the new
 *                      integralities of the columns in the form of
 *                      `kHighsVarType` constants.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsIntegralityByRange(void* highs,
                                            const HighsInt from_col,
                                            const HighsInt to_col,
                                            const HighsInt* integrality);

/**
 * Change the integrality of multiple columns given by an array of indices.
 *
 * @param highs             A pointer to the Highs instance.
 * @param num_set_entries   The number of columns to change.
 * @param set               An array of size [num_set_entries] with the indices
 *                          of the columns to change.
 * @param integrality       An array of length [num_set_entries] with the new
 *                          integralities of the columns in the form of
 *                          `kHighsVarType` constants.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsIntegralityBySet(void* highs,
                                          const HighsInt num_set_entries,
                                          const HighsInt* set,
                                          const HighsInt* integrality);

/**
 * Change the integrality of multiple columns given by a mask.
 *
 * @param highs         A pointer to the Highs instance.
 * @param mask          An array of length [num_col] with 1 if the column
 *                      integrality should be changed and 0 otherwise.
 * @param integrality   An array of length [num_col] with the new
 *                      integralities of the columns in the form of
 *                      `kHighsVarType` constants.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsIntegralityByMask(void* highs, const HighsInt* mask,
                                           const HighsInt* integrality);

/**
 * Clear the integrality of all columns
 *
 * @param highs         A pointer to the Highs instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_clearIntegrality(void* highs);

/**
 * Change the objective coefficient of a column.
 *
 * @param highs     A pointer to the Highs instance.
 * @param col       The index of the column fo change.
 * @param cost      The new objective coefficient.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColCost(void* highs, const HighsInt col,
                             const double cost);

/**
 * Change the cost coefficients of multiple adjacent columns.
 *
 * @param highs     A pointer to the Highs instance.
 * @param from_col  The index of the first column whose cost changes.
 * @param to_col    The index of the last column whose cost changes.
 * @param cost      An array of length [to_col - from_col + 1] with the new
 *                  objective coefficients.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsCostByRange(void* highs, const HighsInt from_col,
                                     const HighsInt to_col, const double* cost);

/**
 * Change the cost of multiple columns given by an array of indices.
 *
 * @param highs             A pointer to the Highs instance.
 * @param num_set_entries   The number of columns to change.
 * @param set               An array of size [num_set_entries] with the indices
 *                          of the columns to change.
 * @param cost              An array of length [num_set_entries] with the new
 *                          costs of the columns.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsCostBySet(void* highs, const HighsInt num_set_entries,
                                   const HighsInt* set, const double* cost);

/**
 * Change the cost of multiple columns given by a mask.
 *
 * @param highs     A pointer to the Highs instance.
 * @param mask      An array of length [num_col] with 1 if the column
 *                  cost should be changed and 0 otherwise.
 * @param cost      An array of length [num_col] with the new costs.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsCostByMask(void* highs, const HighsInt* mask,
                                    const double* cost);

/**
 * Change the variable bounds of a column.
 *
 * @param highs     A pointer to the Highs instance.
 * @param col       The index of the column whose bounds are to change.
 * @param lower     The new lower bound.
 * @param upper     The new upper bound.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColBounds(void* highs, const HighsInt col,
                               const double lower, const double upper);

/**
 * Change the variable bounds of multiple adjacent columns.
 *
 * @param highs     A pointer to the Highs instance.
 * @param from_col  The index of the first column whose bound changes.
 * @param to_col    The index of the last column whose bound changes.
 * @param lower     An array of length [to_col - from_col + 1] with the new
 *                  lower bounds.
 * @param upper     An array of length [to_col - from_col + 1] with the new
 *                  upper bounds.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsBoundsByRange(void* highs, const HighsInt from_col,
                                       const HighsInt to_col,
                                       const double* lower,
                                       const double* upper);

/**
 * Change the bounds of multiple columns given by an array of indices.
 *
 * @param highs             A pointer to the Highs instance.
 * @param num_set_entries   The number of columns to change.
 * @param set               An array of size [num_set_entries] with the indices
 *                          of the columns to change.
 * @param lower             An array of length [num_set_entries] with the new
 *                          lower bounds.
 * @param upper             An array of length [num_set_entries] with the new
 *                          upper bounds.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper);

/**
 * Change the variable bounds of multiple columns given by a mask.
 *
 * @param highs     A pointer to the Highs instance.
 * @param mask      An array of length [num_col] with 1 if the column
 *                  bounds should be changed and 0 otherwise.
 * @param lower     An array of length [num_col] with the new lower bounds.
 * @param upper     An array of length [num_col] with the new upper bounds.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeColsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower, const double* upper);

/**
 * Change the bounds of a row.
 *
 * @param highs     A pointer to the Highs instance.
 * @param row       The index of the row whose bounds are to change.
 * @param lower     The new lower bound.
 * @param upper     The new upper bound.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeRowBounds(void* highs, const HighsInt row,
                               const double lower, const double upper);

/**
 * Change the variable bounds of multiple adjacent rows.
 *
 * @param highs     A pointer to the Highs instance.
 * @param from_row  The index of the first row whose bound changes.
 * @param to_row    The index of the last row whose bound changes.
 * @param lower     An array of length [to_row - from_row + 1] with the new
 *                  lower bounds.
 * @param upper     An array of length [to_row - from_row + 1] with the new
 *                  upper bounds.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeRowsBoundsByRange(void* highs, const HighsInt from_row,
                                       const HighsInt to_row,
                                       const double* lower,
                                       const double* upper);

/**
 * Change the bounds of multiple rows given by an array of indices.
 *
 * @param highs             A pointer to the Highs instance.
 * @param num_set_entries   The number of rows to change.
 * @param set               An array of size [num_set_entries] with the indices
 *                          of the rows to change.
 * @param lower             An array of length [num_set_entries] with the new
 *                          lower bounds.
 * @param upper             An array of length [num_set_entries] with the new
 *                          upper bounds.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeRowsBoundsBySet(void* highs,
                                     const HighsInt num_set_entries,
                                     const HighsInt* set, const double* lower,
                                     const double* upper);

/**
 * Change the bounds of multiple rows given by a mask.
 *
 * @param highs     A pointer to the Highs instance.
 * @param mask      An array of length [num_row] with 1 if the row
 *                  bounds should be changed and 0 otherwise.
 * @param lower     An array of length [num_row] with the new lower bounds.
 * @param upper     An array of length [num_row] with the new upper bounds.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeRowsBoundsByMask(void* highs, const HighsInt* mask,
                                      const double* lower, const double* upper);

/**
 * Change a coefficient in the constraint matrix.
 *
 * @param highs     A pointer to the Highs instance.
 * @param row       The index of the row to change.
 * @param col       The index of the column to change.
 * @param value     The new constraint coefficient.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_changeCoeff(void* highs, const HighsInt row, const HighsInt col,
                           const double value);

/**
 * Get the objective sense.
 *
 * @param highs     A pointer to the Highs instance.
 * @param sense     The location in which the current objective sense should be
 *                  placed. The sense is a `kHighsObjSense` constant.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getObjectiveSense(const void* highs, HighsInt* sense);

/**
 * Get the objective offset.
 *
 * @param highs     A pointer to the Highs instance.
 * @param offset    The location in which the current objective offset should be
 *                  placed.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getObjectiveOffset(const void* highs, double* offset);

/**
 * Get data associated with multiple adjacent columns from the model.
 *
 * To query the constraint coefficients, this function should be called twice.
 *
 * First, call this function with `matrix_start`, `matrix_index`, and
 * `matrix_value` as `NULL`. This call will populate `num_nz` with the number of
 * nonzero elements in the corresponding section of the constraint matrix.
 *
 * Second, allocate new `matrix_index` and `matrix_value` arrays of length
 * `num_nz` and call this function again to populate the new arrays with their
 * contents.
 *
 * @param highs         A pointer to the Highs instance.
 * @param from_col      The first column for which to query data for.
 * @param to_col        The last column (inclusive) for which to query data for.
 * @param num_col       A HighsInt populated with the number of columns got from
 *                      the model (this should equal `to_col - from_col + 1`).
 * @param costs         An array of size [to_col - from_col + 1] for the column
 *                      cost coefficients.
 * @param lower         An array of size [to_col - from_col + 1] for the column
 *                      lower bounds.
 * @param upper         An array of size [to_col - from_col + 1] for the column
 *                      upper bounds.
 * @param num_nz        A HighsInt to be populated with the number of non-zero
 *                      elements in the constraint matrix.
 * @param matrix_start  An array of size [to_col - from_col + 1] with the start
 *                      indices of each column in `matrix_index` and
 *                      `matrix_value`.
 * @param matrix_index  An array of size [num_nz] with the row indices of each
 *                      element in the constraint matrix.
 * @param matrix_value  An array of size [num_nz] with the non-zero elements of
 *                      the constraint matrix.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getColsByRange(const void* highs, const HighsInt from_col,
                              const HighsInt to_col, HighsInt* num_col,
                              double* costs, double* lower, double* upper,
                              HighsInt* num_nz, HighsInt* matrix_start,
                              HighsInt* matrix_index, double* matrix_value);

/**
 * Get data associated with multiple columns given by an array.
 *
 * This function is identical to `Highs_getColsByRange`, except for how the
 * columns are specified.
 *
 * @param num_set_indices   The number of indices in `set`.
 * @param set               An array of size [num_set_entries] with the column
 *                          indices to get.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getColsBySet(const void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_col,
                            double* costs, double* lower, double* upper,
                            HighsInt* num_nz, HighsInt* matrix_start,
                            HighsInt* matrix_index, double* matrix_value);

/**
 * Get data associated with multiple columns given by a mask.
 *
 * This function is identical to `Highs_getColsByRange`, except for how the
 * columns are specified.
 *
 * @param mask  An array of length [num_col] containing a `1` to get the column
 *              and `0` otherwise.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getColsByMask(const void* highs, const HighsInt* mask,
                             HighsInt* num_col, double* costs, double* lower,
                             double* upper, HighsInt* num_nz,
                             HighsInt* matrix_start, HighsInt* matrix_index,
                             double* matrix_value);

/**
 * Get data associated with multiple adjacent rows from the model.
 *
 * To query the constraint coefficients, this function should be called twice.
 *
 * First, call this function with `matrix_start`, `matrix_index`, and
 * `matrix_value` as `NULL`. This call will populate `num_nz` with the number of
 * nonzero elements in the corresponding section of the constraint matrix.
 *
 * Second, allocate new `matrix_index` and `matrix_value` arrays of length
 * `num_nz` and call this function again to populate the new arrays with their
 * contents.
 *
 * @param highs         A pointer to the Highs instance.
 * @param from_row      The first row for which to query data for.
 * @param to_row        The last row (inclusive) for which to query data for.
 * @param num_row       A HighsInt to be populated with the number of rows got
 *                      from the model.
 * @param lower         An array of size [to_row - from_row + 1] for the row
 *                      lower bounds.
 * @param upper         An array of size [to_row - from_row + 1] for the row
 *                      upper bounds.
 * @param num_nz        A HighsInt to be populated with the number of non-zero
 *                      elements in the constraint matrix.
 * @param matrix_start  An array of size [to_row - from_row + 1] with the start
 *                      indices of each row in `matrix_index` and
 *                      `matrix_value`.
 * @param matrix_index  An array of size [num_nz] with the column indices of
 *                      each element in the constraint matrix.
 * @param matrix_value  An array of size [num_nz] with the non-zero elements of
 *                      the constraint matrix.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getRowsByRange(const void* highs, const HighsInt from_row,
                              const HighsInt to_row, HighsInt* num_row,
                              double* lower, double* upper, HighsInt* num_nz,
                              HighsInt* matrix_start, HighsInt* matrix_index,
                              double* matrix_value);

/**
 * Get data associated with multiple rows given by an array.
 *
 * This function is identical to `Highs_getRowsByRange`, except for how the
 * rows are specified.
 *
 * @param num_set_indices   The number of indices in `set`.
 * @param set               An array of size [num_set_entries] containing the
 *                          row indices to get.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getRowsBySet(const void* highs, const HighsInt num_set_entries,
                            const HighsInt* set, HighsInt* num_row,
                            double* lower, double* upper, HighsInt* num_nz,
                            HighsInt* matrix_start, HighsInt* matrix_index,
                            double* matrix_value);

/**
 * Get data associated with multiple rows given by a mask.
 *
 * This function is identical to `Highs_getRowsByRange`, except for how the
 * rows are specified.
 *
 * @param mask  An array of length [num_row] containing a `1` to get the row and
 *              `0` otherwise.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getRowsByMask(const void* highs, const HighsInt* mask,
                             HighsInt* num_row, double* lower, double* upper,
                             HighsInt* num_nz, HighsInt* matrix_start,
                             HighsInt* matrix_index, double* matrix_value);
/**
 * Get the name of a row.
 *
 * @param row   The index of the row to query.
 * @param name  A pointer in which to store the name of the row. This must have
 *              length `kHighsMaximumStringLength`.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getRowName(const void* highs, const HighsInt row, char* name);

/**
 * Get the index of a row from its name.
 *
 * If multiple rows have the same name, or if no row exists with `name`, this
 * function returns `kHighsStatusError`.
 *
 * @param name A pointer of the name of the row to query.
 * @param row  A pointer in which to store the index of the row
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getRowByName(const void* highs, const char* name, HighsInt* row);

/**
 * Get the name of a column.
 *
 * @param col   The index of the column to query.
 * @param name  A pointer in which to store the name of the column. This must
 *              have length `kHighsMaximumStringLength`.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getColName(const void* highs, const HighsInt col, char* name);

/**
 * Get the index of a column from its name.
 *
 * If multiple columns have the same name, or if no column exists with `name`,
 * this function returns `kHighsStatusError`.
 *
 * @param name A pointer of the name of the column to query.
 * @param col  A pointer in which to store the index of the column
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getColByName(const void* highs, const char* name, HighsInt* col);

/**
 * Get the integrality of a column.
 *
 * @param col          The index of the column to query.
 * @param integrality  A HighsInt in which the integrality of the column should
 *                     be placed. The integer is one of the `kHighsVarTypeXXX`
 *                     constants.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getColIntegrality(const void* highs, const HighsInt col,
                                 HighsInt* integrality);

/**
 * Delete multiple adjacent columns.
 *
 * @param highs     A pointer to the Highs instance.
 * @param from_col  The index of the first column to delete.
 * @param to_col    The index of the last column to delete.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_deleteColsByRange(void* highs, const HighsInt from_col,
                                 const HighsInt to_col);

/**
 * Delete multiple columns given by an array of indices.
 *
 * @param highs             A pointer to the Highs instance.
 * @param num_set_entries   The number of columns to delete.
 * @param set               An array of size [num_set_entries] with the indices
 *                          of the columns to delete.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_deleteColsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set);

/**
 * Delete multiple columns given by a mask.
 *
 * @param highs     A pointer to the Highs instance.
 * @param mask      An array of length [num_col] with 1 if the column
 *                  should be deleted and 0 otherwise.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_deleteColsByMask(void* highs, HighsInt* mask);

/**
 * Delete multiple adjacent rows.
 *
 * @param highs     A pointer to the Highs instance.
 * @param from_row  The index of the first row to delete.
 * @param to_row    The index of the last row to delete.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_deleteRowsByRange(void* highs, const HighsInt from_row,
                                 const HighsInt to_row);

/**
 * Delete multiple rows given by an array of indices.
 *
 * @param highs             A pointer to the Highs instance.
 * @param num_set_entries   The number of rows to delete.
 * @param set               An array of size [num_set_entries] with the indices
 *                          of the rows to delete.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_deleteRowsBySet(void* highs, const HighsInt num_set_entries,
                               const HighsInt* set);

/**
 * Delete multiple rows given by a mask.
 *
 * @param highs     A pointer to the Highs instance.
 * @param mask      An array of length [num_row] with `1` if the row should be
 *                  deleted and `0` otherwise. The new index of any column not
 *                  deleted is stored in place of the value `0`.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_deleteRowsByMask(void* highs, HighsInt* mask);

/**
 * Scale a column by a constant.
 *
 * Scaling a column modifies the elements in the constraint matrix, the variable
 * bounds, and the objective coefficient.
 *
 * @param highs     A pointer to the Highs instance.
 * @param col       The index of the column to scale.
 * @param scaleval  The value by which to scale the column. If `scaleval < 0`,
 *                  the variable bounds flipped.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_scaleCol(void* highs, const HighsInt col, const double scaleval);

/**
 * Scale a row by a constant.
 *
 * @param highs     A pointer to the Highs instance.
 * @param row       The index of the row to scale.
 * @param scaleval  The value by which to scale the row. If `scaleval < 0`, the
 *                  row bounds are flipped.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_scaleRow(void* highs, const HighsInt row, const double scaleval);

/**
 * Return the value of infinity used by HiGHS.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The value of infinity used by HiGHS.
 */
double Highs_getInfinity(const void* highs);

/**
 * Return the size of integers used by HiGHS.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The size of integers used by HiGHS.
 */
HighsInt Highs_getSizeofHighsInt(const void* highs);

/**
 * Return the number of columns in the model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The number of columns in the model.
 */
HighsInt Highs_getNumCol(const void* highs);

/**
 * Return the number of rows in the model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The number of rows in the model.
 */
HighsInt Highs_getNumRow(const void* highs);

/**
 * Return the number of nonzeros in the constraint matrix of the model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The number of nonzeros in the constraint matrix of the model.
 */
HighsInt Highs_getNumNz(const void* highs);

/**
 * Return the number of nonzeroes in the Hessian matrix of the model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The number of nonzeroes in the Hessian matrix of the model.
 */
HighsInt Highs_getHessianNumNz(const void* highs);

/**
 * Return the number of columns in the presolved model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The number of columns in the presolved model.
 */
HighsInt Highs_getPresolvedNumCol(const void* highs);

/**
 * Return the number of rows in the presolved model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The number of rows in the presolved model.
 */
HighsInt Highs_getPresolvedNumRow(const void* highs);

/**
 * Return the number of nonzeros in the constraint matrix of the presolved
 * model.
 *
 * @param highs     A pointer to the Highs instance.
 *
 * @returns The number of nonzeros in the constraint matrix of the presolved
 * model.
 */
HighsInt Highs_getPresolvedNumNz(const void* highs);

/**
 * Get the data from a HiGHS model.
 *
 * The input arguments have the same meaning (in a different order) to those
 * used in `Highs_passModel`.
 *
 * Note that all arrays must be pre-allocated to the correct size before calling
 * `Highs_getModel`. Use the following query methods to check the appropriate
 * size:
 *  - `Highs_getNumCol`
 *  - `Highs_getNumRow`
 *  - `Highs_getNumNz`
 *  - `Highs_getHessianNumNz`
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getModel(const void* highs, const HighsInt a_format,
                        const HighsInt q_format, HighsInt* num_col,
                        HighsInt* num_row, HighsInt* num_nz,
                        HighsInt* hessian_num_nz, HighsInt* sense,
                        double* offset, double* col_cost, double* col_lower,
                        double* col_upper, double* row_lower, double* row_upper,
                        HighsInt* a_start, HighsInt* a_index, double* a_value,
                        HighsInt* q_start, HighsInt* q_index, double* q_value,
                        HighsInt* integrality);

/**
 * Get the data from a HiGHS LP.
 *
 * The input arguments have the same meaning (in a different order) to those
 * used in `Highs_passModel`.
 *
 * Note that all arrays must be pre-allocated to the correct size before calling
 * `Highs_getModel`. Use the following query methods to check the appropriate
 * size:
 *  - `Highs_getNumCol`
 *  - `Highs_getNumRow`
 *  - `Highs_getNumNz`
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getLp(const void* highs, const HighsInt a_format,
                     HighsInt* num_col, HighsInt* num_row, HighsInt* num_nz,
                     HighsInt* sense, double* offset, double* col_cost,
                     double* col_lower, double* col_upper, double* row_lower,
                     double* row_upper, HighsInt* a_start, HighsInt* a_index,
                     double* a_value, HighsInt* integrality);

/**
 * Get the data from a HiGHS presolved LP.
 *
 * The input arguments have the same meaning (in a different order) to those
 * used in `Highs_passModel`.
 *
 * Note that all arrays must be pre-allocated to the correct size before calling
 * `Highs_getModel`. Use the following query methods to check the appropriate
 * size:
 *  - `Highs_getPresolvedNumCol`
 *  - `Highs_getPresolvedNumRow`
 *  - `Highs_getPresolvedNumNz`
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getPresolvedLp(const void* highs, const HighsInt a_format,
                              HighsInt* num_col, HighsInt* num_row,
                              HighsInt* num_nz, HighsInt* sense, double* offset,
                              double* col_cost, double* col_lower,
                              double* col_upper, double* row_lower,
                              double* row_upper, HighsInt* a_start,
                              HighsInt* a_index, double* a_value,
                              HighsInt* integrality);
/**
 * Get the name of a column of the presolved LP.
 *
 * @param col   The index of the column to query.
 * @param name  A pointer in which to store the name of the column. This must
 *              have length `kHighsMaximumStringLength`.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getPresolvedColName(const void* highs, const HighsInt col,
                                   char* name);

/**
 * Get the name of a row of the presolved LP.
 *
 * @param row   The index of the row to query.
 * @param name  A pointer in which to store the name of the row. This must
 *              have length `kHighsMaximumStringLength`.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getPresolvedRowName(const void* highs, const HighsInt row,
                                   char* name);

/**
 * Get the data from a HiGHS IIS LP.
 *
 * The input arguments have the same meaning (in a different order) to those
 * used in `Highs_passModel`.
 *
 * Note that all arrays must be pre-allocated to the correct size before calling
 * `Highs_getModel`. Use the following query methods to check the appropriate
 * size:
 *  - `Highs_getPresolvedNumCol`
 *  - `Highs_getPresolvedNumRow`
 *  - `Highs_getPresolvedNumNz`
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getIisLp(const void* highs, const HighsInt a_format,
                        HighsInt* num_col, HighsInt* num_row, HighsInt* num_nz,
                        HighsInt* sense, double* offset, double* col_cost,
                        double* col_lower, double* col_upper, double* row_lower,
                        double* row_upper, HighsInt* a_start, HighsInt* a_index,
                        double* a_value, HighsInt* integrality);

/**
 * Get the LP corresponding to a MIP with non-continuous variables
 * fixed at a MIP solution
 *
 * The input arguments have the same meaning (in a different order) to those
 * used in `Highs_passModel`.
 *
 * Note that all arrays must be pre-allocated to the correct size before calling
 * `Highs_getModel`. Use the following query methods to check the appropriate
 * size:
 *  - `Highs_getNumCol`
 *  - `Highs_getNumRow`
 *  - `Highs_getNumNz`
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getFixedLp(const void* highs, const HighsInt a_format,
                          HighsInt* num_col, HighsInt* num_row,
                          HighsInt* num_nz, HighsInt* sense, double* offset,
                          double* col_cost, double* col_lower,
                          double* col_upper, double* row_lower,
                          double* row_upper, HighsInt* a_start,
                          HighsInt* a_index, double* a_value);

/**
 * Set a primal (and possibly dual) solution as a starting point, then run
 * crossover to compute a basic feasible solution.
 *
 * @param highs      A pointer to the Highs instance.
 * @param num_col    The number of variables.
 * @param num_row    The number of rows.
 * @param col_value  An array of length [num_col] with optimal primal solution
 *                   for each column.
 * @param col_dual   An array of length [num_col] with optimal dual solution for
 *                   each column. May be `NULL`, in which case no dual solution
 *                   is passed.
 * @param row_dual   An array of length [num_row] with optimal dual solution for
 *                   each row. . May be `NULL`, in which case no dual solution
 *                   is passed.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_crossover(void* highs, const HighsInt num_col,
                         const HighsInt num_row, const double* col_value,
                         const double* col_dual, const double* row_dual);

/**
 * Compute the ranging information for all costs and bounds. For
 * nonbasic variables the ranging information is relative to the
 * active bound. For basic variables the ranging information relates
 * to...
 *
 * For any values that are not required, pass NULL.
 *
 * @param highs                  A pointer to the Highs instance.
 * @param col_cost_up_value      The upper range of the cost value
 * @param col_cost_up_objective  The objective at the upper cost range
 * @param col_cost_up_in_var     The variable entering the basis at the upper
 *                               cost range
 * @param col_cost_up_ou_var     The variable leaving the basis at the upper
 *                               cost range
 * @param col_cost_dn_value      The lower range of the cost value
 * @param col_cost_dn_objective  The objective at the lower cost range
 * @param col_cost_dn_in_var     The variable entering the basis at the lower
 *                               cost range
 * @param col_cost_dn_ou_var     The variable leaving the basis at the lower
 *                               cost range
 * @param col_bound_up_value     The upper range of the column bound value
 * @param col_bound_up_objective The objective at the upper column bound range
 * @param col_bound_up_in_var    The variable entering the basis at the upper
 *                               column bound range
 * @param col_bound_up_ou_var    The variable leaving the basis at the upper
 *                               column bound range
 * @param col_bound_dn_value     The lower range of the column bound value
 * @param col_bound_dn_objective The objective at the lower column bound range
 * @param col_bound_dn_in_var    The variable entering the basis at the lower
 *                               column bound range
 * @param col_bound_dn_ou_var    The variable leaving the basis at the lower
 *                               column bound range
 * @param row_bound_up_value     The upper range of the row bound value
 * @param row_bound_up_objective The objective at the upper row bound range
 * @param row_bound_up_in_var    The variable entering the basis at the upper
 *                               row bound range
 * @param row_bound_up_ou_var    The variable leaving the basis at the upper row
 *                               bound range
 * @param row_bound_dn_value     The lower range of the row bound value
 * @param row_bound_dn_objective The objective at the lower row bound range
 * @param row_bound_dn_in_var    The variable entering the basis at the lower
 *                               row bound range
 * @param row_bound_dn_ou_var    The variable leaving the basis at the lower row
 *                               bound range
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getRanging(
    void* highs,
    //
    double* col_cost_up_value, double* col_cost_up_objective,
    HighsInt* col_cost_up_in_var, HighsInt* col_cost_up_ou_var,
    double* col_cost_dn_value, double* col_cost_dn_objective,
    HighsInt* col_cost_dn_in_var, HighsInt* col_cost_dn_ou_var,
    double* col_bound_up_value, double* col_bound_up_objective,
    HighsInt* col_bound_up_in_var, HighsInt* col_bound_up_ou_var,
    double* col_bound_dn_value, double* col_bound_dn_objective,
    HighsInt* col_bound_dn_in_var, HighsInt* col_bound_dn_ou_var,
    double* row_bound_up_value, double* row_bound_up_objective,
    HighsInt* row_bound_up_in_var, HighsInt* row_bound_up_ou_var,
    double* row_bound_dn_value, double* row_bound_dn_objective,
    HighsInt* row_bound_dn_in_var, HighsInt* row_bound_dn_ou_var);

/**
 * Compute the solution corresponding to a (possibly weighted) sum of
 * (allowable) infeasibilities in an LP/MIP.
 *
 * If local penalties are not defined, pass NULL, and the global
 * penalty will be used. Negative penalty values imply that the bound
 * or RHS value cannot be violated
 *
 * @param highs                             A pointer to the Highs instance.
 * @param const double global_lower_penalty The penalty for violating lower
 * bounds on variables
 * @param const double global_upper_penalty The penalty for violating upper
 * bounds on variables
 * @param const double global_rhs_penalty   The penalty for violating constraint
 * RHS values
 * @param const double* local_lower_penalty The penalties for violating specific
 * lower bounds on variables
 * @param const double* local_upper_penalty The penalties for violating specific
 * upper bounds on variables
 * @param const double* local_rhs_penalty   The penalties for violating specific
 * constraint RHS values
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */

HighsInt Highs_feasibilityRelaxation(void* highs,
                                     const double global_lower_penalty,
                                     const double global_upper_penalty,
                                     const double global_rhs_penalty,
                                     const double* local_lower_penalty,
                                     const double* local_upper_penalty,
                                     const double* local_rhs_penalty);

/**
 * Attempt to compute an irreducible infeasibility subsystem (IIS) for
 * an LP, QP, or the relaxation of a MIP. If no IIS is found, then the
 * number of IIS columns and rows will be zero.
 *
 * @param highs                      A pointer to the Highs instance.
 * @param const HighsInt iis_num_col Number of columns in the IIS.
 * @param const HighsInt iis_num_row Number of rows in the IIS.
 * @param const HighsInt* col_index  An array of length [iis_num_col], to be
 *                                   filled with the indices of original
 *                                   variables in the IIS.
 * @param const HighsInt* row_index  An array of length [iis_num_col], to be
 *                                   filled with the indices of original
 *                                   constraints in the IIS.
 * @param const HighsInt* col_bound  An array of length [iis_num_col], to be
 *                                   filled with the bound status of variables
 *                                   in the IIS.
 * @param const HighsInt* row_bound  An array of length [iis_num_col], to be
 *                                   filled with the bound status of constraints
 *                                   in the IIS.
 * @param const HighsInt* col_status An array of length [num_col], to be
 *                                   filled with the IIS status of all original
 *                                   variables.
 * @param const HighsInt* row_status n array of length [num_col], to be
 *                                   filled with the IIS status of all original
 *                                   constraints.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_getIis(void* highs, HighsInt* iis_num_col, HighsInt* iis_num_row,
                      HighsInt* col_index, HighsInt* row_index,
                      HighsInt* col_bound, HighsInt* row_bound,
                      HighsInt* col_status, HighsInt* row_status);
/**
 * Releases all resources held by the scheduler owned by the calling thread, if
 * any.
 *
 * Although the scheduler instance is created internally by a Highs
 * instance, any subsequent Highs instances on the same thread share the
 * scheduler. Hence, calling `Highs_destroy` does not free the
 * scheduler's memory, so `Highs_resetThreadScheduler` must also be
 * called.
 *
 * It is thread-safe to call this function only from the thread that owns the
 * scheduler (not from a worker thread) and when no run() or presolve() is in
 * progress on that thread.
 *
 * After this function has terminated, it is guaranteed that eventually all
 * previously created scheduler threads will terminate and allocated memory will
 * be released.
 *
 * After this function has returned, the option value for the number of threads
 * may be altered to a new value before the next call to `Highs_run` or one of
 * the `Highs_XXXcall` methods.
 *
 * @param blocking   If the `blocking` parameter has a nonzero value, then this
 *                   function will not return until all memory is freed, which
 *                   might be desirable when debugging heap memory, but it
 *                   requires the calling thread to wait for all scheduler
 *                   threads to wake-up which is usually not necessary.
 *
 * @returns No status is returned since the function call cannot fail. Calling
 * this function while on a worker thread or during the execution of a Highs
 * function may cause crashes, but cannot be detected and hence is fully the
 * caller's responsibility.
 */
void Highs_resetThreadScheduler(const HighsInt blocking);

/**
 * Get a void* pointer to a callback data item
 *
 * @param data_out      A pointer to the HighsCallbackDataOut instance.
 * @param item_name     The name of the item.
 *
 * @returns A void* pointer to the callback data item, or NULL if item_name not
 * valid
 */
const void* Highs_getCallbackDataOutItem(const HighsCallbackDataOut* data_out,
                                         const char* item_name);

/**
 * Set a solution within a callback by passing a subset of the values.
 *
 * For any values that are unavailable/unknown, pass kHighsUndefined.
 *
 * @param data_in     A pointer to the callback input data instance.
 * @param num_entries Number of variables in the set
 * @param value       An array of length [num_entries <= num_col] with
 *                    column solution values.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setCallbackSolution(HighsCallbackDataIn* data_in,
                                   const HighsInt num_entries,
                                   const double* value);

/**
 * Set a partial primal solution by passing values for a set of variables,
 * within a valid callback.
 *
 * @param data_in     A pointer to the callback input data instance.
 * @param num_entries Number of variables in the set
 * @param index       Indices of variables in the set
 * @param value       Values of variables in the set
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_setCallbackSparseSolution(HighsCallbackDataIn* data_in,
                                         const HighsInt num_entries,
                                         const HighsInt* index,
                                         const double* value);

/**
 * Finds a feasible solution for a given (partial) primal user solution,
 * within a valid callback.
 *
 * On success, the user solution is updated within the callback input data
 * instance.
 *
 * @returns A `kHighsStatus` constant indicating whether the call succeeded.
 */
HighsInt Highs_repairCallbackSolution(HighsCallbackDataIn* data_in);

/**
 * @brief Initialize optional extras and cache the result.
 *
 * @param directory Directory containing the extras library. NULL means empty.
 * @param search_if_empty If nonzero, search default library directories when
 *                        directory is NULL or empty.
 * @return kHighsStatusOk on success, kHighsStatusError on failure.
 *
 * Only the first call attempts initialization; subsequent calls return the
 * cached result, regardless of arguments, unless extras are explicitly
 * unloaded. Not needed for statically linked extras.
 */
HighsInt Highs_initializeExtras(const char* directory,
                                HighsInt search_if_empty);

// *********************
// * Deprecated methods*
// *********************

/**
 * Return the HiGHS compilation date.
 *
 * @returns The HiGHS compilation date.
 */
static const char* Highs_compilationDate(void);

// These are deprecated because they don't follow the style guide. Constants
// must begin with `k`.
static const HighsInt HighsStatuskError = -1;
static const HighsInt HighsStatuskOk = 0;
static const HighsInt HighsStatuskWarning = 1;

HighsInt Highs_call(const HighsInt num_col, const HighsInt num_row,
                    const HighsInt num_nz, const double* col_cost,
                    const double* col_lower, const double* col_upper,
                    const double* row_lower, const double* row_upper,
                    const HighsInt* a_start, const HighsInt* a_index,
                    const double* a_value, double* col_value, double* col_dual,
                    double* row_value, double* row_dual,
                    HighsInt* col_basis_status, HighsInt* row_basis_status,
                    HighsInt* model_status);

HighsInt Highs_runQuiet(void* highs);

HighsInt Highs_setHighsLogfile(void* highs, const void* logfile);

HighsInt Highs_setHighsOutput(void* highs, const void* outputfile);

HighsInt Highs_getIterationCount(const void* highs);

HighsInt Highs_getSimplexIterationCount(const void* highs);

HighsInt Highs_setHighsBoolOptionValue(void* highs, const char* option,
                                       const HighsInt value);

HighsInt Highs_setHighsIntOptionValue(void* highs, const char* option,
                                      const HighsInt value);

HighsInt Highs_setHighsDoubleOptionValue(void* highs, const char* option,
                                         const double value);

HighsInt Highs_setHighsStringOptionValue(void* highs, const char* option,
                                         const char* value);

HighsInt Highs_setHighsOptionValue(void* highs, const char* option,
                                   const char* value);

HighsInt Highs_getHighsBoolOptionValue(const void* highs, const char* option,
                                       HighsInt* value);

HighsInt Highs_getHighsIntOptionValue(const void* highs, const char* option,
                                      HighsInt* value);

HighsInt Highs_getHighsDoubleOptionValue(const void* highs, const char* option,
                                         double* value);

HighsInt Highs_getHighsStringOptionValue(const void* highs, const char* option,
                                         char* value);

HighsInt Highs_getHighsOptionType(const void* highs, const char* option,
                                  HighsInt* type);

HighsInt Highs_resetHighsOptions(void* highs);

HighsInt Highs_getHighsIntInfoValue(const void* highs, const char* info,
                                    HighsInt* value);

HighsInt Highs_getHighsDoubleInfoValue(const void* highs, const char* info,
                                       double* value);

HighsInt Highs_getNumCols(const void* highs);

HighsInt Highs_getNumRows(const void* highs);

double Highs_getHighsInfinity(const void* highs);

double Highs_getHighsRunTime(const void* highs);

HighsInt Highs_setOptionValue(void* highs, const char* option,
                              const char* value);

HighsInt Highs_getScaledModelStatus(const void* highs);

void Highs_resetGlobalScheduler(const HighsInt blocking);

#ifdef __cplusplus
}
#endif

#endif
