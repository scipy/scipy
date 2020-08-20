/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSIMPLEX_H_
#define SIMPLEX_HSIMPLEX_H_

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsStatus.h"

void setSimplexOptions(
    HighsModelObject& highs_model_object  //!< Model object in which simplex
                                          //!< options are to be set
);

HighsStatus transition(HighsModelObject& highs_model_object  //!< Model object
);

void setNonbasicFlag(const HighsLp& simplex_lp, vector<int>& nonbasicFlag,
                     const HighsBasisStatus* col_status = NULL,
                     const HighsBasisStatus* row_status = NULL);

void setNonbasicMove(const HighsLp& simplex_lp, const HighsScale& scale,
                     const bool have_highs_basis, const HighsBasis& basis,
                     const bool have_highs_solution,
                     const HighsSolution& solution,
                     SimplexBasis& simplex_basis);

void initialiseNonbasicWorkValue(const HighsLp& simplex_lp,
                                 const SimplexBasis& simplex_basis,
                                 HighsSimplexInfo& simplex_info);

bool basisConditionOk(HighsModelObject& highs_model_object);

// Methods not requiring HighsModelObject

bool dual_infeasible(const double value, const double lower, const double upper,
                     const double dual, const double value_tolerance,
                     const double dual_tolerance);

void appendNonbasicColsToBasis(HighsLp& lp, HighsBasis& basis, int XnumNewCol);
void appendNonbasicColsToBasis(HighsLp& lp, SimplexBasis& basis,
                               int XnumNewCol);

void appendBasicRowsToBasis(HighsLp& lp, HighsBasis& basis, int XnumNewRow);
void appendBasicRowsToBasis(HighsLp& lp, SimplexBasis& basis, int XnumNewRow);

void reportBasis(const HighsOptions options, const HighsLp& lp,
                 const HighsBasis& basis);
void reportBasis(const HighsOptions options, const HighsLp& lp,
                 const SimplexBasis& simplex_basis);

void computeDualObjectiveValue(HighsModelObject& highs_model_object,
                               int phase = 2);

void computePrimalObjectiveValue(HighsModelObject& highs_model_object);

int setSourceOutFmBd(const HighsModelObject& highs_model_object,
                     const int column_out);

#ifdef HiGHSDEV
void getPrimalValue(const HighsModelObject& highs_model_object,
                    vector<double>& primal_value);
void analysePrimalObjectiveValue(const HighsModelObject& highs_model_object);
#endif

void initialiseSimplexLpDefinition(HighsModelObject& highs_model);
void initialiseSimplexLpRandomVectors(HighsModelObject& highs_model);

// SCALE:

void scaleHighsModelInit(HighsModelObject& highs_model);

void scaleCosts(HighsModelObject& highs_model);

void scaleFactorRanges(HighsModelObject& highs_model_object,
                       double& min_col_scale, double& max_col_scale,
                       double& min_row_scale, double& max_row_scale);

void scaleSimplexLp(HighsModelObject& highs_model);
bool equilibrationScaleMatrix(HighsModelObject& highs_model);
bool maxValueScaleMatrix(HighsModelObject& highs_model);

HighsStatus deleteScale(const HighsOptions& options, vector<double>& scale,
                        const HighsIndexCollection& index_collection);
// PERMUTE:

void permuteSimplexLp(HighsModelObject& highs_model);

#ifdef HiGHSDEV
// Only used to analyse the row and column status after Crash
void initialise_basic_index(HighsModelObject& highs_model_object);
#endif

void allocateWorkAndBaseArrays(HighsModelObject& highs_model_object);

void initialiseValueAndNonbasicMove(HighsModelObject& highs_model_object);

void initialisePhase2ColBound(HighsModelObject& highs_model_object);

void initialisePhase2RowBound(HighsModelObject& highs_model_object);

void initialiseBound(HighsModelObject& highs_model_object, int phase = 2);

void initialisePhase2ColCost(HighsModelObject& highs_model_object);

void initialisePhase2RowCost(HighsModelObject& highs_model_object);

void initialiseCost(HighsModelObject& highs_model_object, int perturb = 0);

void populateWorkArrays(HighsModelObject& highs_model_object);

#ifdef HiGHSDEV
void reportSimplexProfiling(HighsModelObject& highs_model_object);

#endif
void setRunQuiet(HighsModelObject& highs_model_object);
/**
 * @brief Get the Hager condition number estimate for the basis matrix of a
 * model
 */
double computeBasisCondition(const HighsModelObject& highs_model_object);

void flip_bound(HighsModelObject& highs_model_object, int iCol);

int simplexHandleRankDeficiency(HighsModelObject& highs_model_object);

int computeFactor(HighsModelObject& highs_model_object);

// Compute the primal values (in baseValue) and set the lower and upper bounds
// of basic variables
void computePrimal(HighsModelObject& highs_model_object);

void computeSimplexInfeasible(HighsModelObject& highs_model_object);
void computeSimplexPrimalInfeasible(HighsModelObject& highs_model_object);
void computeSimplexDualInfeasible(HighsModelObject& highs_model_object);

void computeDualInfeasibleWithFlips(HighsModelObject& highs_model_object);

void computeSimplexLpDualInfeasible(HighsModelObject& highs_model_object);

void copySimplexInfeasible(HighsModelObject& highs_model_object);
void copySimplexDualInfeasible(HighsModelObject& highs_model_object);
void copySimplexPrimalInfeasible(HighsModelObject& highs_model_object);

void choosePriceTechnique(const int price_strategy, const double row_ep_density,
                          bool& use_col_price, bool& use_row_price_w_switch);

void computeTableauRowFromPiP(HighsModelObject& highs_model_object,
                              const HVector& row_ep, HVector& row_ap);

void computeDual(HighsModelObject& highs_model_object);

void correctDual(HighsModelObject& highs_model_object,
                 int* free_infeasibility_count);
void correctDual(HighsModelObject& highs_model_object);

// Record the shift in the cost of a particular column
void shift_cost(HighsModelObject& highs_model_object, int iCol, double amount);

// Undo the shift in the cost of a particular column
void shift_back(HighsModelObject& highs_model_object, int iCol);

// The major model updates. Factor calls factor.update; Matrix
// calls matrix.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void update_factor(HighsModelObject& highs_model_object, HVector* column,
                   HVector* row_ep, int* iRow, int* hint);

void update_pivots(HighsModelObject& highs_model_object, int columnIn,
                   int rowOut, int sourceOut);

void update_matrix(HighsModelObject& highs_model_object, int columnIn,
                   int columnOut);

bool reinvertOnNumericalTrouble(const std::string method_name,
                                const HighsModelObject& highs_model_object,
                                double& numerical_trouble_measure,
                                const double alpha_from_col,
                                const double alpha_from_row,
                                const double numerical_trouble_tolerance);

// Analyse the unscaled solution from a Simplex basic solution to get
// suggested feasibility tolerances for resolving the scaled LP
// This sets highs_model_object.unscaled_solution_params_
HighsStatus getNewInfeasibilityTolerancesFromSimplexBasicSolution(
    const HighsModelObject& highs_model_object,
    HighsSolutionParams& get_unscaled_solution_params,
    double& new_scaled_primal_feasibility_tolerance,
    double& new_scaled_dual_feasibility_tolerance);

HighsStatus getInfeasibilitiesAndNewTolerances(
    const HighsOptions& options, const HighsLp& lp, const HighsScale& scale,
    const SimplexBasis& basis, const HighsSimplexInfo& simplex_info,
    const HighsModelStatus scaled_model_status,
    const HighsSolutionParams& unscaled_solution_params,
    const HighsSolutionParams& scaled_solution_params,
    HighsSolutionParams& get_unscaled_solution_params,
    HighsSolutionParams& get_scaled_solution_params,
    double& new_scaled_primal_feasibility_tolerance,
    double& new_scaled_dual_feasibility_tolerance);

void logRebuild(HighsModelObject& highs_model_object, const bool primal,
                const int solve_phase);

void reportSimplexLpStatus(
    HighsSimplexLpStatus&
        simplex_lp_status,    // !< Status of simplex LP to be reported
    const char* message = ""  // !< Message to be written in report
);

void invalidateSimplexLpBasisArtifacts(
    HighsSimplexLpStatus&
        simplex_lp_status  // !< Status of simplex LP whose
                           // basis artifacts are to be invalidated
);

void invalidateSimplexLpBasis(
    HighsSimplexLpStatus& simplex_lp_status  // !< Status of simplex LP whose
                                             // basis is to be invalidated
);

void invalidateSimplexLp(
    HighsSimplexLpStatus&
        simplex_lp_status  // !< Status of simplex LP to be invalidated
);

void updateSimplexLpStatus(
    HighsSimplexLpStatus&
        simplex_lp_status,  // !< Status of simplex LP to be updated
    LpAction action         // !< Action prompting update
);

bool isBasisRightSize(const HighsLp& lp, const SimplexBasis& basis);
#endif  // SIMPLEX_HSIMPLEX_H_
