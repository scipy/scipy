/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsModelUtils.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSMODELUTILS_H_
#define LP_DATA_HIGHSMODELUTILS_H_

#include "lp_data/HighsInfo.h"
#include "model/HighsModel.h"

const std::string kLegalLpFileColRowNameFirstChar =
    "abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "0123456789."
    "!\"#$%&(),.;?@_‘’{}~";
// According to CPLEX (cf fix-2887) names in LP files cannot begin
// with a digit or a full stop. However, HiGHS can read these, and
// models like rgn have row names that are integers, so allow it
const std::string kLegalLpFileColRowNameChar = kLegalLpFileColRowNameFirstChar;

// Analyse lower and upper bounds of a model
void analyseModelBounds(const HighsLogOptions& log_options, const char* message,
                        HighsInt numBd, const std::vector<double>& lower,
                        const std::vector<double>& upper);
bool hasNamesWithSpaces(const HighsLogOptions& log_options, const HighsLp& lp);
bool hasNamesWithSpaces(const HighsLogOptions& log_options, const bool col,
                        const std::vector<std::string>& names);
bool hasIllegalNameForLpFile(const std::vector<std::string>& names);

void writeModelBoundSolution(
    FILE* file, const HighsLogOptions& log_options, const bool columns,
    const HighsInt dim, const std::vector<double>& lower,
    const std::vector<double>& upper, const std::vector<std::string>& names,
    const bool have_primal, const std::vector<double>& primal,
    const bool have_dual, const std::vector<double>& dual,
    const bool have_basis, const std::vector<HighsBasisStatus>& status,
    const HighsVarType* integrality = NULL);

void writeModelObjective(FILE* file, const HighsLogOptions& log_options,
                         const HighsModel& model,
                         const std::vector<double>& primal_solution);

void writeMipBoundGap(FILE* file, const HighsLogOptions& log_options,
                      const double mip_dual_bound, const double mip_gap);

void writeLpObjective(FILE* file, const HighsLogOptions& log_options,
                      const HighsLp& lp,
                      const std::vector<double>& primal_solution);

void writeObjectiveValue(FILE* file, const HighsLogOptions& log_options,
                         const double objective_value);

void writePrimalSolution(FILE* file, const HighsLogOptions& log_options,
                         const HighsLp& lp,
                         const std::vector<double>& primal_solution,
                         const bool sparse, const bool partial = false);

void writeModelSolution(FILE* file, const HighsLogOptions& log_options,
                        const HighsModel& model,
                        const HighsModelStatus model_status,
                        const HighsSolution& solution, const HighsInfo& info,
                        const bool sparse, const bool partial);

bool replaceSpacesByUnderscores(std::string& name);

HighsInt maxNameLength(const HighsLp& lp);
HighsInt maxNameLength(const std::vector<std::string>& names);

HighsStatus normaliseNames(const HighsLogOptions& log_options, HighsLp& lp,
                           HighsFileType type = HighsFileType::kMps);

HighsStatus normaliseNames(const HighsLogOptions& log_options, bool column,
                           HighsInt num_name_required, std::string& name_prefix,
                           HighsInt& name_suffix,
                           std::vector<std::string>& names,
                           HighsNameHash& name_hash,
                           HighsFileType type = HighsFileType::kMps);

HighsFileType getFileType(const std::string& filename);

void writeSolutionFile(FILE* file, const HighsOptions& options,
                       const HighsModel& model, const HighsBasis& basis,
                       const HighsSolution& solution, const HighsInfo& info,
                       const HighsModelStatus model_status,
                       const HighsInt style);

void writeGlpsolCostRow(FILE* file, const HighsLogOptions& log_options,
                        const bool raw, const bool is_mip,
                        const HighsInt row_id,
                        const std::string& objective_name,
                        const double objective_function_value);

void writeGlpsolSolution(FILE* file, const HighsOptions& options,
                         const HighsLp& lp, const HighsBasis& basis,
                         const HighsSolution& solution,
                         const HighsModelStatus model_status,
                         const HighsInfo& info, const bool raw);

void writeOldRawSolution(FILE* file, const HighsLogOptions& log_options,
                         const HighsLp& lp, const HighsBasis& basis,
                         const HighsSolution& solution);

HighsBasisStatus checkedVarHighsNonbasicStatus(
    const HighsBasisStatus ideal_status, const double lower,
    const double upper);

std::string utilModelStatusToString(const HighsModelStatus model_status);

std::string utilSolutionStatusToString(const HighsInt solution_status);

StatusString utilBasisStatusToString(const HighsBasisStatus basis_status);

std::string utilBasisValidityToString(const HighsInt basis_validity);

std::string utilPresolveRuleTypeToString(const HighsInt rule_type);

std::string utilPostsolveStatusToString(
    const HighsPostsolveStatus postsolve_status);

HighsStatus highsStatusFromHighsModelStatus(HighsModelStatus model_status);

std::string statusToString(const HighsBasisStatus status, const double lower,
                           const double upper);
std::string typeToString(const HighsVarType type);

std::string findModelObjectiveName(const HighsLp* lp,
                                   const HighsHessian* hessian = nullptr);

// bool repeatedNames(const std::vector<std::string> name);

#endif
