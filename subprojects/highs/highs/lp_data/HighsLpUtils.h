/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLpUtils.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSLPUTILS_H_
#define LP_DATA_HIGHSLPUTILS_H_

#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsUtils.h"

// class HighsLp;
struct SimplexScale;
struct HighsBasis;
struct HighsSolution;
class HighsOptions;

using std::vector;

void writeBasisFile(FILE*& file, const HighsOptions& options, const HighsLp& lp,
                    const HighsBasis& basis);

HighsStatus getIndexFromName(const HighsLogOptions& log_options,
                             std::string& from_method, const bool is_column,
                             const std::string& name, HighsNameHash& name_hash,
                             HighsInt& index,
                             const std::vector<std::string>& names);

HighsStatus readBasisFile(const HighsLogOptions& log_options, HighsLp& lp,
                          HighsBasis& basis, const std::string& filename);
HighsStatus readBasisStream(const HighsLogOptions& log_options, HighsLp& lp,
                            HighsBasis& basis, std::ifstream& in_file);

// Methods taking HighsLp as an argument
HighsStatus assessLp(HighsLp& lp, const HighsOptions& options);

HighsStatus assessCosts(const HighsOptions& options, const HighsInt ml_col_os,
                        const HighsIndexCollection& index_collection,
                        vector<double>& cost, bool& has_infinite_cost,
                        const double infinite_cost);

HighsStatus assessBounds(const HighsOptions& options, const std::string& type,
                         const HighsInt ml_ix_os,
                         const HighsIndexCollection& index_collection,
                         vector<double>& lower, vector<double>& upper,
                         const double infinite_bound,
                         const HighsVarType* integrality,
                         const std::string* names = nullptr);

HighsStatus cleanBounds(const HighsOptions& options, HighsLp& lp);

HighsStatus userScaleLp(HighsLp& lp, HighsUserScaleData& data,
                        const HighsLogOptions& log_options);

void userScaleLp(HighsLp& lp, HighsUserScaleData& data,
                 const bool apply = true);

void userScaleCosts(const vector<HighsVarType>& integrality, double& offset,
                    vector<double>& cost, HighsUserScaleData& data,
                    const bool apply = true);

void userScaleColBounds(const vector<HighsVarType>& integrality,
                        vector<double>& lower, vector<double>& upper,
                        HighsUserScaleData& data, const bool apply = true);

void userScaleRowBounds(vector<double>& lower, vector<double>& upper,
                        HighsUserScaleData& data, const bool apply = true);

void userScaleMatrix(const vector<HighsVarType>& integrality,
                     HighsSparseMatrix& matrix, HighsUserScaleData& data,
                     const bool apply = true);

HighsStatus userScaleStatus(const HighsLogOptions& log_options,
                            const HighsUserScaleData& data);

HighsStatus assessSemiVariables(HighsLp& lp, const HighsOptions& options,
                                bool& made_semi_variable_mods);
void relaxSemiVariables(HighsLp& lp, bool& made_semi_variable_mods);

bool activeModifiedUpperBounds(const HighsOptions& options, const HighsLp& lp,
                               const std::vector<double>& col_value);

bool considerScaling(const HighsOptions& options, HighsLp& lp);
void scaleLp(const HighsOptions& options, HighsLp& lp,
             const bool force_scaling = false);
bool equilibrationScaleMatrix(const HighsOptions& options, HighsLp& lp,
                              const HighsInt use_scale_strategy);
bool maxValueScaleMatrix(const HighsOptions& options, HighsLp& lp,
                         const HighsInt use_scale_strategy);

HighsStatus applyScalingToLpCol(HighsLp& lp, const HighsInt col,
                                const double colScale);

HighsStatus applyScalingToLpRow(HighsLp& lp, const HighsInt row,
                                const double rowScale);

void unscaleSolution(HighsSolution& solution, const HighsScale& scale);

void appendColsToLpVectors(HighsLp& lp, const HighsInt num_new_col,
                           const vector<double>& colCost,
                           const vector<double>& colLower,
                           const vector<double>& colUpper);

void appendRowsToLpVectors(HighsLp& lp, const HighsInt num_new_row,
                           const vector<double>& rowLower,
                           const vector<double>& rowUpper);

void deleteScale(vector<double>& scale,
                 const HighsIndexCollection& index_collection);

HighsStatus changeLpMatrixCoeff(const HighsLogOptions& log_options, HighsLp& lp,
                                const HighsInt row, const HighsInt col,
                                const double new_value,
                                const double small_matrix_value);

HighsStatus changeLpIntegrality(HighsLp& lp,
                                const HighsIndexCollection& index_collection,
                                const vector<HighsVarType>& new_integrality,
                                const HighsOptions& options);

void changeLpCosts(HighsLp& lp, const HighsIndexCollection& index_collection,
                   const vector<double>& new_col_cost,
                   const double infinite_cost);

void changeLpColBounds(HighsLp& lp,
                       const HighsIndexCollection& index_collection,
                       const vector<double>& new_col_lower,
                       const vector<double>& new_col_upper);

void changeLpRowBounds(HighsLp& lp,
                       const HighsIndexCollection& index_collection,
                       const vector<double>& new_row_lower,
                       const vector<double>& new_row_upper);

void changeBounds(vector<double>& lower, vector<double>& upper,
                  const HighsIndexCollection& index_collection,
                  const vector<double>& new_lower,
                  const vector<double>& new_upper);

/**
 * @brief Report the data of an LP
 */
void reportLp(const HighsLogOptions& log_options,
              const HighsLp& lp,  //!< LP whose data are to be reported
              const HighsLogType report_level = HighsLogType::kInfo
              //!< INFO => scalar [dimensions];
              //!< DETAILED => vector[costs/bounds];
              //!< VERBOSE => vector+matrix
);
/**
 * @brief Report the brief data of an LP
 */
void reportLpBrief(const HighsLogOptions& log_options,
                   const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpDimensions(const HighsLogOptions& log_options,
                        const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpObjSense(const HighsLogOptions& log_options,
                      const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpColVectors(const HighsLogOptions& log_options,
                        const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpRowVectors(const HighsLogOptions& log_options,
                        const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpColMatrix(const HighsLogOptions& log_options,
                       const HighsLp& lp  //!< LP whose data are to be reported
);

void reportMatrix(const HighsLogOptions& log_options,
                  const std::string& message, const HighsInt num_col,
                  const HighsInt num_nz, const HighsInt* start,
                  const HighsInt* index, const double* value);

// Get the number of integer-valued columns in the LP
HighsInt getNumInt(const HighsLp& lp);

// Get the costs for a contiguous set of columns
void getLpCosts(const HighsLp& lp, const HighsInt from_col,
                const HighsInt to_col, double* XcolCost);

// Get the bounds for a contiguous set of columns
void getLpColBounds(const HighsLp& lp, const HighsInt from_col,
                    const HighsInt to_col, double* XcolLower,
                    double* XcolUpper);

// Get the bounds for a contiguous set of rows
void getLpRowBounds(const HighsLp& lp, const HighsInt from_row,
                    const HighsInt to_row, double* XrowLower,
                    double* XrowUpper);

void getLpMatrixCoefficient(const HighsLp& lp, const HighsInt row,
                            const HighsInt col, double* val);
// Analyse the data in an LP problem
void analyseLp(const HighsLogOptions& log_options, const HighsLp& lp);

HighsStatus readSolutionFile(const std::string& filename,
                             const HighsOptions& options, HighsLp& lp,
                             HighsBasis& basis, HighsSolution& solution,
                             const HighsInt style);

HighsStatus readSolutionFileErrorReturn(std::ifstream& in_file);
HighsStatus readSolutionFileReturn(const HighsStatus status,
                                   HighsSolution& solution, HighsBasis& basis,
                                   const HighsSolution& read_solution,
                                   const HighsBasis& read_basis,
                                   std::ifstream& in_file);
bool readSolutionFileIgnoreLineOk(std::ifstream& in_file);
bool readSolutionFileKeywordLineOk(std::string& keyword,
                                   std::ifstream& in_file);
bool readSolutionFileHashKeywordIntLineOk(std::string& hash,
                                          std::string& keyword,
                                          std::string& value_string,
                                          HighsInt& value,
                                          std::ifstream& in_file);
bool readSolutionFileColumnHeaderLineOk(
    std::string& hash, std::string& keyword, std::string& value_string,
    HighsInt& value, std::string& qualifier_string,
    std::stringstream& column_section_line_ss);
bool readSolutionFileColumnLineOk(std::string& id, double& value,
                                  HighsInt& index,
                                  std::stringstream& column_section_line_ss);
bool readSolutionFileIdIgnoreLineOk(std::string& id, std::ifstream& in_file);
bool readSolutionFileIdDoubleLineOk(std::string& id, double& value,
                                    std::ifstream& in_file);

void assessColPrimalSolution(const bool assess_as_mip,
                             const HighsOptions& options, const double primal,
                             const double lower, const double upper,
                             const HighsVarType type, double& col_infeasibility,
                             double& integer_infeasibility);

HighsStatus assessLpPrimalSolution(const std::string& message,
                                   const HighsOptions& options,
                                   const HighsLp& lp,
                                   const HighsSolution& solution, bool& valid,
                                   bool& integral, bool& feasible);

HighsStatus calculateRowValuesQuad(const HighsLp& lp,
                                   const std::vector<double>& col_value,
                                   std::vector<double>& row_value,
                                   const HighsInt report_row = -1);
HighsStatus calculateRowValuesQuad(const HighsLp& lp, HighsSolution& solution,
                                   const HighsInt report_row = -1);

HighsStatus calculateColDualsQuad(const HighsLp& lp, HighsSolution& solution);

bool isColDataNull(const HighsLogOptions& log_options,
                   const double* usr_col_cost, const double* usr_col_lower,
                   const double* usr_col_upper);
bool isRowDataNull(const HighsLogOptions& log_options,
                   const double* usr_row_lower, const double* usr_row_upper);
bool isMatrixDataNull(const HighsLogOptions& log_options,
                      const HighsInt* usr_matrix_start,
                      const HighsInt* usr_matrix_index,
                      const double* usr_matrix_value);

void reportPresolveReductions(const HighsLogOptions& log_options,
                              HighsPresolveStatus presolve_status,
                              const HighsLp& lp, const HighsLp& presolved_lp);

bool isLessInfeasibleDSECandidate(const HighsLogOptions& log_options,
                                  const HighsLp& lp);

HighsLp withoutSemiVariables(const HighsLp& lp, HighsSolution& solution,
                             const double mip_feasibility_tolerance);

void removeRowsOfCountOne(const HighsLogOptions& log_options, HighsLp& lp);

// Get subvectors from data structure of data0, data1, data2 and
// matrix, where the storage of the matrix is compatible with the
// vectors to be extracted from it
//
// Data to be extracted is given by sub_* being nullpointer
//
// * cost, lower and upper bounds for columns, and column-wise LP constraint
// matrix
//
// * lower and upper bounds for rows, and row-wise LP constraint
// * matrix. "cost" is nullptr, and so must be sub_vector_data0
//
void getSubVectors(const HighsIndexCollection& index_collection,
                   const HighsInt data_dim, const double* data0,
                   const double* data1, const double* data2,
                   const HighsSparseMatrix& matrix, HighsInt& num_sub_vector,
                   double* sub_vector_data0, double* sub_vector_data1,
                   double* sub_vector_data2, HighsInt& sub_matrix_num_nz,
                   HighsInt* sub_matrix_start, HighsInt* sub_matrix_index,
                   double* sub_matrix_value);

void getSubVectorsTranspose(const HighsIndexCollection& index_collection,
                            const HighsInt data_dim, const double* data0,
                            const double* data1, const double* data2,
                            const HighsSparseMatrix& matrix,
                            HighsInt& num_sub_vector, double* sub_vector_data0,
                            double* sub_vector_data1, double* sub_vector_data2,
                            HighsInt& sub_matrix_num_nz,
                            HighsInt* sub_matrix_start,
                            HighsInt* sub_matrix_index,
                            double* sub_matrix_value);

std::string highsVarTypeToString(const HighsVarType type);
std::string highsVarTypeToString(const HighsInt type);

void initialiseUserScaleData(const HighsOptions& options,
                             HighsUserScaleData& user_scale_data);

#endif  // LP_DATA_HIGHSLPUTILS_H_
