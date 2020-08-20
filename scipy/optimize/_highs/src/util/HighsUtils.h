/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HUtils.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSUTILS_H_
#define UTIL_HIGHSUTILS_H_

#include <cassert>
#include <string>
#include <vector>

#include "HConfig.h"
#include "lp_data/HighsOptions.h"

struct HighsIndexCollection {
  int dimension_ = -1;
  bool is_interval_ = false;
  int from_ = -1;
  int to_ = -2;
  bool is_set_ = false;
  int set_num_entries_ = -1;
  int* set_ = NULL;
  bool is_mask_ = false;
  int* mask_ = NULL;
};

#ifdef HiGHSDEV
struct HighsValueDistribution {
  std::string distribution_name_;
  std::string value_name_;
  int num_count_;
  int num_zero_;
  int num_one_;
  double min_value_;
  double max_value_;
  std::vector<double> limit_;
  std::vector<int> count_;
  int sum_count_;
};
#endif

struct HighsScatterData {
  int max_num_point_;
  int num_point_;
  int last_point_;
  std::vector<double> value0_;
  std::vector<double> value1_;
  bool have_regression_coeff_;
  double linear_coeff0_;
  double linear_coeff1_;
  double linear_regression_error_;
  double log_coeff0_;
  double log_coeff1_;
  double log_regression_error_;
  int num_error_comparison_;
  int num_awful_linear_;
  int num_awful_log_;
  int num_bad_linear_;
  int num_bad_log_;
  int num_fair_linear_;
  int num_fair_log_;
  int num_better_linear_;
  int num_better_log_;
};

const double awful_regression_error = 2.0;
const double bad_regression_error = 0.2;
const double fair_regression_error = 0.02;

bool assessIndexCollection(const HighsOptions& options,
                           const HighsIndexCollection& index_collection);

bool limitsForIndexCollection(const HighsOptions& options,
                              const HighsIndexCollection& index_collection,
                              int& from_k, int& to_k);

void updateIndexCollectionOutInIndex(
    const HighsIndexCollection& index_collection, int& out_from_ix,
    int& out_to_ix, int& in_from_ix, int& in_to_ix, int& current_set_entry);

int dataSizeOfIndexCollection(const HighsIndexCollection& index_collection);

bool intUserDataNotNull(FILE* logfile, const int* user_data,
                        const std::string name);
bool doubleUserDataNotNull(FILE* logfile, const double* user_data,
                           const std::string name);

double getNorm2(const std::vector<double> values);

/**
 * @brief Logical check of double being +Infinity
 */
bool highs_isInfinity(double val  //!< Value being tested against +Infinity
);
/**
 * @brief Returns the relative difference of two doubles
 */
double highsRelativeDifference(const double v0, const double v1);

#ifdef HiGHSDEV
/**
 * @brief Analyse the values of a vector, assessing how many are in
 * each power of ten, and possibly analyse the distribution of
 * different values
 */
void analyseVectorValues(
    const char* message,             //!< Message to be printed
    int vecDim,                      //!< Dimension of vector
    const std::vector<double>& vec,  //!< Vector of values
    bool analyseValueList = false,   //!< Possibly analyse the distribution of
                                     //!< different values in the vector
    std::string model_name =
        "Unknown"  //!< Model name to report if analysing distribution of
                   //!< different values in the vector
);

void analyseMatrixSparsity(
    const char* message,             //!< Message to be printed
    int numCol,                      //!< Number of columns
    int numRow,                      //!< Number of rows
    const std::vector<int>& Astart,  //!< Matrix column starts
    const std::vector<int>& Aindex   //!< Matrix row indices
);

bool initialiseValueDistribution(const std::string distribution_name,
                                 const std::string value_name,
                                 const double min_value_limit,
                                 const double max_value_limit,
                                 const double base_value_limit,
                                 HighsValueDistribution& value_distribution);

bool updateValueDistribution(const double value,
                             HighsValueDistribution& value_distribution);

int integerPercentage(const int of, const int in);
double doublePercentage(const int of, const int in);

bool printValueDistribution(const HighsValueDistribution& value_distribution,
                            const int mu = 0);
#endif

bool initialiseScatterData(const int max_num_point,
                           HighsScatterData& scatter_data);
bool updateScatterData(const double value0, const double value1,
                       HighsScatterData& scatter_data);
bool regressScatterData(HighsScatterData& scatter_data);
bool predictFromScatterData(const HighsScatterData& scatter_data,
                            const double value0, double& predicted_value1,
                            const bool log_regression = false);
bool printScatterData(std::string name, const HighsScatterData& scatter_data);
void printScatterDataRegressionComparison(std::string name,
                                          const HighsScatterData& scatter_data);
bool computeScatterDataRegressionError(HighsScatterData& scatter_data,
                                       const bool print = false);

#endif  // UTIL_HIGHSUTILS_H_
