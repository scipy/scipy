/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HMpsFF.h
 * @brief
 */
#ifndef IO_HMPSFF_H_
#define IO_HMPSFF_H_

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "io/HighsIO.h"
#include "model/HighsModel.h"
// #include "util/HighsType.h"
#include "util/stringutil.h"

using Triplet = std::tuple<HighsInt, HighsInt, double>;

const std::string mps_comment_chars = "*$";

enum class FreeFormatParserReturnCode {
  kSuccess,
  kParserError,
  kFileNotFound,
  kFixedFormat,
  kTimeout,
};

namespace free_format_parser {

// private:
using wall_clock = std::chrono::high_resolution_clock;
using time_point = wall_clock::time_point;

double getWallTime();

class HMpsFF {
 public:
  HMpsFF() {}
  FreeFormatParserReturnCode loadProblem(const HighsLogOptions& log_options,
                                         const std::string& filename,
                                         HighsModel& model);

  double time_limit_ = kHighsInf;
  bool warning_issued_ = false;

 private:
  double start_time;

  HighsInt num_row;
  HighsInt num_col;
  HighsInt num_nz;
  std::string mps_name;

  ObjSense obj_sense = ObjSense::kMinimize;  // Minimization by default
  double obj_offset = 0;

  std::vector<HighsInt> a_start;
  std::vector<HighsInt> a_index;
  std::vector<double> a_value;
  std::vector<double> col_cost;
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  std::vector<double> row_lower;
  std::vector<double> row_upper;

  std::vector<std::string> row_names;
  std::vector<std::string> col_names;

  std::vector<HighsVarType> col_integrality;

  HighsInt q_dim;
  std::vector<HighsInt> q_start;
  std::vector<HighsInt> q_index;
  std::vector<double> q_value;

  // Keep track of columns that are binary by default, being columns
  // that are defined as integer by markers in the column section, or
  // as binary by having a BV flag in the BOUNDS section, and without
  // any LI or UI flags in the BOUNDS section
  std::vector<HighsBool> col_binary;

  // Record where the cost row is encountered
  HighsInt cost_row_location;

  // Record whether there are duplicate row or column names, and the
  // name and indices of the first duplicates
  bool has_duplicate_row_name_;
  bool has_duplicate_col_name_;
  std::string duplicate_row_name_;
  HighsInt duplicate_row_name_index0_;
  HighsInt duplicate_row_name_index1_;
  std::string duplicate_col_name_;
  HighsInt duplicate_col_name_index0_;
  HighsInt duplicate_col_name_index1_;

  // Record whether there is a data entry in the RHS section of an MPS
  // file for the objective or a row. Have to be class data members so
  // that they can be used by parseName and addRhs in HMpsFF::parseRhs
  bool has_obj_entry_;
  std::vector<HighsBool> has_row_entry_;

  /// load LP from MPS file as transposed triplet matrix
  HighsInt parseFile(std::string filename);
  HighsInt fillMatrix(const HighsLogOptions& log_options);
  HighsInt fillHessian(const HighsLogOptions& log_options);

  /// how to treat variables that appear in COLUMNS section first
  /// assume them to be binary as in the original IBM interpretation
  /// or integer with default bounds
  bool kintegerVarsInColumnsAreBinary = true;

  enum class Parsekey {
    kName,
    kObjsense,
    kMax,
    kMin,
    kRows,
    kCols,
    kRhs,
    kBounds,
    kRanges,
    kQsection,
    kQmatrix,
    kQuadobj,
    kQcmatrix,
    kCsection,
    kDelayedrows,
    kModelcuts,
    kUsercuts,
    kIndicators,
    kSets,
    kSos,
    kGencons,
    kPwlobj,
    kPwlnam,
    kPwlcon,
    kNone,
    kEnd,
    kFail,
    kComment,
    kFixedFormat,
    kTimeout
  };

  enum class Boundtype { kLe, kEq, kGe };  //, kFr };

  // see
  // https://docs.mosek.com/latest_kkt/capi/mps-format.html#csection-optional
  enum class ConeType { kZero, kQuad, kRQuad, kPExp, kPPow, kDExp, kDPow };

  std::string objective_name;
  std::vector<Boundtype> row_type;
  std::vector<HighsInt> integer_column;

  std::vector<Triplet> entries;
  HessianFormat q_format;
  std::vector<Triplet> q_entries;
  std::vector<std::vector<Triplet>> qrows_entries;
  std::vector<std::pair<HighsInt, double>> coeffobj;

  std::vector<std::string> sos_name;
  std::vector<short> sos_type;
  std::vector<std::vector<std::pair<HighsInt, double>>> sos_entries;

  std::vector<std::string> cone_name;
  std::vector<ConeType> cone_type;
  std::vector<double> cone_param;
  std::vector<std::vector<HighsInt>> cone_entries;
  std::unordered_map<std::string, HighsInt> rowname2idx;
  std::unordered_map<std::string, HighsInt> colname2idx;

  mutable std::string section_args;

  bool timeout();
  bool getMpsLine(std::istream& file, std::string& strline, bool& skip);

  FreeFormatParserReturnCode parse(const HighsLogOptions& log_options,
                                   const std::string& filename);
  // Checks first word of strline and wraps it by it_begin and it_end
  HMpsFF::Parsekey checkFirstWord(std::string& strline, size_t& start,
                                  size_t& end, std::string& word) const;

  // Get index of column from column name, possibly adding new column
  // if no index is found
  HighsInt getColIdx(const std::string& colname, const bool add_if_new = true);

  HMpsFF::Parsekey parseDefault(const HighsLogOptions& log_options,
                                std::istream& file);
  HMpsFF::Parsekey parseObjsense(const HighsLogOptions& log_options,
                                 std::istream& file);
  HMpsFF::Parsekey parseRows(const HighsLogOptions& log_options,
                             std::istream& file);
  HMpsFF::Parsekey parseCols(const HighsLogOptions& log_options,
                             std::istream& file);
  HMpsFF::Parsekey parseRhs(const HighsLogOptions& log_options,
                            std::istream& file);
  HMpsFF::Parsekey parseRanges(const HighsLogOptions& log_options,
                               std::istream& file);
  HMpsFF::Parsekey parseBounds(const HighsLogOptions& log_options,
                               std::istream& file);
  HMpsFF::Parsekey parseHessian(const HighsLogOptions& log_options,
                                std::istream& file,
                                const HMpsFF::Parsekey keyword);
  HMpsFF::Parsekey parseQuadRows(const HighsLogOptions& log_options,
                                 std::istream& file,
                                 const HMpsFF::Parsekey keyword);
  HMpsFF::Parsekey parseQuadMatrix(const HighsLogOptions& log_options,
                                   std::istream& file,
                                   const std::string& section_name,
                                   std::vector<Triplet>& q_entries);

  HMpsFF::Parsekey parseCones(const HighsLogOptions& log_options,
                              std::istream& file);
  HMpsFF::Parsekey parseSos(const HighsLogOptions& log_options,
                            std::istream& file, const HMpsFF::Parsekey keyword);

  bool cannotParseSection(const HighsLogOptions& log_options,
                          const HMpsFF::Parsekey keyword);
  bool allZeroed(const std::vector<double>& value);
  double getValue(const std::string& word, bool& is_nan,
                  const HighsInt id = -1) const;
};

}  // namespace free_format_parser
#endif /* IO_HMPSFF_H_ */
