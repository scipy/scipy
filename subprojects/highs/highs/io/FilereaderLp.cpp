/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderLp.cpp
 * @brief
 */

#include "io/FilereaderLp.h"

#include <array>
#include <cstdarg>
#include <cstdio>

#include "io/LpReader.h"
#include "lp_data/HighsLpUtils.h"

FilereaderRetcode FilereaderLp::readModelFromFile(const HighsOptions& options,
                                                  const std::string filename,
                                                  HighsModel& model) {
  return readLpFile(options.log_options, filename, model);
}

void FilereaderLp::writeToFile(FILE* file, const char* format, ...) {
  va_list argptr;
  va_start(argptr, format);
  std::array<char, LP_MAX_LINE_LENGTH + 1> stringbuffer = {};
  HighsInt tokenlength =
      vsnprintf(stringbuffer.data(), stringbuffer.size(), format, argptr);
  va_end(argptr);
  if (this->linelength + tokenlength >= LP_MAX_LINE_LENGTH) {
    fprintf(file, "\n");
    fprintf(file, "%s", stringbuffer.data());
    this->linelength = tokenlength;
  } else {
    fprintf(file, "%s", stringbuffer.data());
    this->linelength += tokenlength;
  }
}

void FilereaderLp::writeToFileLineEnd(FILE* file) {
  fprintf(file, "\n");
  this->linelength = 0;
}

void FilereaderLp::writeToFileValue(FILE* file, const double value,
                                    const bool force_plus) {
  // As for writeModelAsMps
  if (force_plus) {
    this->writeToFile(file, " %+.17g", value);
  } else {
    this->writeToFile(file, " %.17g", value);
  }
}

void FilereaderLp::writeToFileVar(FILE* file, const std::string& var_name) {
  this->writeToFile(file, " %s", var_name.c_str());
}

void FilereaderLp::writeToFileMatrixRow(FILE* file, const HighsInt iRow,
                                        const HighsSparseMatrix& ar_matrix,
                                        const std::vector<string>& col_names) {
  assert(ar_matrix.isRowwise());

  for (HighsInt iEl = ar_matrix.start_[iRow]; iEl < ar_matrix.start_[iRow + 1];
       iEl++) {
    HighsInt iCol = ar_matrix.index_[iEl];
    double coef = ar_matrix.value_[iEl];
    this->writeToFileValue(file, coef);
    this->writeToFileVar(file, col_names[iCol]);
  }
}

HighsStatus FilereaderLp::writeModelToFile(const HighsOptions& options,
                                           const std::string filename,
                                           const HighsModel& model) {
  const HighsLp& lp = model.lp_;

  const bool ok_names = lp.okNames();
  assert(ok_names);
  if (!ok_names) return HighsStatus::kError;

  // Create a row-wise copy of the matrix
  HighsSparseMatrix ar_matrix = lp.a_matrix_.getRowwise();
  assert(ar_matrix.isRowwise());
  // ar_matrix.ensureRowwise(); is now redundant

  FILE* file = fopen(filename.c_str(), "w");
  if (file == 0) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Cannot open file %s\n", filename.c_str());
    return HighsStatus::kError;
  }
  // write comment at the start of the file
  this->writeToFile(file, "\\ %s", LP_COMMENT_FILESTART);
  this->writeToFileLineEnd(file);

  // write objective
  this->writeToFile(file, "%s",
                    lp.sense_ == ObjSense::kMinimize ? "min" : "max");
  this->writeToFileLineEnd(file);
  this->writeToFile(file, " obj:");
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double coef = lp.col_cost_[iCol];
    if (coef != 0.0) {
      this->writeToFileValue(file, coef);
      this->writeToFileVar(file, lp.col_names_[iCol]);
    }
  }
  this->writeToFile(file,
                    " ");  // ToDo Unnecessary, but only to give empty diff
  if (model.isQp()) {
    this->writeToFile(file, "+ [");
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      for (HighsInt iEl = model.hessian_.start_[iCol];
           iEl < model.hessian_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = model.hessian_.index_[iEl];
        if (iCol <= iRow) {
          double coef = model.hessian_.value_[iEl];
          if (iCol != iRow) coef *= 2;
          if (coef != 0.0) {
            this->writeToFileValue(file, coef);
            this->writeToFileVar(file, lp.col_names_[iCol]);
            this->writeToFile(file, " *");
            this->writeToFileVar(file, lp.col_names_[iRow]);
          }
        }
      }
    }
    this->writeToFile(file,
                      "  ]/2 ");  // ToDo Surely needs only to be one space
  }
  double coef = lp.offset_;
  if (coef != 0) this->writeToFileValue(file, coef);
  this->writeToFileLineEnd(file);

  // write constraint section, lower & upper bounds are one constraint
  // each
  this->writeToFile(file, "st");
  this->writeToFileLineEnd(file);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (lp.row_lower_[iRow] == lp.row_upper_[iRow]) {
      // Equality constraint
      this->writeToFileVar(file, lp.row_names_[iRow]);
      this->writeToFile(file, ":");
      this->writeToFileMatrixRow(file, iRow, ar_matrix, lp.col_names_);
      this->writeToFile(file, " =");
      this->writeToFileValue(file, lp.row_lower_[iRow], true);
      this->writeToFileLineEnd(file);
    } else {
      // Need to distinguish the names when writing out boxed
      // constraint row as two single-sided constraints
      const bool boxed =
          lp.row_lower_[iRow] > -kHighsInf && lp.row_upper_[iRow] < kHighsInf;
      if (lp.row_lower_[iRow] > -kHighsInf) {
        // Has a lower bound
        this->writeToFileVar(file, lp.row_names_[iRow]);
        if (boxed) {
          this->writeToFile(file, "lo:");
        } else {
          this->writeToFile(file, ":");
        }
        this->writeToFileMatrixRow(file, iRow, ar_matrix, lp.col_names_);
        this->writeToFile(file, " >=");
        this->writeToFileValue(file, lp.row_lower_[iRow], true);
        this->writeToFileLineEnd(file);
      }
      if (lp.row_upper_[iRow] < kHighsInf) {
        // Has an upper bound
        this->writeToFileVar(file, lp.row_names_[iRow]);
        if (boxed) {
          this->writeToFile(file, "up:");
        } else {
          this->writeToFile(file, ":");
        }
        this->writeToFileMatrixRow(file, iRow, ar_matrix, lp.col_names_);
        this->writeToFile(file, " <=");
        this->writeToFileValue(file, lp.row_upper_[iRow], true);
        this->writeToFileLineEnd(file);
      }
    }
  }

  // write bounds section
  this->writeToFile(file, "bounds");
  this->writeToFileLineEnd(file);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const bool default_bounds =
        lp.col_lower_[iCol] == 0 && lp.col_upper_[iCol] == kHighsInf;
    if (default_bounds) continue;
    if (lp.col_lower_[iCol] <= -kHighsInf && lp.col_upper_[iCol] >= kHighsInf) {
      // Free variable
      this->writeToFileVar(file, lp.col_names_[iCol]);
      this->writeToFile(file, " free");
    } else if (lp.col_lower_[iCol] == lp.col_upper_[iCol]) {
      // Fixed variable
      this->writeToFileVar(file, lp.col_names_[iCol]);
      this->writeToFile(file, " =");
      this->writeToFileValue(file, lp.col_upper_[iCol], false);
    } else {
      assert(!default_bounds);
      // Non-default bound
      if (lp.col_lower_[iCol] != 0) {
        // Nonzero lower bound
        this->writeToFileValue(file, lp.col_lower_[iCol], false);
        this->writeToFile(file, " <=");
      }
      this->writeToFileVar(file, lp.col_names_[iCol]);
      if (lp.col_upper_[iCol] < kHighsInf) {
        // Finite upper bound
        this->writeToFile(file, " <=");
        this->writeToFileValue(file, lp.col_upper_[iCol], false);
      }
    }
    this->writeToFileLineEnd(file);
  }
  if (lp.integrality_.size() > 0) {
    // write binary section
    this->writeToFile(file, "bin");
    this->writeToFileLineEnd(file);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] == HighsVarType::kInteger) {
        if (lp.col_lower_[iCol] == 0.0 && lp.col_upper_[iCol] == 1.0) {
          this->writeToFileVar(file, lp.col_names_[iCol]);
          this->writeToFileLineEnd(file);
        }
      }
    }

    // write general section
    this->writeToFile(file, "gen");
    this->writeToFileLineEnd(file);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] == HighsVarType::kInteger) {
        if (lp.col_lower_[iCol] != 0.0 || lp.col_upper_[iCol] != 1.0) {
          this->writeToFileVar(file, lp.col_names_[iCol]);
          this->writeToFileLineEnd(file);
        }
      } else if (lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
        this->writeToFileVar(file, lp.col_names_[iCol]);
        this->writeToFileLineEnd(file);
      }
    }

    // write semi section
    this->writeToFile(file, "semi");
    this->writeToFileLineEnd(file);
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
          lp.integrality_[iCol] == HighsVarType::kSemiInteger) {
        this->writeToFileVar(file, lp.col_names_[iCol]);
        this->writeToFileLineEnd(file);
      }
    }
  }
  // write end
  this->writeToFile(file, "end");
  this->writeToFileLineEnd(file);

  fclose(file);
  return HighsStatus::kOk;
}
