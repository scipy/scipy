/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HMPSIO.cpp
 * @brief
 */
#include "io/HMPSIO.h"

#include <algorithm>
#include <cstdio>

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsOptions.h"
#include "util/HighsUtils.h"
#include "util/stringutil.h"

#ifdef ZLIB_FOUND
#include "../extern/zstr/zstr.hpp"
#endif

using std::map;

//
// Read file called filename. Returns 0 if OK and 1 if file can't be opened
//
FilereaderRetcode readMps(
    const HighsLogOptions& log_options, const std::string& filename,
    HighsInt mxNumRow, HighsInt mxNumCol, HighsInt& numRow, HighsInt& numCol,
    ObjSense& objSense, double& objOffset, vector<HighsInt>& Astart,
    vector<HighsInt>& Aindex, vector<double>& Avalue, vector<double>& colCost,
    vector<double>& colLower, vector<double>& colUpper,
    vector<double>& rowLower, vector<double>& rowUpper,
    vector<HighsVarType>& integerColumn, std::string& objective_name,
    vector<std::string>& col_names, vector<std::string>& row_names,
    HighsInt& Qdim, vector<HighsInt>& Qstart, vector<HighsInt>& Qindex,
    vector<double>& Qvalue, HighsInt& cost_row_location, bool& warning_issued,
    const HighsInt keep_n_rows) {
  // Keep track of any warnings that are issued so that
  // Highs::readModel can return HighsStatus::kWarning
  warning_issued = false;
  // MPS file buffer
  numRow = 0;
  numCol = 0;
  Qdim = 0;
  cost_row_location = -1;
  objOffset = 0;
  objSense = ObjSense::kMinimize;
  objective_name = "";

  // Astart.clear() added since setting Astart.push_back(0) in
  // setup_clearModel() messes up the MPS read
  Astart.clear();
  highsLogDev(log_options, HighsLogType::kInfo,
              "readMPS: Trying to open file %s\n", filename.c_str());
#ifdef ZLIB_FOUND
  zstr::ifstream file;
  try {
    file.open(filename, std::ios::in);
  } catch (const strict_fstream::Exception& e) {
    highsLogDev(log_options, HighsLogType::kInfo, e.what());
    return FilereaderRetcode::kFileNotFound;
  }
#else
  std::ifstream file;
  file.open(filename, std::ios::in);
#endif
  if (!file.is_open()) {
    highsLogDev(log_options, HighsLogType::kInfo,
                "readMPS: Not opened file OK\n");
    return FilereaderRetcode::kFileNotFound;
  }
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Opened file  OK\n");
  // Input buffer
  const HighsInt lmax = 128;
  std::array<char, lmax> line;
  std::array<char, 2> flag = {0, 0};
  std::array<double, 3> data;

  HighsInt num_alien_entries = 0;
  HighsVarType integerCol = HighsVarType::kContinuous;

  // Load NAME
  load_mpsLine(file, integerCol, lmax, line.data(), flag.data(), data.data());
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read NAME    OK\n");
  // Load OBJSENSE or ROWS
  load_mpsLine(file, integerCol, lmax, line.data(), flag.data(), data.data());
  if (flag[0] == 'O') {
    // Found OBJSENSE
    load_mpsLine(file, integerCol, lmax, line.data(), flag.data(), data.data());
    std::string sense(&line[2], &line[2] + 3);
    // the sense must be "MAX" or "MIN"
    if (sense.compare("MAX") == 0) {
      objSense = ObjSense::kMaximize;
    } else if (sense.compare("MIN") == 0) {
      objSense = ObjSense::kMinimize;
    } else {
      return FilereaderRetcode::kParserError;
    }
    highsLogDev(log_options, HighsLogType::kInfo,
                "readMPS: Read OBJSENSE OK\n");
    // Load ROWS
    load_mpsLine(file, integerCol, lmax, line.data(), flag.data(), data.data());
  }

  row_names.clear();
  col_names.clear();
  vector<char> rowType;
  map<double, int> rowIndex;
  double objName = 0;
  while (load_mpsLine(file, integerCol, lmax, line.data(), flag.data(),
                      data.data())) {
    if (flag[0] == 'N' &&
        (objName == 0 || keep_n_rows == kKeepNRowsDeleteRows)) {
      // N-row: take the first as the objective and possibly ignore any others
      if (objName == 0) {
        objName = data[1];
        std::string name(&line[4], &line[4] + 8);
        objective_name = trim(name);
        cost_row_location = numRow;
      }
    } else {
      if (mxNumRow > 0 && numRow >= mxNumRow)
        return FilereaderRetcode::kParserError;
      rowType.push_back(flag[0]);
      // rowIndex is used to get the row index from a row name in the
      // COLUMNS, RHS and RANGES section. However, if this contains a
      // reference to a row that isn't in the ROWS section the value
      // of rowIndex is zero. Unless the value associated with the
      // name in rowIndex is one more than the index of the row, this
      // return of zero leads to data relating to row 0 being
      // over-written and (generally) corrupted.
      rowIndex[data[1]] = ++numRow;
      std::string name(&line[4], &line[4] + 8);
      name = trim(name);
      row_names.push_back(name);
    }
  }
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read ROWS    OK\n");

  // Load COLUMNS
  map<double, int> colIndex;
  double lastName = 0;
  // flag[1] is used to indicate whether there is more to read on the
  // line - field 5 non-empty. save_flag1 is used to deduce whether
  // the row name and value are from fields 5 and 6, or 3 and 4
  HighsInt save_flag1 = 0;
  while (load_mpsLine(file, integerCol, lmax, line.data(), flag.data(),
                      data.data())) {
    HighsInt iRow = rowIndex[data[2]] - 1;
    std::string name = "";
    if (iRow >= 0) name = row_names[iRow];
    if (lastName != data[1]) {  // New column
      if (mxNumCol > 0 && numCol >= mxNumCol)
        return FilereaderRetcode::kParserError;
      lastName = data[1];
      // colIndex is used to get the column index from a column name
      // in the BOUNDS section. However, if this contains a reference
      // to a column that isn't in the COLUMNS section the value of
      // colIndex is zero. Unless the value associated with the name
      // in colIndex is one more than the index of the column, this
      // return of zero leads to the bounds on column 0 being
      // over-written and (generally) corrupted.
      colIndex[data[1]] = ++numCol;
      colCost.push_back(0);
      Astart.push_back(Aindex.size());
      integerColumn.push_back(integerCol);
      std::string name(&line[field_2_start],
                       &line[field_2_start] + field_2_width);
      name = trim(name);
      col_names.push_back(name);
    }
    if (data[2] == objName)  // Cost
      colCost.back() = data[0];
    else if (data[0] != 0) {
      HighsInt iRow = rowIndex[data[2]] - 1;
      if (iRow >= 0) {
        if (rowType[iRow] != 'N' || keep_n_rows != kKeepNRowsDeleteEntries) {
          Aindex.push_back(iRow);
          Avalue.push_back(data[0]);
        }
      } else {
        // Spurious row name
        std::string name;
        if (!save_flag1) {
          std::string field_3(&line[field_3_start],
                              &line[field_3_start] + field_3_width);
          name = field_3;
        } else {
          std::string field_5(&line[field_5_start],
                              &line[field_5_start] + field_5_width);
          name = field_5;
        }
        num_alien_entries++;
        highsLogDev(
            log_options, HighsLogType::kInfo,
            "COLUMNS section contains row %-8s not in ROWS    section, line: "
            "%s\n",
            name.c_str(), line);
      }
    }
    save_flag1 = flag[1];
  }
  Astart.push_back(Aindex.size());

  if (num_alien_entries) {
    warning_issued = true;
    highsLogUser(log_options, HighsLogType::kWarning,
                 "COLUMNS section entries contain %8" HIGHSINT_FORMAT
                 " with row not in ROWS  "
                 "  section: ignored\n",
                 num_alien_entries);
  }
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read COLUMNS OK\n");

  // Load RHS
  num_alien_entries = 0;
  vector<double> RHS(numRow, 0);
  save_flag1 = 0;
  while (load_mpsLine(file, integerCol, lmax, line.data(), flag.data(),
                      data.data())) {
    if (data[2] != objName) {
      HighsInt iRow = rowIndex[data[2]] - 1;
      if (iRow >= 0) {
        RHS[iRow] = data[0];
      } else {
        // Spurious row name
        std::string name;
        if (!save_flag1) {
          std::string field_3(&line[field_3_start],
                              &line[field_3_start] + field_3_width);
          name = field_3;
        } else {
          std::string field_5(&line[field_5_start],
                              &line[field_5_start] + field_5_width);
          name = field_5;
        }
        num_alien_entries++;
        highsLogUser(log_options, HighsLogType::kInfo,
                     "RHS     section contains row %-8s not in ROWS    "
                     "section, line: %s\n",
                     name.c_str(), line);
      }
    } else {
      // Treat negation of a RHS entry for the N row as an objective
      // offset. Not all MPS readers do this, so give different
      // reported objective values for problems (eg e226)
      highsLogDev(
          log_options, HighsLogType::kInfo,
          "Using RHS value of %g for N-row in MPS file as negated objective "
          "offset\n",
          data[0]);
      objOffset = -data[0];  // Objective offset
    }
    save_flag1 = flag[1];
  }
  if (num_alien_entries) {
    warning_issued = true;
    highsLogUser(log_options, HighsLogType::kWarning,
                 "RHS     section entries contain %8" HIGHSINT_FORMAT
                 " with row not in ROWS  "
                 "  section: ignored\n",
                 num_alien_entries);
  }
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read RHS     OK\n");

  // Load RANGES
  num_alien_entries = 0;
  rowLower.resize(numRow);
  rowUpper.resize(numRow);
  if (flag[0] == 'R') {
    save_flag1 = 0;
    while (load_mpsLine(file, integerCol, lmax, line.data(), flag.data(),
                        data.data())) {
      HighsInt iRow = rowIndex[data[2]] - 1;
      if (iRow >= 0) {
        if (rowType[iRow] == 'L' || (rowType[iRow] == 'E' && data[0] < 0)) {
          rowLower[iRow] = RHS[iRow] - fabs(data[0]);
          rowUpper[iRow] = RHS[iRow];
        } else {
          rowUpper[iRow] = RHS[iRow] + fabs(data[0]);
          rowLower[iRow] = RHS[iRow];
        }
        rowType[iRow] = 'X';
      } else {
        // Spurious row name
        std::string name;
        if (!save_flag1) {
          std::string field_3(&line[field_3_start],
                              &line[field_3_start] + field_3_width);
          name = field_3;
        } else {
          std::string field_5(&line[field_5_start],
                              &line[field_5_start] + field_5_width);
          name = field_5;
        }
        num_alien_entries++;
        highsLogDev(
            log_options, HighsLogType::kInfo,
            "RANGES  section contains row %-8s not in ROWS    section, line: "
            "%s\n",
            name.c_str(), line);
      }
      save_flag1 = flag[1];
    }
  }

  // Setup bounds for row without 'RANGE'
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    switch (rowType[iRow]) {
      case 'L':
        rowLower[iRow] = -kHighsInf;
        rowUpper[iRow] = RHS[iRow];
        break;
      case 'G':
        rowLower[iRow] = RHS[iRow];
        rowUpper[iRow] = +kHighsInf;
        break;
      case 'E':
        rowLower[iRow] = RHS[iRow];
        rowUpper[iRow] = RHS[iRow];
        break;
      case 'N':
        rowLower[iRow] = -kHighsInf;
        rowUpper[iRow] = +kHighsInf;
        break;
      case 'X':
        break;
    }
  }
  if (num_alien_entries) {
    warning_issued = true;
    highsLogUser(log_options, HighsLogType::kWarning,
                 "RANGES  section entries contain %8" HIGHSINT_FORMAT
                 " with row not in ROWS  "
                 "  section: ignored\n",
                 num_alien_entries);
  }
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read RANGES  OK\n");

  // Load BOUNDS
  num_alien_entries = 0;
  colLower.assign(numCol, 0);
  colUpper.assign(numCol, kHighsInf);
  if (flag[0] == 'B') {
    while (load_mpsLine(file, integerCol, lmax, line.data(), flag.data(),
                        data.data())) {
      // Find the column index associated with the name "data[2]". If
      // the name is in colIndex then the value stored is the true
      // column index plus one. Otherwise 0 will be returned.
      HighsInt iCol = colIndex[data[2]] - 1;
      if (iCol >= 0) {
        switch (flag[0]) {
          case 'O': /*LO*/
            colLower[iCol] = data[0];
            break;
          case 'I': /*MI*/
            colLower[iCol] = -kHighsInf;
            break;
          case 'L': /*PL*/
            colUpper[iCol] = kHighsInf;
            break;
          case 'X': /*FX*/
            colLower[iCol] = data[0];
            colUpper[iCol] = data[0];
            break;
          case 'R': /*FR*/
            colLower[iCol] = -kHighsInf;
            colUpper[iCol] = kHighsInf;
            break;
          case 'P': /*UP*/
            colUpper[iCol] = data[0];
            if (colLower[iCol] == 0 && data[0] < 0) colLower[iCol] = -kHighsInf;
            break;
        }
      } else {
        std::string name(&line[field_3_start],
                         &line[field_3_start] + field_3_width);
        num_alien_entries++;
        highsLogDev(
            log_options, HighsLogType::kInfo,
            "BOUNDS  section contains col %-8s not in COLUMNS section, line: "
            "%s\n",
            name.c_str(), line);
      }
    }
  }
  // Load Hessian
  if (flag[0] == 'Q') {
    warning_issued = true;
    highsLogUser(log_options, HighsLogType::kWarning,
                 "Quadratic section: under development. Assumes "
                 "QUADOBJ/QSECTION section\n");
    Qdim = numCol;
    HighsInt hessian_nz = 0;
    HighsInt previous_col = -1;
    bool has_diagonal = false;
    Qstart.clear();
    while (load_mpsLine(file, integerCol, lmax, line.data(), flag.data(),
                        data.data())) {
      HighsInt iCol0 = colIndex[data[1]] - 1;
      std::string name0 = "";
      if (iCol0 >= 0) name0 = col_names[iCol0];
      HighsInt iCol1 = colIndex[data[2]] - 1;
      std::string name1 = "";
      if (iCol1 >= 0) name1 = col_names[iCol1];
      double value = data[0];
      if (iCol0 != previous_col) {
        // Can't handle columns out of order
        if (iCol0 < previous_col) {
          return FilereaderRetcode::kParserError;
        } else if (iCol0 > previous_col) {
          // A previous column has entries, but not on the diagonal
          if (previous_col >= 0 && !has_diagonal)
            return FilereaderRetcode::kParserError;
          // New column, so set up starts of any intermediate (empty) columns
          for (HighsInt iCol = previous_col + 1; iCol < iCol0; iCol++)
            Qstart.push_back(hessian_nz);
          Qstart.push_back(hessian_nz);
          previous_col = iCol0;
          has_diagonal = false;
        }
      }
      // Assumes QUADOBJ, so iCol1 cannot be less than iCol0
      if (iCol1 < iCol0) return FilereaderRetcode::kParserError;
      if (iCol1 == iCol0) has_diagonal = true;
      if (value) {
        Qindex.push_back(iCol1);
        Qvalue.push_back(value);
        hessian_nz++;
      }
    }
    // Hessian entries have been read, so set up starts of any remaining (empty)
    // columns
    for (HighsInt iCol = previous_col + 1; iCol < numCol; iCol++)
      Qstart.push_back(hessian_nz);
    Qstart.push_back(hessian_nz);
    assert((HighsInt)Qstart.size() == Qdim + 1);
    assert((HighsInt)Qindex.size() == hessian_nz);
    assert((HighsInt)Qvalue.size() == hessian_nz);
  }
  // Determine the number of integer variables and set bounds of [0,1]
  // for integer variables without bounds
  HighsInt num_int = 0;
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    if (integerColumn[iCol] == HighsVarType::kInteger) {
      num_int++;
      if (colUpper[iCol] >= kHighsInf) colUpper[iCol] = 1;
    }
  }
  if (num_alien_entries) {
    warning_issued = true;
    highsLogUser(log_options, HighsLogType::kWarning,
                 "BOUNDS  section entries contain %8" HIGHSINT_FORMAT
                 " with col not in "
                 "COLUMNS section: ignored\n",
                 num_alien_entries);
  }
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read BOUNDS  OK\n");
  highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read ENDATA  OK\n");
  highsLogDev(log_options, HighsLogType::kInfo,
              "readMPS: Model has %" HIGHSINT_FORMAT
              " rows and %" HIGHSINT_FORMAT " columns with %" HIGHSINT_FORMAT
              " integer\n",
              numRow, numCol, num_int);
  // Load ENDATA and close file
  file.close();
  // If there are no integer variables then clear the integrality vector
  if (!num_int) integerColumn.clear();
  return FilereaderRetcode::kOk;
}

bool load_mpsLine(std::istream& file, HighsVarType& integerVar, HighsInt lmax,
                  char* line, char* flag, double* data) {
  HighsInt F1 = 1, F2 = 4, F3 = 14, F4 = 24, F5 = 39, F6 = 49;

  // check the buffer
  if (flag[1]) {
    flag[1] = 0;
    memcpy(&data[2], &line[F5], 8);
    data[0] = atof(&line[F6]);
    return true;
  }

  // try to read some to the line
  for (;;) {
    // Line input
    *line = '\0';
    file.get(line, lmax);
    if (*line == '\0' && file.eof())  // nothing read and EOF
      return false;

    // Line trim   -- to delete tailing white spaces
    HighsInt lcnt = strlen(line) - 1;
    // if file.get() did not stop because it reached the lmax-1 limit,
    // then because it reached a newline char (or eof); lets consume this
    // newline (or do nothing if eof)
    if (lcnt + 1 < lmax - 1) file.get();
    while (isspace(line[lcnt]) && lcnt >= 0) lcnt--;
    if (lcnt <= 0 || line[0] == '*') continue;

    // Line expand -- to get data easier
    lcnt++;
    while (lcnt < F4) line[lcnt++] = ' ';  // For row and bound row name
    if (lcnt == F4) line[lcnt++] = '0';    // For bound value
    line[lcnt] = '\0';

    // Done with section symbol
    if (line[0] != ' ') {
      flag[0] = line[0];
      return false;
    }

    if (line[F3] == '\'') {
      if (line[F3 + 1] == 'M' && line[F3 + 2] == 'A' && line[F3 + 3] == 'R' &&
          line[F3 + 4] == 'K' && line[F3 + 5] == 'E' && line[F3 + 6] == 'R') {
        HighsInt cnter = line[F3 + 8];
        while (line[cnter] != '\'') ++cnter;
        if (line[cnter + 1] == 'I' && line[cnter + 2] == 'N' &&
            line[cnter + 3] == 'T' && line[cnter + 4] == 'O' &&
            line[cnter + 5] == 'R' && line[cnter + 6] == 'G')
          integerVar = HighsVarType::kInteger;
        else if (line[cnter + 1] == 'I' && line[cnter + 2] == 'N' &&
                 line[cnter + 3] == 'T' && line[cnter + 4] == 'E' &&
                 line[cnter + 5] == 'N' && line[cnter + 6] == 'D')
          integerVar = HighsVarType::kContinuous;
        continue;
      }
    }

    // Read major symbol & name
    flag[0] = line[F1 + 1] == ' ' ? line[F1] : line[F1 + 1];
    memcpy(&data[1], &line[F2], 8);
    // Read 1st minor name & value to output
    memcpy(&data[2], &line[F3], 8);
    data[0] = atof(&line[F4]);

    // Keep 2nd minor name & value for future
    if (lcnt > F5) flag[1] = 1;
    break;
  }

  return true;
}

HighsStatus writeModelAsMps(const HighsOptions& options,
                            const std::string& filename,
                            const HighsModel& model, const bool free_format) {
  bool warning_found = false;
  const HighsLp& lp = model.lp_;

  const bool ok_names = lp.okNames();
  assert(ok_names);
  if (!ok_names) return HighsStatus::kError;

  const HighsHessian& hessian = model.hessian_;

  HighsInt max_name_length = maxNameLength(lp);
  bool use_free_format = free_format;
  if (!free_format) {
    if (max_name_length > 8) {
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "Maximum name length is %" HIGHSINT_FORMAT
                   " so using free format rather "
                   "than fixed format\n",
                   max_name_length);
      use_free_format = true;
      warning_found = true;
    }
  }
  // Set a local objective name, creating one if necessary
  const std::string local_objective_name =
      findModelObjectiveName(&lp, &hessian);
  // If there is Hessian data to write out, writeMps assumes that hessian is
  // triangular
  if (hessian.dim_) assert(hessian.format_ == HessianFormat::kTriangular);

  HighsStatus write_status = writeMps(
      options.log_options, filename, lp.model_name_, lp.num_row_, lp.num_col_,
      hessian.dim_, lp.sense_, lp.offset_, lp.col_cost_, lp.col_lower_,
      lp.col_upper_, lp.row_lower_, lp.row_upper_, lp.a_matrix_.start_,
      lp.a_matrix_.index_, lp.a_matrix_.value_, hessian.start_, hessian.index_,
      hessian.value_, lp.integrality_, local_objective_name, lp.col_names_,
      lp.row_names_, use_free_format);
  if (write_status == HighsStatus::kOk && warning_found)
    return HighsStatus::kWarning;
  return write_status;
}

HighsStatus writeMps(
    const HighsLogOptions& log_options, const std::string& filename,
    const std::string& model_name, const HighsInt& num_row,
    const HighsInt& num_col, const HighsInt& q_dim, const ObjSense& sense,
    const double& offset, const vector<double>& col_cost,
    const vector<double>& col_lower, const vector<double>& col_upper,
    const vector<double>& row_lower, const vector<double>& row_upper,
    const vector<HighsInt>& a_start, const vector<HighsInt>& a_index,
    const vector<double>& a_value, const vector<HighsInt>& q_start,
    const vector<HighsInt>& q_index, const vector<double>& q_value,
    const vector<HighsVarType>& integrality, const std::string& objective_name,
    const vector<std::string>& col_names, const vector<std::string>& row_names,
    const bool use_free_format) {
  HighsInt num_no_cost_zero_columns = 0;
  HighsInt num_no_cost_zero_columns_in_bounds_section = 0;
  highsLogDev(log_options, HighsLogType::kInfo,
              "writeMPS: Trying to open file %s\n", filename.c_str());
  FILE* file = fopen(filename.c_str(), "w");
  if (file == 0) {
    highsLogUser(log_options, HighsLogType::kError, "Cannot open file %s\n",
                 filename.c_str());
    return HighsStatus::kError;
  }
  highsLogDev(log_options, HighsLogType::kInfo, "writeMPS: Opened file  OK\n");
  // Check that the names are no longer than 8 characters for fixed format write
  HighsInt max_name_length =
      std::max(maxNameLength(col_names), maxNameLength(row_names));
  if (!use_free_format && max_name_length > 8) {
    highsLogUser(
        log_options, HighsLogType::kError,
        "Cannot write fixed MPS with names of length (up to) %" HIGHSINT_FORMAT
        "\n",
        max_name_length);
    fclose(file);
    return HighsStatus::kError;
  }
  assert(objective_name != "");
  vector<HighsInt> r_ty;
  vector<double> rhs, ranges;
  bool have_rhs = false;
  bool have_ranges = false;
  bool have_bounds = false;
  bool have_int = false;
  r_ty.resize(num_row);
  rhs.assign(num_row, 0);
  ranges.assign(num_row, 0);
  for (HighsInt r_n = 0; r_n < num_row; r_n++) {
    if (row_lower[r_n] == row_upper[r_n]) {
      // Equality constraint - Type E - range = 0
      r_ty[r_n] = MPS_ROW_TY_E;
      rhs[r_n] = row_lower[r_n];
    } else if (!highs_isInfinity(row_upper[r_n])) {
      // Upper bounded constraint - Type L
      r_ty[r_n] = MPS_ROW_TY_L;
      rhs[r_n] = row_upper[r_n];
      if (!highs_isInfinity(-row_lower[r_n])) {
        // Boxed constraint - range = u-l
        ranges[r_n] = row_upper[r_n] - row_lower[r_n];
      }
    } else if (!highs_isInfinity(-row_lower[r_n])) {
      // Lower bounded constraint - Type G
      r_ty[r_n] = MPS_ROW_TY_G;
      rhs[r_n] = row_lower[r_n];
    } else {
      // Free constraint - Type N
      r_ty[r_n] = MPS_ROW_TY_N;
      rhs[r_n] = 0;
    }
  }

  for (HighsInt r_n = 0; r_n < num_row; r_n++) {
    if (rhs[r_n]) {
      have_rhs = true;
      break;
    }
  }
  // Check whether there is an objective offset - which will be defines as a RHS
  // on the cost row
  if (offset) have_rhs = true;
  for (HighsInt r_n = 0; r_n < num_row; r_n++) {
    if (ranges[r_n]) {
      have_ranges = true;
      break;
    }
  }
  have_int = false;
  if (integrality.size()) {
    for (HighsInt c_n = 0; c_n < num_col; c_n++) {
      if (integrality[c_n] == HighsVarType::kInteger ||
          integrality[c_n] == HighsVarType::kSemiContinuous ||
          integrality[c_n] == HighsVarType::kSemiInteger) {
        have_int = true;
        break;
      }
    }
  }
  for (HighsInt c_n = 0; c_n < num_col; c_n++) {
    if (col_lower[c_n]) {
      have_bounds = true;
      break;
    }
    bool discrete = false;
    if (have_int)
      discrete = integrality[c_n] == HighsVarType::kInteger ||
                 integrality[c_n] == HighsVarType::kSemiContinuous ||
                 integrality[c_n] == HighsVarType::kSemiInteger;
    if (!highs_isInfinity(col_upper[c_n]) || discrete) {
      // If the upper bound is finite, or the variable is integer or a
      // semi-variable then there is a BOUNDS section. Integer
      // variables with infinite upper bound are indicated as LI. All
      // semi-variables appear in the BOUNDS section.
      have_bounds = true;
      break;
    }
  }
  highsLogDev(log_options, HighsLogType::kInfo,
              "Model: RHS =     %s\n       RANGES =  %s\n       BOUNDS =  %s\n",
              highsBoolToString(have_rhs).c_str(),
              highsBoolToString(have_ranges).c_str(),
              highsBoolToString(have_bounds).c_str());

  // Field:    1           2          3         4         5         6
  // Columns:  2-3        5-12      15-22     25-36     40-47     50-61 Indexed
  // from 1 Columns:  1-2        4-11      14-21     24-35     39-46     49-60
  // Indexed from 0
  //           1         2         3         4         5         6
  // 0123456789012345678901234567890123456789012345678901234567890
  // x11x22222222xx33333333xx444444444444xxx55555555xx666666666666
  // ROWS
  //  N  ENDCAP
  // COLUMNS
  //     CFOOD01   BAGR01          .00756   BFTT01         .150768
  // RHS
  //     RHSIDE    HCAP01            -20.   CBCAP01            -8.
  // RANGES
  //     RANGE1    VILLKOR2            7.   VILLKOR3            7.
  // BOUNDS
  //  LO BOUND     CFOOD01           850.
  //
  // NB d6cube has (eg)
  // COLUMNS
  //        1      1                   1.   4                  -1.
  //        1      5                  -1.   1151                1.
  // Indexed from 0
  //           1         2         3         4         5         6
  // 0123456789012345678901234567890123456789012345678901234567890
  // x11x22222222xx33333333xx444444444444xxx55555555xx666666666666
  //
  // In fixed format the first column name is "      1 ", and its first entry is
  // in row "1       ".
  //
  // The free format reader thought that it had a name of "1      1" containing
  // spaces.

  fprintf(file, "NAME        %s\n", model_name.c_str());
  const bool use_objsense = true;
  const HighsInt use_sense = use_objsense ? 1 : (HighsInt)sense;
  if (sense == ObjSense::kMaximize) fprintf(file, "OBJSENSE\n  MAX\n");
  fprintf(file, "ROWS\n");
  fprintf(file, " N  %-8s\n", objective_name.c_str());
  for (HighsInt r_n = 0; r_n < num_row; r_n++) {
    if (r_ty[r_n] == MPS_ROW_TY_E) {
      fprintf(file, " E  %-8s\n", row_names[r_n].c_str());
    } else if (r_ty[r_n] == MPS_ROW_TY_G) {
      fprintf(file, " G  %-8s\n", row_names[r_n].c_str());
    } else if (r_ty[r_n] == MPS_ROW_TY_L) {
      fprintf(file, " L  %-8s\n", row_names[r_n].c_str());
    } else {
      fprintf(file, " N  %-8s\n", row_names[r_n].c_str());
    }
  }
  bool integerFg = false;
  HighsInt nIntegerMk = 0;
  fprintf(file, "COLUMNS\n");

  auto noCostZeroColumn = [&](const HighsInt c_n) {
    if (col_cost[c_n]) return false;
    if (a_start[c_n] < a_start[c_n + 1]) return false;
    if (q_dim > 0 && q_value[q_start[c_n]]) return false;
    return true;
  };

  for (HighsInt c_n = 0; c_n < num_col; c_n++) {
    if (noCostZeroColumn(c_n)) {
      // Possibly skip this column as it's zero and has no cost
      num_no_cost_zero_columns++;
      // Give the column a presence by writing out a zero cost
      fprintf(file, "    %-8s  %-8s  %.17g\n", col_names[c_n].c_str(),
              objective_name.c_str(), 0.0);
      continue;
    }
    if (have_int) {
      if (integrality[c_n] == HighsVarType::kInteger && !integerFg) {
        // Start an integer section
        fprintf(file,
                "    MARK%04" HIGHSINT_FORMAT
                "  'MARKER'                 'INTORG'\n",
                nIntegerMk++);
        integerFg = true;
      } else if (integrality[c_n] != HighsVarType::kInteger && integerFg) {
        // End an integer section
        fprintf(file,
                "    MARK%04" HIGHSINT_FORMAT
                "  'MARKER'                 'INTEND'\n",
                nIntegerMk++);
        integerFg = false;
      }
    }
    if (col_cost[c_n] != 0) {
      double v = use_sense * col_cost[c_n];
      fprintf(file, "    %-8s  %-8s  %.17g\n", col_names[c_n].c_str(),
              objective_name.c_str(), v);
    }
    for (HighsInt el_n = a_start[c_n]; el_n < a_start[c_n + 1]; el_n++) {
      double v = a_value[el_n];
      HighsInt r_n = a_index[el_n];
      fprintf(file, "    %-8s  %-8s  %.17g\n", col_names[c_n].c_str(),
              row_names[r_n].c_str(), v);
    }
  }
  // End any integer section
  if (integerFg)
    fprintf(file,
            "    MARK%04" HIGHSINT_FORMAT
            "  'MARKER'                 'INTEND'\n",
            nIntegerMk++);
  have_rhs = true;
  if (have_rhs) {
    fprintf(file, "RHS\n");
    if (offset) {
      // Handle the objective offset as a RHS entry for the cost row
      double v = -use_sense * offset;
      fprintf(file, "    RHS_V     %-8s  %.17g\n", objective_name.c_str(), v);
    }
    for (HighsInt r_n = 0; r_n < num_row; r_n++) {
      double v = rhs[r_n];
      if (v) {
        fprintf(file, "    RHS_V     %-8s  %.17g\n", row_names[r_n].c_str(), v);
      }
    }
  }
  if (have_ranges) {
    fprintf(file, "RANGES\n");
    for (HighsInt r_n = 0; r_n < num_row; r_n++) {
      double v = ranges[r_n];
      if (v) {
        fprintf(file, "    RANGE     %-8s  %.17g\n", row_names[r_n].c_str(), v);
      }
    }
  }
  if (have_bounds) {
    fprintf(file, "BOUNDS\n");
    for (HighsInt c_n = 0; c_n < num_col; c_n++) {
      double lb = col_lower[c_n];
      double ub = col_upper[c_n];
      bool discrete = false;
      bool semi = false;
      if (have_int) {
        semi = integrality[c_n] == HighsVarType::kSemiContinuous ||
               integrality[c_n] == HighsVarType::kSemiInteger;
        discrete = integrality[c_n] == HighsVarType::kInteger || semi;
      }
      if (noCostZeroColumn(c_n)) {
        // Possibly skip this column if it's zero and has no cost
        if (!highs_isInfinity(ub) || lb) {
          // Column would have a bound to report
          num_no_cost_zero_columns_in_bounds_section++;
        }
      }
      if (lb == ub && !semi) {
        // Equal lower and upper bounds: Fixed. Not for semi-variables,
        // since FX would lose the zero value, and the variable type
        // is defined by SC/SI
        fprintf(file, " FX BOUND     %-8s  %.17g\n", col_names[c_n].c_str(),
                lb);
      } else if (highs_isInfinity(-lb) && highs_isInfinity(ub) && !semi) {
        // Infinite lower and upper bounds: Free
        fprintf(file, " FR BOUND     %-8s\n", col_names[c_n].c_str());
      } else {
        if (discrete) {
          // Warn if writing non-integer bounds for integer or semi-integer
          // variables
          if (integrality[c_n] == HighsVarType::kInteger ||
              integrality[c_n] == HighsVarType::kSemiInteger) {
            if (lb > -kHighsInf) {
              HighsInt i_lb = static_cast<HighsInt>(lb);
              double dl = lb - i_lb;
              if (dl)
                highsLogUser(log_options, HighsLogType::kWarning,
                             "Lower bound for integer or semi-integer column "
                             "\"%s\" is %g: not integer\n",
                             col_names[c_n].c_str(), lb);
            }
            if (ub < kHighsInf) {
              HighsInt i_ub = static_cast<HighsInt>(ub);
              double dl = ub - i_ub;
              if (dl)
                highsLogUser(log_options, HighsLogType::kWarning,
                             "Upper bound for integer or semi-integer column "
                             "\"%s\" is %g: not integer\n",
                             col_names[c_n].c_str(), ub);
            }
          }
          if (integrality[c_n] == HighsVarType::kInteger) {
            if (lb == 0 && ub == 1) {
              // Binary
              fprintf(file, " BV BOUND     %-8s\n", col_names[c_n].c_str());
            } else {
              if (!highs_isInfinity(-lb)) {
                // Finite lower bound. No need to state this if LB is
                // zero unless UB is infinite
                if (lb || highs_isInfinity(ub))
                  fprintf(file, " LI BOUND     %-8s  %.17g\n",
                          col_names[c_n].c_str(), lb);
              } else {
                // Infinite lower bound
                fprintf(file, " MI BOUND     %-8s\n", col_names[c_n].c_str());
              }
              if (!highs_isInfinity(ub)) {
                // Finite upper bound
                fprintf(file, " UI BOUND     %-8s  %.17g\n",
                        col_names[c_n].c_str(), ub);
              }
            }
          } else if (integrality[c_n] == HighsVarType::kSemiInteger ||
                     integrality[c_n] == HighsVarType::kSemiContinuous) {
            // Have to use lb and ub to define semi-variables: lb is
            // naturally finite, but what if ub is infinite?
            //
            // writing "inf" is fine, except for with mingw, it
            // appears, so use something bigger than the infinite
            // bounds value - since options.infinite_bound is
            // unavailable
            const double infinite_bound = 1e30;
            const bool inf_lb = lb <= -kHighsInf;
            const double use_lb = inf_lb ? -infinite_bound : lb;
            const bool inf_ub = ub >= kHighsInf;
            const double use_ub = inf_ub ? infinite_bound : ub;
            if (inf_lb)
              highsLogUser(
                  log_options, HighsLogType::kWarning,
                  "Lower bound for semi-variable \"%s\" is %g but writing %g\n",
                  col_names[c_n].c_str(), lb, use_lb);
            if (inf_ub)
              highsLogUser(
                  log_options, HighsLogType::kWarning,
                  "Upper bound for semi-variable \"%s\" is %g but writing %g\n",
                  col_names[c_n].c_str(), ub, use_ub);
            fprintf(file, " LO BOUND     %-8s  %.17g\n", col_names[c_n].c_str(),
                    use_lb);
            if (integrality[c_n] == HighsVarType::kSemiInteger) {
              fprintf(file, " SI BOUND     %-8s  %.17g\n",
                      col_names[c_n].c_str(), use_ub);
            } else {
              // Semi-continuous
              fprintf(file, " SC BOUND     %-8s  %.17g\n",
                      col_names[c_n].c_str(), use_ub);
            }
          }
        } else {
          if (!highs_isInfinity(-lb)) {
            // Lower bounded variable - default is 0
            if (lb) {
              fprintf(file, " LO BOUND     %-8s  %.17g\n",
                      col_names[c_n].c_str(), lb);
            }
          } else {
            // Infinite lower bound
            fprintf(file, " MI BOUND     %-8s\n", col_names[c_n].c_str());
          }
          if (!highs_isInfinity(ub)) {
            // Upper bounded variable
            fprintf(file, " UP BOUND     %-8s  %.17g\n", col_names[c_n].c_str(),
                    ub);
          }
        }
      }
    }
  }
  if (q_dim) {
    // Write out Hessian info
    assert((HighsInt)q_start.size() >= q_dim + 1);
    assert((HighsInt)q_index.size() >= q_start[q_dim]);
    assert((HighsInt)q_value.size() >= q_start[q_dim]);

    // Assumes that Hessian entries are the lower triangle column-wise
    fprintf(file, "QUADOBJ\n");
    for (HighsInt col = 0; col < q_dim; col++) {
      for (HighsInt el = q_start[col]; el < q_start[col + 1]; el++) {
        HighsInt row = q_index[el];
        assert(row >= col);
        // May have explicit zeroes on the diagonal
        if (q_value[el])
          fprintf(file, "    %-8s  %-8s  %.17g\n", col_names[col].c_str(),
                  col_names[row].c_str(), use_sense * q_value[el]);
      }
    }
  }
  fprintf(file, "ENDATA\n");
  if (num_no_cost_zero_columns)
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Model has %" HIGHSINT_FORMAT
                 " zero columns with no costs: %" HIGHSINT_FORMAT
                 " have finite upper bounds "
                 "or nonzero lower bounds and are written in MPS file\n",
                 num_no_cost_zero_columns,
                 num_no_cost_zero_columns_in_bounds_section);
  fclose(file);
  return HighsStatus::kOk;
}
