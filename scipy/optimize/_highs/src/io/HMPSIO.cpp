/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HMPSIO.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "io/HMPSIO.h"

#include <algorithm>

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsOptions.h"
#include "util/HighsUtils.h"
#include "util/stringutil.h"

using std::map;

//
// Read file called filename. Returns 0 if OK and 1 if file can't be opened
//
FilereaderRetcode readMPS(FILE* logfile, const std::string filename,
                          int mxNumRow, int mxNumCol, int& numRow, int& numCol,
                          ObjSense& objSense, double& objOffset,
                          vector<int>& Astart, vector<int>& Aindex,
                          vector<double>& Avalue, vector<double>& colCost,
                          vector<double>& colLower, vector<double>& colUpper,
                          vector<double>& rowLower, vector<double>& rowUpper,
                          vector<int>& integerColumn, vector<string>& col_names,
                          vector<string>& row_names, const int keep_n_rows) {
  // MPS file buffer
  numRow = 0;
  numCol = 0;
  objOffset = 0;
  objSense = ObjSense::MINIMIZE;

  // Astart.clear() added since setting Astart.push_back(0) in
  // setup_clearModel() messes up the MPS read
  Astart.clear();
#ifdef HiGHSDEV
  printf("readMPS: Trying to open file %s\n", filename.c_str());
#endif
  FILE* file = fopen(filename.c_str(), "r");
  if (file == 0) {
#ifdef HiGHSDEV
    printf("readMPS: Not opened file OK\n");
#endif
    return FilereaderRetcode::FILENOTFOUND;
  }
#ifdef HiGHSDEV
  printf("readMPS: Opened file  OK\n");
#endif
  // Input buffer
  const int lmax = 128;
  char line[lmax];
  char flag[2] = {0, 0};
  double data[3];

  int num_alien_entries = 0;
  int integerCol = 0;

  // Load NAME
  load_mpsLine(file, integerCol, lmax, line, flag, data);
#ifdef HiGHSDEV
  printf("readMPS: Read NAME    OK\n");
#endif
  // Load OBJSENSE or ROWS
  load_mpsLine(file, integerCol, lmax, line, flag, data);
  if (flag[0] == 'O') {
    // Found OBJSENSE
    load_mpsLine(file, integerCol, lmax, line, flag, data);
    std::string sense(&line[2], &line[2] + 3);
    // the sense must be "MAX" or "MIN"
    if (sense.compare("MAX") == 0) {
      objSense = ObjSense::MAXIMIZE;
    } else if (sense.compare("MIN") == 0) {
      objSense = ObjSense::MINIMIZE;
    } else {
      return FilereaderRetcode::PARSERERROR;
    }
#ifdef HiGHSDEV
    printf("readMPS: Read OBJSENSE OK\n");
#endif
    // Load ROWS
    load_mpsLine(file, integerCol, lmax, line, flag, data);
  }

  row_names.clear();
  col_names.clear();
  vector<char> rowType;
  map<double, int> rowIndex;
  double objName = 0;
  while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
    if (flag[0] == 'N' &&
        (objName == 0 || keep_n_rows == KEEP_N_ROWS_DELETE_ROWS)) {
      // N-row: take the first as the objective and possibly ignore any others
      if (objName == 0) objName = data[1];
    } else {
      if (mxNumRow > 0 && numRow >= mxNumRow)
        return FilereaderRetcode::PARSERERROR;
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
#ifdef HiGHSDEV
  printf("readMPS: Read ROWS    OK\n");
#endif

  // Load COLUMNS
  map<double, int> colIndex;
  double lastName = 0;
  // flag[1] is used to indicate whether there is more to read on the
  // line - field 5 non-empty. save_flag1 is used to deduce whether
  // the row name and value are from fields 5 and 6, or 3 and 4
  int save_flag1 = 0;
  while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
    int iRow = rowIndex[data[2]] - 1;
    std::string name = "";
    if (iRow >= 0) name = row_names[iRow];
    if (lastName != data[1]) {  // New column
      if (mxNumCol > 0 && numCol >= mxNumCol)
        return FilereaderRetcode::PARSERERROR;
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
      int iRow = rowIndex[data[2]] - 1;
      if (iRow >= 0) {
        if (rowType[iRow] != 'N' || keep_n_rows != KEEP_N_ROWS_DELETE_ENTRIES) {
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
#ifdef HiGHSDEV
        printf(
            "COLUMNS section contains row %-8s not in ROWS    section, line: "
            "%s\n",
            name.c_str(), line);
#endif
      }
    }
    save_flag1 = flag[1];
  }
  Astart.push_back(Aindex.size());

  if (num_alien_entries)
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "COLUMNS section entries contain %8d with row not in ROWS  "
                    "  section: ignored",
                    num_alien_entries);
#ifdef HiGHSDEV
  printf("readMPS: Read COLUMNS OK\n");
#endif

  // Load RHS
  num_alien_entries = 0;
  vector<double> RHS(numRow, 0);
  save_flag1 = 0;
  while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
    if (data[2] != objName) {
      int iRow = rowIndex[data[2]] - 1;
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
        HighsLogMessage(logfile, HighsMessageType::INFO,
                        "RHS     section contains row %-8s not in ROWS    "
                        "section, line: %s",
                        name.c_str(), line);
      }
    } else {
      // Treat negation of a RHS entry for the N row as an objective
      // offset. Not all MPS readers do this, so give different
      // reported objective values for problems (eg e226)
#ifdef HiGHSDEV
      printf(
          "Using RHS value of %g for N-row in MPS file as negated objective "
          "offset\n",
          data[0]);
#endif
      objOffset = -data[0];  // Objective offset
    }
    save_flag1 = flag[1];
  }
  if (num_alien_entries)
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "RHS     section entries contain %8d with row not in ROWS  "
                    "  section: ignored",
                    num_alien_entries);
#ifdef HiGHSDEV
  printf("readMPS: Read RHS     OK\n");
#endif

  // Load RANGES
  num_alien_entries = 0;
  rowLower.resize(numRow);
  rowUpper.resize(numRow);
  if (flag[0] == 'R') {
    save_flag1 = 0;
    while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
      int iRow = rowIndex[data[2]] - 1;
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
#ifdef HiGHSDEV
        printf(
            "RANGES  section contains row %-8s not in ROWS    section, line: "
            "%s\n",
            name.c_str(), line);
#endif
      }
      save_flag1 = flag[1];
    }
  }

  // Setup bounds for row without 'RANGE'
  for (int iRow = 0; iRow < numRow; iRow++) {
    switch (rowType[iRow]) {
      case 'L':
        rowLower[iRow] = -HIGHS_CONST_INF;
        rowUpper[iRow] = RHS[iRow];
        break;
      case 'G':
        rowLower[iRow] = RHS[iRow];
        rowUpper[iRow] = +HIGHS_CONST_INF;
        break;
      case 'E':
        rowLower[iRow] = RHS[iRow];
        rowUpper[iRow] = RHS[iRow];
        break;
      case 'N':
        rowLower[iRow] = -HIGHS_CONST_INF;
        rowUpper[iRow] = +HIGHS_CONST_INF;
        break;
      case 'X':
        break;
    }
  }
  if (num_alien_entries)
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "RANGES  section entries contain %8d with row not in ROWS  "
                    "  section: ignored",
                    num_alien_entries);
#ifdef HiGHSDEV
  printf("readMPS: Read RANGES  OK\n");
#endif

  // Load BOUNDS
  num_alien_entries = 0;
  colLower.assign(numCol, 0);
  colUpper.assign(numCol, HIGHS_CONST_INF);
  if (flag[0] == 'B') {
    while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
      // Find the column index associated woith the name "data[2]". If
      // the name is in colIndex then the value stored is the true
      // column index plus one. Otherwise 0 will be returned.
      int iCol = colIndex[data[2]] - 1;
      if (iCol >= 0) {
        switch (flag[0]) {
          case 'O': /*LO*/
            colLower[iCol] = data[0];
            break;
          case 'I': /*MI*/
            colLower[iCol] = -HIGHS_CONST_INF;
            break;
          case 'L': /*PL*/
            colUpper[iCol] = HIGHS_CONST_INF;
            break;
          case 'X': /*FX*/
            colLower[iCol] = data[0];
            colUpper[iCol] = data[0];
            break;
          case 'R': /*FR*/
            colLower[iCol] = -HIGHS_CONST_INF;
            colUpper[iCol] = HIGHS_CONST_INF;
            break;
          case 'P': /*UP*/
            colUpper[iCol] = data[0];
            if (colLower[iCol] == 0 && data[0] < 0)
              colLower[iCol] = -HIGHS_CONST_INF;
            break;
        }
      } else {
        std::string name(&line[field_3_start],
                         &line[field_3_start] + field_3_width);
        num_alien_entries++;
#ifdef HiGHSDEV
        printf(
            "BOUNDS  section contains col %-8s not in COLUMNS section, line: "
            "%s\n",
            name.c_str(), line);
#endif
      }
    }
  }
  // Determine the number of integer variables and set bounds of [0,1]
  // for integer variables without bounds
  int num_int = 0;
  for (int iCol = 0; iCol < numCol; iCol++) {
    if (integerColumn[iCol]) {
      num_int++;
      if (colUpper[iCol] >= HIGHS_CONST_INF) colUpper[iCol] = 1;
    }
  }
  if (num_alien_entries)
    HighsLogMessage(logfile, HighsMessageType::WARNING,
                    "BOUNDS  section entries contain %8d with col not in "
                    "COLUMNS section: ignored",
                    num_alien_entries);
#ifdef HiGHSDEV
  printf("readMPS: Read BOUNDS  OK\n");
  printf("readMPS: Read ENDATA  OK\n");
  printf("readMPS: Model has %d rows and %d columns with %d integer\n", numRow,
         numCol, num_int);
#endif
  // Load ENDATA and close file
  fclose(file);
  return FilereaderRetcode::OK;
}

bool load_mpsLine(FILE* file, int& integerVar, int lmax, char* line, char* flag,
                  double* data) {
  int F1 = 1, F2 = 4, F3 = 14, F4 = 24, F5 = 39, F6 = 49;
  char* fgets_rt;

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
    fgets_rt = fgets(line, lmax, file);
    if (fgets_rt == NULL) {
      return false;
    }
    // Line trim   -- to delete tailing white spaces
    int lcnt = strlen(line) - 1;
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
        int cnter = line[F3 + 8];
        while (line[cnter] != '\'') ++cnter;
        if (line[cnter + 1] == 'I' && line[cnter + 2] == 'N' &&
            line[cnter + 3] == 'T' && line[cnter + 4] == 'O' &&
            line[cnter + 5] == 'R' && line[cnter + 6] == 'G')
          integerVar = 1;
        else if (line[cnter + 1] == 'I' && line[cnter + 2] == 'N' &&
                 line[cnter + 3] == 'T' && line[cnter + 4] == 'E' &&
                 line[cnter + 5] == 'N' && line[cnter + 6] == 'D')
          integerVar = 0;
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

HighsStatus writeLpAsMPS(const HighsOptions& options,
                         const std::string filename, const HighsLp& lp,
                         const bool free_format) {
  bool warning_found = false;
  bool have_col_names = lp.col_names_.size();
  bool have_row_names = lp.row_names_.size();
  std::vector<std::string> local_col_names;
  std::vector<std::string> local_row_names;
  local_col_names.resize(lp.numCol_);
  local_row_names.resize(lp.numRow_);
  //
  // Initialise the local names to any existing names
  if (have_col_names) local_col_names = lp.col_names_;
  if (have_row_names) local_row_names = lp.row_names_;
  //
  // Normalise the column names
  int max_col_name_length = HIGHS_CONST_I_INF;
  if (!free_format) max_col_name_length = 8;
  HighsStatus col_name_status = normaliseNames(
      options, "Column", lp.numCol_, local_col_names, max_col_name_length);
  if (col_name_status == HighsStatus::Error) return col_name_status;
  warning_found = col_name_status == HighsStatus::Warning || warning_found;
  //
  // Normalise the row names
  int max_row_name_length = HIGHS_CONST_I_INF;
  if (!free_format) max_row_name_length = 8;
  HighsStatus row_name_status = normaliseNames(
      options, "Row", lp.numRow_, local_row_names, max_row_name_length);
  if (row_name_status == HighsStatus::Error) return col_name_status;
  warning_found = row_name_status == HighsStatus::Warning || warning_found;

  int max_name_length = std::max(max_col_name_length, max_row_name_length);
  bool use_free_format = free_format;
  if (!free_format) {
    if (max_name_length > 8) {
      HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                      "Maximum name length is %d so using free format rather "
                      "than fixed format",
                      max_name_length);
      use_free_format = true;
      warning_found = true;
    }
  }
  HighsStatus write_status = writeMPS(
      options.logfile, filename, lp.numRow_, lp.numCol_, lp.sense_, lp.offset_,
      lp.Astart_, lp.Aindex_, lp.Avalue_, lp.colCost_, lp.colLower_,
      lp.colUpper_, lp.rowLower_, lp.rowUpper_, lp.integrality_,
      local_col_names, local_row_names, use_free_format);
  if (write_status == HighsStatus::OK && warning_found)
    return HighsStatus::Warning;
  return write_status;
}

HighsStatus writeMPS(
    FILE* logfile, const std::string filename, const int& numRow,
    const int& numCol, const ObjSense& objSense, const double& objOffset,
    const vector<int>& Astart, const vector<int>& Aindex,
    const vector<double>& Avalue, const vector<double>& colCost,
    const vector<double>& colLower, const vector<double>& colUpper,
    const vector<double>& rowLower, const vector<double>& rowUpper,
    const vector<int>& integerColumn, const vector<std::string>& col_names,
    const vector<std::string>& row_names, const bool use_free_format) {
  const bool write_zero_no_cost_columns = true;
  int num_zero_no_cost_columns = 0;
  int num_zero_no_cost_columns_in_bounds_section = 0;
#ifdef HiGHSDEV
  printf("writeMPS: Trying to open file %s\n", filename.c_str());
#endif
  FILE* file = fopen(filename.c_str(), "w");
  if (file == 0) {
    HighsLogMessage(logfile, HighsMessageType::ERROR, "Cannot open file %s",
                    filename.c_str());
    return HighsStatus::Error;
  }
#ifdef HiGHSDEV
  printf("writeMPS: Opened file  OK\n");
#endif
  // Check that the names are no longer than 8 characters for fixed format write
  int max_col_name_length = maxNameLength(numCol, col_names);
  int max_row_name_length = maxNameLength(numRow, row_names);
  int max_name_length = std::max(max_col_name_length, max_row_name_length);
  if (!use_free_format && max_name_length > 8) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "Cannot write fixed MPS with names of length (up to) %d",
                    max_name_length);
    return HighsStatus::Error;
  }
  vector<int> r_ty;
  vector<double> rhs, ranges;
  bool have_rhs = false;
  bool have_ranges = false;
  bool have_bounds = false;
  bool have_int = false;
  r_ty.resize(numRow);
  rhs.assign(numRow, 0);
  ranges.assign(numRow, 0);
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (rowLower[r_n] == rowUpper[r_n]) {
      // Equality constraint - Type E - range = 0
      r_ty[r_n] = MPS_ROW_TY_E;
      rhs[r_n] = rowLower[r_n];
    } else if (!highs_isInfinity(rowUpper[r_n])) {
      // Upper bounded constraint - Type L
      r_ty[r_n] = MPS_ROW_TY_L;
      rhs[r_n] = rowUpper[r_n];
      if (!highs_isInfinity(-rowLower[r_n])) {
        // Boxed constraint - range = u-l
        ranges[r_n] = rowUpper[r_n] - rowLower[r_n];
      }
    } else if (!highs_isInfinity(-rowLower[r_n])) {
      // Lower bounded constraint - Type G
      r_ty[r_n] = MPS_ROW_TY_G;
      rhs[r_n] = rowLower[r_n];
    } else {
      // Free constraint - Type N
      r_ty[r_n] = MPS_ROW_TY_N;
      rhs[r_n] = 0;
    }
  }

  for (int r_n = 0; r_n < numRow; r_n++) {
    if (rhs[r_n]) {
      have_rhs = true;
      break;
    }
  }
  // Check whether there is an objective offset - which will be defines as a RHS
  // on the cost row
  if (objOffset) have_rhs = true;
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (ranges[r_n]) {
      have_ranges = true;
      break;
    }
  }
  have_int = false;
  if (integerColumn.size()) {
    for (int c_n = 0; c_n < numCol; c_n++) {
      if (integerColumn[c_n]) {
        have_int = true;
        break;
      }
    }
  }
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (colLower[c_n]) {
      have_bounds = true;
      break;
    }
    bool discrete = false;
    if (have_int) discrete = integerColumn[c_n];
    if (!highs_isInfinity(colUpper[c_n]) || discrete) {
      // If the upper bound is finite, or the variable is integer then there is
      // a BOUNDS section. Integer variables with infinite upper bound are
      // indicated as LI
      have_bounds = true;
      break;
    }
  }
#ifdef HiGHSDEV
  printf("Model: RHS =     %s\n       RANGES =  %s\n       BOUNDS =  %s\n",
         BoolToCharStar(have_rhs), BoolToCharStar(have_ranges),
         BoolToCharStar(have_bounds));
#endif

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
  fprintf(file, "NAME\n");
  fprintf(file, "ROWS\n");
  fprintf(file, " N  COST\n");
  for (int r_n = 0; r_n < numRow; r_n++) {
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
  int nIntegerMk = 0;
  fprintf(file, "COLUMNS\n");
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (Astart[c_n] == Astart[c_n + 1] && colCost[c_n] == 0) {
      // Possibly skip this column as it's zero and has no cost
      num_zero_no_cost_columns++;
      if (write_zero_no_cost_columns) {
        // Give the column a presence by writing out a zero cost
        double v = 0;
        fprintf(file, "    %-8s  COST      %.15g\n", col_names[c_n].c_str(), v);
      }
      continue;
    }
    if (have_int) {
      if (integerColumn[c_n] && !integerFg) {
        // Start an integer section
        fprintf(file, "    MARK%04d  'MARKER'                 'INTORG'\n",
                nIntegerMk);
        nIntegerMk++;
        integerFg = true;
      } else if (!integerColumn[c_n] && integerFg) {
        // End an integer section
        fprintf(file, "    MARK%04d  'MARKER'                 'INTEND'\n",
                nIntegerMk);
        nIntegerMk++;
        integerFg = false;
      }
    }
    if (colCost[c_n] != 0) {
      double v = (int)objSense * colCost[c_n];
      fprintf(file, "    %-8s  COST      %.15g\n", col_names[c_n].c_str(), v);
    }
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      double v = Avalue[el_n];
      int r_n = Aindex[el_n];
      fprintf(file, "    %-8s  %-8s  %.15g\n", col_names[c_n].c_str(),
              row_names[r_n].c_str(), v);
    }
  }
  have_rhs = true;
  if (have_rhs) {
    fprintf(file, "RHS\n");
    if (objOffset) {
      // Handle the objective offset as a RHS entry for the cost row
      double v = -(int)objSense * objOffset;
      fprintf(file, "    RHS_V     COST      %.15g\n", v);
    }
    for (int r_n = 0; r_n < numRow; r_n++) {
      double v = rhs[r_n];
      if (v) {
        fprintf(file, "    RHS_V     %-8s  %.15g\n", row_names[r_n].c_str(), v);
      }
    }
  }
  if (have_ranges) {
    fprintf(file, "RANGES\n");
    for (int r_n = 0; r_n < numRow; r_n++) {
      double v = ranges[r_n];
      if (v) {
        fprintf(file, "    RANGE     %-8s  %.15g\n", row_names[r_n].c_str(), v);
      }
    }
  }
  if (have_bounds) {
    fprintf(file, "BOUNDS\n");
    for (int c_n = 0; c_n < numCol; c_n++) {
      double lb = colLower[c_n];
      double ub = colUpper[c_n];
      bool discrete = false;
      if (have_int) discrete = integerColumn[c_n];
      if (Astart[c_n] == Astart[c_n + 1] && colCost[c_n] == 0) {
        // Possibly skip this column if it's zero and has no cost
        if (!highs_isInfinity(ub) || lb) {
          // Column would have a bound to report
          num_zero_no_cost_columns_in_bounds_section++;
        }
        if (write_zero_no_cost_columns) continue;
      }
      if (lb == ub) {
        // Equal lower and upper bounds: Fixed
        fprintf(file, " FX BOUND     %-8s  %.15g\n", col_names[c_n].c_str(),
                lb);
      } else if (highs_isInfinity(-lb) && highs_isInfinity(ub)) {
        // Infinite lower and upper bounds: Free
        fprintf(file, " FR BOUND     %-8s\n", col_names[c_n].c_str());
      } else {
        if (discrete) {
          if (lb == 0 && ub == 1) {
            // Binary
            fprintf(file, " BV BOUND     %-8s\n", col_names[c_n].c_str());
          } else {
            if (!highs_isInfinity(-lb)) {
              // Finite lower bound. No need to state this if LB is
              // zero unless UB is infinte
              if (lb || highs_isInfinity(ub))
                fprintf(file, " LI BOUND     %-8s  %.15g\n",
                        col_names[c_n].c_str(), lb);
            }
            if (!highs_isInfinity(ub)) {
              // Finite upper bound
              fprintf(file, " UI BOUND     %-8s  %.15g\n",
                      col_names[c_n].c_str(), ub);
            }
          }
        } else {
          if (!highs_isInfinity(-lb)) {
            // Lower bounded variable - default is 0
            if (lb) {
              fprintf(file, " LO BOUND     %-8s  %.15g\n",
                      col_names[c_n].c_str(), lb);
            }
          } else {
            // Infinite lower bound
            fprintf(file, " MI BOUND     %-8s\n", col_names[c_n].c_str());
          }
          if (!highs_isInfinity(ub)) {
            // Upper bounded variable
            fprintf(file, " UP BOUND     %-8s  %.15g\n", col_names[c_n].c_str(),
                    ub);
          }
        }
      }
    }
  }
  fprintf(file, "ENDATA\n");
  //#ifdef HiGHSDEV
  if (num_zero_no_cost_columns) {
    printf(
        "Model has %d zero columns with no costs: %d have finite upper bounds "
        "or nonzero lower bounds",
        num_zero_no_cost_columns, num_zero_no_cost_columns_in_bounds_section);
    if (write_zero_no_cost_columns) {
      printf(" and are written in MPS file\n");
    } else {
      printf(" and are not written in MPS file\n");
    }
  }
  //#endif
  fclose(file);
  return HighsStatus::OK;
}

inline const char* BoolToCharStar(bool b) { return b ? "True" : "False"; }
