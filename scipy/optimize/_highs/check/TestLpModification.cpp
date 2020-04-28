#include "Avgas.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLpUtils.h"

void HighsStatusReport(FILE* logfile, const char* message, HighsStatus status) {
  HighsLogMessage(logfile, HighsMessageType::INFO,
                  "%s: HighsStatus = %d - %s\n", message, (int)status,
                  HighsStatusToString(status).c_str());
}

bool areLpColEqual(const int num_col0, const double* colCost0,
                   const double* colLower0, const double* colUpper0,
                   const int num_nz0, const int* Astart0, const int* Aindex0,
                   const double* Avalue0, const int num_col1,
                   const double* colCost1, const double* colLower1,
                   const double* colUpper1, const int num_nz1,
                   const int* Astart1, const int* Aindex1,
                   const double* Avalue1, const double infinite_bound) {
  if (num_col0 != num_col1) {
    printf("areLpColEqual: %d = num_col0 != num_col1 = %d\n", num_col0,
           num_col1);
    return false;
  }
  if (!num_col0) return true;
  int num_col = num_col0;
  for (int col = 0; col < num_col; col++) {
    if (colCost0[col] != colCost1[col]) {
      printf("areLpColEqual: %g = colCost0[%d] != colCost1[%d] = %g\n",
             colCost0[col], col, col, colCost1[col]);
      return false;
    }
  }
  for (int col = 0; col < num_col; col++) {
    if (colLower0[col] <= -infinite_bound && colLower1[col] <= -infinite_bound)
      continue;
    if (colLower0[col] != colLower1[col]) {
      printf("areLpColEqual: %g = colLower0[%d] != colLower1[%d] = %g\n",
             colLower0[col], col, col, colLower1[col]);
      return false;
    }
    if (colUpper0[col] >= infinite_bound && colUpper1[col] >= infinite_bound)
      continue;
    if (colUpper0[col] != colUpper1[col]) {
      printf("areLpColEqual: %g = colUpper0[%d] != colUpper1[%d] = %g\n",
             colUpper0[col], col, col, colUpper1[col]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    printf("areLpColEqual: %d = num_nz0 != num_nz1 = %d\n", num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (int col = 0; col < num_col; col++) {
    if (Astart0[col] != Astart1[col]) {
      printf("areLpColEqual: %d = Astart0[%d] != Astart1[%d] = %d\n",
             Astart0[col], col, col, Astart1[col]);
      return false;
    }
  }
  int num_nz = num_nz0;
  for (int nz = 0; nz < num_nz; nz++) {
    if (Aindex0[nz] != Aindex1[nz]) {
      printf("areLpColEqual: %d = Aindex0[%d] != Aindex1[%d] = %d\n",
             Aindex0[nz], nz, nz, Aindex1[nz]);
      return false;
    }
    if (Avalue0[nz] != Avalue1[nz]) {
      printf("areLpColEqual: %g = Avalue0[%d] != Avalue1[%d] = %g\n",
             Avalue0[nz], nz, nz, Avalue1[nz]);
      return false;
    }
  }
  return true;
}

bool areLpRowEqual(const int num_row0, const double* rowLower0,
                   const double* rowUpper0, const int num_nz0,
                   const int* ARstart0, const int* ARindex0,
                   const double* ARvalue0, const int num_row1,
                   const double* rowLower1, const double* rowUpper1,
                   const int num_nz1, const int* ARstart1, const int* ARindex1,
                   const double* ARvalue1, const double infinite_bound) {
  if (num_row0 != num_row1) {
    printf("areLpRowEqual: %d = num_row0 != num_row1 = %d\n", num_row0,
           num_row1);
    return false;
  }
  if (!num_row0) return true;
  int num_row = num_row0;
  for (int row = 0; row < num_row; row++) {
    if (rowLower0[row] <= -infinite_bound && rowLower1[row] <= -infinite_bound)
      continue;
    if (rowLower0[row] != rowLower1[row]) {
      printf("areLpRowEqual: %g = rowLower0[%d] != rowLower1[%d] = %g\n",
             rowLower0[row], row, row, rowLower1[row]);
      return false;
    }
    if (rowUpper0[row] >= infinite_bound && rowUpper1[row] >= infinite_bound)
      continue;
    if (rowUpper0[row] != rowUpper1[row]) {
      printf("areLpRowEqual: %g = rowUpper0[%d] != rowUpper1[%d] = %g\n",
             rowUpper0[row], row, row, rowUpper1[row]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    printf("areLpRowEqual: %d = num_nz0 != num_nz1 = %d\n", num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (int row = 0; row < num_row; row++) {
    if (ARstart0[row] != ARstart1[row]) {
      printf("areLpRowEqual: %d = ARstart0[%d] != ARstart1[%d] = %d\n",
             ARstart0[row], row, row, ARstart1[row]);
      return false;
    }
  }
  int num_nz = num_nz0;
  for (int nz = 0; nz < num_nz; nz++) {
    if (ARindex0[nz] != ARindex1[nz]) {
      printf("areLpRowEqual: %d = ARindex0[%d] != ARindex1[%d] = %d\n",
             ARindex0[nz], nz, nz, ARindex1[nz]);
      return false;
    }
    if (ARvalue0[nz] != ARvalue1[nz]) {
      printf("areLpRowEqual: %g = ARvalue0[%d] != ARvalue1[%d] = %g\n",
             ARvalue0[nz], nz, nz, ARvalue1[nz]);
      return false;
    }
  }
  return true;
}

bool areLpEqual(const HighsLp lp0, const HighsLp lp1,
                const double infinite_bound) {
  bool return_bool;
  if (lp0.numCol_ > 0 && lp1.numCol_ > 0) {
    int lp0_num_nz = lp0.Astart_[lp0.numCol_];
    int lp1_num_nz = lp1.Astart_[lp1.numCol_];
    return_bool = areLpColEqual(
        lp0.numCol_, &lp0.colCost_[0], &lp0.colLower_[0], &lp0.colUpper_[0],
        lp0_num_nz, &lp0.Astart_[0], &lp0.Aindex_[0], &lp0.Avalue_[0],
        lp1.numCol_, &lp1.colCost_[0], &lp1.colLower_[0], &lp1.colUpper_[0],
        lp1_num_nz, &lp1.Astart_[0], &lp1.Aindex_[0], &lp1.Avalue_[0],
        infinite_bound);
    if (!return_bool) return return_bool;
  }
  if (lp0.numRow_ > 0 && lp1.numRow_ > 0) {
    int lp0_num_nz = 0;
    int lp1_num_nz = 0;
    return_bool = areLpRowEqual(
        lp0.numRow_, &lp0.rowLower_[0], &lp0.rowUpper_[0], lp0_num_nz, NULL,
        NULL, NULL, lp1.numRow_, &lp1.rowLower_[0], &lp1.rowUpper_[0],
        lp1_num_nz, NULL, NULL, NULL, infinite_bound);
  }
  return return_bool;
}

void test_delete_keep(const int row_dim, const bool interval,
                      const int from_row, const int to_row, const bool set,
                      const int num_set_entries, const int* row_set,
                      const bool mask, const int* row_mask) {
  int delete_from_row;
  int delete_to_row;
  int keep_from_row;
  int keep_to_row;
  int current_set_entry;
  if (interval) {
    printf("With index interval [%d, %d] in [%d, %d]\n", from_row, to_row, 0,
           row_dim - 1);
  } else if (set) {
    printf("With index set\n");
    for (int set = 0; set < num_set_entries; set++) printf(" %2d", set);
    printf("\n");
    for (int set = 0; set < num_set_entries; set++)
      printf(" %2d", row_set[set]);
    printf("\n");
  } else {
    printf("With index mask\n");
    for (int row = 0; row < row_dim; row++) printf(" %2d", row);
    printf("\n");
    for (int row = 0; row < row_dim; row++) printf(" %2d", row_mask[row]);
    printf("\n");
  }

  keep_from_row = 0;
  if (interval) {
    keep_to_row = from_row - 1;
  } else if (set) {
    current_set_entry = 0;
    keep_to_row = row_set[0] - 1;
  } else {
    keep_to_row = row_dim;
    for (int row = 0; row < row_dim; row++) {
      if (row_mask[row]) {
        keep_to_row = row - 1;
        break;
      }
    }
  }
  printf("Keep   [%2d, %2d]\n", 0, keep_to_row);
  if (keep_to_row >= row_dim - 1) return;
  for (int k = 0; k < row_dim; k++) {
    updateOutInIx(row_dim, interval, from_row, to_row, set, num_set_entries,
                  row_set, mask, row_mask, delete_from_row, delete_to_row,
                  keep_from_row, keep_to_row, current_set_entry);
    printf("Delete [%2d, %2d]; keep [%2d, %2d]\n", delete_from_row,
           delete_to_row, keep_from_row, keep_to_row);
    if (delete_to_row >= row_dim - 1 || keep_to_row >= row_dim - 1) break;
  }
}

bool test_all_delete_keep(int num_row) {
  // Test the extraction of intervals from interval, set and mask
  bool interval = false;
  bool set = false;
  bool mask = false;
  int row_dim = num_row;

  int from_row = 3;
  int to_row = 6;
  int num_set_entries = 4;
  int row_set[] = {1, 4, 5, 8};
  int row_mask[] = {0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
  int save_from_row = from_row;
  int save_row_set_0 = row_set[0];
  int save_row_mask_0 = row_mask[0];

  int to_pass = 2;  // 2
  for (int pass = 0; pass <= to_pass; pass++) {
    printf("\nTesting delete-keep: pass %d\n", pass);
    if (pass == 1) {
      // Mods to test LH limit behaviour
      from_row = 0;
      row_set[0] = 0;
      row_mask[0] = 1;
    } else if (pass == 2) {
      // Mods to test RH limit behaviour
      from_row = save_from_row;
      to_row = 9;
      row_set[0] = save_row_set_0;
      row_set[3] = 9;
      row_mask[0] = save_row_mask_0;
      row_mask[9] = 1;
    }

    interval = true;
    test_delete_keep(row_dim, interval, from_row, to_row, set, num_set_entries,
                     row_set, mask, row_mask);
    interval = false;
    set = true;
    test_delete_keep(row_dim, interval, from_row, to_row, set, num_set_entries,
                     row_set, mask, row_mask);
    set = false;
    mask = true;
    test_delete_keep(row_dim, interval, from_row, to_row, set, num_set_entries,
                     row_set, mask, row_mask);
  }
  return true;
}

void messageReportLp(const char* message, const HighsLp& lp) {
  HighsOptions options;
  options.output = stdout;
  options.message_level = ML_ALWAYS;
  HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                    "\nReporting LP: %s\n", message);
  reportLp(options, lp, 2);
}

void messageReportMatrix(const char* message, const int num_col,
                         const int num_nz, const int* start, const int* index,
                         const double* value) {
  HighsOptions options;
  options.output = stdout;
  options.message_level = ML_ALWAYS;
  HighsPrintMessage(options.output, options.message_level, ML_VERBOSE,
                    "\nReporting Matrix: %s\n", message);
  reportMatrix(options, message, num_col, num_nz, start, index, value);
}

// No commas in test case name.
TEST_CASE("LP-modification", "[highs_data]") {
  test_all_delete_keep(10);

  HighsOptions options;
  options.message_level = ML_ALWAYS;

  Avgas avgas;
  const int avgas_num_col = 8;
  const int avgas_num_row = 10;
  int num_row = 0;
  int num_row_nz = 0;
  vector<double> rowLower;
  vector<double> rowUpper;
  vector<int> ARstart;
  vector<int> ARindex;
  vector<double> ARvalue;

  for (int row = 0; row < avgas_num_row; row++) {
    avgas.row(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex,
              ARvalue);
  }

  int num_col = 0;
  int num_col_nz = 0;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  for (int col = 0; col < avgas_num_col; col++) {
    avgas.col(col, num_col, num_col_nz, colCost, colLower, colUpper, Astart,
              Aindex, Avalue);
  }

  bool return_bool;
  HighsStatus return_status;
  HighsModelStatus model_status;
  std::string message;

  // Create two empty LPs: one to be initialised as AVGAS by adding
  // all the columns and rows separately, the other to be built by
  // adding piecemeal.
  HighsLp avgas_lp;
  HighsLp lp;

  Highs avgas_highs(options);
  return_status = avgas_highs.passModel(avgas_lp);
  //  printf("passModel: return_status = %s\n",
  //  HighsStatusToString(return_status).c_str());
  REQUIRE(return_status == HighsStatus::OK);

  return_bool = avgas_highs.addCols(num_col, &colCost[0], &colLower[0],
                                    &colUpper[0], 0, NULL, NULL, NULL);
  REQUIRE(return_bool);
  return_bool =
      avgas_highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                          &ARstart[0], &ARindex[0], &ARvalue[0]);
  REQUIRE(return_bool);

  return_status = avgas_highs.writeModel("");
  HighsStatusReport(options.logfile, "avgas_highs.writeModel(\"\")",
                    return_status);
  REQUIRE(return_status == HighsStatus::OK);

  Highs highs(options);
  return_status = highs.passModel(lp);
  //  printf("passModel: return_status = %s\n",
  //  HighsStatusToString(return_status).c_str());
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::NOTSET);

  return_status = highs.run();
  HighsStatusReport(options.logfile, "highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::MODEL_EMPTY);

  // Adding column vectors and matrix to model with no rows returns an error
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              num_col_nz, &Astart[0], &Aindex[0], &Avalue[0]);
  REQUIRE(!return_bool);

  // Adding column vectors to model with no rows returns OK
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              0, NULL, NULL, NULL);
  REQUIRE(return_bool);

  return_status = highs.writeModel("");
  HighsStatusReport(options.logfile, "highs.writeModel(\"\")", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  // Adding row vectors and matrix to model with columns returns OK
  return_bool = highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                              &ARstart[0], &ARindex[0], &ARvalue[0]);
  REQUIRE(return_bool);

  return_status = highs.writeModel("");
  HighsStatusReport(options.logfile, "highs.writeModel(\"\")", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  //  const HighsLp &reference_avgas = avgas_highs.getLp();
  //  const HighsLp &reference_lp = highs.getLp();

  return_bool =
      areLpEqual(highs.getLp(), avgas_highs.getLp(), options.infinite_bound);
  REQUIRE(return_bool);

  return_status = highs.run();
  HighsStatusReport(options.logfile, "highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  double avgas_optimal_objective_value;
  highs.getHighsInfoValue("objective_function_value",
                          avgas_optimal_objective_value);
  double optimal_objective_value;

#ifdef HiGHSDEV
  //  const HighsSolution& solution = highs.getSolution();
  //  const HighsBasis& basis = highs.getBasis();
  highs.reportModelStatusSolutionBasis("After avgas solve");
#endif

  // Getting columns from the LP is OK
  int col1357_col_mask[] = {0, 1, 0, 1, 0, 1, 0, 1};
  int col1357_col_set[] = {1, 3, 5, 7};
  int col1357_illegal_col_set[] = {3, 7, 1, 5};
  int col1357_num_ix = 4;
  int col1357_num_col;
  int col1357_num_nz;
  double* col1357_cost = (double*)malloc(sizeof(double) * col1357_num_ix);
  double* col1357_lower = (double*)malloc(sizeof(double) * col1357_num_ix);
  double* col1357_upper = (double*)malloc(sizeof(double) * col1357_num_ix);
  int* col1357_start = (int*)malloc(sizeof(int) * col1357_num_ix);
  int* col1357_index = (int*)malloc(sizeof(int) * num_col_nz);
  double* col1357_value = (double*)malloc(sizeof(double) * num_col_nz);

  return_bool = highs.getCols(3, 6, col1357_num_col, col1357_cost,
                              col1357_lower, col1357_upper, col1357_num_nz,
                              col1357_start, col1357_index, col1357_value);
  REQUIRE(return_bool == true);

  // Calling getCols using an unordered set should be OK but for now HiGHS
  // returns an error
  return_bool =
      highs.getCols(col1357_num_ix, col1357_illegal_col_set, col1357_num_col,
                    col1357_cost, col1357_lower, col1357_upper, col1357_num_nz,
                    col1357_start, col1357_index, col1357_value);
  REQUIRE(!return_bool);

  return_bool =
      highs.getCols(col1357_num_ix, col1357_col_set, col1357_num_col,
                    col1357_cost, col1357_lower, col1357_upper, col1357_num_nz,
                    col1357_start, col1357_index, col1357_value);
  REQUIRE(return_bool);

  return_bool = highs.getCols(col1357_col_mask, col1357_num_col, col1357_cost,
                              col1357_lower, col1357_upper, col1357_num_nz,
                              col1357_start, col1357_index, col1357_value);
  REQUIRE(return_bool);

  // Try to delete an empty range of cols: OK
  return_bool = highs.deleteCols(0, -1);
  REQUIRE(return_bool);

  // Try to delete more cols than there are: ERROR
  return_bool = highs.deleteCols(0, num_col + 1);
  REQUIRE(!return_bool);

  return_bool = highs.deleteCols(col1357_num_ix, col1357_col_set);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleting columns 1, 3, 5, 7";
  //  messageReportLp(message.c_str(), highs.getLp());
  //  printf("%s", message.c_str()); reportLp(highs.getLp(), 2);
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_bool = highs.addCols(col1357_num_col, col1357_cost, col1357_lower,
                              col1357_upper, col1357_num_nz, col1357_start,
                              col1357_index, col1357_value);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After restoring columns 1, 3, 5, 7\n";
  //  printf("%s", message.c_str()); reportLp(highs.getLp(), 2);
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_status = highs.run();
  HighsStatusReport(options.logfile, "highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  highs.getHighsInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After re-solving");
#endif

  // Delete all the columns: OK
  return_bool = highs.deleteCols(0, num_col - 1);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleting all columns";
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Delete all the rows: OK
  return_bool = highs.deleteRows(0, num_row - 1);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleting all rows";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Adding column vectors to model with no rows returns OK
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              0, NULL, NULL, NULL);
  REQUIRE(return_bool);

  message = "With columns but no rows";
  //  messageReportLp(message.c_str(), highs.getLp());

  // Adding row vectors and matrix to model with columns returns OK
  return_bool = highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                              &ARstart[0], &ARindex[0], &ARvalue[0]);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "With columns but and rows";
  //   messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Getting rows from the LP is OK
  int from_row_ix = 0;
  int to_row_ix = 3;
  int row0135789_row_set[] = {0, 1, 3, 5, 7, 8, 9};
  int row0135789_row_mask[] = {1, 1, 0, 1, 0, 1, 0, 1, 1, 1};
  int row0135789_num_ix = 7;
  int row0135789_num_row;
  int row0135789_num_nz;
  double* row0135789_lower =
      (double*)malloc(sizeof(double) * row0135789_num_ix);
  double* row0135789_upper =
      (double*)malloc(sizeof(double) * row0135789_num_ix);
  int* row0135789_start = (int*)malloc(sizeof(int) * row0135789_num_ix);
  int* row0135789_index = (int*)malloc(sizeof(int) * num_row_nz);
  double* row0135789_value = (double*)malloc(sizeof(double) * num_row_nz);

  return_bool =
      highs.getRows(from_row_ix, to_row_ix, row0135789_num_row,
                    row0135789_lower, row0135789_upper, row0135789_num_nz,
                    row0135789_start, row0135789_index, row0135789_value);
  REQUIRE(return_bool);

  //  messageReportMatrix("Get by interval\nRow   ", row0135789_num_row,
  //  row0135789_num_nz, row0135789_start, row0135789_index, row0135789_value);

  return_bool =
      highs.getRows(row0135789_num_ix, row0135789_row_set, row0135789_num_row,
                    row0135789_lower, row0135789_upper, row0135789_num_nz,
                    row0135789_start, row0135789_index, row0135789_value);
  REQUIRE(return_bool);

  //  messageReportMatrix("Get by set\nRow   ", row0135789_num_row,
  //  row0135789_num_nz, row0135789_start, row0135789_index, row0135789_value);

  return_bool =
      highs.getRows(row0135789_row_mask, row0135789_num_row, row0135789_lower,
                    row0135789_upper, row0135789_num_nz, row0135789_start,
                    row0135789_index, row0135789_value);
  REQUIRE(return_bool);

  //  messageReportMatrix("Get by mask\nRow   ", row0135789_num_row,
  //  row0135789_num_nz, row0135789_start, row0135789_index, row0135789_value);

  return_bool =
      highs.getRows(row0135789_num_ix, row0135789_row_set, row0135789_num_row,
                    row0135789_lower, row0135789_upper, row0135789_num_nz,
                    row0135789_start, row0135789_index, row0135789_value);
  REQUIRE(return_bool);

  return_bool = highs.deleteRows(row0135789_num_ix, row0135789_row_set);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleting rows 0-1, 3, 5, 7-9";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  int row012_row_set[] = {0, 1, 2};
  int row012_row_mask[] = {1, 1, 1};
  int row012_num_ix = 3;
  int row012_num_row;
  int row012_num_nz;
  double* row012_lower = (double*)malloc(sizeof(double) * row012_num_ix);
  double* row012_upper = (double*)malloc(sizeof(double) * row012_num_ix);
  int* row012_start = (int*)malloc(sizeof(int) * row012_num_ix);
  int* row012_index = (int*)malloc(sizeof(int) * num_row_nz);
  double* row012_value = (double*)malloc(sizeof(double) * num_row_nz);

  return_bool = highs.getRows(row012_num_ix, row012_row_set, row012_num_row,
                              row012_lower, row012_upper, row012_num_nz,
                              row012_start, row012_index, row012_value);
  REQUIRE(return_bool);

  return_bool = highs.deleteRows(row012_row_mask);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleting rows 0-2";
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Delete all the columns: OK
  return_bool = highs.deleteCols(0, num_col - 1);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleting all columns";
  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Can't add rows with no columns
  return_bool = highs.addRows(
      row0135789_num_row, row0135789_lower, row0135789_upper, row0135789_num_nz,
      row0135789_start, row0135789_index, row0135789_value);
  REQUIRE(!return_bool);

  // Adding column vectors to model with no rows returns OK
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              0, NULL, NULL, NULL);
  REQUIRE(return_bool);

  return_bool = highs.addRows(
      row0135789_num_row, row0135789_lower, row0135789_upper, row0135789_num_nz,
      row0135789_start, row0135789_index, row0135789_value);
  REQUIRE(return_bool);

  return_bool =
      highs.addRows(row012_num_row, row012_lower, row012_upper, row012_num_nz,
                    row012_start, row012_index, row012_value);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After restoring all rows";
  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_status = highs.run();
  HighsStatusReport(options.logfile, "highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  highs.getHighsInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

#ifdef HiGHSDEV
  message = "After resolve";
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Try to delete an empty range of rows: OK
  return_bool = highs.deleteRows(0, -1);
  REQUIRE(return_bool);

  // Try to delete more rows than there are: ERROR
  return_bool = highs.deleteRows(0, num_row);
  REQUIRE(!return_bool);

  return_bool = highs.getCols(col1357_col_mask, col1357_num_col, col1357_cost,
                              col1357_lower, col1357_upper, col1357_num_nz,
                              col1357_start, col1357_index, col1357_value);
  REQUIRE(return_bool);

  return_bool = highs.deleteCols(col1357_num_ix, col1357_col_set);
  REQUIRE(return_bool);

  int col0123_col_mask[] = {1, 1, 1, 1};
  //  int col0123_col_set[] = {0, 1, 2, 3};
  int col0123_num_ix = 4;
  int col0123_num_col;
  int col0123_num_nz;
  double* col0123_cost = (double*)malloc(sizeof(double) * col0123_num_ix);
  double* col0123_lower = (double*)malloc(sizeof(double) * col0123_num_ix);
  double* col0123_upper = (double*)malloc(sizeof(double) * col0123_num_ix);
  int* col0123_start = (int*)malloc(sizeof(int) * col0123_num_ix);
  int* col0123_index = (int*)malloc(sizeof(int) * num_col_nz);
  double* col0123_value = (double*)malloc(sizeof(double) * num_col_nz);

  return_bool = highs.getCols(col0123_col_mask, col0123_num_col, col0123_cost,
                              col0123_lower, col0123_upper, col0123_num_nz,
                              col0123_start, col0123_index, col0123_value);
  REQUIRE(return_bool);
  //  messageReportMatrix("Get col1357 by mask\nRow   ", col1357_num_col,
  //  col1357_num_nz, col1357_start, col1357_index, col1357_value);
  //  messageReportMatrix("Get col0123 by mask\nRow   ", col0123_num_col,
  //  col0123_num_nz, col0123_start, col0123_index, col0123_value);

  return_bool = highs.deleteRows(0, num_row - 1);
  REQUIRE(return_bool);

  return_bool = highs.deleteCols(col0123_col_mask);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleting all rows and columns";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Adding row vectors to model with no columns returns OK
  return_bool = highs.addRows(row0135789_num_row, row0135789_lower,
                              row0135789_upper, 0, NULL, NULL, NULL);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After restoring 7 rows";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_bool = highs.addRows(row012_num_row, row012_lower, row012_upper, 0,
                              row012_start, row012_index, row012_value);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After restoring all rows";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_bool = highs.addCols(col1357_num_col, col1357_cost, col1357_lower,
                              col1357_upper, col1357_num_nz, col1357_start,
                              col1357_index, col1357_value);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After restoring columns 1, 3, 5, 7";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_status = highs.run();
  HighsStatusReport(options.logfile, "highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

#ifdef HiGHSDEV
  message = "After solving after restoring all rows and columns 1, 3, 5, 7";
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_bool = highs.addCols(col0123_num_col, col0123_cost, col0123_lower,
                              col0123_upper, col0123_num_nz, col0123_start,
                              col0123_index, col0123_value);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After restoring columns 0-3";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  return_status = highs.run();
  HighsStatusReport(options.logfile, "highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  highs.getHighsInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

  return_bool = highs.deleteRows(0, num_row - 1);
  REQUIRE(return_bool);

  return_bool = highs.deleteCols(0, num_col - 1);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "After deleteing all rows and columns";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Adding column vectors to model with no rows returns OK
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              0, NULL, NULL, NULL);
  REQUIRE(return_bool);

#ifdef HiGHSDEV
  message = "With columns but no rows";
  //  messageReportLp(message.c_str(), highs.getLp());
  highs.reportModelStatusSolutionBasis(message);
#endif

  // Adding row vectors and matrix to model with columns returns OK
  return_bool = highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                              &ARstart[0], &ARindex[0], &ARvalue[0]);
  REQUIRE(return_bool);

  col1357_cost[0] = 2.01;
  col1357_cost[1] = 2.31;
  col1357_cost[2] = 2.51;
  col1357_cost[3] = 2.71;
  col1357_lower[0] = 0.01;
  col1357_lower[1] = 0.31;
  col1357_lower[2] = 0.51;
  col1357_lower[3] = 0.71;
  col1357_upper[0] = 1.01;
  col1357_upper[1] = 1.31;
  col1357_upper[2] = 1.51;
  col1357_upper[3] = 1.71;

  row0135789_lower[0] = -0.01;
  row0135789_lower[1] = -0.11;
  row0135789_lower[2] = -0.31;
  row0135789_lower[3] = -0.51;
  row0135789_lower[4] = -0.71;
  row0135789_lower[5] = -0.81;
  row0135789_lower[6] = -0.91;
  row0135789_upper[0] = 3.01;
  row0135789_upper[1] = 3.11;
  row0135789_upper[2] = 3.31;
  row0135789_upper[3] = 3.51;
  row0135789_upper[4] = 3.71;
  row0135789_upper[5] = 3.81;
  row0135789_upper[6] = 3.91;

  // Attempting to set a cost to infinity returns error
  return_bool = highs.changeColCost(7, HIGHS_CONST_INF);
  REQUIRE(!return_bool);

  // Attempting to set a cost to a finite value returns OK
  return_bool = highs.changeColCost(7, 77);
  REQUIRE(return_bool);

  return_bool =
      highs.changeColsCost(col1357_num_ix, col1357_col_set, col1357_cost);
  REQUIRE(return_bool);

  // Attempting to set row bounds with infinite lower bound returns error
  return_bool = highs.changeRowBounds(2, HIGHS_CONST_INF, 3.21);
  REQUIRE(!return_bool);

  return_bool = highs.changeRowBounds(2, -HIGHS_CONST_INF, 3.21);
  REQUIRE(return_bool);

  // Attempting to set col bounds with -infinite upper bound returns error
  return_bool = highs.changeColBounds(2, 0.21, -HIGHS_CONST_INF);
  REQUIRE(!return_bool);

  return_bool = highs.changeColBounds(2, 0.21, HIGHS_CONST_INF);
  REQUIRE(return_bool);

  return_bool = highs.changeRowsBounds(row0135789_num_ix, row0135789_row_set,
                                       row0135789_lower, row0135789_upper);
  REQUIRE(return_bool);

  return_bool = highs.changeColsBounds(col1357_num_ix, col1357_col_set,
                                       col1357_lower, col1357_upper);
  REQUIRE(return_bool);

  //  messageReportLp("After changing costs and bounds", highs.getLp());

  // Return the LP to its original state with a mask

  return_bool = highs.changeColsCost(col1357_col_mask, &colCost[0]);
  REQUIRE(return_bool);

  return_bool = highs.changeColBounds(2, colLower[2], colUpper[2]);
  REQUIRE(return_bool);

  return_bool =
      highs.changeColsBounds(col1357_col_mask, &colLower[0], &colUpper[0]);
  REQUIRE(return_bool);

  return_bool =
      highs.changeRowsBounds(row0135789_row_mask, &rowLower[0], &rowUpper[0]);
  REQUIRE(return_bool);

  return_bool = highs.changeRowBounds(2, rowLower[2], rowUpper[2]);
  REQUIRE(return_bool);

  return_bool =
      areLpEqual(avgas_highs.getLp(), highs.getLp(), options.infinite_bound);
  REQUIRE(return_bool);

  int before_num_col;
  int after_num_col;
  int rm_col;
  int before_num_row;
  int after_num_row;
  int rm_row;

  before_num_col = highs.getNumCols();
  rm_col = 0;
  return_bool = highs.deleteCols(rm_col, rm_col);
  REQUIRE(return_bool);
  after_num_col = highs.getNumCols();
  printf("After removing col %d / %d have %d cols\n", rm_col, before_num_col,
         after_num_col);
  REQUIRE(after_num_col == before_num_col - 1);

  before_num_row = highs.getNumRows();
  rm_row = 0;
  return_bool = highs.deleteRows(rm_row, rm_row);
  REQUIRE(return_bool);
  after_num_row = highs.getNumRows();
  printf("After removing row %d / %d have %d rows\n", rm_row, before_num_row,
         after_num_row);
  REQUIRE(after_num_row == before_num_row - 1);

  before_num_col = highs.getNumCols();
  rm_col = before_num_col - 1;
  return_bool = highs.deleteCols(rm_col, rm_col);
  REQUIRE(return_bool);
  after_num_col = highs.getNumCols();
  printf("After removing col %d / %d have %d cols\n", rm_col, before_num_col,
         after_num_col);
  REQUIRE(after_num_col == before_num_col - 1);

  before_num_row = highs.getNumRows();
  rm_row = before_num_row - 1;
  return_bool = highs.deleteRows(rm_row, rm_row);
  REQUIRE(return_bool);
  after_num_row = highs.getNumRows();
  printf("After removing row %d / %d have %d rows\n", rm_row, before_num_row,
         after_num_row);
  REQUIRE(after_num_row == before_num_row - 1);

  //  messageReportLp("After deleting all rows and columns", highs.getLp());

  //  messageReportLp("After restoring costs and bounds", highs.getLp());
  printf("Finished successfully\n");
  fflush(stdout);
}
