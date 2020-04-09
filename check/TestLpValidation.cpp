#include "Avgas.h"
#include "Highs.h"
#include "HighsIO.h"
#include "HighsLp.h"
#include "HighsLpUtils.h"
#include "HighsOptions.h"
#include "HighsSimplexInterface.h"
#include "HighsStatus.h"
#include "HighsTimer.h"
#include "catch.hpp"

// No commas in test case name.
TEST_CASE("LP-validation", "[highs_data]") {
  // Create an empty LP
  HighsLp lp;
  HighsOptions options;
  HighsTimer timer;
  HighsStatus return_status;
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

  return_status = assessLp(lp, options);
  REQUIRE(return_status == HighsStatus::OK);
  //  reportLp(lp, 2);

  const double my_infinity = 1e30;
  HighsModelObject hmo(lp, options, timer);
  HighsSimplexInterface hsi(hmo);

  return_status =
      hsi.addRows(num_row, &rowLower[0], &rowUpper[0], 0, NULL, NULL, NULL);
  //  printf("addRows: return_status = %s\n",
  //  HighsStatusToString(return_status).c_str());
  REQUIRE(return_status == HighsStatus::OK);
  //  reportLp(lp, 2);

  return_status = hsi.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              num_col_nz, &Astart[0], &Aindex[0], &Avalue[0]);
  //  printf("addCols: return_status = %s\n",
  //  HighsStatusToString(return_status).c_str());
  REQUIRE(return_status == HighsStatus::OK);
  //  reportLp(lp, 2);

  // Create an empty column
  int XnumNewCol = 1;
  int XnumNewNZ = 0;
  vector<double> XcolCost;
  XcolCost.resize(XnumNewCol);
  XcolCost[0] = 1;
  vector<double> XcolLower;
  XcolLower.resize(XnumNewCol);
  XcolLower[0] = 0;
  vector<double> XcolUpper;
  XcolUpper.resize(XnumNewCol);
  XcolUpper[0] = 1e25;
  vector<int> XAstart;
  XAstart.resize(XnumNewCol);
  vector<int> XAindex;
  vector<double> XAvalue;
  // Add an empty column
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::OK);
  XcolUpper[0] = my_infinity;
  //  reportLp(lp, 2);

  // Try to add a column with illegal cost
  XcolCost[0] = my_infinity;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::Error);
  XcolCost[0] = -my_infinity;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::Error);
  XcolCost[0] = 1;

  // Add a column with bound inconsistency due to upper
  XcolUpper[0] = -1;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::Warning);
  XcolUpper[0] = 0;

  // Add a column with bound inconsistency due to lower
  XcolLower[0] = 1;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::Warning);
  XcolLower[0] = 0;

  // Add a column with illegal bound due to lower
  XcolLower[0] = my_infinity;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::Error);
  XcolLower[0] = 0;

  // Add a column with illegal bound due to upper
  XcolUpper[0] = -my_infinity;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::Error);
  XcolUpper[0] = 0;

  // Add a legitimate column
  XcolLower[0] = 0;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], NULL, NULL);
  REQUIRE(return_status == HighsStatus::OK);

  //  reportLp(lp, 2);

  // Add a couple of non-empty columns with some small and large values
  XnumNewCol = 2;
  XnumNewNZ = 7;
  XcolCost.resize(XnumNewCol);
  XcolCost[0] = 1;
  XcolCost[1] = 2;
  XcolLower.resize(XnumNewCol);
  XcolLower[0] = 0;
  XcolLower[1] = 0;
  XcolUpper.resize(XnumNewCol);
  XcolUpper[0] = 1;
  XcolUpper[1] = 1;
  XAstart.resize(XnumNewCol + 1);
  XAindex.resize(XnumNewNZ);
  XAstart[1] = 4;
  XAstart[2] = XnumNewNZ;
  XAindex[0] = 0;
  XAindex[1] = 2;
  XAindex[2] = 3;
  XAindex[3] = 9;
  XAindex[4] = 1;
  XAindex[5] = 3;
  XAindex[6] = 8;
  XAvalue.resize(XnumNewNZ);
  XAvalue[0] = 1;
  XAvalue[1] = 1e-12;
  XAvalue[2] = -1e-20;
  XAvalue[3] = -1;
  XAvalue[4] = -1e60;
  XAvalue[5] = 1e100;
  XAvalue[6] = -1;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], &XAindex[0], &XAvalue[0]);
  REQUIRE(return_status == HighsStatus::Error);

  // Legitimise large matrix entries. Small entries now cause warning
  XAvalue[4] = -1;
  XAvalue[5] = 1;
  return_status =
      hsi.addCols(XnumNewCol, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumNewNZ, &XAstart[0], &XAindex[0], &XAvalue[0]);
  REQUIRE(return_status == HighsStatus::Warning);

  //  reportLp(lp, 2);

  Highs highs(options);

  HighsStatus init_status = highs.passModel(lp);
  REQUIRE(init_status == HighsStatus::OK);

  /*
  HighsStatus write_status =  highs.writeModel("");
  REQUIRE(write_status == HighsStatus::Warning);
  */

  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::OK);

  HighsModelStatus model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::PRIMAL_INFEASIBLE);
}
