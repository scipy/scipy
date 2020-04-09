/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/Avgas.cpp
 * @brief Utilities for tests with AVGAS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Avgas.h"

#include <cassert>
#include <cstdio>  // For printf

void Avgas::row(int row, int& num_row, int& num_row_nz,
                std::vector<double>& rowLower, std::vector<double>& rowUpper,
                std::vector<int>& ARstart, std::vector<int>& ARindex,
                std::vector<double>& ARvalue) {
  rowLower.resize(num_row + 1);
  rowUpper.resize(num_row + 1);
  ARstart.resize(num_row + 1);
  ARstart[num_row] = num_row_nz;
  if (row == 0) {
    rowLower[num_row] = -1;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 2;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 0;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 1;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 1) {
    rowLower[num_row] = -1;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 2;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 2;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 3;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 2) {
    rowLower[num_row] = -1;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 2;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 4;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 5;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 3) {
    rowLower[num_row] = -1;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 2;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 6;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 7;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 4) {
    rowLower[num_row] = -2;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 4;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 0;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 2;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 4;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 6;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 5) {
    rowLower[num_row] = -2;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 4;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 1;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 3;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 5;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 7;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 6) {
    rowLower[num_row] = 0;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 3;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 0;
    ARvalue[num_row_nz] = 2;
    num_row_nz++;
    ARindex[num_row_nz] = 2;
    ARvalue[num_row_nz] = 1;
    num_row_nz++;
    ARindex[num_row_nz] = 6;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 7) {
    rowLower[num_row] = 0;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 4;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 0;
    ARvalue[num_row_nz] = 5;
    num_row_nz++;
    ARindex[num_row_nz] = 2;
    ARvalue[num_row_nz] = 3;
    num_row_nz++;
    ARindex[num_row_nz] = 4;
    ARvalue[num_row_nz] = -3;
    num_row_nz++;
    ARindex[num_row_nz] = 6;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
  } else if (row == 8) {
    rowLower[num_row] = 0;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 4;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 1;
    ARvalue[num_row_nz] = 1;
    num_row_nz++;
    ARindex[num_row_nz] = 3;
    ARvalue[num_row_nz] = -1;
    num_row_nz++;
    ARindex[num_row_nz] = 5;
    ARvalue[num_row_nz] = -3;
    num_row_nz++;
    ARindex[num_row_nz] = 7;
    ARvalue[num_row_nz] = -5;
    num_row_nz++;
  } else if (row == 9) {
    rowLower[num_row] = 0;
    rowUpper[num_row] = 1e31;
    int num_new_nz = 3;
    ARindex.resize(num_row_nz + num_new_nz);
    ARvalue.resize(num_row_nz + num_new_nz);
    ARindex[num_row_nz] = 1;
    ARvalue[num_row_nz] = 1;
    num_row_nz++;
    ARindex[num_row_nz] = 5;
    ARvalue[num_row_nz] = -3;
    num_row_nz++;
    ARindex[num_row_nz] = 7;
    ARvalue[num_row_nz] = -2;
    num_row_nz++;
  } else {
    printf("Avgas: row %d out of range\n", row);
  }
  num_row++;
}

void Avgas::col(int col, int& num_col, int& num_col_nz,
                std::vector<double>& colCost, std::vector<double>& colLower,
                std::vector<double>& colUpper, std::vector<int>& Astart,
                std::vector<int>& Aindex, std::vector<double>& Avalue) {
  colCost.resize(num_col + 1);
  colLower.resize(num_col + 1);
  colUpper.resize(num_col + 1);
  Astart.resize(num_col + 1);
  Astart[num_col] = num_col_nz;
  int num_new_nz = 4;
  if (col == 0) {
    colCost[num_col] = 0;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 0;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 4;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 6;
    Avalue[num_col_nz] = 2;
    num_col_nz++;
    Aindex[num_col_nz] = 7;
    Avalue[num_col_nz] = 5;
    num_col_nz++;
  } else if (col == 1) {
    colCost[num_col] = -2;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 0;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 5;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 8;
    Avalue[num_col_nz] = 1;
    num_col_nz++;
    Aindex[num_col_nz] = 9;
    Avalue[num_col_nz] = 1;
    num_col_nz++;
  } else if (col == 2) {
    colCost[num_col] = -1;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 1;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 4;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 6;
    Avalue[num_col_nz] = 1;
    num_col_nz++;
    Aindex[num_col_nz] = 7;
    Avalue[num_col_nz] = 3;
    num_col_nz++;
  } else if (col == 3) {
    num_new_nz = 3;
    colCost[num_col] = -3;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 1;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 5;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 8;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
  } else if (col == 4) {
    num_new_nz = 3;
    colCost[num_col] = -2;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 2;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 4;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 7;
    Avalue[num_col_nz] = -3;
    num_col_nz++;
  } else if (col == 5) {
    colCost[num_col] = -4;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 2;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 5;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 8;
    Avalue[num_col_nz] = -3;
    num_col_nz++;
    Aindex[num_col_nz] = 9;
    Avalue[num_col_nz] = -3;
    num_col_nz++;
  } else if (col == 6) {
    colCost[num_col] = -3;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 3;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 4;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 6;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 7;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
  } else if (col == 7) {
    colCost[num_col] = -5;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_col_nz + num_new_nz);
    Avalue.resize(num_col_nz + num_new_nz);
    Aindex[num_col_nz] = 3;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 5;
    Avalue[num_col_nz] = -1;
    num_col_nz++;
    Aindex[num_col_nz] = 8;
    Avalue[num_col_nz] = -5;
    num_col_nz++;
    Aindex[num_col_nz] = 9;
    Avalue[num_col_nz] = -2;
    num_col_nz++;
  } else {
    printf("Avgas: col %d out of range\n", col);
  }
  num_col++;
}
