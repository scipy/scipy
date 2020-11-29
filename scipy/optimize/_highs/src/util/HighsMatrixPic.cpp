/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsMatrixPic.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "util/HighsMatrixPic.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>

HighsStatus writeLpMatrixPicToFile(const HighsOptions& options,
                                   const std::string fileprefix,
                                   const HighsLp& lp) {
  return writeMatrixPicToFile(options, fileprefix, lp.numRow_, lp.numCol_,
                              lp.Astart_, lp.Aindex_);
}

HighsStatus writeMatrixPicToFile(const HighsOptions& options,
                                 const std::string fileprefix, const int numRow,
                                 const int numCol,
                                 const std::vector<int>& Astart,
                                 const std::vector<int>& Aindex) {
  std::vector<int> ARlength;
  std::vector<int> ARstart;
  std::vector<int> ARindex;
  assert(numRow > 0);
  assert(numCol > 0);
  const int numNz = Astart[numCol];
  ARlength.assign(numRow, 0);
  ARstart.resize(numRow + 1);
  ARindex.resize(numNz);
  for (int iEl = 0; iEl < numNz; iEl++) ARlength[Aindex[iEl]]++;
  ARstart[0] = 0;
  for (int iRow = 0; iRow < numRow; iRow++)
    ARstart[iRow + 1] = ARstart[iRow] + ARlength[iRow];
  for (int iCol = 0; iCol < numCol; iCol++) {
    for (int iEl = Astart[iCol]; iEl < Astart[iCol + 1]; iEl++) {
      int iRow = Aindex[iEl];
      ARindex[ARstart[iRow]++] = iCol;
    }
  }
  ARstart[0] = 0;
  for (int iRow = 0; iRow < numRow; iRow++)
    ARstart[iRow + 1] = ARstart[iRow] + ARlength[iRow];

  return writeRmatrixPicToFile(options, fileprefix, numRow, numCol, ARstart,
                               ARindex);
}

HighsStatus writeRmatrixPicToFile(const HighsOptions& options,
                                  const std::string fileprefix,
                                  const int numRow, const int numCol,
                                  const std::vector<int>& ARstart,
                                  const std::vector<int>& ARindex) {
  if (fileprefix == "") return HighsStatus::Error;
  std::string filename = fileprefix + ".pbm";
  std::ofstream f;
  f.open(filename, std::ios::out);
  const int border_width = 1;
  const int max_num_pixel_wide = 1600;
  const int max_num_pixel_deep = 900;
  const int max_num_matrix_pixel_wide = max_num_pixel_wide - 2 * border_width;
  const int max_num_matrix_pixel_deep = max_num_pixel_deep - 2 * border_width;
  int num_col_per_pixel = 1;
  int num_row_per_pixel = 1;
  if (numCol > max_num_matrix_pixel_wide) {
    num_col_per_pixel = numCol / max_num_matrix_pixel_wide;
    if (num_col_per_pixel * max_num_matrix_pixel_wide < numCol)
      num_col_per_pixel++;
  }
  if (numRow > max_num_matrix_pixel_deep) {
    num_row_per_pixel = numRow / max_num_matrix_pixel_deep;
    if (num_row_per_pixel * max_num_matrix_pixel_deep < numRow)
      num_row_per_pixel++;
  }
  const int dim_per_pixel = std::max(num_col_per_pixel, num_row_per_pixel);
  int num_pixel_wide = numCol / dim_per_pixel;
  if (dim_per_pixel * num_pixel_wide < numCol) num_pixel_wide++;
  int num_pixel_deep = numRow / dim_per_pixel;
  if (dim_per_pixel * num_pixel_deep < numRow) num_pixel_deep++;
  // Account for the borders
  num_pixel_wide += 2;
  num_pixel_deep += 2;
  assert(num_pixel_wide <= max_num_pixel_wide);
  assert(num_pixel_deep <= max_num_pixel_deep);

  HighsLogMessage(
      options.logfile, HighsMessageType::INFO,
      "Representing LP constraint matrix sparsity pattern %dx%d .pbm file,"
      " mapping entries in square of size %d onto one pixel",
      num_pixel_wide, num_pixel_deep, dim_per_pixel);

  std::vector<int> value;
  value.assign(num_pixel_wide, 0);
  f << "P1" << std::endl;
  f << num_pixel_wide << " " << num_pixel_deep << std::endl;
  int pic_num_row = 0;
  // Top border
  for (int pixel = 0; pixel < num_pixel_wide; pixel++) f << "1 ";
  f << std::endl;
  pic_num_row++;
  int from_row = 0;
  for (;;) {
    int to_row = std::min(from_row + dim_per_pixel, numRow);
    for (int iRow = from_row; iRow < to_row; iRow++) {
      for (int iEl = ARstart[iRow]; iEl < ARstart[iRow + 1]; iEl++) {
        int iCol = ARindex[iEl];
        int pixel = iCol / dim_per_pixel;
        assert(pixel < num_pixel_wide - 2);
        value[pixel] = 1;
      }
    }
    // LH border
    f << "1 ";
    for (int pixel = 0; pixel < num_pixel_wide - 2; pixel++)
      f << value[pixel] << " ";
    // LH border
    f << "1 " << std::endl;
    pic_num_row++;
    for (int pixel = 0; pixel < num_pixel_wide - 2; pixel++) value[pixel] = 0;
    if (to_row == numRow) break;
    from_row = to_row;
  }

  // Bottom border
  for (int pixel = 0; pixel < num_pixel_wide; pixel++) f << "1 ";
  f << std::endl;
  pic_num_row++;
  assert(pic_num_row == num_pixel_deep);

  return HighsStatus::OK;
}
