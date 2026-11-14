/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsMatrixPic.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "util/HighsMatrixPic.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>

HighsStatus writeLpMatrixPicToFile(const HighsOptions& options,
                                   const std::string& fileprefix,
                                   const HighsLp& lp) {
  return writeMatrixPicToFile(options, fileprefix, lp.num_row_, lp.num_col_,
                              lp.a_matrix_.start_, lp.a_matrix_.index_);
}

HighsStatus writeHessianPicToFile(const HighsOptions& options,
                                  const std::string& fileprefix,
                                  const HighsHessian& hessian) {
  return writeMatrixPicToFile(options, fileprefix, hessian.dim_, hessian.dim_,
                              hessian.start_, hessian.index_);
}

HighsStatus writeMatrixPicToFile(const HighsOptions& options,
                                 const std::string& fileprefix,
                                 const HighsInt numRow, const HighsInt numCol,
                                 const std::vector<HighsInt>& Astart,
                                 const std::vector<HighsInt>& Aindex) {
  if (numRow == 0 || numCol == 0) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Cannot generate image of matrix with a zero dimension\n");
    return HighsStatus::kError;
  }
  assert(numRow > 0);
  assert(numCol > 0);
  std::vector<HighsInt> ARlength;
  std::vector<HighsInt> ARstart;
  std::vector<HighsInt> ARindex;
  const HighsInt numNz = Astart[numCol];
  ARlength.assign(numRow, 0);
  ARstart.resize(numRow + 1);
  ARindex.resize(numNz);
  for (HighsInt iEl = 0; iEl < numNz; iEl++) ARlength[Aindex[iEl]]++;
  ARstart[0] = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++)
    ARstart[iRow + 1] = ARstart[iRow] + ARlength[iRow];
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    for (HighsInt iEl = Astart[iCol]; iEl < Astart[iCol + 1]; iEl++) {
      HighsInt iRow = Aindex[iEl];
      ARindex[ARstart[iRow]++] = iCol;
    }
  }
  ARstart[0] = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++)
    ARstart[iRow + 1] = ARstart[iRow] + ARlength[iRow];

  return writeRmatrixPicToFile(options, fileprefix, numRow, numCol, ARstart,
                               ARindex);
}

HighsStatus writeRmatrixPicToFile(const HighsOptions& options,
                                  const std::string& fileprefix,
                                  const HighsInt numRow, const HighsInt numCol,
                                  const std::vector<HighsInt>& ARstart,
                                  const std::vector<HighsInt>& ARindex) {
  if (fileprefix == "") return HighsStatus::kError;
  std::string filename = fileprefix + ".pbm";
  std::ofstream f;
  f.open(filename, std::ios::out);
  const HighsInt border_width = 1;
  const HighsInt max_num_pixel_wide = 1600;
  const HighsInt max_num_pixel_deep = 900;
  const HighsInt max_num_matrix_pixel_wide =
      max_num_pixel_wide - 2 * border_width;
  const HighsInt max_num_matrix_pixel_deep =
      max_num_pixel_deep - 2 * border_width;
  HighsInt num_col_per_pixel = 1;
  HighsInt num_row_per_pixel = 1;
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
  const HighsInt dim_per_pixel = std::max(num_col_per_pixel, num_row_per_pixel);
  HighsInt num_pixel_wide = numCol / dim_per_pixel;
  if (dim_per_pixel * num_pixel_wide < numCol) num_pixel_wide++;
  HighsInt num_pixel_deep = numRow / dim_per_pixel;
  if (dim_per_pixel * num_pixel_deep < numRow) num_pixel_deep++;
  // Account for the borders
  num_pixel_wide += 2;
  num_pixel_deep += 2;
  assert(num_pixel_wide <= max_num_pixel_wide);
  assert(num_pixel_deep <= max_num_pixel_deep);

  highsLogUser(options.log_options, HighsLogType::kInfo,
               "Representing matrix sparsity pattern %" HIGHSINT_FORMAT
               "x%" HIGHSINT_FORMAT
               " .pbm file,"
               " mapping entries in square of size %" HIGHSINT_FORMAT
               " onto one pixel\n",
               num_pixel_wide, num_pixel_deep, dim_per_pixel);

  std::vector<HighsInt> value;
  value.assign(num_pixel_wide, 0);
  f << "P1" << std::endl;
  f << num_pixel_wide << " " << num_pixel_deep << std::endl;
  HighsInt pic_num_row = 0;
  // Top border
  for (HighsInt pixel = 0; pixel < num_pixel_wide; pixel++) f << "1 ";
  f << std::endl;
  pic_num_row++;
  HighsInt from_row = 0;
  for (;;) {
    HighsInt to_row = std::min(from_row + dim_per_pixel, numRow);
    for (HighsInt iRow = from_row; iRow < to_row; iRow++) {
      for (HighsInt iEl = ARstart[iRow]; iEl < ARstart[iRow + 1]; iEl++) {
        HighsInt iCol = ARindex[iEl];
        HighsInt pixel = iCol / dim_per_pixel;
        assert(pixel < num_pixel_wide - 2);
        value[pixel] = 1;
      }
    }
    // LH border
    f << "1 ";
    for (HighsInt pixel = 0; pixel < num_pixel_wide - 2; pixel++)
      f << value[pixel] << " ";
    // LH border
    f << "1 " << std::endl;
    pic_num_row++;
    for (HighsInt pixel = 0; pixel < num_pixel_wide - 2; pixel++)
      value[pixel] = 0;
    if (to_row == numRow) break;
    from_row = to_row;
  }

  // Bottom border
  for (HighsInt pixel = 0; pixel < num_pixel_wide; pixel++) f << "1 ";
  f << std::endl;
  pic_num_row++;
  assert(pic_num_row == num_pixel_deep);

  return HighsStatus::kOk;
}
