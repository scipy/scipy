/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsMatrixPic.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSMATRIXPIC_H_
#define UTIL_HIGHSMATRIXPIC_H_

#include <string>
#include <vector>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"

HighsStatus writeLpMatrixPicToFile(const HighsOptions& options,
                                   const std::string fileprefix,
                                   const HighsLp& lp);

HighsStatus writeMatrixPicToFile(const HighsOptions& options,
                                 const std::string fileprefix, const int numRow,
                                 const int numCol,
                                 const std::vector<int>& Astart,
                                 const std::vector<int>& Aindex);

HighsStatus writeRmatrixPicToFile(const HighsOptions& options,
                                  const std::string fileprefix,
                                  const int numRow, const int numCol,
                                  const std::vector<int>& ARstart,
                                  const std::vector<int>& ARindex);

#endif  // UTIL_HIGHSMATRIXPIC_H_
