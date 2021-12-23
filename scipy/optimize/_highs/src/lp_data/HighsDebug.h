/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsDebug.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSDEBUG_H_
#define SIMPLEX_HIGHSDEBUG_H_

#include <string>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsStatus.h"

//#include "lp_data/HighsOptions.h"

HighsStatus debugDebugToHighsStatus(const HighsDebugStatus debug_status);

HighsDebugStatus debugWorseStatus(const HighsDebugStatus status0,
                                  const HighsDebugStatus status1);

bool debugVectorRightSize(const std::vector<double> v, const int right_size);

bool debugVectorRightSize(const std::vector<int> v, const int right_size);

#endif  // SIMPLEX_HIGHSDEBUG_H_
