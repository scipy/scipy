/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsAnalysis.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_ANALYSIS_H_
#define LP_DATA_HIGHS_ANALYSIS_H_

#include <vector>

#include "HConfig.h"
#include "util/HighsTimer.h"

//#ifdef HiGHSDEV
struct HighsTimerClock {
  HighsTimerClock(HighsTimer& timer) : timer_(timer) {}

  HighsTimer& timer_;
  std::vector<int> clock_;
};
//#endif

#endif /* LP_DATA_HIGHS_ANALYSIS_H_ */
