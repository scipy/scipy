/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsAnalysis.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_ANALYSIS_H_
#define LP_DATA_HIGHS_ANALYSIS_H_

#include <vector>

#include "HConfig.h"
#include "util/HighsTimer.h"

struct HighsTimerClock {
  HighsTimer* timer_pointer_;
  std::vector<HighsInt> clock_;
};
#endif /* LP_DATA_HIGHS_ANALYSIS_H_ */
