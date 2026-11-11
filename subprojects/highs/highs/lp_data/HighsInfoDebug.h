/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInfoDebug.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_INFO_DEBUG_H_
#define LP_DATA_HIGHS_INFO_DEBUG_H_

#include "lp_data/HStruct.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
// #include "lp_data/HighsLp.h"

HighsDebugStatus debugInfo(const HighsOptions& options, const HighsLp& lp,
                           const HighsBasis& basis,
                           const HighsSolution& solution, const HighsInfo& info,
                           const HighsModelStatus model_status);

HighsDebugStatus debugNoInfo(const HighsInfo& info);

#endif
