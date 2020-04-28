/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapperEmpty.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IPM_IPX_WRAPPER_EMPTY_H_
#define IPM_IPX_WRAPPER_EMPTY_H_

#include "ipm/IpxStatus.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"

HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, bool& imprecise_solution,
                       HighsBasis& highs_basis, HighsSolution& highs_solution,
                       HighsIterationCounts& iteration_counts,
                       HighsModelStatus& unscaled_model_status,
                       HighsSolutionParams& unscaled_solution_params) {
  unscaled_model_status = HighsModelStatus::NOTSET;
  return HighsStatus::Error;
}

#endif
