/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PRESOLVE_ICRASHX_H_
#define PRESOLVE_ICRASHX_H_

#include <iostream>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"

HighsStatus callCrossover(const HighsOptions& options, const HighsLp& lp,
                          HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& highs_callback);

#endif
