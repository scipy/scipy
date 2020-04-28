/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolve.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSSOLVE_H_
#define LP_DATA_HIGHSSOLVE_H_

#include "lp_data/HighsModelUtils.h"
HighsStatus solveLp(HighsModelObject& highs_model_object, const string message);
HighsStatus solveUnconstrainedLp(HighsModelObject& highs_model_object);
#endif  // LP_DATA_HIGHSSOLVE_H_
