/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexReport.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSIMPLEXREPORT_H_
#define SIMPLEX_HSIMPLEXREPORT_H_

#include "lp_data/HighsModelObject.h"
#include "simplex/SimplexConst.h"

void reportSimplexPhaseIterations(const HighsModelObject& highs_model_object,
                                  const SimplexAlgorithm algorithm,
                                  const bool initialise = false);
#endif  // SIMPLEX_HSIMPLEXREPORT_H_
