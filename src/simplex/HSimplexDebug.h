/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSIMPLEXDEBUG_H_
#define SIMPLEX_HSIMPLEXDEBUG_H_

#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsOptions.h"
#include "simplex/SimplexConst.h"

HighsDebugStatus debugComputePrimal(const HighsModelObject& highs_model_object,
                                    const std::vector<double>& primal_rhs);

HighsDebugStatus debugComputeDual(const HighsModelObject& highs_model_object,
                                  const std::vector<double>& previous_dual,
                                  const std::vector<double>& basic_costs,
                                  const std::vector<double>& row_dual);

HighsDebugStatus debugUpdatedObjectiveValue(
    HighsModelObject& highs_model_object, const SimplexAlgorithm algorithm,
    const int phase, const std::string message);

HighsDebugStatus debugUpdatedObjectiveValue(
    const HighsModelObject& highs_model_object,
    const SimplexAlgorithm algorithm);

HighsDebugStatus debugFixedNonbasicMove(
    const HighsModelObject& highs_model_object);
HighsDebugStatus debugNonbasicMove(const HighsModelObject& highs_model_object);
HighsDebugStatus debugBasisCondition(const HighsModelObject& highs_model_object,
                                     const std::string message);
HighsDebugStatus debugCleanup(HighsModelObject& highs_model_object,
                              const std::vector<double>& original_dual);
#endif  // SIMPLEX_HSIMPLEXDEBUG_H_
