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

#include <set>

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
HighsDebugStatus debugFreeListNumEntries(
    const HighsModelObject& highs_model_object, const std::set<int>& freeList);

HighsDebugStatus debugDualChuzcFail(
    const HighsOptions& options, const int workCount,
    const std::vector<std::pair<int, double>>& workData, const double* workDual,
    const double selectTheta, const double remainTheta);

void debugDualChuzcWorkDataAndGroupReport(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const std::string message,
    const int report_workCount,
    const std::vector<std::pair<int, double>>& report_workData,
    const std::vector<int>& report_workGroup);
HighsDebugStatus debugDualChuzcWorkDataAndGroup(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const int workCount, const int alt_workCount,
    const int breakIndex, const int alt_breakIndex,
    const std::vector<std::pair<int, double>>& workData,
    const std::vector<std::pair<int, double>>& sorted_workData,
    const std::vector<int>& workGroup, const std::vector<int>& alt_workGroup);

HighsDebugStatus debugSimplexBasicSolution(
    const string message, const HighsModelObject& highs_model_object);

HighsDebugStatus debugSimplexInfoBasisConsistent(
    const HighsModelObject& highs_model_object);

HighsDebugStatus debugSimplexHighsSolutionDifferences(
    const HighsModelObject& highs_model_object);

HighsDebugStatus debugAssessSolutionNormDifference(const HighsOptions& options,
                                                   const std::string type,
                                                   const double difference);
#endif  // SIMPLEX_HSIMPLEXDEBUG_H_
