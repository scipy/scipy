/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mip/HighsObjectiveFunction.h"

#include <numeric>

#include "../extern/pdqsort/pdqsort.h"
#include "lp_data/HighsLp.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsIntegers.h"

HighsObjectiveFunction::HighsObjectiveFunction(const HighsMipSolver& mipsolver)
    : model(mipsolver.model_) {
  objectiveNonzeros.reserve(model->num_col_);

  for (HighsInt i = 0; i < model->num_col_; ++i) {
    if (model->col_cost_[i] != 0.0) {
      objectiveNonzeros.push_back(i);
    }
  }

  colToPartition.resize(model->num_col_, -1);
  cliquePartitionStart.resize(1);

  if (objectiveNonzeros.empty()) {
    numBinary = 0;
    numIntegral = 0;
    objIntScale = 1.0;
    return;
  }

  numIntegral = std::partition(
                    objectiveNonzeros.begin(), objectiveNonzeros.end(),
                    [&](HighsInt i) { return !mipsolver.isColContinuous(i); }) -
                objectiveNonzeros.begin();

  if (numIntegral == 0)
    numBinary = 0;
  else
    numBinary =
        std::partition(objectiveNonzeros.begin(),
                       objectiveNonzeros.begin() + numIntegral,
                       [&](HighsInt i) {
                         return mipsolver.model_->col_lower_[i] == 0.0 &&
                                mipsolver.model_->col_upper_[i] == 1.0;
                       }) -
        objectiveNonzeros.begin();

  objectiveVals.reserve(objectiveNonzeros.size());
  for (HighsInt i : objectiveNonzeros)
    objectiveVals.push_back(model->col_cost_[i]);

  objIntScale = 0.0;
}

void HighsObjectiveFunction::setupCliquePartition(
    const HighsDomain& globaldom, HighsCliqueTable& cliqueTable) {
  // skip if not more than 1 binary column is present
  if (numBinary <= 1) return;

  std::vector<HighsCliqueTable::CliqueVar> clqvars;

  auto binaryEnd = objectiveNonzeros.begin() + numBinary;

  for (auto it = objectiveNonzeros.begin(); it != binaryEnd; ++it) {
    HighsInt col = *it;
    clqvars.emplace_back(col, model->col_cost_[col] < 0.0);
  }

  cliqueTable.cliquePartition(model->col_cost_, clqvars, cliquePartitionStart);
  HighsInt numPartitions = cliquePartitionStart.size() - 1;
  if (numPartitions == numBinary)
    cliquePartitionStart.resize(1);
  else {
    HighsInt p = 0;
    HighsInt k = 0;

    for (HighsInt i = 0; i < numPartitions; ++i) {
      if (cliquePartitionStart[i + 1] - cliquePartitionStart[i] == 1) continue;

      cliquePartitionStart[p] = k;

      for (HighsInt j = cliquePartitionStart[i];
           j < cliquePartitionStart[i + 1]; ++j) {
        colToPartition[clqvars[j].col] = k++;
      }

      ++p;
    }

    cliquePartitionStart[p] = k;
    cliquePartitionStart.resize(p + 1);

    pdqsort(objectiveNonzeros.begin(), objectiveNonzeros.begin() + numBinary,
            [&](HighsInt i, HighsInt j) {
              return std::make_pair((HighsUInt)colToPartition[i],
                                    HighsHashHelpers::hash(i)) <
                     std::make_pair((HighsUInt)colToPartition[j],
                                    HighsHashHelpers::hash(j));
            });

    for (HighsInt i = 0; i < numBinary; ++i)
      objectiveVals[i] = model->col_cost_[objectiveNonzeros[i]];
  }
}

void HighsObjectiveFunction::checkIntegrality(double epsilon) {
  if (numIntegral == (HighsInt)objectiveNonzeros.size()) {
    if (numIntegral) {
      objIntScale =
          HighsIntegers::integralScale(objectiveVals, epsilon, epsilon);
      if (objIntScale * kHighsTiny > epsilon) objIntScale = 0.0;
    } else
      objIntScale = 1.0;
  }
}
