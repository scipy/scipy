/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_OBJECTIVE_FUNCTION_H_
#define HIGHS_OBJECTIVE_FUNCTION_H_

#include <algorithm>
#include <cassert>
#include <vector>

#include "util/HighsType.h"

class HighsCliqueTable;
class HighsDomain;
class HighsLp;
class HighsMipSolver;

class HighsObjectiveFunction {
  const HighsLp* model;
  double objIntScale;
  HighsInt numIntegral;
  HighsInt numBinary;
  std::vector<HighsInt> objectiveNonzeros;
  std::vector<double> objectiveVals;
  std::vector<HighsInt> cliquePartitionStart;
  std::vector<HighsInt> colToPartition;

 public:
  HighsObjectiveFunction(const HighsMipSolver& mipsolver);

  void setupCliquePartition(const HighsDomain& globaldom,
                            HighsCliqueTable& cliqueTable);

  void checkIntegrality(double epsilon);

  /// returns the vector of column indices with nonzero objective value
  /// They will be ordered so that binary columns come first
  const std::vector<HighsInt>& getObjectiveNonzeros() const {
    return objectiveNonzeros;
  }

  const std::vector<double>& getObjectiveValuesPacked() const {
    return objectiveVals;
  }

  HighsInt getNumBinariesInObjective() const { return numBinary; }

  const std::vector<HighsInt>& getCliquePartitionStarts() const {
    return cliquePartitionStart;
  }

  HighsInt getNumCliquePartitions() const {
    return cliquePartitionStart.size() - 1;
  }

  HighsInt getColCliquePartition(HighsInt col) const {
    return colToPartition[col];
  }

  double integralScale() const { return objIntScale; }

  bool isIntegral() const { return objIntScale != 0.0; }

  bool isZero() const { return objectiveNonzeros.empty(); }
};

#endif
