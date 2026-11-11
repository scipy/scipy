/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsCutGeneration.h
 * @brief Class that generates cuts from single row relaxations
 *
 *
 */

#ifndef MIP_HIGHS_CUT_GENERATION_H_
#define MIP_HIGHS_CUT_GENERATION_H_

#include <cstdint>
#include <vector>

#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"
#include "util/HighsType.h"

class HighsLpRelaxation;
class HighsTransformedLp;
class HighsCutPool;
class HighsDomain;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsCutGeneration {
 private:
  const HighsLpRelaxation& lpRelaxation;
  HighsCutPool& cutpool;
  HighsRandom randgen;
  std::vector<HighsInt> cover;
  HighsCDouble coverweight;
  HighsCDouble lambda;
  std::vector<double> upper;
  std::vector<double> solval;
  std::vector<HighsBool> complementation;
  std::vector<HighsBool> isintegral;
  const double feastol;
  const double epsilon;

  double* vals;
  HighsInt* inds;
  HighsCDouble rhs;
  bool integralSupport;
  bool integralCoefficients;
  HighsInt rowlen;
  double initialScale;

  std::vector<HighsInt> integerinds;
  std::vector<double> deltas;

  std::vector<double> tmpVals;
  std::vector<HighsInt> tmpInds;
  std::vector<HighsBool> tmpComplementation;
  std::vector<double> tmpSolval;

  bool determineCover(bool lpSol = true);

  void separateLiftedKnapsackCover();

  bool separateLiftedMixedBinaryCover();

  bool separateLiftedMixedIntegerCover();

  bool cmirCutGenerationHeuristic(double minEfficacy,
                                  bool onlyInitialCMIRScale = false);

  double scale(double val);

  bool postprocessCut(const HighsDomain& globaldom);

  bool preprocessBaseInequality(bool& hasUnboundedInts, bool& hasGeneralInts,
                                bool& hasContinuous);

  void flipComplementation(HighsInt index);

  void removeComplementation();

  void updateViolationAndNorm(HighsInt index, double aj, double& violation,
                              double& norm) const;

  bool tryGenerateCut(std::vector<HighsInt>& inds, std::vector<double>& vals,
                      bool hasUnboundedInts, bool hasGeneralInts,
                      bool hasContinuous, double minEfficacy,
                      bool onlyInitialCMIRScale = false,
                      bool allowRejectCut = true, bool lpSol = true);

 public:
  HighsCutGeneration(const HighsLpRelaxation& lpRelaxation,
                     HighsCutPool& cutpool);

  /// separates the LP solution for the given single row relaxation
  bool generateCut(HighsTransformedLp& transLp, std::vector<HighsInt>& inds,
                   std::vector<double>& vals, double& rhs,
                   bool onlyInitialCMIRScale = false);

  /// generate a conflict from the given proof constraint which cuts of the
  /// given local domain
  bool generateConflict(const HighsDomain& localdom,
                        const HighsDomain& globaldom,
                        std::vector<HighsInt>& proofinds,
                        std::vector<double>& proofvals, double& proofrhs);

  /// applies postprocessing to an externally generated cut and adds it to the
  /// cutpool if it is violated enough
  bool finalizeAndAddCut(const HighsDomain& globaldom,
                         std::vector<HighsInt>& inds, std::vector<double>& vals,
                         double& rhs);
};

#endif
