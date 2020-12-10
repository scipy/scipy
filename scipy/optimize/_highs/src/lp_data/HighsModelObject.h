/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef LP_DATA_HIGHS_MODEL_OBJECT_H_
#define LP_DATA_HIGHS_MODEL_OBJECT_H_

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "simplex/HFactor.h"
#include "simplex/HMatrix.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "util/HighsRandom.h"
#include "util/HighsTimer.h"
// include Sensitivity(or Ranging) header

// Class to communicate data between the simplex solver and the class
// Highs below. Sensitivity data structure would be added here. Only
// include essential data.
class HighsModelObject {
 public:
  HighsModelObject(HighsLp& lp, HighsOptions& options, HighsTimer& timer)
      : lp_(lp), options_(options), timer_(timer), simplex_analysis_(timer) {}

  HighsLp& lp_;
  HighsOptions& options_;
  HighsTimer& timer_;

  HighsModelStatus unscaled_model_status_ = HighsModelStatus::NOTSET;
  HighsModelStatus scaled_model_status_ = HighsModelStatus::NOTSET;

  HighsSolutionParams unscaled_solution_params_;
  HighsSolutionParams scaled_solution_params_;
  HighsIterationCounts iteration_counts_;
  HighsBasis basis_;
  HighsSolution solution_;

  HighsLp simplex_lp_;
  SimplexBasis simplex_basis_;
  HighsSimplexInfo simplex_info_;
  HighsSimplexLpStatus simplex_lp_status_;
  //  HighsRanging ranging_;
  HighsScale scale_;
  HMatrix matrix_;
  HFactor factor_;
  HighsSimplexAnalysis simplex_analysis_;
  HighsRandom random_;

  bool report_model_operations_clock = false;
};

#endif  // LP_DATA_HIGHS_MODEL_OBJECT_H_
