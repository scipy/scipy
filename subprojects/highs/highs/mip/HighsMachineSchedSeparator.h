/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsMachineSchedSeparator.h
 * @brief Class for separating release-date cuts.
 *
 * Given a set of jobs N with start times s_j, processing times p_ij,
 * release times r_j, and binary dependency variables y_ij,
 * where y = 1 -> s_j >= s_i + p_ij.
 * To simplify the problem: We set p_j = min{p_ji : i \in N / j}
 * A valid inequality is then:
 * s_j >= min{r_{i} : i \in N} + \sum_{i \in N / j} p_i y_ij
 *
 * General idea: Each pair of jobs has some non-zero wait time between them
 * that depends on which is ordered first, e.g., if job 1 comes after job 2,
 * then job_1 >= job_2 + 3, however if job 2 comes after job 1,
 * then job_2 >= job_1 + 4. We simplify this by assuming the processing time
 * of a job is the minimum across all pairs where the job comes first, e.g.,
 * job_1 >= job_2 + 3 && job_1 >= job_3 + 5, then the "processing time" of
 * job 1 is 3. If all jobs have greater than zero pair-wise distance, then
 * no jobs can be processed simultaneously. Our derived inequality
 * exploits this.
 *
 * For more info on the cuts, see Section 4.1 in
 * Recognition and Exploitation of Single-Machine Scheduling Subproblems in
 * Mixed Integer Programs, Reinout Lambertus Henricus Wijfjes, Msc Thesis, 2022,
 * or Cutting Plane Algorithm for the Single Machine Scheduling Problem
 * with Release Times, G. L. Nemhauser and M. W. P. Savelsbergh, 1992.
 */

#ifndef MIP_HIGHS_MACHINE_SCHED_SEPARATOR_H_
#define MIP_HIGHS_MACHINE_SCHED_SEPARATOR_H_

#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsSeparator.h"
#include "util/HighsRandom.h"

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsMachineSchedSeparator : public HighsSeparator {
 private:
  HighsRandom randgen;
  bool has_single_machine_schedule = false;
  bool already_tried = false;

  bool findSingleMachineScheduleClique(std::vector<std::vector<double>>& vals,
                                       std::vector<std::vector<HighsInt>>& inds,
                                       std::vector<double>& rhss,
                                       const HighsDomain& globaldom,
                                       const HighsMipSolver& mipsolver);

 public:
  void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                          HighsLpAggregator& lpAggregator,
                          HighsTransformedLp& transLp,
                          HighsCutPool& cutpool) override;

  HighsMachineSchedSeparator(const HighsMipSolver& mipsolver)
      : HighsSeparator(mipsolver, kMachineSchedSepaString) {
    randgen.initialise(mipsolver.options_mip_->random_seed);
    has_single_machine_schedule = false;
    already_tried = false;
  }
};

#endif
