/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/MipTimer.h
 * @brief Indices of mip iClocks
 */
#ifndef MIP_MIPTIMER_H_
#define MIP_MIPTIMER_H_

// Clocks for profiling the MIP solver
enum iClockMip : int {
  kMipClockTotal = 0,
  kMipClockPresolve = kToSubSolver,
  kMipClockSolve,
  kMipClockPostsolve,
  // Level 1
  kFromMipClock,
  kMipClockInit = kFromMipClock,
  kMipClockRunPresolve,
  kMipClockRunSetup,
  kMipClockFeasibilityJump,
  kMipClockTrivialHeuristics,
  kMipClockEvaluateRootNode,
  kMipClockPerformAging0,
  kMipClockSearch,
  // Search
  kMipClockProbingPresolve,
  kMipClockEnumerationPresolve,
  kMipClockPerformAging1,
  kMipClockDive,
  kMipClockOpenNodesToQueue0,
  kMipClockDomainPropagate,
  kMipClockPruneInfeasibleNodes,
  kMipClockUpdateLocalDomain,
  kMipClockNodeSearch,

  // Evaluate root node
  kMipClockStartSymmetryDetection,
  kMipClockStartAnalyticCentreComputation,
  kMipClockEvaluateRootLp,
  kMipClockSeparateLpCuts,
  kMipClockRandomizedRounding,
  kMipClockPerformRestart,
  kMipClockRootSeparation,
  kMipClockFinishAnalyticCentreComputation,
  kMipClockRootCentralRounding,
  kMipClockRootSeparationRound0,
  kMipClockRootHeuristicsReducedCost,
  kMipClockRootSeparationRound1,
  kMipClockRootHeuristicsRens,
  kMipClockRootSeparationRound2,
  kMipClockRootFeasibilityPump,
  kMipClockRootSeparationRound3,
  //  kMipClock@,
  //  kMipClock@,
  //  kMipClock@,
  //  kMipClock@,

  kMipClockEvaluateRootNode0,
  kMipClockEvaluateRootNode1,
  kMipClockEvaluateRootNode2,

  // Dive
  kMipClockDiveEvaluateNode,
  kMipClockDivePrimalHeuristics,
  kMipClockTheDive,
  kMipClockBacktrackPlunge,
  kMipClockPerformAging2,

  // Dive primal heuristics
  kMipClockDiveRandomizedRounding,
  kMipClockDiveRens,
  kMipClockDiveRins,

  // NodeSearch
  kMipClockCurrentNodeToQueue,
  kMipClockSearchBacktrack,
  kMipClockNodePrunedLoop,
  kMipClockOpenNodesToQueue1,
  kMipClockEvaluateNode1,
  kMipClockNodeSearchSeparation,
  kMipClockStoreBasis,
  //  kMipClock@,

  // Separation
  kMipClockRootSeparationRound,
  kMipClockRootSeparationFinishAnalyticCentreComputation,
  kMipClockRootSeparationCentralRounding,
  kMipClockRootSeparationEvaluateRootLp,
  kMipClockImplboundSepa,
  kMipClockCliqueSepa,
  kMipClockTableauSepa,
  kMipClockPathAggrSepa,
  kMipClockModKSepa,
  kMipClockMachineSchedSepa,

  // LP solves
  kMipClockDuSimplexBasisSolveLp,
  kMipClockDuSimplexNoBasisSolveLp,
  kMipClockHipoSolveAnalyticCentreLp,
  kMipClockIpxSolveAnalyticCentreLp,
  kMipClockHipoSolveLp,
  kMipClockIpxSolveLp,
  kMipClockPrSimplexBasisSolveLp,
  kMipClockPrSimplexNoBasisSolveLp,

  // Sub-MIP solves
  kMipClockSubMipSolve,

  kMipClockProbingImplications,

  kLastMipClock = kMipClockProbingImplications,
  kToMipClock = kLastMipClock + 1
};

const HighsInt kNumThreadMipClock = kLastMipClock;

static const double kMipClockTolerancePercentReport = 0.1;

inline void initialiseMipProfilingNames(std::vector<std::string>& name) {
  assert(name.size() == static_cast<size_t>(kToMipClock));
  // Level 1 - Should correspond to kMipClockTotal
  name[kMipClockInit] = "Initialise";
  name[kMipClockRunPresolve] = "Run presolve";
  name[kMipClockRunSetup] = "Run setup";
  name[kMipClockFeasibilityJump] = "Feasibility jump";
  name[kMipClockTrivialHeuristics] = "Trivial heuristics";
  name[kMipClockEvaluateRootNode] = "Evaluate root node";
  name[kMipClockPerformAging0] = "Perform aging 0";
  name[kMipClockSearch] = "Search";
  // kMipClockPostsolve

  // Evaluate root node
  name[kMipClockStartSymmetryDetection] = "Start symmetry detection";
  name[kMipClockStartAnalyticCentreComputation] = "A-centre - start";
  name[kMipClockEvaluateRootLp] = "Evaluate root LP";
  name[kMipClockSeparateLpCuts] = "Separate LP cuts";
  name[kMipClockRandomizedRounding] = "Randomized rounding";
  name[kMipClockPerformRestart] = "Perform restart";
  name[kMipClockRootSeparation] = "Root separation";
  name[kMipClockFinishAnalyticCentreComputation] = "A-centre - finish";
  name[kMipClockRootCentralRounding] = "Root central rounding";
  name[kMipClockRootSeparationRound0] = "Root separation round 0";
  name[kMipClockRootHeuristicsReducedCost] = "Root heuristics reduced cost";
  name[kMipClockRootSeparationRound1] = "Root separation round 1";
  name[kMipClockRootHeuristicsRens] = "Root heuristics RENS";
  name[kMipClockRootSeparationRound2] = "Root separation round 2";
  name[kMipClockRootFeasibilityPump] = "Root feasibility pump";
  name[kMipClockRootSeparationRound3] = "Root separation round 3";
  name[kMipClockEvaluateRootNode0] = "kMipClockEvaluateRootNode0";
  name[kMipClockEvaluateRootNode1] = "kMipClockEvaluateRootNode1";
  name[kMipClockEvaluateRootNode2] = "kMipClockEvaluateRootNode2";

  // Separation
  name[kMipClockRootSeparationRound] = "Separation";
  name[kMipClockRootSeparationFinishAnalyticCentreComputation] =
      "A-centre - finish";
  name[kMipClockRootSeparationCentralRounding] = "Central rounding";
  name[kMipClockRootSeparationEvaluateRootLp] = "Evaluate root LP";

  /*
    clock[kMipClockImplboundSepa] =
    timer_pointer->clock_def(kImplboundSepaString.c_str());
    clock[kMipClockCliqueSepa] =
    timer_pointer->clock_def(kCliqueSepaString.c_str());
    clock[kMipClockTableauSepa] =
    timer_pointer->clock_def(kTableauSepaString.c_str());
    clock[kMipClockPathAggrSepa] =
    timer_pointer->clock_def(kPathAggrSepaString.c_str());
    clock[kMipClockModKSepa] =
    timer_pointer->clock_def(kModKSepaString.c_str());
  */
  // Presolve - Should correspond to kMipClockRunPresolve
  name[kMipClockProbingPresolve] = "Probing - presolve";
  name[kMipClockEnumerationPresolve] = "Enumeration - presolve";

  // Search - Should correspond to kMipClockSearch
  name[kMipClockPerformAging1] = "Perform aging 1";
  name[kMipClockDive] = "Dive";
  name[kMipClockOpenNodesToQueue0] = "Open nodes to queue 0";
  name[kMipClockDomainPropagate] = "Domain propagate";
  name[kMipClockPruneInfeasibleNodes] = "Prune infeasible nodes";
  name[kMipClockUpdateLocalDomain] = "Update local domain";
  name[kMipClockNodeSearch] = "Node search";

  // Dive - Should correspond to kMipClockDive
  name[kMipClockDiveEvaluateNode] = "Evaluate node";
  name[kMipClockDivePrimalHeuristics] = "Dive primal heuristics";
  name[kMipClockTheDive] = "The dive";
  name[kMipClockBacktrackPlunge] = "Backtrack plunge";
  name[kMipClockPerformAging2] = "Perform aging 2";

  // Primal heuristics - Should correspond to kMipDiveClockPrimalHeuristics
  name[kMipClockDiveRandomizedRounding] = "Dive Randomized rounding";
  name[kMipClockDiveRens] = "Dive RENS";
  name[kMipClockDiveRins] = "Dive RINS";

  // Node search
  name[kMipClockCurrentNodeToQueue] = "Current node to queue";
  name[kMipClockSearchBacktrack] = "Search backtrack";
  name[kMipClockNodePrunedLoop] = "Pruned loop search";
  name[kMipClockOpenNodesToQueue1] = "Open nodes to queue 1";
  name[kMipClockEvaluateNode1] = "Evaluate node 1";
  name[kMipClockNodeSearchSeparation] = "Node search separation";
  name[kMipClockStoreBasis] = "Store basis";

  name[kMipClockProbingImplications] = "Probing - implications";
};

class MipTimer {
 public:
  void initialiseMipClocks(HighsTimerClock& mip_timer_clock,
                           const HighsInt thread_mip_clock_offset) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    std::vector<HighsInt>& clock = mip_timer_clock.clock_;

    clock.resize(kLastMipClock);
    clock[kMipClockTotal] = 0;
    clock[kMipClockPresolve] = timer_pointer->clock_def("MIP presolve");
    clock[kMipClockSolve] = timer_pointer->clock_def("MIP solve");
    clock[kMipClockPostsolve] = timer_pointer->clock_def("MIP postsolve");

    // Sometimes the analytic centre clock isn't stopped - because it
    // runs on a separate thread. Although it would be good to
    // understand this better, for now don't assert that this clock
    // has stopped in HighsTimer.h. This is done with a hard-coded
    // clock IDs that need to equal clock[kMipClockHipoSolveAnalyticCentreLp]
    // and clock[kMipClockIpxSolveAnalyticCentreLp]
    //
    // Define the clocks for evaluating the LPs first, so that they
    // aren't changed by inserting new clocks
    clock[kMipClockDuSimplexBasisSolveLp] =
        timer_pointer->clock_def("Solve LP - du simplex basis");
    assert(clock[kMipClockDuSimplexBasisSolveLp] ==
           thread_mip_clock_offset + 7);
    clock[kMipClockDuSimplexNoBasisSolveLp] =
        timer_pointer->clock_def("Solve LP - du simplex no basis");
    clock[kMipClockHipoSolveAnalyticCentreLp] =
        timer_pointer->clock_def("Solve LP: HiPO analytic centre");
    clock[kMipClockIpxSolveAnalyticCentreLp] =
        timer_pointer->clock_def("Solve LP: IPX analytic centre");
    clock[kMipClockHipoSolveLp] = timer_pointer->clock_def("Solve LP: HiPO");
    clock[kMipClockIpxSolveLp] = timer_pointer->clock_def("Solve LP: IPX");
    clock[kMipClockPrSimplexBasisSolveLp] =
        timer_pointer->clock_def("Solve LP - pr simplex basis");
    clock[kMipClockPrSimplexNoBasisSolveLp] =
        timer_pointer->clock_def("Solve LP - pr simplex no basis");

    // Level 1 - Should correspond to kMipClockTotal
    clock[kMipClockInit] = timer_pointer->clock_def("Initialise");
    clock[kMipClockRunPresolve] = timer_pointer->clock_def("Run presolve");
    clock[kMipClockRunSetup] = timer_pointer->clock_def("Run setup");
    clock[kMipClockFeasibilityJump] =
        timer_pointer->clock_def("Feasibility jump");
    clock[kMipClockTrivialHeuristics] =
        timer_pointer->clock_def("Trivial heuristics");
    clock[kMipClockEvaluateRootNode] =
        timer_pointer->clock_def("Evaluate root node");
    clock[kMipClockPerformAging0] = timer_pointer->clock_def("Perform aging 0");
    clock[kMipClockSearch] = timer_pointer->clock_def("Search");
    // kMipClockPostsolve

    // Evaluate root node
    clock[kMipClockStartSymmetryDetection] =
        timer_pointer->clock_def("Start symmetry detection");
    clock[kMipClockStartAnalyticCentreComputation] =
        timer_pointer->clock_def("A-centre - start");
    clock[kMipClockEvaluateRootLp] =
        timer_pointer->clock_def("Evaluate root LP");
    clock[kMipClockSeparateLpCuts] =
        timer_pointer->clock_def("Separate LP cuts");
    clock[kMipClockRandomizedRounding] =
        timer_pointer->clock_def("Randomized rounding");
    clock[kMipClockPerformRestart] =
        timer_pointer->clock_def("Perform restart");
    clock[kMipClockRootSeparation] =
        timer_pointer->clock_def("Root separation");
    clock[kMipClockFinishAnalyticCentreComputation] =
        timer_pointer->clock_def("A-centre - finish");
    clock[kMipClockRootCentralRounding] =
        timer_pointer->clock_def("Root central rounding");
    clock[kMipClockRootSeparationRound0] =
        timer_pointer->clock_def("Root separation round 0");
    clock[kMipClockRootHeuristicsReducedCost] =
        timer_pointer->clock_def("Root heuristics reduced cost");
    clock[kMipClockRootSeparationRound1] =
        timer_pointer->clock_def("Root separation round 1");
    clock[kMipClockRootHeuristicsRens] =
        timer_pointer->clock_def("Root heuristics RENS");
    clock[kMipClockRootSeparationRound2] =
        timer_pointer->clock_def("Root separation round 2");
    clock[kMipClockRootFeasibilityPump] =
        timer_pointer->clock_def("Root feasibility pump");
    clock[kMipClockRootSeparationRound3] =
        timer_pointer->clock_def("Root separation round 3");
    //    clock[kMipClock@] = timer_pointer->clock_def("@");

    clock[kMipClockEvaluateRootNode0] =
        timer_pointer->clock_def("kMipClockEvaluateRootNode0");
    clock[kMipClockEvaluateRootNode1] =
        timer_pointer->clock_def("kMipClockEvaluateRootNode1");
    clock[kMipClockEvaluateRootNode2] =
        timer_pointer->clock_def("kMipClockEvaluateRootNode2");

    // Separation
    clock[kMipClockRootSeparationRound] =
        timer_pointer->clock_def("Separation");
    clock[kMipClockRootSeparationFinishAnalyticCentreComputation] =
        timer_pointer->clock_def("A-centre - finish");
    clock[kMipClockRootSeparationCentralRounding] =
        timer_pointer->clock_def("Central rounding");
    clock[kMipClockRootSeparationEvaluateRootLp] =
        timer_pointer->clock_def("Evaluate root LP");

    /*
    clock[kMipClockImplboundSepa] =
        timer_pointer->clock_def(kImplboundSepaString.c_str());
    clock[kMipClockCliqueSepa] =
        timer_pointer->clock_def(kCliqueSepaString.c_str());
    clock[kMipClockTableauSepa] =
        timer_pointer->clock_def(kTableauSepaString.c_str());
    clock[kMipClockPathAggrSepa] =
        timer_pointer->clock_def(kPathAggrSepaString.c_str());
    clock[kMipClockModKSepa] =
        timer_pointer->clock_def(kModKSepaString.c_str());
    clock[kMipClockMachineSchedSepa] =
        timer_pointer->clock_def(kMachineSchedSepaString.c_str());
    */
    // Presolve - Should correspond to kMipClockRunPresolve
    clock[kMipClockProbingPresolve] =
        timer_pointer->clock_def("Probing - presolve");

    clock[kMipClockEnumerationPresolve] =
        timer_pointer->clock_def("Enumeration - presolve");

    // Search - Should correspond to kMipClockSearch
    clock[kMipClockPerformAging1] = timer_pointer->clock_def("Perform aging 1");
    clock[kMipClockDive] = timer_pointer->clock_def("Dive");
    clock[kMipClockOpenNodesToQueue0] =
        timer_pointer->clock_def("Open nodes to queue 0");
    clock[kMipClockDomainPropagate] =
        timer_pointer->clock_def("Domain propagate");
    clock[kMipClockPruneInfeasibleNodes] =
        timer_pointer->clock_def("Prune infeasible nodes");
    clock[kMipClockUpdateLocalDomain] =
        timer_pointer->clock_def("Update local domain");
    clock[kMipClockNodeSearch] = timer_pointer->clock_def("Node search");
    //    clock[kMipClock@] = timer_pointer->clock_def("@");

    // Dive - Should correspond to kMipClockDive
    clock[kMipClockDiveEvaluateNode] =
        timer_pointer->clock_def("Evaluate node");
    clock[kMipClockDivePrimalHeuristics] =
        timer_pointer->clock_def("Dive primal heuristics");
    clock[kMipClockTheDive] = timer_pointer->clock_def("The dive");
    clock[kMipClockBacktrackPlunge] =
        timer_pointer->clock_def("Backtrack plunge");
    clock[kMipClockPerformAging2] = timer_pointer->clock_def("Perform aging 2");

    // Primal heuristics - Should correspond to kMipDiveClockPrimalHeuristics
    clock[kMipClockDiveRandomizedRounding] =
        timer_pointer->clock_def("Dive Randomized rounding");
    clock[kMipClockDiveRens] = timer_pointer->clock_def("Dive RENS");
    clock[kMipClockDiveRins] = timer_pointer->clock_def("Dive RINS");

    // Node search
    clock[kMipClockCurrentNodeToQueue] =
        timer_pointer->clock_def("Current node to queue");
    clock[kMipClockSearchBacktrack] =
        timer_pointer->clock_def("Search backtrack");
    clock[kMipClockNodePrunedLoop] =
        timer_pointer->clock_def("Pruned loop search");
    clock[kMipClockOpenNodesToQueue1] =
        timer_pointer->clock_def("Open nodes to queue 1");
    clock[kMipClockEvaluateNode1] = timer_pointer->clock_def("Evaluate node 1");
    clock[kMipClockNodeSearchSeparation] =
        timer_pointer->clock_def("Node search separation");
    clock[kMipClockStoreBasis] = timer_pointer->clock_def("Store basis");
    //    clock[] = timer_pointer->clock_def("");

    // Sub-MIP clock
    clock[kMipClockSubMipSolve] = timer_pointer->clock_def("Sub-MIP solves");

    clock[kMipClockProbingImplications] =
        timer_pointer->clock_def("Probing - implications");
    //    clock[] = timer_pointer->clock_def("");
  };

  bool reportMipClockList(const char* grepStamp,
                          const std::vector<HighsInt>& mip_clock_list,
                          const HighsTimerClock& mip_timer_clock,
                          const HighsInt kMipClockIdeal = kMipClockTotal,
                          const double tolerance_percent_report_ = -1) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    if (!timer_pointer->printf_flag) return false;
    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    HighsInt mip_clock_list_size = mip_clock_list.size();
    std::vector<HighsInt> clockList;
    clockList.resize(mip_clock_list_size);
    for (HighsInt en = 0; en < mip_clock_list_size; en++) {
      clockList[en] = clock[mip_clock_list[en]];
    }
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kMipClockIdeal]];
    const double tolerance_percent_report =
        tolerance_percent_report_ >= 0 ? tolerance_percent_report_ : 1e-8;
    return timer_pointer->reportOnTolerance(
        grepStamp, clockList, ideal_sum_time, tolerance_percent_report);
  };

  void csvMipClockList(const std::string& grep_query,
                       const std::string& model_name,
                       const std::vector<HighsInt>& mip_clock_list,
                       const HighsTimerClock& mip_timer_clock,
                       const HighsInt kMipClockIdeal, const bool header,
                       const bool end_line) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    if (!timer_pointer->printf_flag) return;
    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kMipClockIdeal]];
    if (ideal_sum_time < 1e-2) return;
    const HighsInt num_clock = mip_clock_list.size();
    if (header) {
      printf("grep_%s,model,ideal", grep_query.c_str());
      for (HighsInt iX = 0; iX < num_clock; iX++) {
        HighsInt iclock = clock[mip_clock_list[iX]];
        printf(",%s", timer_pointer->clock_names[iclock].c_str());
      }
      printf(",Unaccounted");
      if (end_line) printf("\n");
      return;
    }
    double sum_time = 0;
    printf("grep_%s,%s,%11.4g", grep_query.c_str(), model_name.c_str(),
           ideal_sum_time);
    for (HighsInt iX = 0; iX < num_clock; iX++) {
      HighsInt iclock = clock[mip_clock_list[iX]];
      double time = timer_pointer->read(iclock);
      sum_time += time;
      printf(",%11.4g", time);
    }
    printf(",%11.4g", ideal_sum_time - sum_time);
    if (end_line) printf("\n");
  }

  void reportMipCoreClock(const HighsTimerClock& mip_timer_clock) {
    //    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    const std::vector<HighsInt> mip_clock_list{
        kMipClockPresolve, kMipClockSolve, kMipClockPostsolve};
    reportMipClockList("MipCore_", mip_clock_list, mip_timer_clock,
                       kMipClockTotal);
  };

  void reportMipLevel1Clock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockInit,
                                               kMipClockRunPresolve,
                                               kMipClockRunSetup,
                                               kMipClockFeasibilityJump,
                                               kMipClockTrivialHeuristics,
                                               kMipClockEvaluateRootNode,
                                               kMipClockPerformAging0,
                                               kMipClockSearch,
                                               kMipClockPostsolve};
    reportMipClockList("MipLevl1", mip_clock_list, mip_timer_clock,
                       kMipClockTotal, kMipClockTolerancePercentReport);
  };

  void reportMipSolveLpClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockDuSimplexBasisSolveLp,
        kMipClockDuSimplexNoBasisSolveLp,
        kMipClockPrSimplexBasisSolveLp,
        kMipClockPrSimplexNoBasisSolveLp,
        kMipClockHipoSolveAnalyticCentreLp,
        kMipClockIpxSolveAnalyticCentreLp,
        kMipClockHipoSolveLp,
        kMipClockIpxSolveLp};
    reportMipClockList("MipSlvLp", mip_clock_list, mip_timer_clock,
                       kMipClockTotal);  //, kMipClockTolerancePercentReport);
  };

  void reportMipSubMipSolveClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockSubMipSolve};
    reportMipClockList("MipSlvLp", mip_clock_list, mip_timer_clock,
                       kMipClockTotal);  //, kMipClockTolerancePercentReport);
  };

  void reportMipPresolveClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockProbingPresolve,
                                               kMipClockEnumerationPresolve};
    reportMipClockList("MipPrslv", mip_clock_list, mip_timer_clock,
                       kMipClockRunPresolve, kMipClockTolerancePercentReport);
  };

  void reportAltEvaluateRootNodeClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{kMipClockEvaluateRootNode0,
                                               kMipClockEvaluateRootNode1,
                                               kMipClockEvaluateRootNode2};
    reportMipClockList(
        "AltEvaluateRootNode", mip_clock_list, mip_timer_clock,
        kMipClockEvaluateRootNode);  //, kMipClockTolerancePercentReport);
  };

  void reportMipEvaluateRootNodeClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockStartSymmetryDetection,
        kMipClockStartAnalyticCentreComputation,
        kMipClockEvaluateRootLp,
        kMipClockSeparateLpCuts,
        kMipClockRandomizedRounding,
        kMipClockPerformRestart,
        kMipClockRootSeparation,
        kMipClockFinishAnalyticCentreComputation,
        kMipClockRootCentralRounding,
        kMipClockRootSeparationRound0,
        kMipClockRootHeuristicsReducedCost,
        kMipClockRootSeparationRound1,
        kMipClockRootHeuristicsRens,
        kMipClockRootSeparationRound2,
        kMipClockRootFeasibilityPump,
        kMipClockRootSeparationRound3
        //	kMipClock@,
        //	kMipClock@
    };
    reportMipClockList(
        "MipEvaluateRootNode", mip_clock_list, mip_timer_clock,
        kMipClockEvaluateRootNode);  //, kMipClockTolerancePercentReport);
  };

  void reportMipRootSeparationClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockRootSeparationRound,
        kMipClockRootSeparationFinishAnalyticCentreComputation,
        kMipClockRootSeparationCentralRounding,
        kMipClockRootSeparationEvaluateRootLp};
    reportMipClockList(
        "MipRootSeparation", mip_clock_list, mip_timer_clock,
        kMipClockRootSeparation);  //, kMipClockTolerancePercentReport);
  };

  void reportMipSearchClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockPerformAging1,        kMipClockDive,
        kMipClockOpenNodesToQueue0,    kMipClockDomainPropagate,
        kMipClockPruneInfeasibleNodes, kMipClockUpdateLocalDomain,
        kMipClockNodeSearch,
        //	kMipClock@
    };
    reportMipClockList("MipSerch", mip_clock_list, mip_timer_clock,
                       kMipClockSearch, kMipClockTolerancePercentReport);
  };

  void reportMipDiveClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockDiveEvaluateNode, kMipClockDivePrimalHeuristics,
        kMipClockTheDive, kMipClockBacktrackPlunge, kMipClockPerformAging2};
    reportMipClockList("MipDive_", mip_clock_list, mip_timer_clock,
                       kMipClockDive, kMipClockTolerancePercentReport);
  };

  void reportMipDivePrimalHeuristicsClock(
      const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockDiveRandomizedRounding, kMipClockDiveRens, kMipClockDiveRins};
    reportMipClockList("MipDivePrimalHeuristics", mip_clock_list,
                       mip_timer_clock, kMipClockDivePrimalHeuristics,
                       kMipClockTolerancePercentReport);
  };

  void reportMipNodeSearchClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockCurrentNodeToQueue, kMipClockNodePrunedLoop,
        //      kMipClockSearchBacktrack,
        kMipClockOpenNodesToQueue1, kMipClockEvaluateNode1,
        kMipClockNodeSearchSeparation};  //, kMipClockStoreBasis};
    reportMipClockList(
        "MipNodeSearch", mip_clock_list, mip_timer_clock,
        kMipClockNodeSearch);  //, kMipClockTolerancePercentReport);
  };

  void reportMipSeparationClock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockImplboundSepa, kMipClockCliqueSepa, kMipClockTableauSepa,
        kMipClockPathAggrSepa,  kMipClockModKSepa,   kMipClockMachineSchedSepa};
    reportMipClockList("MipSeparation", mip_clock_list, mip_timer_clock,
                       kMipClockTotal);  //, kMipClockTolerancePercentReport);
  };

  void csvMipClock(const std::string& model_name,
                   const HighsTimerClock& mip_timer_clock, const bool header,
                   const bool end_line) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockRunPresolve, kMipClockEvaluateRootNode,
        kMipClockDivePrimalHeuristics, kMipClockTheDive, kMipClockNodeSearch};
    csvMipClockList("csvMIP", model_name, mip_clock_list, mip_timer_clock,
                    kMipClockTotal, header, end_line);
  };

  void csvEvaluateRootNodeClock(const std::string& model_name,
                                const HighsTimerClock& mip_timer_clock,
                                const bool header, const bool end_line) {
    const std::vector<HighsInt> mip_clock_list{
        kMipClockStartSymmetryDetection,
        kMipClockStartAnalyticCentreComputation,
        kMipClockEvaluateRootLp,
        kMipClockSeparateLpCuts,
        kMipClockRandomizedRounding,
        kMipClockPerformRestart,
        kMipClockRootSeparation,
        kMipClockFinishAnalyticCentreComputation,
        kMipClockRootCentralRounding,
        kMipClockRootSeparationRound0,
        kMipClockRootHeuristicsReducedCost,
        kMipClockRootSeparationRound1,
        kMipClockRootHeuristicsRens,
        kMipClockRootSeparationRound2,
        kMipClockRootFeasibilityPump,
        kMipClockRootSeparationRound3};
    csvMipClockList("csvRootNode", model_name, mip_clock_list, mip_timer_clock,
                    kMipClockEvaluateRootNode, header, end_line);
  };

  void reportFjClock(std::string& model,
                     const HighsTimerClock& mip_timer_clock) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    if (!timer_pointer->printf_flag) return;
    HighsInt iClock = mip_timer_clock.clock_[kMipClockFeasibilityJump];
    const double fj_time = timer_pointer->read(iClock);
    const double total_time = timer_pointer->read();
    const double pct = total_time > 0 ? 1e2 * fj_time / total_time : 0;
    printf("grepFK,%s,%s,%d,%g,%g\n", model.c_str(),
           timer_pointer->clock_names[iClock].c_str(),
           int(timer_pointer->clock_num_call[iClock]), fj_time, pct);
  }
};

#endif /* MIP_MIPTIMER_H_ */
