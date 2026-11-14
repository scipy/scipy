#include "pdlp/cupdlp/cupdlp_restart.h"

PDHG_restart_choice PDHG_Check_Restart_GPU(CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPtimers *timers = work->timers;

  // Is it first time called?
  if (timers->nIter == iterates->iLastRestartIter) {
    resobj->dPrimalFeasLastRestart = resobj->dPrimalFeas;
    resobj->dDualFeasLastRestart = resobj->dDualFeas;
    resobj->dDualityGapLastRestart = resobj->dDualityGap;

    resobj->dPrimalFeasLastCandidate = resobj->dPrimalFeas;
    resobj->dDualFeasLastCandidate = resobj->dDualFeas;
    resobj->dDualityGapLastCandidate = resobj->dDualityGap;

    return PDHG_NO_RESTART;
  }

  cupdlp_float muCurrent = PDHG_Restart_Score_GPU(
      work->stepsize->dBeta, work->resobj->dPrimalFeas, work->resobj->dDualFeas,
      work->resobj->dDualityGap);
  cupdlp_float muAverage = PDHG_Restart_Score_GPU(
      work->stepsize->dBeta, work->resobj->dPrimalFeasAverage,
      work->resobj->dDualFeasAverage, work->resobj->dDualityGapAverage);
#if PDLP_DEBUG_LOG
  debugPdlpRestartLog(work->debug_pdlp_log_file_, timers->nIter, muCurrent, muAverage);
#endif
  cupdlp_float muCandidate = 0.0;
  PDHG_restart_choice restart_choice = PDHG_RESTART_TO_AVERAGE;
  if (muCurrent < muAverage) {
    restart_choice = PDHG_RESTART_TO_CURRENT;
    muCandidate = muCurrent;
  } else {
    // restart_choice = PDHG_RESTART_TO_AVERAGE;
    muCandidate = muAverage;
  }

  // Or should we do artificial restart based on iteration count?
  // if (2 * (timers->nIter - iterates->iLastRestartIter) >= timers->nIter) {
  if ((timers->nIter - iterates->iLastRestartIter) >= 0.36 * timers->nIter) {
#if CUPDLP_DEBUG
    cupdlp_printf("Doing artificial restart.");
#endif
  } else {
    cupdlp_float muLastRestart = PDHG_Restart_Score_GPU(
        work->stepsize->dBeta, work->resobj->dPrimalFeasLastRestart,
        work->resobj->dDualFeasLastRestart,
        work->resobj->dDualityGapLastRestart);

    // Sufficient decay
    if (muCandidate < 0.2 * muLastRestart) {
#if CUPDLP_DEBUG
      cupdlp_printf("Doing sufficient restart.");
#endif
    } else {
      cupdlp_float muLastCandidate = PDHG_Restart_Score_GPU(
          work->stepsize->dBeta, work->resobj->dPrimalFeasLastCandidate,
          work->resobj->dDualFeasLastCandidate,
          work->resobj->dDualityGapLastCandidate);

      // Necessary decay
      if (muCandidate < 0.8 * muLastRestart && muCandidate > muLastCandidate) {
#if CUPDLP_DEBUG
        cupdlp_printf("Doing necessary restart.");
#endif
      } else {
        restart_choice = PDHG_NO_RESTART;
      }
    }
  }
  // record candidate
  if (muCurrent < muAverage) {
    resobj->dPrimalFeasLastCandidate = resobj->dPrimalFeas;
    resobj->dDualFeasLastCandidate = resobj->dDualFeas;
    resobj->dDualityGapLastCandidate = resobj->dDualityGap;
  } else {
    resobj->dPrimalFeasLastCandidate = resobj->dPrimalFeasAverage;
    resobj->dDualFeasLastCandidate = resobj->dDualFeasAverage;
    resobj->dDualityGapLastCandidate = resobj->dDualityGapAverage;
  }
  //  work->settings->nLogLevel = 2;
  if (restart_choice != PDHG_NO_RESTART) {
    if (muCurrent < muAverage) {
        if (work->settings->nLogLevel > 1)
	  cupdlp_printf("Last restart was iter %d: %s", iterates->iLastRestartIter,
                    "current\n");
    } else {
        if (work->settings->nLogLevel > 1)
	  cupdlp_printf("Last restart was iter %d: %s", iterates->iLastRestartIter,
                    "average\n");
    }
  }
  return restart_choice;
}

cupdlp_bool PDHG_Check_Restart_Merit_Function(CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPresobj *resobj = work->resobj;

  return (
      (fabs(resobj->dDualityGap) > 2.0 * fabs(resobj->dDualityGapAverage)) &&
      (resobj->dPrimalFeas > 2.0 * resobj->dPrimalFeasAverage));
}

cupdlp_float PDHG_Restart_Score_GPU(cupdlp_float weightSquared,
                                    cupdlp_float dPrimalFeas,
                                    cupdlp_float dDualFeas,
                                    cupdlp_float dDualityGap) {
  cupdlp_float dScoreGPU = 0.0;

  dScoreGPU =
      sqrt(weightSquared * dPrimalFeas * dPrimalFeas +
           dDualFeas * dDualFeas / weightSquared + dDualityGap * dDualityGap);

  return dScoreGPU;
}
