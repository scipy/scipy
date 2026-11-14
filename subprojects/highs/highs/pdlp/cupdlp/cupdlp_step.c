//
// Created by chuwen on 23-11-28.
//

#include "pdlp/cupdlp/cupdlp_step.h"

#include "pdlp/cupdlp/cupdlp_defs.h"
#include "pdlp/cupdlp/cupdlp_linalg.h"
#include "pdlp/cupdlp/cupdlp_proj.h"
// #include "cupdlp_scaling.h"
#include "pdlp/cupdlp/cupdlp_solver.h"
#include "pdlp/cupdlp/cupdlp_utils.h"
#include "pdlp/cupdlp/glbopts.h"

// x^{k+1} = proj_{X}(x^k - dPrimalStep * (c - A'y^k))
void PDHG_primalGradientStep(CUPDLPwork *work, CUPDLPvec *xUpdate,
                             const CUPDLPvec *x, const CUPDLPvec *ATy,
                             cupdlp_float dPrimalStepSize) {
  CUPDLPproblem *problem = work->problem;

#if !defined(CUPDLP_CPU) && USE_KERNELS
  cupdlp_pgrad_cuda(xUpdate->data, x->data, problem->cost,
                    ATy->data, problem->lower, problem->upper, dPrimalStepSize,
                    (int)problem->nCols);
#else
  // cupdlp_copy(xUpdate, x, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(xUpdate->data, x->data, cupdlp_float, problem->nCols);

  // AddToVector(xUpdate, -dPrimalStepSize, problem->cost, problem->nCols);
  // AddToVector(xUpdate, dPrimalStepSize, ATy, problem->nCols);

  cupdlp_float alpha = -dPrimalStepSize;
  cupdlp_axpy(work, problem->nCols, &alpha, problem->cost, xUpdate->data);
  alpha = dPrimalStepSize;
  cupdlp_axpy(work, problem->nCols, &alpha, ATy->data, xUpdate->data);

  cupdlp_projub(xUpdate->data, problem->upper, problem->nCols);
  cupdlp_projlb(xUpdate->data, problem->lower, problem->nCols);
#endif
}

// y^{k+1} = proj_{Y}(y^k + dDualStep * (b - A * (2x^{k+1} - x^{k})))
void PDHG_dualGradientStep(CUPDLPwork *work, CUPDLPvec *yUpdate,
                           const CUPDLPvec *y, const CUPDLPvec *Ax,
                           const CUPDLPvec *AxUpdate, cupdlp_float dDualStepSize) {
  CUPDLPproblem *problem = work->problem;

#if !defined(CUPDLP_CPU) && USE_KERNELS
  cupdlp_dgrad_cuda(yUpdate->data, y->data, problem->rhs, Ax->data,
                    AxUpdate->data, dDualStepSize, (int)problem->nRows,
                    (int)problem->nEqs);
#else
  // cupdlp_copy(yUpdate, y, cupdlp_float, problem->nRows);
  CUPDLP_COPY_VEC(yUpdate->data, y->data, cupdlp_float, problem->nRows);

  // AddToVector(yUpdate, dDualStepSize, problem->rhs, problem->nRows);
  // AddToVector(yUpdate, -2.0 * dDualStepSize, AxUpdate, problem->nRows);
  // AddToVector(yUpdate, dDualStepSize, Ax, problem->nRows);

  cupdlp_float alpha = dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, problem->rhs, yUpdate->data);
  alpha = -2.0 * dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, AxUpdate->data, yUpdate->data);
  alpha = dDualStepSize;
  cupdlp_axpy(work, problem->nRows, &alpha, Ax->data, yUpdate->data);

  cupdlp_projPos(yUpdate->data + problem->nEqs, problem->nRows - problem->nEqs);
#endif
}

cupdlp_retcode PDHG_Power_Method(CUPDLPwork *work, cupdlp_float *lambda) {
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPiterates *iterates = work->iterates;

  if (work->settings->nLogLevel>0) 
    cupdlp_printf("Power Method:\n");

  // work->buffer is a CUPDLPvec*, but a pointer to its values is
  // needed for some linalg calls
  cupdlp_float *q = work->buffer->data;

  cupdlp_initvec(q, 1.0, lp->nRows);

  cupdlp_int iter = work->timers->nIter;
  CUPDLPvec *ax = iterates->ax[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];

  double res = 0.0;
  double previous_lambda = 0.0;
  // q is stored in work->buffer->data
  //
  // z is stored in ax->data
  //
  // y is stored in aty->data
  //
  int log_iters = work->settings->nLogLevel > 1;

  if (log_iters)
    cupdlp_printf("It       lambda   dl_lambda    residual\n");
  for (cupdlp_int iter = 0; iter < 20; ++iter) {
    // z = A*A'*q
    //
    // as y = A'q...
    ATy(work, aty, work->buffer);
    // ... then z = Ay
    Ax(work, ax, aty);

    // q = z / norm(z)
    CUPDLP_COPY_VEC(q, ax->data, cupdlp_float, lp->nRows);
    cupdlp_float qNorm = 0.0;
    cupdlp_twoNorm(work, lp->nRows, q, &qNorm);
    cupdlp_scaleVector(work, 1.0 / qNorm, q, lp->nRows);

    // Now compute the Rayleigh quotient of q which, since q'q=1 is
    //
    // lambda = q'AA'q = w^Tw = ||w||_2
    //
    // where w = A'q
    //
    // aty is no longer needed, so w is stored in aty->data
    ATy(work, aty, work->buffer);

    cupdlp_twoNormSquared(work, lp->nCols, aty->data, lambda);

    cupdlp_float alpha = -(*lambda);
    // Now compute the residual between z = A*A'*q (old q) and
    // lambda.q (new q) to just to log progress
    //
    // z := z - lambda.q
    cupdlp_axpy(work, lp->nRows, &alpha, q, ax->data);

    cupdlp_twoNormSquared(work, lp->nCols, ax->data, &res);

    double dl_lambda = fabs(*lambda - previous_lambda);
    previous_lambda = *lambda;

    if (log_iters)
      cupdlp_printf("%2d %12.6g %11.4g %11.4g\n", iter, *lambda, dl_lambda, res);
  }

exit_cleanup:
  return retcode;
}

void PDHG_Compute_Step_Size_Ratio(CUPDLPwork *pdhg) {
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  cupdlp_float dMeanStepSize =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  cupdlp_int iter = pdhg->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];

  // cupdlp_float dDiffPrimal = cupdlp_diffTwoNorm(x, iterates->xLastRestart, problem->nCols);
  // cupdlp_float dDiffDual = cupdlp_diffTwoNorm(y, iterates->yLastRestart, problem->nRows);

  cupdlp_float dDiffPrimal = 0.0;
  cupdlp_diffTwoNorm(pdhg, x->data, iterates->xLastRestart, problem->nCols, &dDiffPrimal);
  cupdlp_float dDiffDual = 0.0;
  cupdlp_diffTwoNorm(pdhg, y->data, iterates->yLastRestart, problem->nRows, &dDiffDual);

  if (fmin(dDiffPrimal, dDiffDual) > 1e-10) {
    cupdlp_float dBetaUpdate = dDiffDual / dDiffPrimal;
    cupdlp_float dLogBetaUpdate =
        0.5 * log(dBetaUpdate) + 0.5 * log(sqrt(stepsize->dBeta));
    stepsize->dBeta = exp(dLogBetaUpdate) * exp(dLogBetaUpdate);
  }

  stepsize->dPrimalStep = dMeanStepSize / sqrt(stepsize->dBeta);
  stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
  stepsize->dTheta = 1.0;
}

void PDHG_Update_Iterate_Constant_Step_Size(CUPDLPwork *pdhg) {
  //            CUPDLP_ASSERT(0);
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;

  cupdlp_int iter = pdhg->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];
  CUPDLPvec *ax = iterates->ax[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];
  CUPDLPvec *xUpdate = iterates->x[(iter + 1) % 2];
  CUPDLPvec *yUpdate = iterates->y[(iter + 1) % 2];
  CUPDLPvec *axUpdate = iterates->ax[(iter + 1) % 2];
  CUPDLPvec *atyUpdate = iterates->aty[(iter + 1) % 2];

  Ax(pdhg, ax, x);
  ATy(pdhg, aty, y);

  // x^{k+1} = proj_{X}(x^k - dPrimalStep * (c - A'y^k))
  PDHG_primalGradientStep(pdhg, xUpdate, x, aty, stepsize->dPrimalStep);

  Ax(pdhg, axUpdate, xUpdate);

  // y^{k+1} = proj_{Y}(y^k + dDualStep * (b - A * (2 x^{k+1} - x^{k})))
  PDHG_dualGradientStep(pdhg, yUpdate, y, ax, axUpdate, stepsize->dDualStep);

  ATy(pdhg, atyUpdate, yUpdate);
}

void PDHG_Update_Iterate_Malitsky_Pock(CUPDLPwork *pdhg) {
  cupdlp_printf("Malitsky-Pock is not implemented\n");
  cupdlp_printf(" - use %d and %d instead", PDHG_FIXED_LINESEARCH,
                PDHG_ADAPTIVE_LINESEARCH);
  exit(-1);
}

cupdlp_retcode PDHG_Update_Iterate_Adaptive_Step_Size(CUPDLPwork *pdhg) {
  cupdlp_retcode retcode = RETCODE_OK;
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;

  cupdlp_int iter = pdhg->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];
  CUPDLPvec *ax = iterates->ax[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];
  CUPDLPvec *xUpdate = iterates->x[(iter + 1) % 2];
  CUPDLPvec *yUpdate = iterates->y[(iter + 1) % 2];
  CUPDLPvec *axUpdate = iterates->ax[(iter + 1) % 2];
  CUPDLPvec *atyUpdate = iterates->aty[(iter + 1) % 2];

  cupdlp_float dStepSizeUpdate =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);

  cupdlp_bool isDone = false;
  // number of steps this round
  int stepIterThis = 0;
  while (!isDone) {
    ++stepsize->nStepSizeIter;
    ++stepIterThis;

    cupdlp_float dPrimalStepUpdate = dStepSizeUpdate / sqrt(stepsize->dBeta);
    cupdlp_float dDualStepUpdate = dStepSizeUpdate * sqrt(stepsize->dBeta);

    // x^{k+1} = proj_{X}(x^k - dPrimalStep * (c - A'y^k))
    PDHG_primalGradientStep(pdhg, xUpdate, x, aty, dPrimalStepUpdate);

    Ax(pdhg, axUpdate, xUpdate);

    // y^{k+1} = proj_{Y}(y^k + dDualStep * (b - A * (2 x^{k+1} - x^{k})))
    PDHG_dualGradientStep(pdhg, yUpdate, y, ax, axUpdate, dDualStepUpdate);

    ATy(pdhg, atyUpdate, yUpdate);

    cupdlp_float dMovement = 0.0;
    cupdlp_float dInteraction = 0.0;

    cupdlp_compute_interaction_and_movement(pdhg, &dMovement, &dInteraction);

#if CUPDLP_DUMP_LINESEARCH_STATS && CUPDLP_DEBUG
    cupdlp_float dInteractiony = 0.0;
    //      Δy' (AΔx)
    cupdlp_diffDotDiff(pdhg, y->data, yUpdate->data, ax->data, axUpdate->data,
                       problem->nRows, &dInteractiony);
#endif

    cupdlp_float dStepSizeLimit;
    if (dInteraction != 0.0) {
      dStepSizeLimit = dMovement / fabs(dInteraction);
    } else {
      dStepSizeLimit = INFINITY;
    }
    if (dStepSizeUpdate <= dStepSizeLimit) {
      isDone = true;
      // break;
    } else {
      CUPDLP_CHECK_TIMEOUT(pdhg);
    }

    cupdlp_float dFirstTerm = (1.0 - pow(stepsize->nStepSizeIter + 1.0,
                                         -PDHG_STEPSIZE_REDUCTION_EXP)) *
                              dStepSizeLimit;
    cupdlp_float dSecondTerm =
        (1.0 + pow(stepsize->nStepSizeIter + 1.0, -PDHG_STEPSIZE_GROWTH_EXP)) *
        dStepSizeUpdate;
    dStepSizeUpdate = fmin(dFirstTerm, dSecondTerm);
#if CUPDLP_DUMP_LINESEARCH_STATS && CUPDLP_DEBUG
    cupdlp_printf(" -- stepsize iteration %d: %f %f\n", stepIterThis,
                  dStepSizeUpdate, dStepSizeLimit);

    cupdlp_printf(" -- PrimalStep DualStep: %f %f\n", stepsize->dPrimalStep,
                  stepsize->dDualStep);
    cupdlp_printf(" -- FirstTerm SecondTerm: %f %f\n", dFirstTerm, dSecondTerm);
    cupdlp_printf(" -- nStepSizeIter: %d\n", stepsize->nStepSizeIter);
    cupdlp_printf(" -- RED_EXP GRO_EXP: %f %f\n", PDHG_STEPSIZE_REDUCTION_EXP,
                  PDHG_STEPSIZE_GROWTH_EXP);

    cupdlp_printf("     -- iteraction(x) interaction(y): %f %f\n", dInteraction,
                  dInteractiony);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    cupdlp_printf("     -- movement (scaled norm)  : %f\n", dMovement);
    if (stepIterThis > 200) break;  // avoid unlimited runs due to bugs.
#endif
  }

  stepsize->dPrimalStep = dStepSizeUpdate / sqrt(stepsize->dBeta);
  stepsize->dDualStep = dStepSizeUpdate * sqrt(stepsize->dBeta);

exit_cleanup:
  return retcode;
}

cupdlp_retcode PDHG_Init_Step_Sizes(CUPDLPwork *pdhg) {
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPstepsize *stepsize = pdhg->stepsize;

  if (stepsize->eLineSearchMethod == PDHG_FIXED_LINESEARCH) {
    CUPDLP_CALL(PDHG_Power_Method(pdhg, &stepsize->dPrimalStep));
    cupdlp_float power_method_lambda = stepsize->dPrimalStep;
    // PDLP Intial primal weight = norm(cost) / norm(rhs) = sqrt(beta)
    // cupdlp_float a = twoNormSquared(problem->cost, problem->nCols);
    // cupdlp_float b = twoNormSquared(problem->rhs, problem->nRows);
    cupdlp_float a = 0.0;
    cupdlp_float b = 0.0;
    cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
    cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);

    if (fmin(a, b) > 1e-6) {
      stepsize->dBeta = a / b;
    } else {
      stepsize->dBeta = 1.0;
    }

    stepsize->dPrimalStep = 0.8 / sqrt(stepsize->dPrimalStep);
    stepsize->dDualStep = stepsize->dPrimalStep;
    stepsize->dPrimalStep /= sqrt(stepsize->dBeta);
    stepsize->dDualStep *= sqrt(stepsize->dBeta);
    if (pdhg->settings->nLogLevel > 1)
      cupdlp_printf("Initial step sizes from power method lambda = %g: primal = %g; dual = %g\n",
		    power_method_lambda, stepsize->dPrimalStep, stepsize->dDualStep);	  
  } else {
    stepsize->dTheta = 1.0;

    // PDLP Intial primal weight = norm(cost) / norm(rhs) = sqrt(beta)
    // cupdlp_float a = twoNormSquared(problem->cost, problem->nCols);
    // cupdlp_float b = twoNormSquared(problem->rhs, problem->nRows);
    cupdlp_float a = 0.0;
    cupdlp_float b = 0.0;
    cupdlp_twoNormSquared(pdhg, problem->nCols, problem->cost, &a);
    cupdlp_twoNormSquared(pdhg, problem->nRows, problem->rhs, &b);

    if (fmin(a, b) > 1e-6) {
      stepsize->dBeta = a / b;
    } else {
      stepsize->dBeta = 1.0;
    }
    // infNorm can be avoid by previously calculated infNorm of csc matrix
    stepsize->dPrimalStep =
        // (1.0 / infNorm(problem->data->csc_matrix->colMatElem,
        // problem->data->csc_matrix->nMatElem)) /
        (1.0 / problem->data->csc_matrix->MatElemNormInf) /
        sqrt(stepsize->dBeta);
    stepsize->dDualStep = stepsize->dPrimalStep * stepsize->dBeta;
    iterates->dLastRestartBeta = stepsize->dBeta;
  }

  iterates->iLastRestartIter = 0;
  stepsize->dSumPrimalStep = 0;
  stepsize->dSumDualStep = 0;

exit_cleanup:
  return retcode;
}

void PDHG_Compute_Average_Iterate(CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  cupdlp_float dPrimalScale =
      stepsize->dSumPrimalStep > 0.0 ? 1.0 / stepsize->dSumPrimalStep : 1.0;
  cupdlp_float dDualScale =
      stepsize->dSumDualStep > 0.0 ? 1.0 / stepsize->dSumDualStep : 1.0;

  // cupdlp_scaleVector(iterates->xAverage, iterates->xSum, dPrimalScale,
  // lp->nCols); cupdlp_scaleVector(iterates->yAverage, iterates->ySum,
  // dDualScale, lp->nRows);

  CUPDLP_COPY_VEC(iterates->xAverage->data, iterates->xSum, cupdlp_float,
                  lp->nCols);
  CUPDLP_COPY_VEC(iterates->yAverage->data, iterates->ySum, cupdlp_float,
                  lp->nRows);
  cupdlp_scaleVector(work, dPrimalScale, iterates->xAverage->data, lp->nCols);
  cupdlp_scaleVector(work, dDualScale, iterates->yAverage->data, lp->nRows);

  // Ax(work, iterates->axAverage, iterates->xAverage);
  // ATyCPU(work, iterates->atyAverage, iterates->yAverage);
  Ax(work, iterates->axAverage, iterates->xAverage);
  ATy(work, iterates->atyAverage, iterates->yAverage);

#if PDLP_DEBUG_LOG
  // print norm of x_average
  cupdlp_float debug_pdlp_data_x_average_norm = 0.0;
  cupdlp_twoNormSquared(work, lp->nCols, iterates->xAverage
->data, &debug_pdlp_data_x_average_norm);
  work->debug_pdlp_data_.x_average_norm = debug_pdlp_data_x_average_norm;

  cupdlp_float debug_pdlp_data_ax_average_norm = 0.0;
  cupdlp_twoNormSquared(work, lp->nRows, iterates->axAverage->data, &debug_pdlp_data_ax_average_norm);
  work->debug_pdlp_data_.ax_average_norm = debug_pdlp_data_ax_average_norm;

  // Uncomment this once Yanyu is computing aty_average_norm 
  cupdlp_float debug_pdlp_data_aty_average_norm = 0.0;
  cupdlp_twoNormSquared(work, lp->nCols, iterates->atyAverage->data, &debug_pdlp_data_aty_average_norm);
  work->debug_pdlp_data_.aty_average_norm = debug_pdlp_data_aty_average_norm;
#endif
}

void PDHG_Update_Average(CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  cupdlp_int iter = work->timers->nIter;
  CUPDLPvec *xUpdate = iterates->x[(iter + 1) % 2];
  CUPDLPvec *yUpdate = iterates->y[(iter + 1) % 2];

  // PDLP weighs average iterates in this way
  cupdlp_float dMeanStepSize =
      sqrt(stepsize->dPrimalStep * stepsize->dDualStep);
  // AddToVector(iterates->xSum, dMeanStepSize, xUpdate, lp->nCols);
  // AddToVector(iterates->ySum, dMeanStepSize, yUpdate, lp->nRows);
  cupdlp_axpy(work, lp->nCols, &dMeanStepSize, xUpdate->data, iterates->xSum);
  cupdlp_axpy(work, lp->nRows, &dMeanStepSize, yUpdate->data, iterates->ySum);

  stepsize->dSumPrimalStep += dMeanStepSize;
  stepsize->dSumDualStep += dMeanStepSize;
}

cupdlp_retcode PDHG_Update_Iterate(CUPDLPwork *pdhg) {
  cupdlp_retcode retcode = RETCODE_OK;

#if PDHG_USE_TIMERS
  CUPDLPtimers *timers = pdhg->timers;
  ++timers->nUpdateIterateCalls;
  cupdlp_float dStartTime = getTimeStamp();
#endif

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPiterates *iterates = pdhg->iterates;

  switch (stepsize->eLineSearchMethod) {
    case PDHG_FIXED_LINESEARCH:
      PDHG_Update_Iterate_Constant_Step_Size(pdhg);
      break;
    case PDHG_MALITSKY_POCK_LINESEARCH:
      PDHG_Update_Iterate_Malitsky_Pock(pdhg);
      break;
    case PDHG_ADAPTIVE_LINESEARCH:
      CUPDLP_CALL(PDHG_Update_Iterate_Adaptive_Step_Size(pdhg));
      break;
  }

  PDHG_Update_Average(pdhg);


#if PDHG_USE_TIMERS
  timers->dUpdateIterateTime += getTimeStamp() - dStartTime;
#endif

exit_cleanup:
  return RETCODE_OK;
}
