
#include "pdlp/cupdlp/cupdlp_solver.h"

#include "pdlp/cupdlp/cupdlp_defs.h"
#include "pdlp/cupdlp/cupdlp_linalg.h"
#include "pdlp/cupdlp/cupdlp_proj.h"
#include "pdlp/cupdlp/cupdlp_restart.h"
#include "pdlp/cupdlp/cupdlp_step.h"
#include "pdlp/cupdlp/cupdlp_utils.h"
#include "pdlp/cupdlp/glbopts.h"

void PDHG_Compute_Primal_Feasibility(CUPDLPwork *work, cupdlp_float *primalResidual,
                                    const cupdlp_float *ax, const cupdlp_float *x,
                                    cupdlp_float *dPrimalFeasibility,
                                    cupdlp_float *dPrimalObj) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPscaling *scaling = work->scaling;

  // primal variable violation

  // todo, add this
  //    *dPrimalObj = Dotprod_Neumaier(problem->cost, x, lp->nCols);
  cupdlp_dot(work, lp->nCols, x, problem->cost, dPrimalObj);
  *dPrimalObj = *dPrimalObj * problem->sense_origin + problem->offset;

#if !defined(CUPDLP_CPU) && USE_KERNELS
  primalResidual = work->buffer2;

  cupdlp_primal_feasibility_kernel_cuda(primalResidual, ax, problem->rhs,
                                        work->rowScale, scaling->ifScaled,
                                        problem->nEqs, lp->nRows);

  cupdlp_twoNorm(work, lp->nRows, primalResidual, dPrimalFeasibility);

#else

  // cupdlp_copy(primalResidual, ax, cupdlp_float, lp->nRows);
  CUPDLP_COPY_VEC(primalResidual, ax, cupdlp_float, lp->nRows);

  // AddToVector(primalResidual, -1.0, problem->rhs, lp->nRows);
  cupdlp_float alpha = -1.0;
  cupdlp_axpy(work, lp->nRows, &alpha, problem->rhs, primalResidual);

  //  double dPrimalFeas = 0.0; // Redundant

  // todo, check
  //  cupdlp_projNegative(primalResidual + problem->nEqs, primalResidual +
  //  problem->nEqs, lp->nRows - problem->nEqs);
  //

  cupdlp_projNeg(primalResidual + problem->nEqs, lp->nRows - problem->nEqs);

  if (scaling->ifScaled) {
    // cupdlp_edot(primalResidual, scaling->rowScale, lp->nRows);
    // cupdlp_edot(primalResidual, scaling->rowScale_gpu, lp->nRows);

    // cupdlp_copy_vec(work->buffer3, scaling->rowScale, cupdlp_float,
    // lp->nRows); cupdlp_edot(primalResidual, work->buffer3, lp->nRows);

    cupdlp_edot(primalResidual, work->rowScale, lp->nRows);
  }

  cupdlp_twoNorm(work, lp->nRows, primalResidual, dPrimalFeasibility);

#endif
}

void PDHG_Compute_Dual_Feasibility(CUPDLPwork *work, cupdlp_float *dualResidual,
                                   const cupdlp_float *aty, const cupdlp_float *x,
                                   const cupdlp_float *y, cupdlp_float *dDualFeasibility,
                                   cupdlp_float *dDualObj, cupdlp_float *dComplementarity,
                                   cupdlp_float *dSlackPos, cupdlp_float *dSlackNeg) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPscaling *scaling = work->scaling;
  // todo, compute Neumaier
  //    *dDualObj = Dotprod_Neumaier(problem->rhs, y, lp->nRows);
  cupdlp_dot(work, lp->nRows, y, problem->rhs, dDualObj);

  #if !defined(CUPDLP_CPU) && USE_KERNELS
  dualResidual = work->buffer2;

  cupdlp_float alpha;
  cupdlp_dual_feasibility_kernel_1_cuda(dualResidual, aty,
                                        problem->cost, lp->nCols);

  cupdlp_dual_feasibility_kernel_2_cuda(dSlackPos, dualResidual,
                                        problem->hasLower, lp->nCols);

  cupdlp_float temp = 0.0;
  cupdlp_dot(work, lp->nCols, dSlackPos, resobj->dLowerFiltered, &temp);
  *dDualObj += temp;

  cupdlp_dual_feasibility_kernel_3_cuda(dSlackNeg, dualResidual,
                                        problem->hasUpper, lp->nCols);

  cupdlp_dot(work, lp->nCols, dSlackNeg, resobj->dUpperFiltered, &temp);
  *dDualObj -= temp;

  *dDualObj = *dDualObj * problem->sense_origin + problem->offset;

  alpha = -1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, dSlackPos, dualResidual);
  alpha = 1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, dSlackNeg, dualResidual);

  if (scaling->ifScaled) {
    cupdlp_edot(dualResidual, work->colScale, lp->nCols);
  }
  cupdlp_twoNorm(work, lp->nCols, dualResidual, dDualFeasibility);

#else
  // *dComplementarity = 0.0;

  // @note:
  // original dual residual in pdlp:
  // they compute:
  //    violation   +  reduced cost
  //  |max(-y, 0)|  + |(I-\Pi)(c-Α'\nu)|
  // compute c - A'y

  CUPDLP_COPY_VEC(dualResidual, aty, cupdlp_float, lp->nCols);
  cupdlp_float alpha = -1.0;
  cupdlp_scaleVector(work, alpha, dualResidual, lp->nCols);

  alpha = 1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, problem->cost, dualResidual);

  //    julia version
  //        function compute_reduced_costs_from_primal_gradient_kernel!(
  //            primal_gradient::CuDeviceVector{Float64},
  //            isfinite_variable_lower_bound::CuDeviceVector{Bool},
  //            isfinite_variable_upper_bound::CuDeviceVector{Bool},
  //            num_variables::Int64,
  //            reduced_costs::CuDeviceVector{Float64},
  //            reduced_costs_violation::CuDeviceVector{Float64},
  //        )
  //            tx = threadIdx().x + (blockDim().x * (blockIdx().x - 0x1))
  //            if tx <= num_variables
  //                @inbounds begin
  //                    reduced_costs[tx] = max(primal_gradient[tx], 0.0) *
  //                    isfinite_variable_lower_bound[tx] +
  //                    min(primal_gradient[tx], 0.0) *
  //                    isfinite_variable_upper_bound[tx]
  //
  //                    reduced_costs_violation[tx] = primal_gradient[tx] -
  //                    reduced_costs[tx]
  //                end
  //            end
  //            return
  //        end

  CUPDLP_COPY_VEC(dSlackPos, dualResidual, cupdlp_float, lp->nCols);

  cupdlp_projPos(dSlackPos, lp->nCols);

  cupdlp_edot(dSlackPos, problem->hasLower, lp->nCols);

  cupdlp_float temp = 0.0;
  // cupdlp_dot(work, lp->nCols, x, dSlackPos, &temp);
  // *dComplementarity += temp;
  // cupdlp_dot(work, lp->nCols, dSlackPos, resobj->dLowerFiltered, &temp);
  // *dComplementarity -= temp;
  cupdlp_dot(work, lp->nCols, dSlackPos, resobj->dLowerFiltered, &temp);
  *dDualObj += temp;

  CUPDLP_COPY_VEC(dSlackNeg, dualResidual, cupdlp_float, lp->nCols);

  cupdlp_projNeg(dSlackNeg, lp->nCols);

  cupdlp_scaleVector(work, -1.0, dSlackNeg, lp->nCols);

  cupdlp_edot(dSlackNeg, problem->hasUpper, lp->nCols);

  // cupdlp_dot(work, lp->nCols, x, dSlackNeg, &temp);
  // *dComplementarity -= temp;
  // cupdlp_dot(work, lp->nCols, dSlackNeg, resobj->dUpperFiltered, &temp);
  // *dComplementarity += temp;
  cupdlp_dot(work, lp->nCols, dSlackNeg, resobj->dUpperFiltered, &temp);
  *dDualObj -= temp;

  *dDualObj = *dDualObj * problem->sense_origin + problem->offset;

  alpha = -1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, dSlackPos, dualResidual);
  alpha = 1.0;
  cupdlp_axpy(work, lp->nCols, &alpha, dSlackNeg, dualResidual);

  if (scaling->ifScaled) {
    // cupdlp_edot(dualResidual, scaling->colScale, lp->nCols);
    // cupdlp_edot(dualResidual, scaling->colScale_gpu, lp->nCols);

    // cupdlp_copy_vec(work->buffer3, scaling->colScale, cupdlp_float,
    // lp->nCols); cupdlp_edot(dualResidual, work->buffer3, lp->nCols);

    cupdlp_edot(dualResidual, work->colScale, lp->nCols);
  }

  cupdlp_twoNorm(work, lp->nCols, dualResidual, dDualFeasibility);

#endif
}

void PDHG_Compute_Primal_Infeasibility(CUPDLPwork *work, const cupdlp_float *y,
                                       const cupdlp_float *dSlackPos,
                                       const cupdlp_float *dSlackNeg,
                                       const cupdlp_float *aty,
                                       const cupdlp_float dualObj,
                                       cupdlp_float *dPrimalInfeasObj,
                                       cupdlp_float *dPrimalInfeasRes) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPscaling *scaling = work->scaling;

  cupdlp_float alpha;

  cupdlp_float yNrmSq = 1.0;
  cupdlp_float slackPosNrmSq = 1.0;
  cupdlp_float slackNegSq = 1.0;
  cupdlp_float dScale = 1.0;

  // cupdlp_float dConstrResSq = 0.0;
  // y and lambda must be feasible, no need to check bound
  // cupdlp_float dBoundLbResSq = 0.0;
  // cupdlp_float dBoundUbResSq = 0.0;

#if !defined(CUPDLP_CPU) && USE_KERNELS
  cupdlp_twoNormSquared(work, problem->data->nRows, y, &yNrmSq);
  cupdlp_twoNormSquared(work, problem->data->nCols, dSlackPos, &slackPosNrmSq);
  cupdlp_twoNormSquared(work, problem->data->nCols, dSlackNeg, &slackNegSq);
  dScale = sqrt(yNrmSq + slackPosNrmSq + slackNegSq);
  // dScale /= sqrt(problem->data->nRows + 2 * problem->data->nCols);
  if (dScale < 1e-8) {
    dScale = 1.0;
  }

  // dual obj
  *dPrimalInfeasObj =
      (dualObj - problem->offset) / problem->sense_origin / dScale;

  // dual constraints [ATy1 + GTy2 + lambda]

  cupdlp_primal_infeasibility_kernel_cuda(work->buffer2, aty, dSlackPos,
                                          dSlackNeg, work->colScale, 1.0 / dScale,
                                          scaling->ifScaled, problem->data->nCols);

  // cupdlp_twoNormSquared(work, problem->data->nCols, resobj->dualInfeasConstr,
  //                       &dConstrResSq);
  cupdlp_twoNorm(work, problem->data->nCols, work->buffer2,
                 dPrimalInfeasRes);

  // dual bound
  // always satisfied, no need to check

#else

  // y, ldb ray
  CUPDLP_COPY_VEC(resobj->dualInfeasRay, y, cupdlp_float, problem->data->nRows);
  CUPDLP_COPY_VEC(resobj->dualInfeasLbRay, dSlackPos, cupdlp_float,
                  problem->data->nCols);
  CUPDLP_COPY_VEC(resobj->dualInfeasUbRay, dSlackNeg, cupdlp_float,
                  problem->data->nCols);
  cupdlp_twoNormSquared(work, problem->data->nRows, resobj->dualInfeasRay,
                        &yNrmSq);
  cupdlp_twoNormSquared(work, problem->data->nCols, resobj->dualInfeasLbRay,
                        &slackPosNrmSq);
  cupdlp_twoNormSquared(work, problem->data->nCols, resobj->dualInfeasUbRay,
                        &slackNegSq);
  dScale = sqrt(yNrmSq + slackPosNrmSq + slackNegSq);
  // dScale /= sqrt(problem->data->nRows + 2 * problem->data->nCols);
  if (dScale < 1e-8) {
    dScale = 1.0;
  }
  cupdlp_scaleVector(work, 1 / dScale, resobj->dualInfeasRay,
                     problem->data->nRows);
  cupdlp_scaleVector(work, 1 / dScale, resobj->dualInfeasLbRay,
                     problem->data->nCols);
  cupdlp_scaleVector(work, 1 / dScale, resobj->dualInfeasUbRay,
                     problem->data->nCols);

  // dual obj
  *dPrimalInfeasObj =
      (dualObj - problem->offset) / problem->sense_origin / dScale;

  // dual constraints [ATy1 + GTy2 + lambda]
  CUPDLP_COPY_VEC(resobj->dualInfeasConstr, aty, cupdlp_float,
                  problem->data->nCols);
  cupdlp_scaleVector(work, 1.0 / dScale, resobj->dualInfeasConstr,
                     problem->data->nCols);
  alpha = 1.0;
  cupdlp_axpy(work, problem->data->nCols, &alpha, resobj->dualInfeasLbRay,
              resobj->dualInfeasConstr);
  alpha = -1.0;
  cupdlp_axpy(work, problem->data->nCols, &alpha, resobj->dualInfeasUbRay,
              resobj->dualInfeasConstr);
  if (scaling->ifScaled) {
    cupdlp_edot(resobj->dualInfeasConstr, work->colScale, problem->data->nCols);
  }
  // cupdlp_twoNormSquared(work, problem->data->nCols, resobj->dualInfeasConstr,
  //                       &dConstrResSq);
  cupdlp_twoNorm(work, problem->data->nCols, resobj->dualInfeasConstr,
                 dPrimalInfeasRes);

  // dual bound
  // always satisfied, no need to check

#endif

}

void PDHG_Compute_Dual_Infeasibility(CUPDLPwork *work, const cupdlp_float *x,
                                     const cupdlp_float *ax,
                                     const cupdlp_float primalObj,
                                     cupdlp_float *dDualInfeasObj,
                                     cupdlp_float *dDualInfeasRes) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPscaling *scaling = work->scaling;
  cupdlp_float pScale = 1.0;
  cupdlp_float pConstrResSq = 0.0;
  cupdlp_float pBoundLbResSq = 0.0;
  cupdlp_float pBoundUbResSq = 0.0;

  #if !defined(CUPDLP_CPU) && USE_KERNELS
  // x ray
  cupdlp_twoNorm(work, problem->data->nCols, x, &pScale);
  // pScale /= sqrt(problem->data->nCols);
  if (pScale < 1e-8) {
    pScale = 1.0;
  }

  // primal obj
  *dDualInfeasObj =
      (primalObj - problem->offset) / problem->sense_origin / pScale;

  // primal constraints [Ax, min(Gx, 0)]
  cupdlp_dual_infeasibility_kernel_constr_cuda(work->buffer2, ax,
                                              work->rowScale, 1.0 / pScale,
                                              scaling->ifScaled,
                                              problem->nEqs,
                                              problem->data->nRows);
  cupdlp_twoNormSquared(work, problem->data->nRows, work->buffer2,
                        &pConstrResSq);

  // primal bound
  // lb
  cupdlp_dual_infeasibility_kernel_lb_cuda(work->buffer2, x,
                                           problem->hasLower,
                                           work->colScale,
                                           1.0 / pScale, scaling->ifScaled,
                                           problem->data->nCols);
  cupdlp_twoNormSquared(work, problem->data->nCols, work->buffer2,
                        &pBoundLbResSq);
  // ub
  cupdlp_dual_infeasibility_kernel_ub_cuda(work->buffer2, x,
                                           problem->hasUpper,
                                           work->colScale,
                                           1.0 / pScale, scaling->ifScaled,
                                           problem->data->nCols);
  cupdlp_twoNormSquared(work, problem->data->nCols, work->buffer2,
                        &pBoundUbResSq);

  // sum up
  *dDualInfeasRes = sqrt(pConstrResSq + pBoundLbResSq + pBoundUbResSq);

#else

  // x ray
  CUPDLP_COPY_VEC(resobj->primalInfeasRay, x, cupdlp_float,
                  problem->data->nCols);
  cupdlp_twoNorm(work, problem->data->nCols, resobj->primalInfeasRay, &pScale);
  // pScale /= sqrt(problem->data->nCols);
  if (pScale < 1e-8) {
    pScale = 1.0;
  }
  cupdlp_scaleVector(work, 1.0 / pScale, resobj->primalInfeasRay,
                     problem->data->nCols);

  // primal obj
  *dDualInfeasObj =
      (primalObj - problem->offset) / problem->sense_origin / pScale;

  // primal constraints [Ax, min(Gx, 0)]
  CUPDLP_COPY_VEC(resobj->primalInfeasConstr, ax, cupdlp_float,
                  problem->data->nRows);
  cupdlp_scaleVector(work, 1.0 / pScale, resobj->primalInfeasConstr,
                     problem->data->nRows);
  cupdlp_projNeg(resobj->primalInfeasConstr + problem->nEqs,
                 problem->data->nRows - problem->nEqs);
  if (scaling->ifScaled) {
    cupdlp_edot(resobj->primalInfeasConstr, work->rowScale,
                problem->data->nRows);
  }
  cupdlp_twoNormSquared(work, problem->data->nRows, resobj->primalInfeasConstr,
                        &pConstrResSq);

  // primal bound
  // lb
  CUPDLP_COPY_VEC(resobj->primalInfeasBound, resobj->primalInfeasRay,
                  cupdlp_float, problem->data->nCols);
  cupdlp_projNeg(resobj->primalInfeasBound, problem->data->nCols);
  cupdlp_edot(resobj->primalInfeasBound, problem->hasLower,
              problem->data->nCols);
  if (scaling->ifScaled) {
    cupdlp_ediv(resobj->primalInfeasBound, work->colScale,
                problem->data->nCols);
  }
  cupdlp_twoNormSquared(work, problem->data->nCols, resobj->primalInfeasBound,
                        &pBoundLbResSq);
  // ub
  CUPDLP_COPY_VEC(resobj->primalInfeasBound, resobj->primalInfeasRay,
                  cupdlp_float, problem->data->nCols);
  cupdlp_projPos(resobj->primalInfeasBound, problem->data->nCols);
  cupdlp_edot(resobj->primalInfeasBound, problem->hasUpper,
              problem->data->nCols);
  if (scaling->ifScaled) {
    cupdlp_ediv(resobj->primalInfeasBound, work->colScale,
                problem->data->nCols);
  }
  cupdlp_twoNormSquared(work, problem->data->nCols, resobj->primalInfeasBound,
                        &pBoundUbResSq);

  // sum up
  *dDualInfeasRes = sqrt(pConstrResSq + pBoundLbResSq + pBoundUbResSq);

#endif
}

// must be called after PDHG_Compute_Residuals(CUPDLPwork *work)
// because it needs resobj->dSlackPos and resobj->dSlackNeg
void PDHG_Compute_Infeas_Residuals(CUPDLPwork *work) {
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPresobj *resobj = work->resobj;


  cupdlp_int iter = work->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];
  CUPDLPvec *ax = iterates->ax[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];

  // current solution
  PDHG_Compute_Primal_Infeasibility(work, y->data, resobj->dSlackPos,
                                    resobj->dSlackNeg, aty->data,
                                    resobj->dDualObj, &resobj->dPrimalInfeasObj,
                                    &resobj->dPrimalInfeasRes);
  PDHG_Compute_Dual_Infeasibility(work, x->data, ax->data,
                                  resobj->dPrimalObj, &resobj->dDualInfeasObj,
                                  &resobj->dDualInfeasRes);

  // average solution
  PDHG_Compute_Primal_Infeasibility(
      work, iterates->yAverage->data, resobj->dSlackPosAverage,
      resobj->dSlackNegAverage, iterates->atyAverage->data,
      resobj->dDualObjAverage, &resobj->dPrimalInfeasObjAverage,
      &resobj->dPrimalInfeasResAverage);
  PDHG_Compute_Dual_Infeasibility(
      work, iterates->xAverage->data, iterates->axAverage->data,
      resobj->dPrimalObjAverage, &resobj->dDualInfeasObjAverage,
      &resobj->dDualInfeasResAverage);

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDHG_Compute_Residuals(CUPDLPwork *work) {
#if problem_USE_TIMERS
  ++problem->nComputeResidualsCalls;
  double dStartTime = getTimeStamp();
#endif
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPresobj *resobj = work->resobj;
  CUPDLPiterates *iterates = work->iterates;
  CUPDLPscaling *scaling = work->scaling;
  CUPDLPsettings *settings = work->settings;

  cupdlp_int iter = work->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];
  CUPDLPvec *ax = iterates->ax[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];


  PDHG_Compute_Primal_Feasibility(work, resobj->primalResidual,
                                  ax->data, x->data,
                                  &resobj->dPrimalFeas, &resobj->dPrimalObj);
  PDHG_Compute_Dual_Feasibility(
      work, resobj->dualResidual, aty->data, x->data,
      y->data, &resobj->dDualFeas, &resobj->dDualObj,
      &resobj->dComplementarity, resobj->dSlackPos, resobj->dSlackNeg);

  PDHG_Compute_Primal_Feasibility(
      work, resobj->primalResidualAverage, iterates->axAverage->data,
      iterates->xAverage->data, &resobj->dPrimalFeasAverage,
      &resobj->dPrimalObjAverage);
  PDHG_Compute_Dual_Feasibility(
      work, resobj->dualResidualAverage, iterates->atyAverage->data,
      iterates->xAverage->data, iterates->yAverage->data,
      &resobj->dDualFeasAverage, &resobj->dDualObjAverage,
      &resobj->dComplementarityAverage, resobj->dSlackPosAverage,
      resobj->dSlackNegAverage);

  // resobj->dPrimalObj /= (scaling->dObjScale * scaling->dObjScale);
  // resobj->dDualObj /= (scaling->dObjScale * scaling->dObjScale);
  resobj->dDualityGap = resobj->dPrimalObj - resobj->dDualObj;
  resobj->dRelObjGap =
      fabs(resobj->dPrimalObj - resobj->dDualObj) /
      (1.0 + fabs(resobj->dPrimalObj) + fabs(resobj->dDualObj));

  // resobj->dPrimalObjAverage /= scaling->dObjScale * scaling->dObjScale;
  // resobj->dDualObjAverage /= scaling->dObjScale * scaling->dObjScale;
  resobj->dDualityGapAverage =
      resobj->dPrimalObjAverage - resobj->dDualObjAverage;
  resobj->dRelObjGapAverage =
      fabs(resobj->dPrimalObjAverage - resobj->dDualObjAverage) /
      (1.0 + fabs(resobj->dPrimalObjAverage) + fabs(resobj->dDualObjAverage));

#if problem_USE_TIMERS
  problem->dComputeResidualsTime += getTimeStamp() - dStartTime;
#endif
}

void PDHG_Init_Variables(const cupdlp_int* has_variables, CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;

  cupdlp_int iter = work->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];
  CUPDLPvec *ax = iterates->ax[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];

  if (!*has_variables) 
    CUPDLP_ZERO_VEC(x->data, cupdlp_float, lp->nCols);

  // XXX: PDLP Does not project x0,  so we uncomment for 1-1 comparison

  PDHG_Project_Bounds(work, x->data);
  double x_norm;
  cupdlp_twoNorm(work, lp->nCols, x->data, &x_norm);
  if (work->settings->nLogLevel>0)
    cupdlp_printf("||x0||_2 = %e\n", x_norm);
  
  // cupdlp_zero(iterates->y, cupdlp_float, lp->nRows);
  if (!*has_variables)
    CUPDLP_ZERO_VEC(y->data, cupdlp_float, lp->nRows);

  // Ax(work, iterates->ax, iterates->x);
  // ATyCPU(work, iterates->aty, iterates->y);
  Ax(work, ax, x);
  ATy(work, aty, y);

#if PDLP_DEBUG_LOG
  double ax_norm;
  cupdlp_twoNorm(work, lp->nRows, ax->data, &ax_norm);
  if (work->settings->nLogLevel>0)
    cupdlp_printf("||Ax0||_2 = %e\n", ax_norm);
  work->debug_pdlp_data_.ax_norm = ax_norm;
  double aty_norm;
  cupdlp_twoNorm(work, lp->nCols, aty->data, &aty_norm);
  work->debug_pdlp_data_.aty_norm = aty_norm;
#endif

  // cupdlp_zero(iterates->xSum, cupdlp_float, lp->nCols);
  // cupdlp_zero(iterates->ySum, cupdlp_float, lp->nRows);
  // cupdlp_zero(iterates->xAverage, cupdlp_float, lp->nCols);
  // cupdlp_zero(iterates->yAverage, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->xSum, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->ySum, cupdlp_float, lp->nRows);
  CUPDLP_ZERO_VEC(iterates->xAverage->data, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yAverage->data, cupdlp_float, lp->nRows);

  PDHG_Project_Bounds(work, iterates->xSum);
  PDHG_Project_Bounds(work, iterates->xAverage->data);

  stepsize->dSumPrimalStep = 0.0;
  stepsize->dSumDualStep = 0.0;

  CUPDLP_ZERO_VEC(iterates->xLastRestart, cupdlp_float, lp->nCols);
  CUPDLP_ZERO_VEC(iterates->yLastRestart, cupdlp_float, lp->nRows);
}

/* TODO: this function seems considering
 *       l1 <= Ax <= u1
 *       l2 <=  x <= u2
 *       needs rewritten for current formulation
 */
void PDHG_Check_Data(CUPDLPwork *work) {
  CUPDLPproblem *problem = work->problem;
  CUPDLPdata *lp = problem->data;
  CUPDLPstepsize *stepsize = work->stepsize;
  CUPDLPiterates *iterates = work->iterates;
  cupdlp_int nFreeCol = 0;
  cupdlp_int nFixedCol = 0;
  cupdlp_int nUpperCol = 0;
  cupdlp_int nLowerCol = 0;
  cupdlp_int nRangedCol = 0;
  cupdlp_int nFreeRow = 0;
  cupdlp_int nFixedRow = 0;
  cupdlp_int nUpperRow = 0;
  cupdlp_int nLowerRow = 0;
  cupdlp_int nRangedRow = 0;


  cupdlp_int iter = work->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];

  for (cupdlp_int iSeq = 0; iSeq < lp->nCols; ++iSeq) {
    cupdlp_bool hasLower = problem->lower[iSeq] > -INFINITY;
    cupdlp_bool hasUpper = problem->upper[iSeq] < +INFINITY;

    if (!hasLower && !hasUpper) {
      ++nFreeCol;
      if (work->settings->nLogLevel>0) 
	cupdlp_printf("Warning: variable %d is free.", iSeq);
    }

    if (hasLower && hasUpper) {
      if (problem->lower[iSeq] == problem->upper[iSeq]) {
        ++nFixedCol;
        // cupdlp_printf( "Warning: variable %d is fixed.", iSeq);
      } else
        ++nRangedCol;
    }

    if (hasLower) {
      // XXX: uncommented for PDLP comparison
      // CUPDLP_ASSERT(x[iSeq] >= problem->lower[iSeq]);
      nLowerCol += !hasUpper;
    }

    if (hasUpper) {
      // XXX: uncommented for PDLP comparison
      // CUPDLP_ASSERT(x[iSeq] <= problem->upper[iSeq]);
      nUpperCol += !hasLower;
    }
  }

  for (cupdlp_int iSeq = lp->nCols; iSeq < lp->nCols; ++iSeq) {
    cupdlp_bool hasLower = problem->lower[iSeq] > -INFINITY;
    cupdlp_bool hasUpper = problem->upper[iSeq] < +INFINITY;

    if (!hasLower && !hasUpper) {
      ++nFreeRow;
      if (work->settings->nLogLevel>0) 
	cupdlp_printf("Warning: row %d is free.", iSeq - lp->nCols);
    }

    if (hasLower && hasUpper) {
      if (problem->lower[iSeq] == problem->upper[iSeq])
        ++nFixedRow;
      else
        ++nRangedRow;
    }

    if (hasLower) {
      // CUPDLP_ASSERT(x[iSeq] >= problem->lower[iSeq]);
      nLowerRow += !hasUpper;
    }

    if (hasUpper) {
      // CUPDLP_ASSERT(x[iSeq] <= problem->upper[iSeq]);
      nUpperRow += !hasLower;
    }
  }

  for (cupdlp_int iRow = 0; iRow < lp->nRows; ++iRow) {
    CUPDLP_ASSERT(y->data[iRow] < +INFINITY);
    CUPDLP_ASSERT(y->data[iRow] > -INFINITY);
  }

  for (cupdlp_int iRow = 0; iRow < lp->nRows; ++iRow) {
    if (problem->data->csr_matrix->rowMatBeg[iRow + 1] -
            problem->data->csr_matrix->rowMatBeg[iRow] ==
        1) {
      if (work->settings->nLogLevel>0) 
	cupdlp_printf("Warning: row %d is a singleton row.", iRow);
    }
  }

  CUPDLP_ASSERT(nRangedRow == 0);
  if (work->settings->nLogLevel>0) {
    cupdlp_printf("nFreeCol  : %d\n", nFreeCol);
    cupdlp_printf("nFixedCol : %d\n", nFixedCol);
    cupdlp_printf("nRangedCol: %d\n", nRangedCol);
    cupdlp_printf("nLowerCol : %d\n", nLowerCol);
    cupdlp_printf("nUpperCol : %d\n", nUpperCol);
    cupdlp_printf("nFreeRow  : %d\n", nFreeRow);
    cupdlp_printf("nFixedRow : %d\n", nFixedRow);
    cupdlp_printf("nRangedRow: %d\n", nRangedRow);
    cupdlp_printf("nLowerRow : %d\n", nLowerRow);
    cupdlp_printf("nUpperRow : %d\n", nUpperRow);
  }

  // We need to test problems ranged row-bounds more carefully.
  CUPDLP_ASSERT(nRangedRow == 0);
}

termination_code PDHG_Check_Primal_Infeasibility(
    CUPDLPwork *pdhg, cupdlp_float dPrimalInfeasObj,
    cupdlp_float dPrimalInfeasRes) {
  CUPDLPresobj *resobj = pdhg->resobj;

  termination_code primalCode = FEASIBLE;

  if (dPrimalInfeasObj > 0.0) {
    if (dPrimalInfeasRes < resobj->dFeasTol * dPrimalInfeasObj)
      primalCode = INFEASIBLE;
  }

  return primalCode;
}

termination_code PDHG_Check_Dual_Infeasibility(CUPDLPwork *pdhg,
                                               cupdlp_float dDualInfeasObj,
                                               cupdlp_float dDualInfeasRes) {
  CUPDLPresobj *resobj = pdhg->resobj;

  termination_code dualCode = FEASIBLE;

  if (dDualInfeasObj < 0.0) {
    if (dDualInfeasRes < -resobj->dFeasTol * dDualInfeasObj)
      dualCode = INFEASIBLE;
  }

  return dualCode;
}

termination_code PDHG_Check_Infeasibility(CUPDLPwork *pdhg, int bool_print) {
  CUPDLPresobj *resobj = pdhg->resobj;
  termination_code t_code = FEASIBLE;

  // current solution

  // primal infeasibility
  if (PDHG_Check_Primal_Infeasibility(pdhg, resobj->dPrimalInfeasObj,
                                      resobj->dPrimalInfeasRes) == INFEASIBLE) {
    resobj->primalCode = INFEASIBLE;
    resobj->termInfeasIterate = LAST_ITERATE;
    t_code = INFEASIBLE_OR_UNBOUNDED;
  }

  // dual infeasibility
  if (PDHG_Check_Dual_Infeasibility(pdhg, resobj->dDualInfeasObj,
                                    resobj->dDualInfeasRes) == INFEASIBLE) {
    resobj->dualCode = INFEASIBLE;
    resobj->termInfeasIterate = LAST_ITERATE;
    t_code = INFEASIBLE_OR_UNBOUNDED;
  }

  // average solution
  // primal infeasibility
  if (PDHG_Check_Primal_Infeasibility(pdhg, resobj->dPrimalInfeasObjAverage,
                                      resobj->dPrimalInfeasResAverage) ==
      INFEASIBLE) {
    resobj->primalCode = INFEASIBLE;
    resobj->termInfeasIterate = AVERAGE_ITERATE;
    t_code = INFEASIBLE_OR_UNBOUNDED;
  }

  // dual infeasibility
  if (PDHG_Check_Dual_Infeasibility(pdhg, resobj->dDualInfeasObjAverage,
                                    resobj->dDualInfeasResAverage) ==
      INFEASIBLE) {
    resobj->dualCode = INFEASIBLE;
    resobj->termInfeasIterate = AVERAGE_ITERATE;
    t_code = INFEASIBLE_OR_UNBOUNDED;
  }

  if (bool_print) {
    printf("Last iter:\n");
    printf("  Primal obj = %+.4e, res = %+.4e\n", resobj->dPrimalInfeasObj,
           resobj->dPrimalInfeasRes);
    printf("  Dual   obj = %+.4e, res = %+.4e\n", resobj->dDualInfeasObj,
           resobj->dDualInfeasRes);
    printf("Average iter:\n");
    printf("  Primal obj = %+.4e, res = %+.4e\n",
           resobj->dPrimalInfeasObjAverage, resobj->dPrimalInfeasResAverage);
    printf("  Dual   obj = %+.4e, res = %+.4e\n", resobj->dDualInfeasObjAverage,
           resobj->dDualInfeasResAverage);
  }

  return t_code;
}

cupdlp_bool PDHG_Check_Termination(CUPDLPwork *pdhg, int bool_print) {
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPscaling *scaling = pdhg->scaling;
#if PDHG_DISPLAY_TERMINATION_CHECK
  // todo, check, is it correct
  if (bool_print) {
    cupdlp_printf(
        "Termination check: %e|%e  %e|%e  %e|%e\n", resobj->dPrimalFeas,
        settings->dPrimalTol * (1.0 + scaling->dNormRhs), resobj->dDualFeas,
        settings->dDualTol * (1.0 + scaling->dNormCost), resobj->dRelObjGap,
        settings->dGapTol);
  }

#endif
  int bool_pass =
    (resobj->dPrimalFeas < settings->dPrimalTol * (1.0 + scaling->dNormRhs)) &&
    (resobj->dDualFeas < settings->dDualTol * (1.0 + scaling->dNormCost));
  bool_pass = bool_pass && (resobj->dRelObjGap < settings->dGapTol);
  return bool_pass;
}

cupdlp_bool PDHG_Check_Termination_Average(CUPDLPwork *pdhg, int bool_print) {
  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPscaling *scaling = pdhg->scaling;
#if PDHG_DISPLAY_TERMINATION_CHECK
  if (bool_print) {
    cupdlp_printf("Termination check: %e|%e  %e|%e  %e|%e\n",
                  resobj->dPrimalFeasAverage,
                  settings->dPrimalTol * (1.0 + scaling->dNormRhs),
                  resobj->dDualFeasAverage,
                  settings->dDualTol * (1.0 + scaling->dNormCost),
                  resobj->dRelObjGapAverage, settings->dGapTol);
  }
#endif
  int bool_pass = ((resobj->dPrimalFeasAverage <
                    settings->dPrimalTol * (1.0 + scaling->dNormRhs)) &&
                   (resobj->dDualFeasAverage <
                    settings->dDualTol * (1.0 + scaling->dNormCost)) &&
                   (resobj->dRelObjGapAverage < settings->dGapTol));
  return bool_pass;
}

void PDHG_Print_Header(CUPDLPwork *pdhg) {
  // cupdlp_printf("%9s  %15s  %15s   %8s  %8s  %10s  %8s %7s\n", "Iter",
  //               "Primal.Obj", "Dual.Obj", "Gap", "Compl", "Primal.Inf",
  //               "Dual.Inf", "Time");
  cupdlp_printf("%9s  %15s  %15s   %8s  %10s  %8s %7s\n", "Iter", "Primal.Obj",
                "Dual.Obj", "Gap", "Primal.Inf", "Dual.Inf", "Time");
}

void PDHG_Print_Iter(CUPDLPwork *pdhg) {
  /* Format time as xxx.yy for < 1000s and as integer afterwards. */
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;
  char timeString[8];
  if (timers->dSolvingTime < 100.0)
    cupdlp_snprintf(timeString, 8, "%6.2fs", timers->dSolvingTime);
  else
    cupdlp_snprintf(timeString, 8, "%6ds", (cupdlp_int)timers->dSolvingTime);

  CUPDLPscaling *scaling = pdhg->scaling;
  const double gap_log_value = resobj->dRelObjGap;
  const double primal_feas_log_value = resobj->dPrimalFeas / (1.0 + scaling->dNormRhs);
  const double dual_feas_log_value = resobj->dDualFeas / (1.0 + scaling->dNormCost);
  cupdlp_printf("%9d  %+15.8e  %+15.8e  %+8.2e  %10.2e  %8.2e %7s [L]\n",
                timers->nIter, resobj->dPrimalObj, resobj->dDualObj,
		//                resobj->dDualityGap, resobj->dPrimalFeas, resobj->dDualFeas,
		gap_log_value, primal_feas_log_value, dual_feas_log_value,
                timeString);
}

void PDHG_Print_Iter_Average(CUPDLPwork *pdhg) {
  /* Format time as xxx.yy for < 1000s and as integer afterwards. */
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPtimers *timers = pdhg->timers;
  char timeString[8];
  if (timers->dSolvingTime < 100.0)
    cupdlp_snprintf(timeString, 8, "%6.2fs", timers->dSolvingTime);
  else
    cupdlp_snprintf(timeString, 8, "%6ds", (cupdlp_int)timers->dSolvingTime);

  CUPDLPscaling *scaling = pdhg->scaling;
  const double gap_log_value = resobj->dRelObjGapAverage;
  const double primal_feas_log_value = resobj->dPrimalFeasAverage / (1.0 + scaling->dNormRhs);
  const double dual_feas_log_value = resobj->dDualFeasAverage / (1.0 + scaling->dNormCost);
  cupdlp_printf("%9d  %+15.8e  %+15.8e  %+8.2e  %10.2e  %8.2e %7s [A]\n",
                timers->nIter, resobj->dPrimalObjAverage,
                resobj->dDualObjAverage,
		//		resobj->dDualityGapAverage, resobj->dPrimalFeasAverage, resobj->dDualFeasAverage,
		gap_log_value, primal_feas_log_value, dual_feas_log_value,
                timeString);
}

void PDHG_Compute_SolvingTime(CUPDLPwork *pdhg) {
  CUPDLPtimers *timers = pdhg->timers;
  timers->dSolvingTime = getTimeStamp() - timers->dSolvingBeg;
}

cupdlp_retcode PDHG_Solve(const cupdlp_int* has_variables, CUPDLPwork *pdhg) {
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPstepsize *stepsize = pdhg->stepsize;
  CUPDLPsettings *settings = pdhg->settings;
  CUPDLPresobj *resobj = pdhg->resobj;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPtimers *timers = pdhg->timers;
  CUPDLPscaling *scaling = pdhg->scaling;

  timers->nIter = 0;
  timers->dSolvingBeg = getTimeStamp();

#if PDLP_DEBUG_LOG
  pdhg->debug_pdlp_log_file_ = fopen("cuPDLP.log", "w");
  assert(pdhg->debug_pdlp_log_file_);
  debugPdlpDataInitialise(&pdhg->debug_pdlp_data_);
#endif
  // PDHG_Init_Data does nothing!
  PDHG_Init_Data(pdhg);

  CUPDLP_CALL(PDHG_Init_Step_Sizes(pdhg));

  PDHG_Init_Variables(has_variables, pdhg);

  // todo: translate check_data into cuda or do it on cpu
  // PDHG_Check_Data(pdhg);

  // PDHG_Print_Header(pdhg);

  // Repeat the iteration logging header periodically if logging style
  // is minimal (pdhg->settings->nLogLevel=1), initialising
  // iter_log_since_header so that an initial header is printed
  const int iter_log_between_header = 50;
  int iter_log_since_header = iter_log_between_header;
#if PDLP_DEBUG_LOG
  debugPdlpIterHeaderLog(pdhg->debug_pdlp_log_file_);
#endif
  // Main loop!
  for (timers->nIter = 0; timers->nIter < settings->nIterLim; ++timers->nIter) {
#if PDLP_DEBUG_LOG
    debugPdlpIterLog(pdhg->debug_pdlp_log_file_,
		     timers->nIter,
		     &pdhg->debug_pdlp_data_,
		     pdhg->stepsize->dBeta,
		     pdhg->stepsize->dPrimalStep,
		     pdhg->stepsize->dDualStep);
#endif
    PDHG_Compute_SolvingTime(pdhg);
    
#if CUPDLP_DUMP_ITERATES_STATS && CUPDLP_DEBUG
    PDHG_Dump_Stats(pdhg);
#endif
    int bool_checking = (timers->nIter < 10) ||
                        (timers->nIter == (settings->nIterLim - 1)) ||
                        (timers->dSolvingTime > settings->dTimeLim);
    int bool_print = 0;
#if CUPDLP_DEBUG
    bool_checking = (bool_checking || !(timers->nIter % CUPDLP_DEBUG_INTERVAL));
    bool_print = bool_checking;
#else
    bool_checking =
        (bool_checking || !(timers->nIter % CUPDLP_RELEASE_INTERVAL));
    bool_print =
        (bool_checking && !(timers->nIter % (CUPDLP_RELEASE_INTERVAL *
                                             settings->nLogInterval))) ||
        (timers->nIter == (settings->nIterLim - 1)) ||
        (timers->dSolvingTime > settings->dTimeLim);
    // Ensure that bool_print is false if the logging level has been
    // set to 0 via the HiGHS option
    bool_print = pdhg->settings->nLogLevel>0 && bool_print;
#endif
    // Full printing is false only if the logging level has been set
    // to 0 or 1 via the HiGHS option
    int full_print = pdhg->settings->nLogLevel >= 2;
    if (bool_checking) {
      PDHG_Compute_Average_Iterate(pdhg);
      PDHG_Compute_Residuals(pdhg);
      PDHG_Compute_Infeas_Residuals(pdhg);

      if (bool_print) {
	// With reduced printing, the header is only needed for the
	// first iteration since only average iteration printing is
	// carried out
        if (full_print ||
	    iter_log_since_header == iter_log_between_header) {
	  PDHG_Print_Header(pdhg);
	  iter_log_since_header = 0;
	}
	int print_iter = timers->nIter == 0 || full_print;
	int print_iter_average = timers->nIter > 0 || full_print;
	if (print_iter) PDHG_Print_Iter(pdhg);
	if (print_iter_average) PDHG_Print_Iter_Average(pdhg);
	iter_log_since_header++;
      }

#if PDLP_DEBUG_LOG
      debugPdlpFeasOptLog(pdhg->debug_pdlp_log_file_,
			  pdhg->timers->nIter, 
			  pdhg->resobj->dPrimalObj,
			  pdhg->resobj->dDualObj,
			  pdhg->resobj->dRelObjGap,
			  pdhg->resobj->dPrimalFeas / (1.0 + pdhg->scaling->dNormRhs),
			  pdhg->resobj->dDualFeas / (1.0 + pdhg->scaling->dNormCost), "[L]");
      debugPdlpFeasOptLog(pdhg->debug_pdlp_log_file_,
			  pdhg->timers->nIter, 
			  pdhg->resobj->dPrimalObjAverage,
			  pdhg->resobj->dDualObjAverage,
			  pdhg->resobj->dRelObjGapAverage,
			  pdhg->resobj->dPrimalFeasAverage / (1.0 + pdhg->scaling->dNormRhs),
			  pdhg->resobj->dDualFeasAverage / (1.0 + pdhg->scaling->dNormCost), "[A]");
      debugPdlpIterHeaderLog(pdhg->debug_pdlp_log_file_);
#endif
      // Termination check printing is only done when printing is full
      int termination_print = bool_print && full_print;
      if (PDHG_Check_Termination(pdhg, termination_print)) {
        // cupdlp_printf("Optimal current solution.\n");
        resobj->termIterate = LAST_ITERATE;
        resobj->termCode = OPTIMAL;
        break;
      }

      if (PDHG_Check_Termination_Average(pdhg, termination_print)) {
	// cupdlp_printf("Optimal average solution.\n");
	
    cupdlp_int iter = pdhg->timers->nIter;
    CUPDLPvec *x = iterates->x[iter % 2];
    CUPDLPvec *y = iterates->y[iter % 2];
    CUPDLPvec *ax = iterates->ax[iter % 2];
    CUPDLPvec *aty = iterates->aty[iter % 2];

    CUPDLP_COPY_VEC(x->data, iterates->xAverage->data, cupdlp_float, problem->nCols);
    CUPDLP_COPY_VEC(y->data, iterates->yAverage->data, cupdlp_float, problem->nRows);
    CUPDLP_COPY_VEC(ax->data, iterates->axAverage->data, cupdlp_float, problem->nRows);
    CUPDLP_COPY_VEC(aty->data, iterates->atyAverage->data, cupdlp_float, problem->nCols);
    CUPDLP_COPY_VEC(resobj->dSlackPos, resobj->dSlackPosAverage, cupdlp_float, problem->nCols);
    CUPDLP_COPY_VEC(resobj->dSlackNeg, resobj->dSlackNegAverage, cupdlp_float, problem->nCols);

	resobj->termIterate = AVERAGE_ITERATE;
	resobj->termCode = OPTIMAL;
	break;
      }

      if (PDHG_Check_Infeasibility(pdhg, 0) == INFEASIBLE_OR_UNBOUNDED) {
        // cupdlp_printf("Infeasible or unbounded.\n");
        // if (resobj->primalCode == INFEASIBLE && resobj->dualCode == FEASIBLE)
        // {
        //   resobj->dualCode = UNBOUNDED;
        // } else if (resobj->primalCode == FEASIBLE &&
        //            resobj->dualCode == INFEASIBLE) {
        //   resobj->primalCode = UNBOUNDED;
        // }
        // resobj->termCode = resobj->primalCode;
        resobj->termCode = INFEASIBLE_OR_UNBOUNDED;
        break;
      }

      if (timers->dSolvingTime > settings->dTimeLim) {
        // cupdlp_printf("Time limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      if (timers->nIter >= (settings->nIterLim - 1)) {
        // cupdlp_printf("Iteration limit reached.\n");
        resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
        break;
      }

      PDHG_Restart_Iterate(pdhg);
    }

    CUPDLPvec *ax = iterates->ax[timers->nIter % 2];
    CUPDLPvec *aty = iterates->aty[timers->nIter % 2];

#if PDLP_DEBUG_LOG
    double debug_pdlp_data_ax_norm = 0.0;
    cupdlp_twoNorm(pdhg, problem->nRows, ax->data,
                    &debug_pdlp_data_ax_norm);
    pdhg->debug_pdlp_data_.ax_norm = debug_pdlp_data_ax_norm;
    double debug_pdlp_data_aty_norm = 0.0;
    cupdlp_twoNorm(pdhg, problem->nCols, aty->data, &debug_pdlp_data_aty_norm);
    pdhg->debug_pdlp_data_.aty_norm = debug_pdlp_data_aty_norm;
#endif
    // CUPDLP_CALL(PDHG_Update_Iterate(pdhg));
    if (PDHG_Update_Iterate(pdhg) == RETCODE_FAILED) {
      // cupdlp_printf("Time limit reached.\n");
      resobj->termCode = TIMELIMIT_OR_ITERLIMIT;
      break;
    }
  }

  // print at last
  if (pdhg->settings->nLogLevel>0 && timers->nIter > 0) {
    int full_print = pdhg->settings->nLogLevel >= 2;
    if (full_print) PDHG_Print_Header(pdhg);
    int print_iter = resobj->termIterate == LAST_ITERATE || full_print;
    int print_iter_average = resobj->termIterate == AVERAGE_ITERATE || full_print;
    if (print_iter) PDHG_Print_Iter(pdhg);
    if (print_iter_average) PDHG_Print_Iter_Average(pdhg);
  }

  if (pdhg->settings->nLogLevel>0) {
    cupdlp_printf("\n");
    cupdlp_printf("%-27s ", "Solving information:");

    switch (resobj->termCode) {
    case OPTIMAL:
      if (resobj->termIterate == LAST_ITERATE) {
        cupdlp_printf("Optimal current solution.\n");
      } else if (resobj->termIterate == AVERAGE_ITERATE) {
        cupdlp_printf("Optimal average solution.\n");
      }
      break;
    case TIMELIMIT_OR_ITERLIMIT:
      if (timers->dSolvingTime > settings->dTimeLim) {
        cupdlp_printf("Time limit reached.\n");
      } else if (timers->nIter >= (settings->nIterLim - 1)) {
        cupdlp_printf("Iteration limit reached.\n");
      }
      break;
    case INFEASIBLE_OR_UNBOUNDED:
      if (resobj->primalCode == INFEASIBLE && resobj->dualCode == FEASIBLE) {
        cupdlp_printf("Infeasible or unbounded: primal infeasible.");
      } else if (resobj->primalCode == FEASIBLE &&
                 resobj->dualCode == INFEASIBLE) {
        cupdlp_printf("Infeasible or unbounded: dual infeasible.");
      } else {
        cupdlp_printf(
		      "Infeasible or unbounded: both primal and dual infeasible.");
      }
      
      if (resobj->termInfeasIterate == LAST_ITERATE) {
        cupdlp_printf(" [L]\n");
      } else if (resobj->termInfeasIterate == AVERAGE_ITERATE) {
        cupdlp_printf(" [A]\n");
      }
      break;
    default:
      cupdlp_printf("Unexpected.\n");
      break;
    }

    if (resobj->termCode == OPTIMAL && resobj->termIterate == AVERAGE_ITERATE) {
      cupdlp_printf("%27s %+15.8e\n",
		    "Primal objective:", resobj->dPrimalObjAverage);
      cupdlp_printf("%27s %+15.8e\n", "Dual objective:", resobj->dDualObjAverage);
      cupdlp_printf("%27s %8.2e / %8.2e\n",
		    "Primal infeas (abs/rel):", resobj->dPrimalFeasAverage,
		    resobj->dPrimalFeasAverage / (1.0 + scaling->dNormRhs));
      cupdlp_printf("%27s %8.2e / %8.2e\n",
		    "Dual infeas (abs/rel):", resobj->dDualFeasAverage,
		    resobj->dDualFeasAverage / (1.0 + scaling->dNormCost));
      cupdlp_printf("%27s %8.2e / %8.2e\n",
		    "Duality gap (abs/rel):", fabs(resobj->dDualityGapAverage),
		    resobj->dRelObjGapAverage);
    } else {
      cupdlp_printf("%27s %+15.8e\n", "Primal objective:", resobj->dPrimalObj);
      cupdlp_printf("%27s %+15.8e\n", "Dual objective:", resobj->dDualObj);
      cupdlp_printf("%27s %8.2e / %8.2e\n",
		    "Primal infeas (abs/rel):", resobj->dPrimalFeas,
		    resobj->dPrimalFeas / (1.0 + scaling->dNormRhs));
      cupdlp_printf("%27s %8.2e / %8.2e\n",
		    "Dual infeas (abs/rel):", resobj->dDualFeas,
		    resobj->dDualFeas / (1.0 + scaling->dNormCost));
      cupdlp_printf("%27s %8.2e / %8.2e\n",
		    "Duality gap (abs/rel):", fabs(resobj->dDualityGap),
		    resobj->dRelObjGap);
    }
    cupdlp_printf("%27s %d\n", "Number of iterations:", timers->nIter);
    cupdlp_printf("\n");
  }

#if PDHG_USE_TIMERS
  if (pdhg->settings->nLogLevel>1) {
    cupdlp_printf("Timing information:\n");
    // cupdlp_printf("%21s %e in %d iterations\n", "Total solver time",
    //               timers->dSolvingTime, timers->nIter);
    cupdlp_printf(
		  "%21s %e in %d iterations\n", "Total solver time",
		  timers->dSolvingTime + timers->dScalingTime + timers->dPresolveTime,
		  timers->nIter);
    cupdlp_printf("%21s %e in %d iterations\n", "Solve time",
		  timers->dSolvingTime, timers->nIter);
    cupdlp_printf("%21s %e \n", "Iters per sec",
		  timers->nIter / timers->dSolvingTime);
    cupdlp_printf("%21s %e\n", "Scaling time", timers->dScalingTime);
    cupdlp_printf("%21s %e\n", "Presolve time", timers->dPresolveTime);
    cupdlp_printf("%21s %e in %d calls\n", "Ax", timers->dAxTime,
		  timers->nAxCalls);
    cupdlp_printf("%21s %e in %d calls\n", "Aty", timers->dAtyTime,
		  timers->nAtyCalls);
    cupdlp_printf("%21s %e in %d calls\n", "ComputeResiduals",
		  timers->dComputeResidualsTime, timers->nComputeResidualsCalls);
    cupdlp_printf("%21s %e in %d calls\n", "UpdateIterates",
		  timers->dUpdateIterateTime, timers->nUpdateIterateCalls);
  }
#endif

#ifndef CUPDLP_CPU
  if (pdhg->settings->nLogLevel>1) {
    cupdlp_printf("\n");
    cupdlp_printf("GPU Timing information:\n");
    cupdlp_printf("%21s %e\n", "CudaPrepare", timers->CudaPrepareTime);
    cupdlp_printf("%21s %e\n", "Alloc&CopyMatToDevice",
		  timers->AllocMem_CopyMatToDeviceTime);
    cupdlp_printf("%21s %e\n", "CopyVecToDevice", timers->CopyVecToDeviceTime);
    cupdlp_printf("%21s %e\n", "DeviceMatVecProd", timers->DeviceMatVecProdTime);
    cupdlp_printf("%21s %e\n\n", "CopyVecToHost", timers->CopyVecToHostTime);
  }
#endif

exit_cleanup:
  //fclose(pdhg->debug_pdlp_log_file_);
  return retcode;
}

cupdlp_retcode PDHG_PreSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                              cupdlp_int *constraint_new_idx,
                              cupdlp_int *constraint_type,
                              cupdlp_float *col_value, cupdlp_float *col_dual,
                              cupdlp_float *row_value, cupdlp_float *row_dual,
                              cupdlp_int *value_valid, cupdlp_int *dual_valid) {
  cupdlp_retcode retcode = RETCODE_OK;
  if (!*value_valid || !*dual_valid) {
    // Nothing to presolve
    return retcode;
  }

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPscaling *scaling = pdhg->scaling;
  CUPDLPresobj *resobj = pdhg->resobj;
  cupdlp_float sense = problem->sense_origin;

  CUPDLPvec *x = iterates->x[0];
  CUPDLPvec *y = iterates->y[0];

  assert(constraint_new_idx);
  assert(constraint_type);

  // Have to build cuPDLP-C solution in CPU row and column buffers, so
  // that only copy operations are required to put the cuPDLP-C
  // solution into GPU buffers
  cupdlp_float *col_buffer = NULL;
  cupdlp_float *row_buffer = NULL;
  CUPDLP_INIT_DOUBLE(col_buffer, problem->nCols);
  CUPDLP_INIT_DOUBLE(row_buffer, problem->nRows);


  CUPDLP_ZERO_VEC(x->data, cupdlp_float, problem->nCols);
  CUPDLP_ZERO_VEC(y->data, cupdlp_float, problem->nRows);


  int iCol = 0;
  for (; iCol < nCols_origin; iCol++)
    col_buffer[iCol] = col_value[iCol];

  for (int iRow = 0; iRow < problem->nRows; iRow++) {
    const double mu = constraint_type[iRow] == 1 ? -1 : 1;
    row_buffer[constraint_new_idx[iRow]] = sense * mu * row_dual[iRow];
    if (constraint_type[iRow] == 3) 
      col_buffer[iCol++] = row_value[iRow];
  }

  // Copy into GPU buffers
  CUPDLP_COPY_VEC(x->data, col_buffer, cupdlp_float, problem->nCols);
  CUPDLP_COPY_VEC(y->data, row_buffer, cupdlp_float, problem->nRows);

  // Scale
  if (scaling->ifScaled) {
    cupdlp_edot(x->data, pdhg->colScale, problem->nCols);
    cupdlp_edot(y->data, pdhg->rowScale, problem->nRows);
  }
exit_cleanup:
  // free buffer
  CUPDLP_FREE(col_buffer);
  CUPDLP_FREE(row_buffer);
  return retcode;
}

cupdlp_retcode PDHG_PostSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                              cupdlp_int *constraint_new_idx,
                              cupdlp_int *constraint_type,
                              cupdlp_float *col_value, cupdlp_float *col_dual,
                              cupdlp_float *row_value, cupdlp_float *row_dual,
                              cupdlp_int *value_valid, cupdlp_int *dual_valid) {
  cupdlp_retcode retcode = RETCODE_OK;

  CUPDLPproblem *problem = pdhg->problem;
  CUPDLPiterates *iterates = pdhg->iterates;
  CUPDLPscaling *scaling = pdhg->scaling;
  CUPDLPresobj *resobj = pdhg->resobj;
  cupdlp_float sense = problem->sense_origin;

  // flag
  cupdlp_int col_value_flag = 0;
  cupdlp_int col_dual_flag = 0;
  cupdlp_int row_value_flag = 0;
  cupdlp_int row_dual_flag = 0;

  // allocate buffer
  cupdlp_float *col_buffer = NULL;
  cupdlp_float *row_buffer = NULL;
  cupdlp_float *col_buffer2 = NULL;
  // no need for row_buffer2
  // cupdlp_float *row_buffer2 = NULL;
  CUPDLP_INIT_DOUBLE(col_buffer, problem->nCols);
  CUPDLP_INIT_DOUBLE(row_buffer, problem->nRows);
  CUPDLP_INIT_DOUBLE(col_buffer2, problem->nCols);
  // CUPDLP_INIT_DOUBLE(row_buffer2, problem->nRows);

  
  cupdlp_int iter = pdhg->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];
  CUPDLPvec *ax = iterates->ax[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];

  /*
    for (int iCol = 0; iCol < problem->nCols; iCol++)
    printf("PDHG_PostSolve: Col %d primal value = %9.3g\n", iCol, x->data[iCol]);
    for (int iRow = 0; iRow < problem->nRows; iRow++)
    printf("PDHG_PostSolve: Row %d   dual value = %9.3g\n", iRow, y->data[iRow]);
  */
  
  // unscale
  if (scaling->ifScaled) {
    cupdlp_ediv(x->data, pdhg->colScale, problem->nCols);
    cupdlp_ediv(y->data, pdhg->rowScale, problem->nRows);
    cupdlp_edot(resobj->dSlackPos, pdhg->colScale, problem->nCols);
    cupdlp_edot(resobj->dSlackNeg, pdhg->colScale, problem->nCols);
    cupdlp_edot(ax->data, pdhg->rowScale, problem->nRows);
    cupdlp_edot(aty->data, pdhg->colScale, problem->nCols);
  }

  // col value: extract x from (x, z)
  if (col_value) {
    CUPDLP_COPY_VEC(col_value, x->data, cupdlp_float, nCols_origin);

    col_value_flag = 1;
  }

  // row value
  if (row_value) {
    if (constraint_new_idx) {
      CUPDLP_COPY_VEC(row_buffer, ax->data, cupdlp_float,
                      problem->nRows);

      // un-permute row value
      for (int i = 0; i < problem->nRows; i++) {
        row_value[i] = row_buffer[constraint_new_idx[i]];
      }
    } else {
      CUPDLP_COPY_VEC(row_value, ax->data, cupdlp_float,
                      problem->nRows);
    }

    if (constraint_type) {
      CUPDLP_COPY_VEC(col_buffer, x->data, cupdlp_float,
                      problem->nCols);

      // EQ = 0, LEQ = 1, GEQ = 2, BOUND = 3
      for (int i = 0, j = 0; i < problem->nRows; i++) {
        if (constraint_type[i] == 1) {  // LEQ: multiply -1
          row_value[i] = -row_value[i];
        } else if (constraint_type[i] == 3) {  // BOUND: get Ax from Ax - z
          row_value[i] = row_value[i] + col_buffer[nCols_origin + j];
          j++;
        }
      }
    }

    row_value_flag = 1;
  }

  // col duals of l <= x <= u
  if (col_dual) {
    CUPDLP_COPY_VEC(col_buffer, resobj->dSlackPos, cupdlp_float, nCols_origin);
    CUPDLP_COPY_VEC(col_buffer2, resobj->dSlackNeg, cupdlp_float, nCols_origin);

    for (int i = 0; i < nCols_origin; i++) {
      col_dual[i] = col_buffer[i] - col_buffer2[i];
    }

    ScaleVector(sense, col_dual, nCols_origin);

    col_dual_flag = 1;
  }

  // row dual: recover y
  if (row_dual) {
    if (constraint_new_idx) {
      CUPDLP_COPY_VEC(row_buffer, y->data, cupdlp_float,
                      problem->nRows);
      // un-permute row dual
      for (int i = 0; i < problem->nRows; i++) {
        row_dual[i] = row_buffer[constraint_new_idx[i]];
      }
    } else {
      CUPDLP_COPY_VEC(row_dual, y->data, cupdlp_float,
                      problem->nRows);
    }

    ScaleVector(sense, row_dual, problem->nRows);

    if (constraint_type) {
      // EQ = 0, LEQ = 1, GEQ = 2, BOUND = 3
      for (int i = 0; i < problem->nRows; i++) {
        if (constraint_type[i] == 1) {  // LEQ: multiply -1
          row_dual[i] = -row_dual[i];
        }
      }
    }

    row_dual_flag = 1;
  }

  // valid
  if (value_valid) {
    *value_valid = col_value_flag && row_value_flag;
  }

  if (dual_valid) {
    *dual_valid = col_dual_flag && row_dual_flag;
  }

exit_cleanup:
  // free buffer
  CUPDLP_FREE(col_buffer);
  CUPDLP_FREE(row_buffer);
  CUPDLP_FREE(col_buffer2);
  // CUPDLP_FREE(row_buffer2);

  return retcode;
}

cupdlp_retcode LP_SolvePDHG(
    CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam,
    cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp,
    cupdlp_int nCols_origin, cupdlp_float *col_value, cupdlp_float *col_dual,
    cupdlp_float *row_value, cupdlp_float *row_dual, cupdlp_int *value_valid,
    cupdlp_int *dual_valid, cupdlp_bool ifSaveSol, char *fp_sol,
    cupdlp_int *constraint_new_idx, cupdlp_int *constraint_type,
    cupdlp_int *model_status, cupdlp_int* num_iter) {
  cupdlp_retcode retcode = RETCODE_OK;

 // Set the parameters first - which is silent
  CUPDLP_CALL(PDHG_SetUserParam(pdhg, ifChangeIntParam, intParam,
                                ifChangeFloatParam, floatParam));

 // Call PDHG_PrintHugeCUPDHG() if logging level (set in
 // PDHG_SetUserParam) is verbose
  if (pdhg->settings->nLogLevel > 1) 
    PDHG_PrintHugeCUPDHG();

#if !defined(CUPDLP_CPU)
  if (pdhg->settings->nLogLevel > 0)
    print_cuda_info(pdhg->cusparsehandle);
#endif

 CUPDLP_CALL(PDHG_PreSolve(pdhg, nCols_origin, constraint_new_idx,
 	    constraint_type, col_value, col_dual, row_value,
 	    row_dual, value_valid, dual_valid));

  cupdlp_int has_variables = (*value_valid + *dual_valid) != 0;

 if (pdhg->settings->nLogLevel>0 && has_variables)
   cupdlp_printf("Hot starting with given column primal values and row dual values\n");

  CUPDLP_CALL(PDHG_Solve(&has_variables, pdhg));

  *model_status = (cupdlp_int)pdhg->resobj->termCode;
  *num_iter = (cupdlp_int)pdhg->timers->nIter;

  CUPDLP_CALL(PDHG_PostSolve(pdhg, nCols_origin, constraint_new_idx,
                             constraint_type, col_value, col_dual, row_value,
                             row_dual, value_valid, dual_valid));

  if (fp)
    writeJson(fp, pdhg);

  if (ifSaveSol && fp_sol) {
    if (strcmp(fp, fp_sol) != 0) {
      writeSol(fp_sol, nCols_origin, pdhg->problem->nRows, col_value, col_dual,
               row_value, row_dual);
    } else {
      if (pdhg->settings->nLogLevel>0) 
	cupdlp_printf("Warning: fp and fp_sol are the same, stop saving solution.\n");
    }
  }

exit_cleanup:
  // Move this to CupdlpWrapper so pdhg->debug_pdlp_log_file_ can
  // still be used
  //
  //  PDHG_Destroy(&pdhg);
  return retcode;
}
