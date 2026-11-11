//
// Created by LJS on 23-11-30.
// Same as the JULIA CPU version
//

#include "cupdlp_scaling.h"

#include "cupdlp_linalg.h"
// #include "cupdlp_scaling.h"
#include "cupdlp_utils.h"

// This version disable dScalingTarget, which is the target of scaled matrix
// elements
// cupdlp_retcode scale_problem(CUPDLPwork *w, cupdlp_float *col_scaling,
// cupdlp_float *row_scaling)

static cupdlp_retcode scale_problem(CUPDLPcsc *csc, cupdlp_float *cost,
                             cupdlp_float *lower, cupdlp_float *upper,
                             cupdlp_float *rhs, cupdlp_float *col_scaling,
                             cupdlp_float *row_scaling) {
  cupdlp_retcode retcode = RETCODE_OK;
  cupdlp_int nRows = csc->nRows;
  cupdlp_int nCols = csc->nCols;

  cupdlp_cdiv(cost, col_scaling, nCols);
  cupdlp_cdot(lower, col_scaling, nCols);
  cupdlp_cdot(upper, col_scaling, nCols);
  cupdlp_cdiv(rhs, row_scaling, nRows);

  // row scaling
  for (int i = 0; i < csc->colMatBeg[nCols]; i++) {
    csc->colMatElem[i] /= row_scaling[csc->colMatIdx[i]];
  }
  // col scaling
  for (int i = 0; i < nCols; i++) {
    for (int j = csc->colMatBeg[i]; j < csc->colMatBeg[i + 1]; j++) {
      csc->colMatElem[j] /= col_scaling[i];
    }
  }

  // w->scaling->ifScaled = 1;

exit_cleanup:
  return retcode;
}

static cupdlp_retcode cupdlp_ruiz_scaling(CUPDLPcsc *csc, cupdlp_float *cost,
                                   cupdlp_float *lower,
                                   cupdlp_float *upper, cupdlp_float *rhs,
                                   CUPDLPscaling *scaling)
// cupdlp_retcode cupdlp_ruiz_scaling(CUPDLPwork *work, cupdlp_int max_iter,
// cupdlp_float norm)
{
  cupdlp_retcode retcode = RETCODE_OK;

  cupdlp_int nRows = csc->nRows;
  cupdlp_int nCols = csc->nCols;

  cupdlp_float *current_col_scaling = NULL;  // for variable
  cupdlp_float *current_row_scaling = NULL;  // for constraint
  CUPDLP_INIT_ZERO_DOUBLE(current_col_scaling, nCols);
  CUPDLP_INIT_ZERO_DOUBLE(current_row_scaling, nRows);

  for (cupdlp_int i = 0; i < scaling->RuizTimes; i++) {
    cupdlp_zero(current_col_scaling, cupdlp_float, nCols);
    cupdlp_zero(current_row_scaling, cupdlp_float, nRows);

    if (csc != NULL) {
      for (int j = 0; j < nCols; j++) {
        if (csc->colMatBeg[j] == csc->colMatBeg[j + 1]) {
          current_col_scaling[j] = 0;
        } else {
          current_col_scaling[j] = SQRTF(GenNorm(
              &csc->colMatElem[csc->colMatBeg[j]],
              csc->colMatBeg[j + 1] - csc->colMatBeg[j], scaling->RuizNorm));
        }
      }
    }
    for (int j = 0; j < nCols; j++) {
      if (current_col_scaling[j] == 0.0) {
        current_col_scaling[j] = 1.0;
      }
    }

    if (scaling->RuizNorm == INFINITY) {
      if (nRows > 0 && csc != NULL) {
        // inf norm of rows of csc
        for (int j = 0; j < csc->colMatBeg[nCols]; j++) {
          if (current_row_scaling[csc->colMatIdx[j]] <
              ABS(csc->colMatElem[j])) {
            current_row_scaling[csc->colMatIdx[j]] = ABS(csc->colMatElem[j]);
          }
        }
        for (int j = 0; j < nRows; j++) {
          if (current_row_scaling[j] == 0.0) {
            current_row_scaling[j] = 1.0;
          } else {
            current_row_scaling[j] = SQRTF(current_row_scaling[j]);
          }
        }
      }
    } else {
      cupdlp_printf("Currently only support infinity norm for Ruiz scaling\n");
      exit(1);
    }

    // apply scaling
    // scale_problem(work, current_col_scaling, current_row_scaling);
    scale_problem(csc, cost, lower, upper, rhs, current_col_scaling,
                  current_row_scaling);

    // update scaling
    cupdlp_cdot(scaling->colScale, current_col_scaling, nCols);
    cupdlp_cdot(scaling->rowScale, current_row_scaling, nRows);
  }
exit_cleanup:
  cupdlp_free(current_col_scaling);
  cupdlp_free(current_row_scaling);
  return retcode;
}

static cupdlp_retcode cupdlp_l2norm_scaling(CUPDLPcsc *csc, cupdlp_float *cost,
                                     cupdlp_float *lower,
                                     cupdlp_float *upper,
                                     cupdlp_float *rhs,
                                     CUPDLPscaling *scaling)
// cupdlp_retcode cupdlp_l2norm_scaling(CUPDLPwork *work)
{
  cupdlp_retcode retcode = RETCODE_OK;
  cupdlp_int nRows = csc->nRows;
  cupdlp_int nCols = csc->nCols;

  cupdlp_float *current_col_scaling = NULL;  // for variable
  cupdlp_float *current_row_scaling = NULL;  // for constraint
  CUPDLP_INIT_ZERO_DOUBLE(current_col_scaling, nCols);
  CUPDLP_INIT_ZERO_DOUBLE(current_row_scaling, nRows);

  if (nRows > 0) {
    for (int j = 0; j < nCols; j++) {
      if (csc->colMatBeg[j] == csc->colMatBeg[j + 1]) {
        current_col_scaling[j] = 1.0;
      } else {
        current_col_scaling[j] =
            SQRTF(GenNorm(&csc->colMatElem[csc->colMatBeg[j]],
                          csc->colMatBeg[j + 1] - csc->colMatBeg[j], 2.0));
      }
    }

    for (int i = 0; i < csc->colMatBeg[nCols]; i++) {
      current_row_scaling[csc->colMatIdx[i]] += pow(csc->colMatElem[i], 2.0);
    }
    for (int i = 0; i < nRows; i++) {
      current_row_scaling[i] = SQRTF(SQRTF(current_row_scaling[i]));
      if (current_row_scaling[i] == 0.0) {
        current_row_scaling[i] = 1.0;
      }
    }
  }
  // apply scaling
  // scale_problem(work, current_col_scaling, current_row_scaling);
  scale_problem(csc, cost, lower, upper, rhs, current_col_scaling,
                current_row_scaling);

  // update scaling
  cupdlp_cdot(scaling->colScale, current_col_scaling, nCols);
  cupdlp_cdot(scaling->rowScale, current_row_scaling, nRows);

exit_cleanup:
  cupdlp_free(current_col_scaling);
  cupdlp_free(current_row_scaling);
  return retcode;
}

static cupdlp_retcode cupdlp_pc_scaling(CUPDLPcsc *csc, cupdlp_float *cost,
                                 cupdlp_float *lower, cupdlp_float *upper,
                                 cupdlp_float *rhs, CUPDLPscaling *scaling)
// cupdlp_retcode cupdlp_pc_scaling(CUPDLPwork *work, cupdlp_float alpha)
{
  cupdlp_retcode retcode = RETCODE_OK;
  cupdlp_int nRows = csc->nRows;
  cupdlp_int nCols = csc->nCols;
  cupdlp_float alpha = scaling->PcAlpha;

  cupdlp_float *current_col_scaling = NULL;  // for variable
  cupdlp_float *current_row_scaling = NULL;  // for constraint
  CUPDLP_INIT_ZERO_DOUBLE(current_col_scaling, nCols);
  CUPDLP_INIT_ZERO_DOUBLE(current_row_scaling, nRows);

  if (alpha > 2.0 || alpha < 0.0) {
    cupdlp_printf("alpha should be in [0, 2]\n");
    exit(1);
  }

  if (nRows > 0) {
    for (int i = 0; i < nCols; i++) {
      for (int j = csc->colMatBeg[i]; j < csc->colMatBeg[i + 1]; j++) {
        current_col_scaling[i] += POWF(ABS(csc->colMatElem[j]), alpha);
      }
      current_col_scaling[i] = SQRTF(POWF(current_col_scaling[i], 1.0 / alpha));
      if (current_col_scaling[i] == 0.0) {
        current_col_scaling[i] = 1.0;
      }
    }

    for (int i = 0; i < csc->colMatBeg[nCols]; i++) {
      current_row_scaling[csc->colMatIdx[i]] +=
          POWF(ABS(csc->colMatElem[i]), 2.0 - alpha);
    }
    for (int i = 0; i < nRows; i++) {
      current_row_scaling[i] =
          SQRTF(POWF(current_row_scaling[i], 1.0 / (2.0 - alpha)));
      if (current_row_scaling[i] == 0.0) {
        current_row_scaling[i] = 1.0;
      }
    }
  }

  // apply scaling
  // scale_problem(work, current_col_scaling, current_row_scaling);
  scale_problem(csc, cost, lower, upper, rhs, current_col_scaling,
                current_row_scaling);

  // update scaling
  cupdlp_cdot(scaling->colScale, current_col_scaling, nCols);
  cupdlp_cdot(scaling->rowScale, current_row_scaling, nRows);

exit_cleanup:
  cupdlp_free(current_col_scaling);
  cupdlp_free(current_row_scaling);
  return retcode;
}

cupdlp_retcode PDHG_Scale_Data(cupdlp_int log_level, CUPDLPcsc *csc,
                               cupdlp_int ifScaling, CUPDLPscaling *scaling,
                               cupdlp_float *cost, cupdlp_float *lower,
                               cupdlp_float *upper, cupdlp_float *rhs) {
  cupdlp_retcode retcode = RETCODE_OK;
  // scaling->dObjScale = 1.0;

#if CUPDLP_DEBUG
  //------------------- for debug ------------------
  cupdlp_float dMinElem = OUR_DBL_MAX;
  cupdlp_float dMaxElem = 0.0;
  cupdlp_float dAvgElem = 0.0;
  cupdlp_int nRows = csc->nRows;
  cupdlp_int nCols = csc->nCols;

  for (cupdlp_int iMatElem = csc->colMatBeg[0];
       iMatElem < csc->colMatBeg[nCols]; iMatElem++) {
    cupdlp_float dAbsElem = fabs(csc->colMatElem[iMatElem]);
    if (dAbsElem != 0.0) {
      dMaxElem = fmax(dMaxElem, dAbsElem);
      dMinElem = fmin(dMinElem, dAbsElem);
    }
    dAvgElem += dAbsElem;
  }
  dAvgElem /= csc->colMatBeg[nCols];

  if (log_level) {
    cupdlp_printf("Problem before rescaling:\n");
    cupdlp_printf(
        "Absolute value of nonzero constraint matrix elements: largest=%f, "
        "smallest=%f, avg=%f\n",
        dMaxElem, dMinElem, dAvgElem);
  }

  // calculate the three statistics of objective vector
  dMinElem = OUR_DBL_MAX;
  dMaxElem = 0.0;
  dAvgElem = 0.0;
  for (cupdlp_int iCol = 0; iCol < nCols; iCol++) {
    cupdlp_float dAbsElem = fabs(cost[iCol]);
    if (dAbsElem != 0.0) {
      dMaxElem = fmax(dMaxElem, dAbsElem);
      dMinElem = fmin(dMinElem, dAbsElem);
    }
    dAvgElem += dAbsElem;
  }
  dAvgElem /= nCols;

  // if (log_level) {}
  cupdlp_printf(
      "Absolute value of objective vector elements: largest=%f, smallest=%f, "
      "avg=%f\n",
      dMaxElem, dMinElem, dAvgElem);
  // calculate the three statistics of rhs vector
  dMinElem = OUR_DBL_MAX;
  dMaxElem = 0.0;
  dAvgElem = 0.0;
  for (cupdlp_int iRow = 0; iRow < nRows; iRow++) {
    cupdlp_float dAbsElem = fabs(rhs[iRow]);
    if (dAbsElem != 0.0) {
      dMaxElem = fmax(dMaxElem, dAbsElem);
      dMinElem = fmin(dMinElem, dAbsElem);
    }
    dAvgElem += dAbsElem;
  }
  dAvgElem /= nRows;
  cupdlp_printf(
      "Absolute value of rhs vector elements: largest=%f, smallest=%f, "
      "avg=%f\n",
      dMaxElem, dMinElem, dAvgElem);
//------------------- for debug ------------------
#endif

  if (ifScaling) {
    if (log_level) cupdlp_printf("--------------------------------------------------\n");
    if (log_level) cupdlp_printf("running scaling\n");

    if (scaling->ifRuizScaling) {
      if (log_level) cupdlp_printf("- use Ruiz scaling\n");
      CUPDLP_CALL(cupdlp_ruiz_scaling(csc, cost, lower, upper, rhs, scaling))
      scaling->ifScaled = 1;
    }
    if (scaling->ifL2Scaling) {
      if (log_level) cupdlp_printf("- use L2 scaling\n");
      CUPDLP_CALL(cupdlp_l2norm_scaling(csc, cost, lower, upper, rhs, scaling))
      scaling->ifScaled = 1;
    }
    if (scaling->ifPcScaling) {
      if (log_level) cupdlp_printf("- use PC scaling\n");
      CUPDLP_CALL(cupdlp_pc_scaling(csc, cost, lower, upper, rhs, scaling))
      scaling->ifScaled = 1;
    }

    if (log_level) cupdlp_printf("--------------------------------------------------\n");
  }

  /* make sure the csr matrix is also scaled*/
  // csc2csr(data->csr_matrix, csc);
#if CUPDLP_DEBUG
  //------------------- for debug ------------------
  dMinElem = OUR_DBL_MAX;
  dMaxElem = 0.0;
  dAvgElem = 0.0;
  for (cupdlp_int iMatElem = csc->colMatBeg[0];
       iMatElem < csc->colMatBeg[nCols]; iMatElem++) {
    cupdlp_float dAbsElem = fabs(csc->colMatElem[iMatElem]);
    if (dAbsElem != 0.0) {
      dMaxElem = fmax(dMaxElem, dAbsElem);
      dMinElem = fmin(dMinElem, dAbsElem);
    }
    dAvgElem += dAbsElem;
  }
  dAvgElem /= csc->colMatBeg[nCols];

  if (log_level) cupdlp_printf("Problem after rescaling:\n");
  if (log_level) cupdlp_printf(
      "Absolute value of nonzero constraint matrix elements: largest=%f, "
      "smallest=%f, avg=%f\n",
      dMaxElem, dMinElem, dAvgElem);

  // calculate the three statistics of objective vector
  dMinElem = OUR_DBL_MAX;
  dMaxElem = 0.0;
  dAvgElem = 0.0;
  for (cupdlp_int iCol = 0; iCol < nCols; iCol++) {
    cupdlp_float dAbsElem = fabs(cost[iCol]);
    if (dAbsElem != 0.0) {
      dMaxElem = fmax(dMaxElem, dAbsElem);
      dMinElem = fmin(dMinElem, dAbsElem);
    }
    dAvgElem += dAbsElem;
  }
  dAvgElem /= nCols;
  if (log_level) cupdlp_printf(
      "Absolute value of objective vector elements: largest=%f, smallest=%f, "
      "avg=%f\n",
      dMaxElem, dMinElem, dAvgElem);
  // calculate the three statistics of rhs vector
  dMinElem = OUR_DBL_MAX;
  dMaxElem = 0.0;
  dAvgElem = 0.0;
  for (cupdlp_int iRow = 0; iRow < nRows; iRow++) {
    cupdlp_float dAbsElem = fabs(rhs[iRow]);
    if (dAbsElem != 0.0) {
      dMaxElem = fmax(dMaxElem, dAbsElem);
      dMinElem = fmin(dMinElem, dAbsElem);
    }
    dAvgElem += dAbsElem;
  }
  dAvgElem /= nRows;
  if (log_level) cupdlp_printf(
      "Absolute value of rhs vector elements: largest=%f, smallest=%f, "
      "avg=%f\n",
      dMaxElem, dMinElem, dAvgElem);
//------------------- for debug ------------------
#endif

exit_cleanup:

  return retcode;
}

cupdlp_retcode Init_Scaling(cupdlp_int log_level, CUPDLPscaling* scaling,
                            cupdlp_int ncols, cupdlp_int nrows,
                            cupdlp_float* cost, cupdlp_float* rhs) {
  cupdlp_retcode retcode = RETCODE_OK;

  scaling->rowScale = NULL;
  scaling->colScale = NULL;

  scaling->ifRuizScaling = 1;
  scaling->ifL2Scaling = 0;
  scaling->ifPcScaling = 1;

  // todo, read these paras
  scaling->RuizTimes = 10;
  scaling->RuizNorm = INFINITY;
  scaling->PcAlpha = 1.0;
  CUPDLP_INIT_DOUBLE(scaling->colScale, ncols);
  CUPDLP_INIT_DOUBLE(scaling->rowScale, nrows);

  for (cupdlp_int iCol = 0; iCol < ncols; iCol++) scaling->colScale[iCol] = 1.0;
  for (cupdlp_int iRow = 0; iRow < nrows; iRow++) scaling->rowScale[iRow] = 1.0;

  scaling->dNormCost = twoNorm(cost, ncols);
  scaling->dNormRhs = twoNorm(rhs, nrows);
  if (log_level)
    cupdlp_printf("Using cost norm = %9.3g and RHS norm = %9.3g\n",
                  scaling->dNormCost, scaling->dNormRhs);

exit_cleanup:
  return retcode;
}
