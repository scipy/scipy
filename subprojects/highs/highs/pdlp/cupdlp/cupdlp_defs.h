#ifndef CUPDLP_H_GUARD
#define CUPDLP_H_GUARD

// The below already defined in CMake.
// #define CUPDLP_CPU
// #define CUPDLP_DEBUG (1)
#define CUPDLP_TIMER

#include "HConfig.h"

#ifdef CUPDLP_GPU
#include "cuda/cupdlp_cuda_kernels.cuh"
#include "cuda/cupdlp_cudalinalg.cuh"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <string.h>

#include "glbopts.h"

#define PDHG_USE_TIMERS (1)
#define USE_MY_BLAS (1)
#define USE_KERNELS (1)

#define PDHG_STEPSIZE_REDUCTION_EXP \
  (0.3)  // Parameters in PDLP adaptive linesearch
#define PDHG_STEPSIZE_GROWTH_EXP (0.6)
#define PDHG_USE_TIMERS (1)
#define PDHG_DISPLAY_TERMINATION_CHECK (1)

#define PDHG_PROJECT_INITIAL_X (0)
#define SHOW_DEFINE(x) printf("%s=%s\n", #x, STR(x))

#define CUPDLP_DEBUG_INTERVAL (40)
#define CUPDLP_RELEASE_INTERVAL (40)
#define CUPDLP_DUMP_ITERATES_STATS (1)
#define CUPDLP_DUMP_LINESEARCH_STATS (1)
#define CUPDLP_INEXACT_EPS (1e-4)

#define PDLP_DEBUG_LOG (0)

typedef struct CUPDLP_CUDA_DENSE_VEC CUPDLPvec;
typedef struct CUPDLP_DENSE_MATRIX CUPDLPdense;
typedef struct CUPDLP_CSR_MATRIX CUPDLPcsr;
typedef struct CUPDLP_CSC_MATRIX CUPDLPcsc;
typedef struct CUPDLP_DATA CUPDLPdata;
typedef struct CUPDLP_SETTINGS CUPDLPsettings;
typedef struct CUPDLP_PROBLEM CUPDLPproblem;
typedef struct CUPDLP_RES_OBJ CUPDLPresobj;
typedef struct CUPDLP_ITERATES CUPDLPiterates;
typedef struct CUPDLP_STEPSIZE CUPDLPstepsize;
typedef struct CUPDLP_SCALING CUPDLPscaling;
typedef struct CUPDLP_TIMERS CUPDLPtimers;
typedef struct CUPDLP_WORK CUPDLPwork;
typedef cupdlp_int cupdlp_retcode;

typedef enum {
  OPTIMAL = 0,
  INFEASIBLE,
  UNBOUNDED,
  INFEASIBLE_OR_UNBOUNDED,
  TIMELIMIT_OR_ITERLIMIT,
  FEASIBLE,
} termination_code;

typedef enum {
  LAST_ITERATE = 0,
  AVERAGE_ITERATE,
} termination_iterate;

typedef enum {
  PDHG_FIXED_LINESEARCH = 0,
  PDHG_MALITSKY_POCK_LINESEARCH,
  PDHG_ADAPTIVE_LINESEARCH
} pdhg_linesearch;

typedef enum {
  PDHG_WITHOUT_RESTART = 0,
  PDHG_GPU_RESTART,
  PDHG_CPU_RESTART,
} pdhg_restart;

typedef enum {
  CPU = 0,
  SINGLE_GPU,
  MULTI_GPU,
} CUPDLP_DEVICE;

typedef enum {
  DENSE = 0,
  CSR,
  CSC,
  CSR_CSC,
} CUPDLP_MATRIX_FORMAT;

typedef enum {
  N_ITER_LIM = 0,
  IF_SCALING,
  I_SCALING_METHOD,
  E_LINE_SEARCH_METHOD,
  E_RESTART_METHOD,
  IF_RUIZ_SCALING,
  IF_L2_SCALING,
  IF_PC_SCALING,
  N_LOG_LEVEL,
  N_LOG_INTERVAL,
  IF_PRESOLVE,
  I_INF_NORM_ABS_LOCAL_TERMINATION,
  N_INT_USER_PARAM
} CUPDLP_INT_USER_PARAM_INDEX;
  //#define N_INT_USER_PARAM 12
typedef enum {
  D_SCALING_LIMIT = 0,
  D_PRIMAL_TOL,
  D_DUAL_TOL,
  D_GAP_TOL,
  D_FEAS_TOL,
  D_TIME_LIM,
  N_FLOAT_USER_PARAM
} CUPDLP_FLOAT_USER_PARAM_INDEX;
  //#define N_FLOAT_USER_PARAM 6

// used in sparse matrix-dense vector multiplication
struct CUPDLP_CUDA_DENSE_VEC {
  cupdlp_int len;
  cupdlp_float *data;
#ifndef CUPDLP_CPU
  cusparseDnVecDescr_t cuda_vec;
#endif
};

struct CUPDLP_DENSE_MATRIX {
  cupdlp_int nRows;
  cupdlp_int nCols;
  cupdlp_float *data;
};

struct CUPDLP_CSR_MATRIX {
  cupdlp_int nRows;
  cupdlp_int nCols;
  cupdlp_int nMatElem;
  cupdlp_int *rowMatBeg;
  cupdlp_int *rowMatIdx;
  cupdlp_float *rowMatElem;
#ifndef CUPDLP_CPU
  // Pointers to GPU vectors
  cusparseSpMatDescr_t cuda_csr;
#endif
};

struct CUPDLP_CSC_MATRIX {
  cupdlp_int nRows;
  cupdlp_int nCols;
  cupdlp_int nMatElem;
  cupdlp_int *colMatBeg;
  cupdlp_int *colMatIdx;
  cupdlp_float *colMatElem;

  // Used to avoid implementing NormInf on cuda
  cupdlp_float MatElemNormInf;
#ifndef CUPDLP_CPU
  // Pointers to GPU vectors
  cusparseSpMatDescr_t cuda_csc;
#endif
};

struct CUPDLP_DATA {
  cupdlp_int nRows;
  cupdlp_int nCols;
  CUPDLP_MATRIX_FORMAT matrix_format;
  CUPDLPdense *dense_matrix;
  CUPDLPcsr *csr_matrix;
  CUPDLPcsc *csc_matrix;
  CUPDLP_DEVICE device;
};

struct CUPDLP_SETTINGS {
  // scaling
  cupdlp_int ifScaling;
  cupdlp_int iScalingMethod;
  cupdlp_float dScalingLimit;

  // termination criteria
  cupdlp_float dPrimalTol;
  cupdlp_float dDualTol;
  cupdlp_float dGapTol;

  // max iter and time
  cupdlp_int nIterLim;
  cupdlp_float dTimeLim;

  // Logging
  cupdlp_int nLogLevel;
  cupdlp_int nLogInterval;

  // restart
  pdhg_restart eRestartMethod;
};

// some elements are duplicated from CUPDLP_DATA
struct CUPDLP_PROBLEM {
  /* Copy of LP problem with permuted columns. */
  CUPDLPdata *data;
  // cupdlp_int nMatElem;
  // cupdlp_int *colMatBeg;
  // cupdlp_int *colMatIdx;
  // cupdlp_float *colMatElem;
  // cupdlp_int *rowMatBeg;
  // cupdlp_int *rowMatIdx;
  // cupdlp_float *rowMatElem;
  cupdlp_float *lower;
  cupdlp_float *upper;
  cupdlp_float *cost;  // cost for minimization
  cupdlp_float *rhs;
  cupdlp_float dMaxCost;
  cupdlp_float dMaxRhs;
  cupdlp_float dMaxRowBound;
  cupdlp_int nRows;
  cupdlp_int nCols;
  cupdlp_int nEqs;
  cupdlp_float *hasLower;
  cupdlp_float *hasUpper;
  cupdlp_float
      offset;  // true objVal = c'x * sig + offset, sig = 1 (min) or -1 (max)
  cupdlp_float sense_origin;  // sig = 1 (min) or -1 (max)
};

struct CUPDLP_RES_OBJ {
  /* residuals and objectives */
  cupdlp_float dFeasTol;
  cupdlp_float dPrimalObj;
  cupdlp_float dDualObj;
  cupdlp_float dDualityGap;
  cupdlp_float dComplementarity;
  cupdlp_float dPrimalFeas;
  cupdlp_float dDualFeas;
  cupdlp_float dRelObjGap;
  cupdlp_float *primalResidual;
  cupdlp_float *dualResidual;
  cupdlp_float *dSlackPos;
  cupdlp_float *dSlackNeg;
  cupdlp_float *dSlackPosAverage;
  cupdlp_float *dSlackNegAverage;
  cupdlp_float *dLowerFiltered;
  cupdlp_float *dUpperFiltered;

  /* for infeasibility detection */
  termination_code primalCode;
  termination_code dualCode;
  termination_iterate termInfeasIterate;

  cupdlp_float dPrimalInfeasObj;
  cupdlp_float dDualInfeasObj;
  cupdlp_float dPrimalInfeasRes;
  cupdlp_float dDualInfeasRes;

  cupdlp_float dPrimalInfeasObjAverage;
  cupdlp_float dDualInfeasObjAverage;
  cupdlp_float dPrimalInfeasResAverage;
  cupdlp_float dDualInfeasResAverage;

  // buffers
  cupdlp_float *primalInfeasRay;     // x / norm(x)
  cupdlp_float *primalInfeasConstr;  // [Ax, min(Gx, 0)]
  cupdlp_float *primalInfeasBound;   // primal bound violation
  cupdlp_float *dualInfeasRay;       // y / norm(y, lbd)
  cupdlp_float *dualInfeasLbRay;     // lbd^+ / norm(y, lbd)
  cupdlp_float *dualInfeasUbRay;     // lbd^- / norm(y, lbd)
  cupdlp_float *dualInfeasConstr;    // ATy1 + GTy2 + lambda
  // cupdlp_float *dualInfeasBound;     // dual bound violation

  cupdlp_float dPrimalObjAverage;
  cupdlp_float dDualObjAverage;
  cupdlp_float dDualityGapAverage;
  cupdlp_float dComplementarityAverage;
  cupdlp_float dPrimalFeasAverage;
  cupdlp_float dDualFeasAverage;
  cupdlp_float dRelObjGapAverage;
  cupdlp_float *primalResidualAverage;
  cupdlp_float *dualResidualAverage;

  cupdlp_float dPrimalFeasLastRestart;
  cupdlp_float dDualFeasLastRestart;
  cupdlp_float dDualityGapLastRestart;

  cupdlp_float dPrimalFeasLastCandidate;
  cupdlp_float dDualFeasLastCandidate;
  cupdlp_float dDualityGapLastCandidate;

  termination_code termCode;
  termination_iterate termIterate;
};

struct CUPDLP_ITERATES {
  /* iterates */
  cupdlp_int nRows;
  cupdlp_int nCols;
  //  todo, CPU VERSION, check
  //        cupdlp_float *x;
  //        cupdlp_float *y;
  //        cupdlp_float *xUpdate;
  //        cupdlp_float *yUpdate;
  //
  //        cupdlp_int iLastRestartIter;
  //        cupdlp_float dLastRestartDualityGap;
  //        cupdlp_float dLastRestartBeta;
  //        cupdlp_float *xSum;
  //        cupdlp_float *ySum;
  //        cupdlp_float *xAverage;
  //        cupdlp_float *yAverage;
  //        cupdlp_float *xLastRestart;
  //        cupdlp_float *yLastRestart;
  //
  //        cupdlp_float *ax;
  //        cupdlp_float *axUpdate;
  //        cupdlp_float *axAverage;
  //        cupdlp_float *aty;
  //        cupdlp_float *atyUpdate;
  //        cupdlp_float *atyAverage;

  cupdlp_int iLastRestartIter;
  cupdlp_float dLastRestartDualityGap;
  cupdlp_float dLastRestartBeta;
  cupdlp_float *xSum;
  cupdlp_float *ySum;

  cupdlp_float *xLastRestart;
  cupdlp_float *yLastRestart;

  CUPDLPvec *x[2]; // in iteration k, x^k stored in x[k % 2], x^{k+1} created in x[k + 1 % 2]
  CUPDLPvec *y[2];
  CUPDLPvec *ax[2];
  CUPDLPvec *aty[2];
  CUPDLPvec *xAverage, *yAverage, *axAverage, *atyAverage;

};

struct CUPDLP_STEPSIZE {
  /* stepsize */
  pdhg_linesearch eLineSearchMethod;  // 0 = FixedStep
  cupdlp_float dPrimalStep;
  cupdlp_float dDualStep;
  cupdlp_float dSumPrimalStep;
  cupdlp_float dSumDualStep;
  // Stepsize ratio,
  //  \beta = dBeta = dDualStep / dPrimalStep,
  //    in the paper, primal weight is the \omega:
  //    \omega = \sqrt\beta
  cupdlp_float dBeta;
  cupdlp_float dTheta;  // Used in Malitsky-Pock stepsize
  cupdlp_int nStepSizeIter;
};

struct CUPDLP_SCALING {
  /* scaling */
  cupdlp_int ifScaled;
  cupdlp_float *rowScale;
  cupdlp_float *colScale;

  /*new scaling*/
  cupdlp_int ifRuizScaling;
  cupdlp_int ifL2Scaling;
  cupdlp_int ifPcScaling;
  cupdlp_int RuizTimes;
  cupdlp_float RuizNorm;
  cupdlp_float PcAlpha;

  /* original 2 norm */
  cupdlp_float dNormCost;
  cupdlp_float dNormRhs;
};

struct CUPDLP_TIMERS {
  /* timers */
  cupdlp_int nIter;
  cupdlp_float dSolvingTime;
  cupdlp_float dSolvingBeg;
  cupdlp_float dScalingTime;
  cupdlp_float dPresolveTime;
#if PDHG_USE_TIMERS
  cupdlp_float dAtyTime;
  cupdlp_float dAxTime;
  cupdlp_float dComputeResidualsTime;
  cupdlp_float dUpdateIterateTime;
  cupdlp_int nAtyCalls;
  cupdlp_int nAxCalls;
  cupdlp_int nComputeResidualsCalls;
  cupdlp_int nUpdateIterateCalls;
#endif
#ifndef CUPDLP_CPU
  // GPU timers
  cupdlp_float AllocMem_CopyMatToDeviceTime;
  cupdlp_float CopyVecToDeviceTime;
  cupdlp_float DeviceMatVecProdTime;
  cupdlp_float CopyVecToHostTime;
  cupdlp_float FreeDeviceMemTime;
  cupdlp_float CudaPrepareTime;
#endif
};

struct DebugPdlpData {
  double ax_norm;
  double aty_norm;
  double ax_average_norm;
  double aty_average_norm;
  double x_average_norm;
};
    
struct CUPDLP_WORK {
  CUPDLPproblem *problem;
  CUPDLPsettings *settings;
  CUPDLPresobj *resobj;
  CUPDLPiterates *iterates;
  CUPDLPstepsize *stepsize;
  CUPDLPscaling *scaling;
  CUPDLPtimers *timers;
  // cupdlp_float *buffer;
  CUPDLPvec *buffer;
  cupdlp_float *buffer2;
  cupdlp_float *buffer3;

  cupdlp_float *rowScale;
  cupdlp_float *colScale;
#ifndef CUPDLP_CPU
  // CUDAmv *MV;
  cusparseHandle_t cusparsehandle;
  void *dBuffer_csc_ATy;
  void *dBuffer_csr_Ax;

  // cusparseDnVecDescr_t vecbuffer;
  cublasHandle_t cublashandle;
#endif

#if PDLP_DEBUG_LOG 
  FILE* debug_pdlp_log_file_;
  struct DebugPdlpData debug_pdlp_data_;
#endif
};

#ifdef __cplusplus
}
#endif
#endif
