/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* HiGHS link code */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/* GAMS API */
#include "gevmcc.h"
#include "gmomcc.h"

typedef struct optRec* optHandle_t;

/* HiGHS API */
#include "Highs.h"
#include "io/FilereaderLp.h"
#include "io/FilereaderMps.h"
#include "io/LoadOptions.h" /* for loadOptionsFromFile */

#if defined(_WIN32)
#if !defined(STDCALL)
#define STDCALL __stdcall
#endif
#if !defined(DllExport)
#define DllExport __declspec(dllexport)
#endif
#else
#if !defined(STDCALL)
#define STDCALL
#endif
#if !defined(DllExport)
#define DllExport
#endif
#endif

struct gamshighs_s {
  gmoHandle_t gmo;
  gevHandle_t gev;
  int debug;

  Highs* highs;
  HighsLp* lp;
  HighsOptions* options;
};
typedef struct gamshighs_s gamshighs_t;

static void gevprint(int level, const char* msg, void* msgcb_data) {
  gevHandle_t gev = (gevHandle_t)msgcb_data;
  gevLogPChar(gev, msg);
}

static void gevlog(HighsMessageType type, const char* msg, void* msgcb_data) {
  gevHandle_t gev = (gevHandle_t)msgcb_data;
  if (type == HighsMessageType::INFO)
    gevLogPChar(gev, msg);
  else
    gevLogStatPChar(gev, msg);
}

static enum gmoVarEquBasisStatus translateBasisStatus(HighsBasisStatus status) {
  switch (status) {
    case HighsBasisStatus::BASIC:
      return gmoBstat_Basic;
    case HighsBasisStatus::LOWER:
      return gmoBstat_Lower;
    case HighsBasisStatus::NONBASIC:
    case HighsBasisStatus::SUPER:
    case HighsBasisStatus::ZERO:
      return gmoBstat_Super;
    case HighsBasisStatus::UPPER:
      return gmoBstat_Upper;
  }
  // this should never happen
  return gmoBstat_Super;
}

static HighsBasisStatus translateBasisStatus(enum gmoVarEquBasisStatus status) {
  switch (status) {
    case gmoBstat_Basic:
      return HighsBasisStatus::BASIC;
    case gmoBstat_Lower:
      return HighsBasisStatus::LOWER;
    case gmoBstat_Super:
      return HighsBasisStatus::SUPER;
    case gmoBstat_Upper:
      return HighsBasisStatus::UPPER;
  }
  // this should never happen
  return HighsBasisStatus::SUPER;
}

static int setupOptions(gamshighs_t* gh) {
  assert(gh != NULL);
  assert(gh->options == NULL);

  gh->options = new HighsOptions;

  gh->options->time_limit = gevGetDblOpt(gh->gev, gevResLim);
  if (gevGetIntOpt(gh->gev, gevIterLim) != ITERLIM_INFINITY)
    gh->options->simplex_iteration_limit = gevGetIntOpt(gh->gev, gevIterLim);

  if (gevGetIntOpt(gh->gev, gevUseCutOff))
    gh->options->dual_objective_value_upper_bound =
        gevGetDblOpt(gh->gev, gevCutOff);

  if (gmoOptFile(gh->gmo) > 0) {
    char optfilename[GMS_SSSIZE];
    gmoNameOptFile(gh->gmo, optfilename);
    gh->options->options_file = optfilename;
    if (!loadOptionsFromFile(*gh->options)) return 1;
  }

  gh->options->printmsgcb = gevprint;
  gh->options->logmsgcb = gevlog;
  gh->options->msgcb_data = (void*)gh->gev;

  return 0;
}

static int setupProblem(gamshighs_t* gh) {
  int numCol;
  int numRow;
  int numNz;
  int i;
  int rc = 1;
  HighsSolution sol;

  assert(gh != NULL);
  assert(gh->options != NULL);
  assert(gh->highs == NULL);
  assert(gh->lp == NULL);

  gh->highs = new Highs(*gh->options);

  numCol = gmoN(gh->gmo);
  numRow = gmoM(gh->gmo);
  numNz = gmoNZ(gh->gmo);

  gh->lp = new HighsLp();

  gh->lp->numRow_ = numRow;
  gh->lp->numCol_ = numCol;
  gh->lp->nnz_ = numNz;

  /* columns */
  gh->lp->colUpper_.resize(numCol);
  gh->lp->colLower_.resize(numCol);
  gmoGetVarLower(gh->gmo, &gh->lp->colLower_[0]);
  gmoGetVarUpper(gh->gmo, &gh->lp->colUpper_[0]);

  /* objective */
  gh->lp->colCost_.resize(numCol);
  gmoGetObjVector(gh->gmo, &gh->lp->colCost_[0], NULL);
  if (gmoSense(gh->gmo) == gmoObj_Min)
    gh->lp->sense_ = OBJSENSE_MINIMIZE;
  else
    gh->lp->sense_ = OBJSENSE_MAXIMIZE;
  gh->lp->offset_ = gmoObjConst(gh->gmo);

  /* row left- and right-hand-side */
  gh->lp->rowLower_.resize(numRow);
  gh->lp->rowUpper_.resize(numRow);
  for (i = 0; i < numRow; ++i) {
    switch (gmoGetEquTypeOne(gh->gmo, i)) {
      case gmoequ_E:
        gh->lp->rowLower_[i] = gh->lp->rowUpper_[i] = gmoGetRhsOne(gh->gmo, i);
        break;

      case gmoequ_G:
        gh->lp->rowLower_[i] = gmoGetRhsOne(gh->gmo, i);
        gh->lp->rowUpper_[i] = HIGHS_CONST_INF;
        break;

      case gmoequ_L:
        gh->lp->rowLower_[i] = -HIGHS_CONST_INF;
        gh->lp->rowUpper_[i] = gmoGetRhsOne(gh->gmo, i);
        break;

      case gmoequ_N:
      case gmoequ_X:
      case gmoequ_C:
      case gmoequ_B:
        /* these should not occur */
        goto TERMINATE;
    }
  }

  /* coefficients matrix */
  gh->lp->Astart_.resize(numCol + 1);
  gh->lp->Aindex_.resize(numNz);
  gh->lp->Avalue_.resize(numNz);
  gmoGetMatrixCol(gh->gmo, &gh->lp->Astart_[0], &gh->lp->Aindex_[0],
                  &gh->lp->Avalue_[0], NULL);

  gh->highs->passModel(*gh->lp);

  // FilereaderLp().writeModelToFile("highs.lp", *gh->lp);
  // FilereaderMps().writeModelToFile("highs.mps", *gh->lp);

  // pass initial solution
  sol.col_value.resize(numCol);
  sol.col_dual.resize(numCol);
  sol.row_value.resize(numRow);
  sol.row_dual.resize(numRow);
  gmoGetVarL(gh->gmo, &sol.col_value[0]);
  gmoGetVarM(gh->gmo, &sol.col_dual[0]);
  gmoGetEquL(gh->gmo, &sol.row_value[0]);
  gmoGetEquM(gh->gmo, &sol.row_dual[0]);
  gh->highs->setSolution(sol);

  if (gmoHaveBasis(gh->gmo)) {
    HighsBasis basis;
    basis.col_status.resize(numCol);
    basis.row_status.resize(numRow);
    int nbasic = 0;

    for (int i = 0; i < numCol; ++i) {
      basis.col_status[i] = translateBasisStatus(
          (enum gmoVarEquBasisStatus)gmoGetVarStatOne(gh->gmo, i));
      if (basis.col_status[i] == HighsBasisStatus::BASIC) ++nbasic;
    }

    for (int i = 0; i < numRow; ++i) {
      basis.row_status[i] = translateBasisStatus(
          (enum gmoVarEquBasisStatus)gmoGetEquStatOne(gh->gmo, i));
      if (basis.row_status[i] == HighsBasisStatus::BASIC) ++nbasic;
    }

    basis.valid_ = nbasic == numRow;

    gh->highs->setBasis(basis);
  }

  rc = 0;
TERMINATE:

  return rc;
}

static int processSolve(gamshighs_t* gh) {
  assert(gh != NULL);
  assert(gh->highs != NULL);
  assert(gh->lp != NULL);

  gmoHandle_t gmo = gh->gmo;
  Highs* highs = gh->highs;

  gmoSetHeadnTail(gmo, gmoHresused, gevTimeDiffStart(gh->gev));
  gmoSetHeadnTail(gmo, gmoHiterused, highs->getIterationCount());

  // figure out model and solution status and whether we should have a solution
  // to be written
  bool writesol = false;
  switch (highs->getModelStatus()) {
    case HighsModelStatus::NOTSET:
    case HighsModelStatus::LOAD_ERROR:
    case HighsModelStatus::MODEL_ERROR:
    case HighsModelStatus::MODEL_EMPTY:
    case HighsModelStatus::PRESOLVE_ERROR:
    case HighsModelStatus::SOLVE_ERROR:
    case HighsModelStatus::POSTSOLVE_ERROR:
      gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
      gmoSolveStatSet(gmo, gmoSolveStat_SolverErr);
      break;

    case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      // TODO is there a solution to write and is it feasible?
      // gmoModelStatSet(gmo, havesol ? gmoModelStat_InfeasibleIntermed :
      // gmoModelStat_NoSolutionReturned);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      gmoSolveStatSet(gmo, gmoSolveStat_Solver);
      break;

    case HighsModelStatus::PRIMAL_UNBOUNDED:
      // TODO is there a (feasible) solution to write?
      // gmoModelStatSet(gmo, havesol ? gmoModelStat_Unbounded :
      // gmoModelStat_UnboundedNoSolution);
      gmoModelStatSet(gmo, gmoModelStat_UnboundedNoSolution);
      gmoSolveStatSet(gmo, gmoSolveStat_Normal);
      break;

    case HighsModelStatus::PRIMAL_INFEASIBLE:
      // TODO is there an infeasible solution to write?
      // gmoModelStatSet(gmo, havesol ? gmoModelStat_InfeasibleGlobal :
      // gmoModelStat_InfeasibleNoSolution);
      gmoModelStatSet(gmo, gmoModelStat_InfeasibleNoSolution);
      gmoSolveStatSet(gmo, gmoSolveStat_Normal);
      break;

      // case HighsModelStatus::PRIMAL_FEASIBLE:
      /*
      gmoModelStatSet(gmo, gmoModelStat_Feasible);
      gmoSolveStatSet(gmo, gmoSolveStat_Solver);
      writesol = true;
      break;
      */

      // case HighsModelStatus::DUAL_FEASIBLE:
      /*
      gmoModelStatSet(gmo, gmoModelStat_InfeasibleIntermed);
      gmoSolveStatSet(gmo, gmoSolveStat_Solver);
      writesol = true;
      break;
      */

    case HighsModelStatus::OPTIMAL:
      gmoModelStatSet(gmo, gmoModelStat_OptimalGlobal);
      gmoSolveStatSet(gmo, gmoSolveStat_Normal);
      writesol = true;
      break;

    case HighsModelStatus::REACHED_TIME_LIMIT:
      // TODO is there an (feasible) solution to write?
      // gmoModelStatSet(gmo, havesol ? gmoModelStat_InfeasibleIntermed :
      // gmoModelStat_NoSolutionReturned);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      gmoSolveStatSet(gmo, gmoSolveStat_Resource);
      break;

    case HighsModelStatus::REACHED_ITERATION_LIMIT:
      // TODO is there an (feasible) solution to write?
      // gmoModelStatSet(gmo, havesol ? gmoModelStat_InfeasibleIntermed :
      // gmoModelStat_NoSolutionReturned);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      gmoSolveStatSet(gmo, gmoSolveStat_Iteration);
      break;
  }

  if (writesol) {
    const HighsSolution& sol = highs->getSolution();
    assert((int)sol.col_value.size() == gmoN(gmo));
    assert((int)sol.col_dual.size() == gmoN(gmo));
    assert((int)sol.row_value.size() == gmoM(gmo));
    assert((int)sol.row_dual.size() == gmoM(gmo));

    const HighsBasis& basis = highs->getBasis();
    assert(!basis.valid_ || (int)basis.col_status.size() == gmoN(gmo));
    assert(!basis.valid_ || (int)basis.row_status.size() == gmoM(gmo));

    for (int i = 0; i < gmoN(gmo); ++i) {
      gmoVarEquBasisStatus basisstat;
      if (basis.valid_)
        basisstat = translateBasisStatus(basis.col_status[i]);
      else
        basisstat = gmoBstat_Super;

      // TODO change when we can process infeasible or unbounded solutions
      gmoVarEquStatus stat = gmoCstat_OK;

      gmoSetSolutionVarRec(gmo, i, sol.col_value[i], sol.col_dual[i], basisstat,
                           stat);
    }

    for (int i = 0; i < gmoM(gmo); ++i) {
      gmoVarEquBasisStatus basisstat;
      if (basis.valid_)
        basisstat = translateBasisStatus(basis.row_status[i]);
      else
        basisstat = gmoBstat_Super;

      // TODO change when we can process infeasible or unbounded solutions
      gmoVarEquStatus stat = gmoCstat_OK;

      gmoSetSolutionEquRec(gmo, i, sol.row_value[i], sol.row_dual[i], basisstat,
                           stat);
    }

    gmoCompleteObjective(gmo, highs->getObjectiveValue());
  }

  return 0;
}

extern "C" {

void his_Initialize(void) {
  gmoInitMutexes();
  gevInitMutexes();
}

void his_Finalize(void) {
  gmoFiniMutexes();
  gevFiniMutexes();
}

DllExport void STDCALL hisXCreate(void** Cptr) {
  assert(Cptr != NULL);

  *Cptr = calloc(1, sizeof(gamshighs_t));
}

DllExport int STDCALL hiscreate(void** Cptr, char* msgBuf, int msgBufLen) {
  assert(Cptr != NULL);
  assert(msgBufLen > 0);
  assert(msgBuf != NULL);

  *Cptr = calloc(1, sizeof(gamshighs_t));

  msgBuf[0] = 0;

  return 1;
}

DllExport void STDCALL hisXFree(void** Cptr) {
  assert(Cptr != NULL);
  assert(*Cptr != NULL);

  free(*Cptr);
  *Cptr = NULL;

  gmoLibraryUnload();
  gevLibraryUnload();
}

DllExport int STDCALL hisfree(void** Cptr) {
  hisXFree(Cptr);

  return 1;
}

/* comp returns the compatibility mode:
           0: client is too old for the DLL, no compatibility
           1: client version and DLL version are the same, full compatibility
           2: client is older than DLL, but defined as compatible, backward
   compatibility 3: client is newer than DLL, forward compatibility
           FIXME: for now, we just claim full compatibility
 */
DllExport int STDCALL C__hisXAPIVersion(int api, char* Msg, int* comp) {
  *comp = 1;
  return 1;
}

DllExport int STDCALL D__hisXAPIVersion(int api, char* Msg, int* comp) {
  *comp = 1;
  return 1;
}

DllExport int STDCALL C__hisXCheck(const char* funcn, int ClNrArg, int Clsign[],
                                   char* Msg) {
  return 1;
}

DllExport int STDCALL D__hisXCheck(const char* funcn, int ClNrArg, int Clsign[],
                                   char* Msg) {
  return 1;
}

DllExport int STDCALL C__hisReadyAPI(void* Cptr, gmoHandle_t Gptr,
                                     optHandle_t Optr) {
  gamshighs_t* gh;

  assert(Cptr != NULL);
  assert(Gptr != NULL);
  assert(Optr == NULL);

  char msg[256];
  if (!gmoGetReady(msg, sizeof(msg))) return 1;
  if (!gevGetReady(msg, sizeof(msg))) return 1;

  gh = (gamshighs_t*)Cptr;
  gh->gmo = Gptr;
  gh->gev = (gevHandle_t)gmoEnvironment(gh->gmo);

  return 0;
}

#define XQUOTE(x) QUOTE(x)
#define QUOTE(x) #x

DllExport int STDCALL C__hisCallSolver(void* Cptr) {
  int rc = 1;
  gamshighs_t* gh;
  HighsStatus status;

  gh = (gamshighs_t*)Cptr;
  assert(gh->gmo != NULL);
  assert(gh->gev != NULL);

  gevLogStatPChar(gh->gev, "HiGHS " XQUOTE(HIGHS_VERSION_MAJOR) "." XQUOTE(HIGHS_VERSION_MINOR) "." XQUOTE(HIGHS_VERSION_PATCH) " [date: " HIGHS_COMPILATION_DATE ", git hash: " HIGHS_GITHASH "]\n");
  gevLogStatPChar(gh->gev,
                  "Copyright (c) 2019 ERGO-Code under MIT license terms.\n");

  gmoModelStatSet(gh->gmo, gmoModelStat_NoSolutionReturned);
  gmoSolveStatSet(gh->gmo, gmoSolveStat_SystemErr);

  /* get the problem into a normal form */
  gmoObjStyleSet(gh->gmo, gmoObjType_Fun);
  gmoObjReformSet(gh->gmo, 1);
  gmoIndexBaseSet(gh->gmo, 0);
  gmoSetNRowPerm(gh->gmo); /* hide =N= rows */
  gmoMinfSet(gh->gmo, -HIGHS_CONST_INF);
  gmoPinfSet(gh->gmo, HIGHS_CONST_INF);

  if (setupOptions(gh)) goto TERMINATE;

  if (setupProblem(gh)) goto TERMINATE;

  gevTimeSetStart(gh->gev);

  /* solve the problem */
  status = gh->highs->run();
  if (status != HighsStatus::OK) goto TERMINATE;

  /* pass solution, status, etc back to GMO */
  if (processSolve(gh)) goto TERMINATE;

  rc = 0;
TERMINATE:

  delete gh->lp;
  gh->lp = NULL;

  delete gh->highs;
  gh->highs = NULL;

  delete gh->options;
  gh->options = NULL;

  return rc;
}

DllExport int STDCALL C__hisHaveModifyProblem(void* Cptr) { return 1; }

DllExport int STDCALL C__hisModifyProblem(void* Cptr) {
  gamshighs_t* gh = (gamshighs_t*)Cptr;
  assert(gh != NULL);

  /* set the GMO styles, in case someone changed this */
  gmoObjStyleSet(gh->gmo, gmoObjType_Fun);
  gmoObjReformSet(gh->gmo, 1);
  gmoIndexBaseSet(gh->gmo, 0);
  gmoSetNRowPerm(gh->gmo); /* hide =N= rows */
  gmoMinfSet(gh->gmo, -HIGHS_CONST_INF);
  gmoPinfSet(gh->gmo, HIGHS_CONST_INF);

  int maxsize = std::max(gmoN(gh->gmo), gmoM(gh->gmo));

  int jacnz;
  gmoGetJacUpdate(gh->gmo, NULL, NULL, NULL, &jacnz);
  if (jacnz + 1 > maxsize) maxsize = jacnz + 1;

  int* colidx = new int[maxsize];
  int* rowidx = new int[maxsize];
  double* array1 = new double[maxsize];
  double* array2 = new double[maxsize];

  Highs* highs = gh->highs;
  assert(highs != NULL);

  // update objective coefficients
  int nz;
  int nlnz;
  gmoGetObjSparse(gh->gmo, colidx, array1, NULL, &nz, &nlnz);
  assert(nlnz == gmoObjNZ(gh->gmo));
  highs->changeColsCost(nz, colidx, array1);

  // TODO update objective offset

  // update variable bounds
  gmoGetVarLower(gh->gmo, array1);
  gmoGetVarUpper(gh->gmo, array2);
  highs->changeColsBounds(0, gmoN(gh->gmo), array1, array2);

  // update constraint sides
  for (int i = 0; i < gmoM(gh->gmo); ++i) {
    double rhs = gmoGetRhsOne(gh->gmo, i);
    rowidx[i] = 1;
    switch (gmoGetEquTypeOne(gh->gmo, i)) {
      case gmoequ_E:
        array1[i] = rhs;
        array2[i] = rhs;
        break;

      case gmoequ_G:
        array1[i] = rhs;
        array2[i] = HIGHS_CONST_INF;
        break;

      case gmoequ_L:
        array1[i] = -HIGHS_CONST_INF;
        array2[i] = rhs;
        break;

      case gmoequ_N:
      case gmoequ_X:
      case gmoequ_C:
      case gmoequ_B:
        /* these should not occur */
        rowidx[i] = 0;
        break;
    }
  }
  highs->changeRowsBounds(rowidx, array1, array2);

  // update constraint matrix
  gmoGetJacUpdate(gh->gmo, rowidx, colidx, array1, &jacnz);
  for (int i = 0; i < nz; ++i)
    highs->changeCoeff(rowidx[i], colidx[i], array1[i]);

  delete[] array2;
  delete[] array1;
  delete[] rowidx;
  delete[] colidx;

  return 0;
}

void his_Initialize(void);
void his_Finalize(void);

}  // extern "C"
