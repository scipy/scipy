/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file interfaces/OsiHiGHSInterface.cpp
 * @brief Osi/HiGHS interface implementation
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "OsiHiGHSSolverInterface.hpp"

#include <cmath>

#include "CoinWarmStartBasis.hpp"
#include "Highs.h"
#include "HighsLp.h"
#include "HighsOptions.h"
#include "HighsStatus.h"
#include "io/FilereaderMps.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"

static void printtomessagehandler(int level, const char* msg,
                                  void* msgcb_data) {
  assert(msgcb_data != NULL);

  CoinMessageHandler* handler = (CoinMessageHandler*)msgcb_data;

  int len = strlen(msg);
  if (len > 0 && msg[len - 1] == '\n') {
    const_cast<char*>(msg)[len - 1] = '\0';
    handler->message(0, "HiGHS", msg, ' ') << CoinMessageEol;
    const_cast<char*>(msg)[len - 1] = '\n';
  } else
    handler->message(0, "HiGHS", msg, ' ');
}

static void logtomessagehandler(HighsMessageType type, const char* msg,
                                void* msgcb_data) {
  assert(msgcb_data != NULL);

  CoinMessageHandler* handler = (CoinMessageHandler*)msgcb_data;

  // we know log message end with a newline, replace by coin-eol
  int len = strlen(msg);
  assert(len > 0);
  assert(msg[len - 1] == '\n');
  const_cast<char*>(msg)[len - 1] = '\0';

  handler->message(0, "HiGHS", msg, ' ') << CoinMessageEol;

  const_cast<char*>(msg)[len - 1] = '\n';
}

OsiHiGHSSolverInterface::OsiHiGHSSolverInterface()
    //  : status(HighsStatus::Init) {
    : status(HighsStatus::OK) {
  HighsSetMessageCallback(printtomessagehandler, logtomessagehandler,
                          (void*)handler_);

  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::OsiHiGHSSolverInterface()\n");
  this->highs = new Highs();
  this->dummy_solution = new HighsSolution;

  // because HiGHS calls HiGHSSetIO with the Options, which overwrites
  // the previous setting
  this->highs->options_.printmsgcb = printtomessagehandler;
  this->highs->options_.logmsgcb = logtomessagehandler;
  this->highs->options_.msgcb_data = (void*)handler_;

  setStrParam(OsiSolverName, "HiGHS");
}

OsiHiGHSSolverInterface::OsiHiGHSSolverInterface(
    const OsiHiGHSSolverInterface& original)
    : OsiSolverInterface(original),
      //      status(HighsStatus::Init)
      status(HighsStatus::OK) {
  HighsSetMessageCallback(printtomessagehandler, logtomessagehandler,
                          (void*)handler_);

  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::OsiHiGHSSolverInterface()\n");
  this->highs = new Highs();
  this->dummy_solution = new HighsSolution;

  // because HiGHS calls HiGHSSetIO with the Options, whichoverwrites the
  // previous setting
  this->highs->options_.printmsgcb = printtomessagehandler;
  this->highs->options_.logmsgcb = logtomessagehandler;
  this->highs->options_.msgcb_data = (void*)handler_;

  this->highs->passModel(original.highs->getLp());
  setStrParam(OsiSolverName, "HiGHS");
}

OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface() {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface()\n");

  HighsSetMessageCallback(NULL, NULL, NULL);

  delete this->highs;

  if (this->rowRange != NULL) {
    delete[] this->rowRange;
  }

  if (this->rhs != NULL) {
    delete[] this->rhs;
  }

  if (this->rowSense != NULL) {
    delete[] this->rowSense;
  }

  if (this->matrixByCol != NULL) {
    delete this->matrixByCol;
  }
}

OsiSolverInterface* OsiHiGHSSolverInterface::clone(bool copyData) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::clone()\n");
  if (!copyData) {
    OsiHiGHSSolverInterface* cln = new OsiHiGHSSolverInterface();
    return cln;

  } else {
    OsiHiGHSSolverInterface* cln = new OsiHiGHSSolverInterface(*this);
    cln->objOffset = this->objOffset;
    return cln;
  }
}

bool OsiHiGHSSolverInterface::setIntParam(OsiIntParam key, int value) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setIntParam()\n");
  switch (key) {
    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
      this->highs->options_.simplex_iteration_limit = value;
      return true;
    case OsiNameDiscipline:
      // TODO
      return false;
    case OsiLastIntParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::setDblParam(OsiDblParam key, double value) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setDblParam()\n");
  switch (key) {
    case OsiDualObjectiveLimit:
      this->highs->options_.dual_objective_value_upper_bound = value;
      return true;
    case OsiPrimalObjectiveLimit:
      return false;
    case OsiDualTolerance:
      this->highs->options_.dual_feasibility_tolerance = value;
      return true;
    case OsiPrimalTolerance:
      this->highs->options_.primal_feasibility_tolerance = value;
      return true;
    case OsiObjOffset:
      this->objOffset = value;
      return true;
    case OsiLastDblParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::setStrParam(OsiStrParam key,
                                          const std::string& value) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setStrParam(%d, %s)\n",
                    key, value.c_str());
  switch (key) {
    case OsiProbName:
      return OsiSolverInterface::setStrParam(key, value);
    case OsiSolverName:
      return OsiSolverInterface::setStrParam(key, value);
    case OsiLastStrParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::getIntParam(OsiIntParam key, int& value) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getIntParam()\n");
  switch (key) {
    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
      value = this->highs->options_.simplex_iteration_limit;
      return true;
    case OsiNameDiscipline:
      // TODO
      return false;
    case OsiLastIntParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::getDblParam(OsiDblParam key,
                                          double& value) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getDblParam()\n");
  switch (key) {
    case OsiDualObjectiveLimit:
      value = this->highs->options_.dual_objective_value_upper_bound;
      return true;
    case OsiPrimalObjectiveLimit:
      return false;
    case OsiDualTolerance:
      value = this->highs->options_.dual_feasibility_tolerance;
      return true;
    case OsiPrimalTolerance:
      value = this->highs->options_.primal_feasibility_tolerance;
      return true;
    case OsiObjOffset:
      value = this->objOffset;
      return true;
    case OsiLastDblParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::getStrParam(OsiStrParam key,
                                          std::string& value) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getStrParam(%d, %s)\n",
                    key, value.c_str());
  switch (key) {
    case OsiProbName:
      return OsiSolverInterface::getStrParam(key, value);
    case OsiSolverName:
      return OsiSolverInterface::getStrParam(key, value);
    case OsiLastStrParam:
    default:
      return false;
  }
}

void OsiHiGHSSolverInterface::initialSolve() {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::initialSolve()\n");
  this->status = this->highs->run();
}

bool OsiHiGHSSolverInterface::isAbandoned() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::isAbandoned()\n");
  //  return this->status == HighsStatus::NumericalDifficulties;
  return false;
}

bool OsiHiGHSSolverInterface::isProvenOptimal() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::isProvenOptimal()\n");
  //  return (this->status == HighsStatus::Optimal) ||
  //         (this->status == HighsStatus::OK);
  return false;
}

bool OsiHiGHSSolverInterface::isProvenPrimalInfeasible() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isProvenPrimalInfeasible()\n");
  //  return this->status == HighsStatus::Infeasible;
  return false;
}

bool OsiHiGHSSolverInterface::isProvenDualInfeasible() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isProvenDualInfeasible()\n");
  //  return this->status == HighsStatus::Unbounded;
  return false;
}

bool OsiHiGHSSolverInterface::isPrimalObjectiveLimitReached() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isPrimalObjectiveLimitReached()\n");
  return false;
}

bool OsiHiGHSSolverInterface::isDualObjectiveLimitReached() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isDualObjectiveLimitReached()\n");
  //  return this->status == HighsStatus::ReachedDualObjectiveUpperBound;
  return false;
}

bool OsiHiGHSSolverInterface::isIterationLimitReached() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isIterationLimitReached()\n");
  //  return this->status == HighsStatus::ReachedIterationLimit;
  return false;
}

int OsiHiGHSSolverInterface::getNumCols() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getNumCols()\n");
  return this->highs->lp_.numCol_;
}

int OsiHiGHSSolverInterface::getNumRows() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getNumRows()\n");
  return this->highs->lp_.numRow_;
}

int OsiHiGHSSolverInterface::getNumElements() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getNumElements()\n");
  return this->highs->lp_.nnz_;
}

const double* OsiHiGHSSolverInterface::getColLower() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getColLower()\n");
  return &(this->highs->lp_.colLower_[0]);
}

const double* OsiHiGHSSolverInterface::getColUpper() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getColUpper()\n");
  return &(this->highs->lp_.colUpper_[0]);
}

const double* OsiHiGHSSolverInterface::getRowLower() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowLower()\n");
  return &(this->highs->lp_.rowLower_[0]);
}

const double* OsiHiGHSSolverInterface::getRowUpper() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowUpper()\n");
  return &(this->highs->lp_.rowUpper_[0]);
}

const double* OsiHiGHSSolverInterface::getObjCoefficients() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getObjCoefficients()\n");
  return &(this->highs->lp_.colCost_[0]);
}

// TODO: review: 10^20?
double OsiHiGHSSolverInterface::getInfinity() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_NONE,
                    "Calling OsiHiGHSSolverInterface::getInfinity()\n");
  return HIGHS_CONST_INF;
}

const double* OsiHiGHSSolverInterface::getRowRange() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowRange()\n");
  if (this->rowRange != NULL) {
    delete[] this->rowRange;
  }

  int nrows = this->getNumRows();

  if (nrows == 0) {
    return this->rowRange;
  }

  this->rowRange = new double[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute range for row i
    double lo = this->highs->lp_.rowLower_[i];
    double hi = this->highs->lp_.rowUpper_[i];
    double t1;
    char t2;
    this->convertBoundToSense(lo, hi, t2, t1, this->rowRange[i]);
  }

  return this->rowRange;
}

const double* OsiHiGHSSolverInterface::getRightHandSide() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRightHandSide()\n");
  if (this->rhs != NULL) {
    delete[] this->rhs;
  }

  int nrows = this->getNumRows();

  if (nrows == 0) {
    return this->rhs;
  }

  this->rhs = new double[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute rhs for row i
    double lo = this->highs->lp_.rowLower_[i];
    double hi = this->highs->lp_.rowUpper_[i];
    double t1;
    char t2;
    this->convertBoundToSense(lo, hi, t2, this->rhs[i], t1);
  }

  return this->rhs;
}

const char* OsiHiGHSSolverInterface::getRowSense() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowSense()\n");
  if (this->rowSense != NULL) {
    delete[] this->rowSense;
  }

  int nrows = this->getNumRows();

  if (nrows == 0) {
    return this->rowSense;
  }

  this->rowSense = new char[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute sense for row i
    double lo = this->highs->lp_.rowLower_[i];
    double hi = this->highs->lp_.rowUpper_[i];
    double t1, t2;
    this->convertBoundToSense(lo, hi, this->rowSense[i], t1, t2);
  }

  return this->rowSense;
}

const CoinPackedMatrix* OsiHiGHSSolverInterface::getMatrixByCol() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getMatrixByCol()\n");
  if (this->matrixByCol != NULL) {
    delete this->matrixByCol;
  }

  int nrows = this->getNumRows();
  int ncols = this->getNumCols();
  int nelements = this->getNumElements();

  int* len = new int[ncols];
  int* start = new int[ncols + 1];
  int* index = new int[nelements];
  double* value = new double[nelements];

  // copy data
  memcpy(start, &(this->highs->lp_.Astart_[0]), (ncols + 1) * sizeof(int));
  memcpy(index, &(this->highs->lp_.Aindex_[0]), nelements * sizeof(int));
  memcpy(value, &(this->highs->lp_.Avalue_[0]), nelements * sizeof(double));

  for (int i = 0; i < ncols; i++) {
    len[i] = start[i + 1] - start[i];
  }

  this->matrixByCol = new CoinPackedMatrix();

  this->matrixByCol->assignMatrix(true, nrows, ncols, nelements, value, index,
                                  start, len);
  assert(this->matrixByCol->getNumCols() == ncols);
  assert(this->matrixByCol->getNumRows() == nrows);

  return this->matrixByCol;
}

const CoinPackedMatrix* OsiHiGHSSolverInterface::getMatrixByRow() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getMatrixByRow()\n");
  if (this->matrixByRow != NULL) {
    delete this->matrixByRow;
  }
  this->matrixByRow = new CoinPackedMatrix();
  this->matrixByRow->reverseOrderedCopyOf(*this->getMatrixByCol());

  return this->matrixByRow;
}

double OsiHiGHSSolverInterface::getObjSense() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getObjSense()\n");
  return (double)this->highs->lp_.sense_;
}

void OsiHiGHSSolverInterface::setObjSense(double s) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setObjSense()\n");
  ObjSense pass_sense = ObjSense::MINIMIZE;
  if (s == (double)ObjSense::MAXIMIZE) pass_sense = ObjSense::MAXIMIZE;
  this->highs->changeObjectiveSense(pass_sense);
}

void OsiHiGHSSolverInterface::addRow(const CoinPackedVectorBase& vec,
                                     const double rowlb, const double rowub) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::addRow()\n");
  bool success = this->highs->addRow(rowlb, rowub, vec.getNumElements(),
                                     vec.getIndices(), vec.getElements());
  assert(success);
  if (!success) {
    HighsPrintMessage(
        options.output, options.message_level, ML_ALWAYS,
        "Return from OsiHiGHSSolverInterface::addRow() is false\n");
  }
}

void OsiHiGHSSolverInterface::addRow(const CoinPackedVectorBase& vec,
                                     const char rowsen, const double rowrhs,
                                     const double rowrng) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::addRow()\n");
  // Assign arbitrary values so that compilation is clean
  double lb = 0;
  double ub = 1e200;
  this->convertSenseToBound(rowsen, rowrhs, rowrng, lb, ub);
  this->addRow(vec, lb, ub);
}

void OsiHiGHSSolverInterface::addCol(const CoinPackedVectorBase& vec,
                                     const double collb, const double colub,
                                     const double obj) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::addCol()\n");
  bool success = this->highs->addCol(obj, collb, colub, vec.getNumElements(),
                                     vec.getIndices(), vec.getElements());
  assert(success);
  if (!success) {
    HighsPrintMessage(
        options.output, options.message_level, ML_ALWAYS,
        "Return from OsiHiGHSSolverInterface::addCol() is false\n");
  }
}

void OsiHiGHSSolverInterface::deleteCols(const int num, const int* colIndices) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::deleteCols()\n");
  this->highs->deleteCols(num, colIndices);
}

void OsiHiGHSSolverInterface::deleteRows(const int num, const int* rowIndices) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::deleteRows()\n");
  this->highs->deleteRows(num, rowIndices);
}

void OsiHiGHSSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
                                            double*& collb, double*& colub,
                                            double*& obj, double*& rowlb,
                                            double*& rowub) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::assignProblem()\n");
  loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
  delete matrix;
  matrix = 0;
  delete[] collb;
  collb = 0;
  delete[] colub;
  colub = 0;
  delete[] obj;
  obj = 0;
  delete[] rowlb;
  rowlb = 0;
  delete[] rowub;
  rowub = 0;
}

void OsiHiGHSSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
                                          const double* collb,
                                          const double* colub,
                                          const double* obj, const char* rowsen,
                                          const double* rowrhs,
                                          const double* rowrng) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
  int numRow = matrix.getNumRows();

  double* rowlb = new double[numRow];
  double* rowub = new double[numRow];

  char* myrowsen = (char*)rowsen;
  bool rowsennull = false;
  double* myrowrhs = (double*)rowrhs;
  bool rowrhsnull = false;
  double* myrowrng = (double*)rowrng;
  bool rowrngnull = false;

  if (rowsen == NULL) {
    rowsennull = true;
    myrowsen = new char[numRow];
    for (int i = 0; i < numRow; i++) {
      myrowsen[i] = 'G';
    }
  }

  if (rowrhs == NULL) {
    rowsennull = true;
    myrowrhs = new double[numRow];
    for (int i = 0; i < numRow; i++) {
      myrowrhs[i] = 0.0;
    }
  }

  if (rowrng == NULL) {
    rowrngnull = true;
    myrowrng = new double[numRow];
    for (int i = 0; i < numRow; i++) {
      myrowrng[i] = 0.0;
    }
  }

  for (int i = 0; i < numRow; i++) {
    this->convertSenseToBound(myrowsen[i], myrowrhs[i], myrowrng[i], rowlb[i],
                              rowub[i]);
  }

  this->loadProblem(matrix, collb, colub, obj, rowlb, rowub);

  delete[] rowlb;
  delete[] rowub;

  if (rowsennull) {
    delete[] myrowsen;
  }

  if (rowrhsnull) {
    delete[] myrowrhs;
  }

  if (rowrngnull) {
    delete[] myrowrng;
  }
}

void OsiHiGHSSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
                                            double*& collb, double*& colub,
                                            double*& obj, char*& rowsen,
                                            double*& rowrhs, double*& rowrng) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::assignProblem()\n");
  loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
  delete matrix;
  matrix = 0;
  delete[] collb;
  collb = 0;
  delete[] colub;
  colub = 0;
  delete[] obj;
  obj = 0;
  delete[] rowsen;
  rowsen = 0;
  delete[] rowrhs;
  rowrhs = 0;
  delete[] rowrng;
  rowrng = 0;
}

void OsiHiGHSSolverInterface::loadProblem(
    const int numcols, const int numrows, const CoinBigIndex* start,
    const int* index, const double* value, const double* collb,
    const double* colub, const double* obj, const double* rowlb,
    const double* rowub) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
  double oldObjSense = this->getObjSense();

  HighsLp lp;

  lp.numRow_ = numrows;
  lp.numCol_ = numcols;
  lp.nnz_ = start[numcols];

  // setup HighsLp data structures
  lp.colCost_.resize(numcols);
  lp.colUpper_.resize(numcols);
  lp.colLower_.resize(numcols);

  lp.rowLower_.resize(numrows);
  lp.rowUpper_.resize(numrows);

  lp.Astart_.resize(numcols + 1);
  lp.Aindex_.resize(start[numcols]);
  lp.Avalue_.resize(start[numcols]);

  // copy data
  if (obj != NULL) {
    lp.colCost_.assign(obj, obj + numcols);
  } else {
    lp.colCost_.assign(numcols, 0.0);
  }

  if (collb != NULL) {
    lp.colLower_.assign(collb, collb + numcols);
  } else {
    lp.colLower_.assign(numcols, 0.0);
  }

  if (colub != NULL) {
    lp.colUpper_.assign(colub, colub + numcols);
  } else {
    lp.colUpper_.assign(numcols, HIGHS_CONST_INF);
  }

  if (rowlb != NULL) {
    lp.rowLower_.assign(rowlb, rowlb + numrows);
  } else {
    lp.rowLower_.assign(numrows, -HIGHS_CONST_INF);
  }

  if (rowub != NULL) {
    lp.rowUpper_.assign(rowub, rowub + numrows);
  } else {
    lp.rowUpper_.assign(numrows, HIGHS_CONST_INF);
  }

  lp.Astart_.assign(start, start + numcols + 1);
  lp.Aindex_.assign(index, index + start[numcols]);
  lp.Avalue_.assign(value, value + start[numcols]);
  this->highs->passModel(lp);
  this->setObjSense(oldObjSense);
}

void OsiHiGHSSolverInterface::loadProblem(
    const int numcols, const int numrows, const CoinBigIndex* start,
    const int* index, const double* value, const double* collb,
    const double* colub, const double* obj, const char* rowsen,
    const double* rowrhs, const double* rowrng) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
  double* rowlb = new double[numrows];
  double* rowub = new double[numrows];

  for (int i = 0; i < numrows; i++) {
    this->convertSenseToBound(rowsen[i], rowrhs[i], rowrng[i], rowlb[i],
                              rowub[i]);
  }

  this->loadProblem(numcols, numrows, start, index, value, collb, colub, obj,
                    rowlb, rowub);

  delete[] rowlb;
  delete[] rowub;
}

void OsiHiGHSSolverInterface::loadProblem(
    const CoinPackedMatrix& matrix, const double* collb, const double* colub,
    const double* obj, const double* rowlb, const double* rowub) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
  bool transpose = false;
  if (!matrix.isColOrdered()) {
    transpose = true;
    // ToDo: remove this hack
    //((CoinPackedMatrix *)&matrix)->transpose();
    ((CoinPackedMatrix*)&matrix)->reverseOrdering();
  }

  int numCol = matrix.getNumCols();
  int numRow = matrix.getNumRows();
  int nnz = matrix.getNumElements();

  int* start = new int[numCol + 1];
  int* index = new int[nnz];
  double* value = new double[nnz];

  // get matrix data
  // const CoinBigIndex *vectorStarts = matrix.getVectorStarts();
  const int* vectorLengths = matrix.getVectorLengths();
  const double* elements = matrix.getElements();
  const int* indices = matrix.getIndices();

  // set matrix in HighsLp
  start[0] = 0;
  int nz = 0;
  for (int i = 0; i < numCol; i++) {
    start[i + 1] = start[i] + vectorLengths[i];
    CoinBigIndex first = matrix.getVectorFirst(i);
    for (int j = 0; j < vectorLengths[i]; j++) {
      index[nz] = indices[first + j];
      value[nz] = elements[first + j];
      nz++;
    }
  }
  assert(nnz == nz);

  this->loadProblem(numCol, numRow, start, index, value, collb, colub, obj,
                    rowlb, rowub);

  if (transpose) {
    //((CoinPackedMatrix)matrix).transpose();
    ((CoinPackedMatrix*)&matrix)->reverseOrdering();
  }

  delete[] start;
  delete[] index;
  delete[] value;
}

/// Read a problem in MPS format from the given filename.
// int OsiHiGHSSolverInterface::readMps(const char *filename,
//   const char *extension)
// {
//   HighsOptions& options = this->highs->options_;
//   HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
//                     "Calling OsiHiGHSSolverInterface::readMps()\n");

//   HighsLp lp;

//   highs->options_.filename = std::string(filename) + "." +
//   std::string(extension);

//   FilereaderRetcode rc = FilereaderMps().readModelFromFile(highs->options_,
//   lp); if (rc != FilereaderRetcode::OK)
// 	  return (int)rc;
//   this->setDblParam(OsiDblParam::OsiObjOffset, lp.offset_);
//   highs->passModel(lp);

//   return 0;
// }

/// Write the problem into an mps file of the given filename.
void OsiHiGHSSolverInterface::writeMps(const char* filename,
                                       const char* extension,
                                       double objSense) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::writeMps()\n");

  std::string fullname = std::string(filename) + "." + std::string(extension);

  if (objSense != 0.0) {
    // HiGHS doesn't do funny stuff with the objective sense, so use Osi's
    // method if something strange is requested
    OsiSolverInterface::writeMpsNative(fullname.c_str(), NULL, NULL, 0, 2,
                                       objSense);
    return;
  }

  FilereaderMps frmps;
  HighsStatus rc =
      frmps.writeModelToFile(highs->options_, fullname.c_str(), highs->lp_);

  if (rc != HighsStatus::OK)
    throw CoinError("Creating MPS file failed", "writeMps",
                    "OsiHiGHSSolverInterface", __FILE__, __LINE__);
}

void OsiHiGHSSolverInterface::passInMessageHandler(
    CoinMessageHandler* handler) {
  OsiSolverInterface::passInMessageHandler(handler);

  HighsSetMessageCallback(printtomessagehandler, logtomessagehandler,
                          (void*)handler);

  // because HiGHS calls HiGHSSetIO with the Options, whichoverwrites the
  // previous setting
  this->highs->options_.printmsgcb = printtomessagehandler;
  this->highs->options_.logmsgcb = logtomessagehandler;
  this->highs->options_.msgcb_data = (void*)handler_;
}

const double* OsiHiGHSSolverInterface::getColSolution() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getColSolution()\n");
  if (!highs) {
    return nullptr;
  } else {
    if (highs->solution_.col_value.size() == 0) {
      double num_cols = highs->lp_.numCol_;
      this->dummy_solution->col_value.resize(num_cols);
      for (int col = 0; col < highs->lp_.numCol_; col++) {
        if (highs->lp_.colLower_[col] <= 0 && highs->lp_.colUpper_[col] >= 0)
          dummy_solution->col_value[col] = 0;
        else if (std::fabs(highs->lp_.colLower_[col] <
                           std::fabs(highs->lp_.colUpper_[col])))
          dummy_solution->col_value[col] = highs->lp_.colLower_[col];
        else
          dummy_solution->col_value[col] = highs->lp_.colUpper_[col];
      }
      return &dummy_solution->col_value[0];
    }
  }

  return &highs->solution_.col_value[0];
}

const double* OsiHiGHSSolverInterface::getRowPrice() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowPrice()\n");
  if (!highs)
    return nullptr;
  else {
    if (highs->solution_.row_dual.size() == 0) {
      double num_cols = highs->lp_.numCol_;
      this->dummy_solution->row_dual.resize(num_cols);
      for (int col = 0; col < highs->lp_.numCol_; col++) {
        if (highs->lp_.colLower_[col] <= 0 && highs->lp_.colUpper_[col] >= 0)
          dummy_solution->row_dual[col] = 0;
        else if (std::fabs(highs->lp_.colLower_[col] <
                           std::fabs(highs->lp_.colUpper_[col])))
          dummy_solution->row_dual[col] = highs->lp_.colLower_[col];
        else
          dummy_solution->row_dual[col] = highs->lp_.colUpper_[col];
      }
      return &dummy_solution->row_dual[0];
    }
  }

  return &highs->solution_.row_dual[0];
}

const double* OsiHiGHSSolverInterface::getReducedCost() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getReducedCost()\n");
  if (!highs)
    return nullptr;
  else {
    if (highs->solution_.col_dual.size() == 0) {
      const HighsLp& lp = highs->lp_;
      double num_cols = lp.numCol_;
      this->dummy_solution->col_dual.resize(num_cols);
      for (int col = 0; col < num_cols; col++) {
        dummy_solution->col_dual[col] = lp.colCost_[col];
        for (int i = lp.Astart_[col]; i < lp.Astart_[col + 1]; i++) {
          const int row = lp.Aindex_[i];
          assert(row >= 0);
          assert(row < lp.numRow_);

          dummy_solution->col_dual[col] +=
              dummy_solution->row_dual[row] * lp.Avalue_[i];
        }
      }
      return &dummy_solution->col_dual[0];
    }
  }

  return &highs->solution_.col_dual[0];
}

const double* OsiHiGHSSolverInterface::getRowActivity() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowActivity()\n");
  if (!highs)
    return nullptr;
  else {
    if (highs->solution_.row_value.size() == 0) {
      double num_cols = highs->lp_.numCol_;
      this->dummy_solution->row_value.resize(num_cols);
      for (int col = 0; col < highs->lp_.numCol_; col++) {
        if (highs->lp_.colLower_[col] <= 0 && highs->lp_.colUpper_[col] >= 0)
          dummy_solution->row_value[col] = 0;
        else if (std::fabs(highs->lp_.colLower_[col] <
                           std::fabs(highs->lp_.colUpper_[col])))
          dummy_solution->row_value[col] = highs->lp_.colLower_[col];
        else
          dummy_solution->row_value[col] = highs->lp_.colUpper_[col];
      }
      return &dummy_solution->row_value[0];
    }
  }

  return &highs->solution_.row_value[0];
}

double OsiHiGHSSolverInterface::getObjValue() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getObjValue()\n");
  double objVal = 0.0;
  if (true || highs->solution_.col_value.size() == 0) {
    const double* sol = this->getColSolution();
    const double* cost = this->getObjCoefficients();
    int ncols = this->getNumCols();

    objVal = -this->objOffset;
    for (int i = 0; i < ncols; i++) {
      objVal += sol[i] * cost[i];
    }
  } else {
    this->highs->getHighsInfoValue("objective_function_value", objVal);
  }

  return objVal;
}

int OsiHiGHSSolverInterface::getIterationCount() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getIterationCount()\n");
  if (!highs) {
    return 0;
  }
  int iteration_count;
  this->highs->getHighsInfoValue("simplex_iteration_count", iteration_count);
  return iteration_count;
}

void OsiHiGHSSolverInterface::setRowPrice(const double* rowprice) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,

                    "Calling OsiHiGHSSolverInterface::setRowPrice()\n");
  if (!rowprice) return;
  HighsSolution solution;
  solution.row_dual.resize(highs->lp_.numRow_);
  for (int row = 0; row < highs->lp_.numRow_; row++)
    solution.row_dual[row] = rowprice[row];

  /*HighsStatus result =*/highs->setSolution(solution);
}

void OsiHiGHSSolverInterface::setColSolution(const double* colsol) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setColSolution()\n");
  if (!colsol) return;
  HighsSolution solution;
  solution.col_value.resize(highs->lp_.numCol_);
  for (int col = 0; col < highs->lp_.numCol_; col++)
    solution.col_value[col] = colsol[col];

  /*HighsStatus result =*/highs->setSolution(solution);
}

void OsiHiGHSSolverInterface::applyRowCut(const OsiRowCut& rc) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::applyRowCut()\n");
}

void OsiHiGHSSolverInterface::applyColCut(const OsiColCut& cc) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::applyColCut()\n");
}

void OsiHiGHSSolverInterface::setContinuous(int index) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setContinuous()\n");
}

void OsiHiGHSSolverInterface::setInteger(int index) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setInteger()\n");
}

bool OsiHiGHSSolverInterface::isContinuous(int colNumber) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::isContinuous()\n");
  return true;
}

void OsiHiGHSSolverInterface::setRowType(int index, char sense,
                                         double rightHandSide, double range) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setRowType()\n");
  // Assign arbitrary values so that compilation is clean
  double lo = 0;
  double hi = 1e200;
  this->convertSenseToBound(sense, rightHandSide, range, lo, hi);
  this->setRowBounds(index, lo, hi);
}

void OsiHiGHSSolverInterface::setRowLower(int elementIndex,
                                          double elementValue) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setRowLower()\n");

  double upper = this->getRowUpper()[elementIndex];

  this->highs->changeRowBounds(elementIndex, elementValue, upper);
}

void OsiHiGHSSolverInterface::setRowUpper(int elementIndex,
                                          double elementValue) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setRowUpper()\n");
  double lower = this->getRowLower()[elementIndex];
  this->highs->changeRowBounds(elementIndex, lower, elementValue);
}

void OsiHiGHSSolverInterface::setColLower(int elementIndex,
                                          double elementValue) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setColLower()\n");
  double upper = this->getColUpper()[elementIndex];
  this->highs->changeColBounds(elementIndex, elementValue, upper);
}

void OsiHiGHSSolverInterface::setColUpper(int elementIndex,
                                          double elementValue) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setColUpper()\n");
  double lower = this->getColLower()[elementIndex];
  this->highs->changeColBounds(elementIndex, lower, elementValue);
}

void OsiHiGHSSolverInterface::setObjCoeff(int elementIndex,
                                          double elementValue) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setObjCoeff()\n");
  this->highs->changeColCost(elementIndex, elementValue);
}

std::vector<double*> OsiHiGHSSolverInterface::getDualRays(int maxNumRays,
                                                          bool fullRay) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getDualRays()\n");
  return std::vector<double*>(0);
}

std::vector<double*> OsiHiGHSSolverInterface::getPrimalRays(
    int maxNumRays) const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getPrimalRays()\n");
  return std::vector<double*>(0);
}

CoinWarmStart* OsiHiGHSSolverInterface::getEmptyWarmStart() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getEmptyWarmStart()\n");
  return (dynamic_cast<CoinWarmStart*>(new CoinWarmStartBasis()));
}

CoinWarmStart* OsiHiGHSSolverInterface::getWarmStart() const {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getWarmStart()\n");
  if (!highs) return NULL;

  if (highs->basis_.col_status.size() == 0 ||
      highs->basis_.row_status.size() == 0)
    return NULL;

  int num_cols = highs->lp_.numCol_;
  int num_rows = highs->lp_.numRow_;

  int* cstat = new int[num_cols];
  int* rstat = new int[num_rows];

  getBasisStatus(cstat, rstat);

  CoinWarmStartBasis* warm_start = new CoinWarmStartBasis();
  warm_start->setSize(num_cols, num_rows);

  for (int i = 0; i < num_rows; ++i)
    warm_start->setArtifStatus(i, CoinWarmStartBasis::Status(rstat[i]));
  for (int i = 0; i < num_cols; ++i)
    warm_start->setStructStatus(i, CoinWarmStartBasis::Status(cstat[i]));

  return warm_start;
}

bool OsiHiGHSSolverInterface::setWarmStart(const CoinWarmStart* warmstart) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setWarmStart()\n");
  return false;
}

void OsiHiGHSSolverInterface::resolve() {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::resolve()\n");
  this->status = this->highs->run();
}

void OsiHiGHSSolverInterface::setRowBounds(int elementIndex, double lower,
                                           double upper) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setRowBounds()\n");

  this->highs->changeRowBounds(elementIndex, lower, upper);
}

void OsiHiGHSSolverInterface::setColBounds(int elementIndex, double lower,
                                           double upper) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setColBounds()\n");

  this->highs->changeColBounds(elementIndex, lower, upper);
}

void OsiHiGHSSolverInterface::setRowSetBounds(const int* indexFirst,
                                              const int* indexLast,
                                              const double* boundList) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setRowSetBounds()\n");
  OsiSolverInterface::setRowSetBounds(indexFirst, indexLast - 1, boundList);
}

void OsiHiGHSSolverInterface::setColSetBounds(const int* indexFirst,
                                              const int* indexLast,
                                              const double* boundList) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setColSetBounds()\n");
  OsiSolverInterface::setColSetBounds(indexFirst, indexLast - 1, boundList);
}

void OsiHiGHSSolverInterface::branchAndBound() {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::branchAndBound()\n");
  // TODO
}

void OsiHiGHSSolverInterface::setObjCoeffSet(const int* indexFirst,
                                             const int* indexLast,
                                             const double* coeffList) {
  HighsOptions& options = this->highs->options_;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setObjCoeffSet()\n");
  OsiSolverInterface::setObjCoeffSet(indexFirst, indexLast - 1, coeffList);
}

int OsiHiGHSSolverInterface::canDoSimplexInterface() const { return 0; }

/* Osi return codes:
0: free
1: basic
2: upper
3: lower
*/
void OsiHiGHSSolverInterface::getBasisStatus(int* cstat, int* rstat) const {
  if (!highs) return;

  if (highs->basis_.col_status.size() == 0 ||
      highs->basis_.row_status.size() == 0)
    return;

  for (size_t i = 0; i < highs->basis_.col_status.size(); ++i)
    switch (highs->basis_.col_status[i]) {
      case HighsBasisStatus::BASIC:
        cstat[i] = 1;
        break;
      case HighsBasisStatus::LOWER:
        cstat[i] = 3;
        break;
      case HighsBasisStatus::UPPER:
        cstat[i] = 2;
        break;
      case HighsBasisStatus::SUPER:
        cstat[i] = 0;
        break;
      case HighsBasisStatus::ZERO:
        cstat[i] = 0;
        break;
      case HighsBasisStatus::NONBASIC:
        cstat[i] = 3;
        break;
    }

  for (size_t i = 0; i < highs->basis_.row_status.size(); ++i)
    switch (highs->basis_.row_status[i]) {
      case HighsBasisStatus::BASIC:
        rstat[i] = 1;
        break;
      case HighsBasisStatus::LOWER:
        rstat[i] = 3;
        break;
      case HighsBasisStatus::UPPER:
        rstat[i] = 2;
        break;
      case HighsBasisStatus::SUPER:
        rstat[i] = 0;
        break;
      case HighsBasisStatus::ZERO:
        rstat[i] = 0;
        break;
      case HighsBasisStatus::NONBASIC:
        rstat[i] = 3;
        break;
    }
}

void OsiHiGHSSolverInterface ::setRowNames(OsiNameVec& srcNames, int srcStart,
                                           int len, int tgtStart) {}

void OsiHiGHSSolverInterface ::setColNames(OsiNameVec& srcNames, int srcStart,
                                           int len, int tgtStart) {}

void OsiSolverInterfaceMpsUnitTest(
    const std::vector<OsiSolverInterface*>& vecSiP, const std::string& mpsDir) {
}
