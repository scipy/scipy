#ifndef HIPO_KKT_MATRIX_H
#define HIPO_KKT_MATRIX_H

#include <atomic>
#include <vector>

#include "Info.h"
#include "Model.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/factorhighs/Symbolic.h"

namespace hipo {

struct KktMatrix {
  std::vector<Int> ptrNE;
  std::vector<Int> rowsNE;
  std::vector<double> valNE;
  std::vector<Int> ptrA_rw, idxA_rw;
  std::vector<Int> corr_A;
  std::atomic<Int64> NE_nz_limit{kHighsIInf};

  bool NE_ready = false;
  bool AS_ready = false;

  std::vector<Int> ptrAS;
  std::vector<Int> rowsAS;
  std::vector<double> valAS;

  Symbolic S;

  const Model& model;
  const Regularisation& regul;
  Info& info;
  const Logger& logger;
  const Options& options;

  KktMatrix(const Model& model, const Regularisation& regul, Info& info,
            const Logger& logger, const Options& options);

  Int buildASstructure();
  Int buildASvalues(const std::vector<double>& scaling);
  void freeASmemory();

  Int buildNEstructure();
  Int buildNEvalues(const std::vector<double>& scaling);
  void freeNEmemory();

  Int n() const;
  Int nz() const;
  std::string nla() const;
  bool isNE() const;
  bool isAS() const;
  const std::vector<Int>& iperm() const;
};

}  // namespace hipo

#endif