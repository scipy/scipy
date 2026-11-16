/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma once
#ifndef __SRC_LIB_BASIS_HPP__
#define __SRC_LIB_BASIS_HPP__

#include <map>
#include <vector>

#include "qpsolver/instance.hpp"
#include "qpsolver/pricing.hpp"
#include "qpsolver/qpconst.hpp"
#include "qpsolver/runtime.hpp"
#include "qpsolver/snippets.hpp"
#include "util/HFactor.h"
#include "util/HVector.h"
#include "util/HVectorBase.h"

class Basis {
  HVector buffer_vec2hvec;
  QpVector Ztprod_res;
  QpVector buffer_Zprod;

  HVector& vec2hvec(const QpVector& vec) {
    buffer_vec2hvec.clear();
    for (HighsInt i = 0; i < vec.num_nz; i++) {
      buffer_vec2hvec.index[i] = vec.index[i];
      buffer_vec2hvec.array[vec.index[i]] = vec.value[vec.index[i]];
    }
    buffer_vec2hvec.count = vec.num_nz;
    buffer_vec2hvec.packFlag = true;
    return buffer_vec2hvec;
  }

  QpVector& hvec2vec(const HVector& hvec, QpVector& target) {
    target.reset();
    for (HighsInt i = 0; i < hvec.count; i++) {
      target.index[i] = hvec.index[i];
      target.value[target.index[i]] = hvec.array[hvec.index[i]];
    }
    // for (HighsInt i = 0; i < hvec.size; i++) {
    //   target.index[i] = hvec.index[i];
    //   target.value[i] = hvec.array[i];
    // }
    target.num_nz = hvec.count;
    return target;
  }

  QpVector hvec2vec(const HVector& hvec) {
    QpVector vec(hvec.size);

    return hvec2vec(hvec, vec);
  }

  Runtime& runtime;
  HFactor basisfactor;
  HighsInt updatessinceinvert = 0;

  MatrixBase Atran;

  // indices of active constraints in basis
  std::vector<HighsInt> active_constraint_index;

  // ids of constraints that are in the basis but not active
  // I need to extract those columns to get Z
  std::vector<HighsInt> non_active_constraint_index;

  // ids of constraints that are in the basis
  std::vector<HighsInt> baseindex;

  std::map<int, BasisStatus> basisstatus;

  // index i: -1 if constraint not in basis, [0, num_var] if
  // constraint in basis (active or not)
  std::vector<HighsInt> constraintindexinbasisfactor;

  void build();

  // buffer to avoid recreating vectors
  QpVector buffer_column_aq;
  QpVector buffer_row_ep;

  // buffers to prevent multiple btran/ftran
  HighsInt buffered_q = -1;
  HighsInt buffered_p = -1;
  HVector row_ep;
  HVector col_aq;

  bool reinversion_hint = false;

 public:
  Basis(Runtime& rt, std::vector<HighsInt> active,
        std::vector<BasisStatus> atlower, std::vector<HighsInt> inactive);

  bool getreinversionhint() { return reinversion_hint; }

  void rebuild();

  HighsInt getnupdatessinceinvert() { return updatessinceinvert; }

  HighsInt getnumactive() const { return active_constraint_index.size(); };

  HighsInt getnuminactive() const {
    return non_active_constraint_index.size();
  };

  const std::vector<HighsInt>& getactive() const {
    return active_constraint_index;
  };

  const std::vector<HighsInt>& getinactive() const {
    return non_active_constraint_index;
  };

  const std::vector<HighsInt>& getindexinfactor() const {
    return constraintindexinbasisfactor;
  };

  BasisStatus getstatus(HighsInt conid) { return basisstatus[conid]; };

  void report();

  // move that constraint into V section basis (will correspond to
  // Nullspace from now on)
  void deactivate(HighsInt conid);

  QpSolverStatus activate(const Settings& settings, HighsInt conid,
                          BasisStatus atlower, HighsInt nonactivetoremove,
                          Pricing* pricing);

  void updatebasis(const Settings& settings, HighsInt newactivecon,
                   HighsInt droppedcon, Pricing* pricing);

  QpVector btran(const QpVector& rhs, bool buffer = false, HighsInt p = -1);

  QpVector ftran(const QpVector& rhs, bool buffer = false, HighsInt q = -1);

  QpVector& btran(const QpVector& rhs, QpVector& target, bool buffer = false,
                  HighsInt p = -1);

  QpVector& ftran(const QpVector& rhs, QpVector& target, bool buffer = false,
                  HighsInt q = -1);

  QpVector recomputex(const Instance& inst);

  void write(std::string filename);

  QpVector& Ztprod(const QpVector& rhs, QpVector& target, bool buffer = false,
                   HighsInt q = -1);

  QpVector& Zprod(const QpVector& rhs, QpVector& target);
};

#endif
