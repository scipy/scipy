/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "basis.hpp"

#include <cassert>
#include <memory>

Basis::Basis(Runtime& rt, std::vector<HighsInt> active,
             std::vector<BasisStatus> status, std::vector<HighsInt> inactive)
    : Ztprod_res(rt.instance.num_var),
      buffer_Zprod(rt.instance.num_var),
      runtime(rt),
      buffer_column_aq(rt.instance.num_var),
      buffer_row_ep(rt.instance.num_var) {
  buffer_vec2hvec.setup(rt.instance.num_var);

  for (HighsInt i = 0; i < runtime.instance.num_var + runtime.instance.num_con;
       i++) {
    basisstatus[i] = BasisStatus::kInactive;
  }

  for (size_t i = 0; i < active.size(); i++) {
    active_constraint_index.push_back(active[i]);
    basisstatus[active_constraint_index[i]] = status[i];
  }
  for (size_t i = 0; i < inactive.size(); i++) {
    non_active_constraint_index.push_back(inactive[i]);
    basisstatus[non_active_constraint_index[i]] = BasisStatus::kInactiveInBasis;
  }

  Atran = rt.instance.A.t();

  col_aq.setup(rt.instance.num_var);
  row_ep.setup(rt.instance.num_var);

  build();
}

void Basis::build() {
  //  report();

  updatessinceinvert = 0;

  baseindex.resize(active_constraint_index.size() +
                   non_active_constraint_index.size());
  constraintindexinbasisfactor.clear();

  basisfactor = HFactor();

  constraintindexinbasisfactor.assign(Atran.num_row + Atran.num_col, -1);
  assert((HighsInt)(non_active_constraint_index.size() +
                    active_constraint_index.size()) == Atran.num_row);

  HighsInt counter = 0;
  for (HighsInt i : non_active_constraint_index) {
    baseindex[counter++] = i;
  }
  for (HighsInt i : active_constraint_index) {
    baseindex[counter++] = i;
  }

  basisfactor.setup(Atran.num_col, Atran.num_row, Atran.start.data(),
                    Atran.index.data(), Atran.value.data(), baseindex.data());
  basisfactor.build();

  for (size_t i = 0;
       i < active_constraint_index.size() + non_active_constraint_index.size();
       i++) {
    constraintindexinbasisfactor[baseindex[i]] = i;
  }
}

void Basis::rebuild() {
  //  report();

  updatessinceinvert = 0;
  constraintindexinbasisfactor.clear();

  constraintindexinbasisfactor.assign(Atran.num_row + Atran.num_col, -1);
  assert((HighsInt)(non_active_constraint_index.size() +
                    active_constraint_index.size()) == Atran.num_row);

  basisfactor.build();

  for (size_t i = 0;
       i < active_constraint_index.size() + non_active_constraint_index.size();
       i++) {
    constraintindexinbasisfactor[baseindex[i]] = i;
  }
  reinversion_hint = false;
}

void Basis::report() {
  //
  // Basis of dimension qp_num_var, analogous to primal simplex
  // nonbasic variables, partitioned into
  //
  // * Indices of active variables/constraints, so index values in {0,
  // * ..., qp_num_con-1}. These are listed in active_constraint_index
  // * and, for each the value of basisstatus will be ActiveAtLower or
  // * ActiveAtUpper
  //
  // * Indices of inactive variables, so index values in {qp_num_con,
  // * ..., qp_num_con+qp_num_var-1} (or possibly also from {0, ...,
  // * qp_num_con-1}?) used to complete the basis. The number of
  // * inactive variables defines the dimension of the null space.
  //
  // Remaining qp_num_con indices may be degenerate, otherwise they
  // are off their bounds. They are analogous to primal simplex basic
  // variables, in that their values are solved for.
  //
  // Hence the correspondence between the QP basis and a HiGHS
  // (simplex) basis is as follows
  //
  // For variables or constraints
  //
  // BasisStatus::kInactive: HighsBasisStatus::kBasic
  //
  // BasisStatus::kActiveAtLower: HighsBasisStatus::kLower
  //
  // BasisStatus::kActiveAtUpper: HighsBasisStatus::kUpper
  //
  // BasisStatus::kInactiveInBasis: HighsBasisStatus::kNonbasic
  //
  //
  const HighsInt qp_num_var = Atran.num_row;
  const HighsInt qp_num_con = Atran.num_col;
  const HighsInt num_active_in_basis = active_constraint_index.size();
  const HighsInt num_inactive_in_basis = non_active_constraint_index.size();

  HighsInt num_var_inactive = 0;
  HighsInt num_var_active_at_lower = 0;
  HighsInt num_var_active_at_upper = 0;
  HighsInt num_var_inactive_in_basis = 0;
  HighsInt num_con_inactive = 0;
  HighsInt num_con_active_at_lower = 0;
  HighsInt num_con_active_at_upper = 0;
  HighsInt num_con_inactive_in_basis = 0;

  for (HighsInt i = 0; i < qp_num_var; i++) {
    switch (basisstatus[qp_num_con + i]) {
      case BasisStatus::kInactive:
        num_var_inactive++;
        continue;
      case BasisStatus::kActiveAtLower:
        num_var_active_at_lower++;
        continue;
      case BasisStatus::kActiveAtUpper:
        num_var_active_at_upper++;
        continue;
      case BasisStatus::kInactiveInBasis:
        num_var_inactive_in_basis++;
        continue;
      default:
        assert(111 == 123);
    }
  }

  for (HighsInt i = 0; i < qp_num_con; i++) {
    switch (basisstatus[i]) {
      case BasisStatus::kInactive:
        num_con_inactive++;
        continue;
      case BasisStatus::kActiveAtLower:
        num_con_active_at_lower++;
        continue;
      case BasisStatus::kActiveAtUpper:
        num_con_active_at_upper++;
        continue;
      case BasisStatus::kInactiveInBasis:
        num_con_inactive_in_basis++;
        continue;
      default:
        assert(111 == 123);
    }
  }

  const bool print_basis = num_inactive_in_basis + num_active_in_basis < 100;
  if (print_basis) {
    printf("basis: ");
    for (HighsInt a_ : active_constraint_index) {
      if (a_ < qp_num_con) {
        printf("c%-3d ", int(a_));
      } else {
        printf("v%-3d ", int(a_ - qp_num_con));
      }
    }
    printf(" - ");
    for (HighsInt n_ : non_active_constraint_index) {
      if (n_ < qp_num_con) {
        printf("c%-3d ", int(n_));
      } else {
        printf("v%-3d ", int(n_ - qp_num_con));
      }
    }
    printf("\n");
  }

  printf("Basis::report: QP(%6d [inact %6d; act %6d], %6d)", int(qp_num_var),
         int(num_inactive_in_basis), int(num_active_in_basis), int(qp_num_con));
  printf(
      " (inact / lo / up / basis) for var (%6d / %6d / %6d / %6d) and con (%6d "
      "/ %6d / %6d / %6d)\n",
      int(num_var_inactive), int(num_var_active_at_lower),
      int(num_var_active_at_upper), int(num_var_inactive_in_basis),
      int(num_con_inactive), int(num_con_active_at_lower),
      int(num_con_active_at_upper), int(num_con_inactive_in_basis));
  assert(qp_num_var == num_inactive_in_basis + num_active_in_basis);
  assert(qp_num_con == num_var_inactive + num_con_inactive);
  assert(num_inactive_in_basis ==
         num_var_inactive_in_basis + num_con_inactive_in_basis);
  assert(num_active_in_basis ==
         num_var_active_at_lower + num_var_active_at_upper +
             num_con_active_at_lower + num_con_active_at_upper);
}

// move that constraint into V section basis (will correspond to Nullspace
// from now on)
void Basis::deactivate(HighsInt conid) {
  // printf("deact %" HIGHSINT_FORMAT "\n", conid);
  assert(contains(active_constraint_index, conid));
  basisstatus[conid] = BasisStatus::kInactiveInBasis;
  remove(active_constraint_index, conid);
  non_active_constraint_index.push_back(conid);
}

QpSolverStatus Basis::activate(const Settings& settings, HighsInt conid,
                               BasisStatus newstatus,
                               HighsInt nonactivetoremove, Pricing* pricing) {
  // printf("activ %" HIGHSINT_FORMAT "\n", conid);
  if (!contains(active_constraint_index, (HighsInt)conid)) {
    basisstatus[nonactivetoremove] = BasisStatus::kInactive;
    basisstatus[conid] = newstatus;
    active_constraint_index.push_back(conid);
  } else {
    printf("Degeneracy? constraint %" HIGHSINT_FORMAT " already in basis\n",
           conid);
    return QpSolverStatus::DEGENERATE;
  }

  // printf("drop %d\n", nonactivetoremove);
  // remove non-active row from basis
  HighsInt rowtoremove = constraintindexinbasisfactor[nonactivetoremove];

  baseindex[rowtoremove] = conid;
  remove(non_active_constraint_index, nonactivetoremove);
  updatebasis(settings, conid, nonactivetoremove, pricing);

  if (updatessinceinvert != 0) {
    constraintindexinbasisfactor[nonactivetoremove] = -1;
    constraintindexinbasisfactor[conid] = rowtoremove;
  }
  return QpSolverStatus::OK;
}

void Basis::updatebasis(const Settings& settings, HighsInt newactivecon,
                        HighsInt droppedcon, Pricing* pricing) {
  if (newactivecon == droppedcon) {
    return;
  }

  const HighsInt kHintNotChanged = 99999;
  HighsInt hint = kHintNotChanged;

  HighsInt droppedcon_rowindex = constraintindexinbasisfactor[droppedcon];
  if (buffered_p != droppedcon) {
    row_ep.clear();
    row_ep.packFlag = true;
    row_ep.index[0] = droppedcon_rowindex;
    row_ep.array[droppedcon_rowindex] = 1.0;
    row_ep.count = 1;
    basisfactor.btranCall(row_ep, 1.0);
  }

  pricing->update_weights(hvec2vec(col_aq), hvec2vec(row_ep), droppedcon,
                          newactivecon);
  HighsInt row_out = droppedcon_rowindex;

  basisfactor.update(&col_aq, &row_ep, &row_out, &hint);

  updatessinceinvert++;
  if (updatessinceinvert >= settings.reinvertfrequency ||
      hint != kHintNotChanged) {
    reinversion_hint = true;
  }
  // since basis changed, buffered values are no longer valid
  buffered_p = -1;
  buffered_q = -1;
}

QpVector& Basis::btran(const QpVector& rhs, QpVector& target, bool buffer,
                       HighsInt p) {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.btranCall(rhs_hvec, 1.0);
  if (buffer) {
    row_ep.copy(&rhs_hvec);
    for (HighsInt i = 0; i < rhs_hvec.packCount; i++) {
      row_ep.packIndex[i] = rhs_hvec.packIndex[i];
      row_ep.packValue[i] = rhs_hvec.packValue[i];
    }
    row_ep.packCount = rhs_hvec.packCount;
    row_ep.packFlag = rhs_hvec.packFlag;
    buffered_q = p;
  }
  return hvec2vec(rhs_hvec, target);
}

QpVector Basis::btran(const QpVector& rhs, bool buffer, HighsInt p) {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.btranCall(rhs_hvec, 1.0);
  if (buffer) {
    row_ep.copy(&rhs_hvec);
    for (HighsInt i = 0; i < rhs_hvec.packCount; i++) {
      row_ep.packIndex[i] = rhs_hvec.packIndex[i];
      row_ep.packValue[i] = rhs_hvec.packValue[i];
    }
    row_ep.packCount = rhs_hvec.packCount;
    row_ep.packFlag = rhs_hvec.packFlag;
    buffered_q = p;
  }
  return hvec2vec(rhs_hvec);
}

QpVector& Basis::ftran(const QpVector& rhs, QpVector& target, bool buffer,
                       HighsInt q) {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.ftranCall(rhs_hvec, 1.0);
  if (buffer) {
    col_aq.copy(&rhs_hvec);
    for (HighsInt i = 0; i < rhs_hvec.packCount; i++) {
      col_aq.packIndex[i] = rhs_hvec.packIndex[i];
      col_aq.packValue[i] = rhs_hvec.packValue[i];
    }
    col_aq.packCount = rhs_hvec.packCount;
    col_aq.packFlag = rhs_hvec.packFlag;
    buffered_q = q;
  }

  return hvec2vec(rhs_hvec, target);
}

QpVector Basis::ftran(const QpVector& rhs, bool buffer, HighsInt q) {
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.ftranCall(rhs_hvec, 1.0);
  if (buffer) {
    col_aq.copy(&rhs_hvec);
    for (HighsInt i = 0; i < rhs_hvec.packCount; i++) {
      col_aq.packIndex[i] = rhs_hvec.packIndex[i];
      col_aq.packValue[i] = rhs_hvec.packValue[i];
    }
    col_aq.packCount = rhs_hvec.packCount;
    col_aq.packFlag = rhs_hvec.packFlag;
    buffered_q = q;
  }
  return hvec2vec(rhs_hvec);
}

QpVector Basis::recomputex(const Instance& inst) {
  assert((HighsInt)active_constraint_index.size() == inst.num_var);
  QpVector rhs(inst.num_var);

  for (HighsInt i = 0; i < inst.num_var; i++) {
    HighsInt con = active_constraint_index[i];
    if (constraintindexinbasisfactor[con] == -1) {
      printf("error\n");
    }
    if (basisstatus[con] == BasisStatus::kActiveAtLower) {
      if (con < inst.num_con) {
        rhs.value[constraintindexinbasisfactor[con]] = inst.con_lo[con];
      } else {
        rhs.value[constraintindexinbasisfactor[con]] =
            inst.var_lo[con - inst.num_con];
      }
    } else {
      if (con < inst.num_con) {
        rhs.value[constraintindexinbasisfactor[con]] = inst.con_up[con];
        // rhs.value[i] = inst.con_up[con];
      } else {
        rhs.value[constraintindexinbasisfactor[con]] =
            inst.var_up[con - inst.num_con];
        // rhs.value[i] = inst.var_up[con - inst.num_con];
      }
    }

    rhs.index[i] = i;
    rhs.num_nz++;
  }
  HVector rhs_hvec = vec2hvec(rhs);
  basisfactor.btranCall(rhs_hvec, 1.0);
  return hvec2vec(rhs_hvec);
}

QpVector& Basis::Ztprod(const QpVector& rhs, QpVector& target, bool buffer,
                        HighsInt q) {
  ftran(rhs, Ztprod_res, buffer, q);

  target.reset();
  for (size_t i = 0; i < non_active_constraint_index.size(); i++) {
    HighsInt nonactive = non_active_constraint_index[i];
    HighsInt idx = constraintindexinbasisfactor[nonactive];
    target.index[i] = static_cast<HighsInt>(i);
    target.value[i] = Ztprod_res.value[idx];
  }
  target.resparsify();
  return target;
}

QpVector& Basis::Zprod(const QpVector& rhs, QpVector& target) {
  buffer_Zprod.reset();
  buffer_Zprod.dim = target.dim;
  for (HighsInt i = 0; i < rhs.num_nz; i++) {
    HighsInt nz = rhs.index[i];
    HighsInt nonactive = non_active_constraint_index[nz];
    HighsInt idx = constraintindexinbasisfactor[nonactive];
    buffer_Zprod.index[i] = idx;
    buffer_Zprod.value[idx] = rhs.value[nz];
  }
  buffer_Zprod.num_nz = rhs.num_nz;
  return btran(buffer_Zprod, target);
}

// void Basis::write(std::string filename) {
//    FILE* file = fopen(filename.c_str(), "w");

//    fprintf(file, "%lu %lu\n", active_constraint_index.size(),
//    non_active_constraint_index.size()); for (HighsInt i=0;
//    i<active_constraint_index.size(); i++) {
//       fprintf(file, "%" HIGHSINT_FORMAT " %" HIGHSINT_FORMAT "\n",
//       active_constraint_index[i], (HighsInt)rowstatus[i]);
//    }
//    for (HighsInt i=0; i<non_active_constraint_index.size(); i++) {
//       fprintf(file, "%" HIGHSINT_FORMAT " %" HIGHSINT_FORMAT "\n",
//       non_active_constraint_index[i], (HighsInt)rowstatus[i]);
//    }
//    // TODO

//    fclose(file);
// }
