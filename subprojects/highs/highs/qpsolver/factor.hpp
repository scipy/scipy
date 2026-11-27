/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_NEWFACTOR_HPP__
#define __SRC_LIB_NEWFACTOR_HPP__

#include <cassert>
#include <vector>

#include "matrix.hpp"
#include "qpconst.hpp"
#include "runtime.hpp"

using std::min;

class CholeskyFactor {
 private:
  bool uptodate = false;
  HighsInt numberofreduces = 0;

  Runtime& runtime;

  Basis& basis;

  HighsInt current_k = 0;
  HighsInt current_k_max;
  std::vector<double> L;

  bool has_negative_eigenvalue = false;
  std::vector<double> a;

  void resize(HighsInt new_k_max) {
    // If new_k_max is ever zero, then current_k_max = 0, so when
    // 2*current_k_max is passed into resize from expand, new_k_max is
    // still zero. Observed in #3453, but could not reproduce. No harm
    // in ensuring that the value used for new_k_max is at least 1.
    new_k_max = std::max(new_k_max, HighsInt(1));
    std::vector<double> L_old = L;
    L.clear();
    L.resize((new_k_max) * (new_k_max));
    const HighsInt l_size = L.size();
    // Driven by #958, changes made in following lines to avoid array
    // bound error when new_k_max < current_k_max
    HighsInt min_k_max = min(new_k_max, current_k_max);
    for (HighsInt i = 0; i < min_k_max; i++) {
      for (HighsInt j = 0; j < min_k_max; j++) {
        assert(i * (new_k_max) + j < l_size);
        L[i * (new_k_max) + j] = L_old[i * current_k_max + j];
      }
    }
    current_k_max = new_k_max;
  }

 public:
  CholeskyFactor(Runtime& rt, Basis& bas) : runtime(rt), basis(bas) {
    uptodate = false;
    current_k_max =
        max(min((HighsInt)ceil(rt.instance.num_var / 16.0), (HighsInt)1000),
            basis.getnuminactive());
    L.resize(current_k_max * current_k_max);
  }

  // the factor no longer matches the null space: recompute it when next used
  void invalidate() { uptodate = false; }

  QpSolverStatus recompute() {
    std::vector<std::vector<double>> orig;
    HighsInt dim_ns = basis.getinactive().size();
    numberofreduces = 0;

    orig.assign(dim_ns, std::vector<double>(dim_ns, 0.0));
    resize(dim_ns);

    Matrix temp(dim_ns, 0);

    QpVector buffer_Qcol(runtime.instance.num_var);
    QpVector buffer_ZtQi(dim_ns);
    for (HighsInt i = 0; i < runtime.instance.num_var; i++) {
      runtime.instance.Q.mat.extractcol(i, buffer_Qcol);
      basis.Ztprod(buffer_Qcol, buffer_ZtQi);
      temp.append(buffer_ZtQi);
    }
    MatrixBase& temp_t = temp.t();
    for (HighsInt i = 0; i < dim_ns; i++) {
      basis.Ztprod(temp_t.extractcol(i, buffer_Qcol), buffer_ZtQi);
      for (HighsInt j = 0; j < buffer_ZtQi.num_nz; j++) {
        orig[i][buffer_ZtQi.index[j]] = buffer_ZtQi.value[buffer_ZtQi.index[j]];
      }
    }

    for (size_t col = 0; col < orig.size(); col++) {
      for (size_t row = 0; row <= col; row++) {
        double sum = 0;
        if (row == col) {
          for (size_t k = 0; k < row; k++)
            sum += L[k * current_k_max + row] * L[k * current_k_max + row];
          double d_value = orig[row][row] - sum;
          if (d_value <= 0) {
            return QpSolverStatus::NOTPOSITIVDEFINITE;
          }
          L[row * current_k_max + row] = sqrt(d_value);
        } else {
          for (size_t k = 0; k < row; k++)
            sum += (L[k * current_k_max + col] * L[k * current_k_max + row]);
          L[row * current_k_max + col] =
              (orig[col][row] - sum) / L[row * current_k_max + row];
        }
      }
    }
    current_k = dim_ns;
    uptodate = true;

    return QpSolverStatus::OK;
  }

  // Appends the null space column yp to the factor. The new diagonal entry
  // is the Schur complement mu - l'l, with mu = yp'Q yp and l = L^{-1}Z'Q yp.
  // When yp is nearly Q-dependent on the current null space, mu and l'l agree
  // to all digits and their difference is rounding error, possibly negative.
  // For a nonempty null space, the caller passes curvature = p'Qp, with p the
  // Q-conjugate search direction yp - Z(Z'QZ)^{-1}Z'Q yp: this is the same
  // Schur complement, evaluated as a quadratic form without cancellation.
  QpSolverStatus expand(const QpVector& yp, QpVector& gyp, QpVector& l,
                        QpVector& m, const double curvature) {
    if (!uptodate) {
      return QpSolverStatus::OK;
    }
    double mu = gyp * yp;
    l.resparsify();
    double lambda = l.dim > 0 ? curvature : mu;
    if (lambda > 0.0) {
      if (current_k_max <= current_k + 1) {
        resize(current_k_max * 2);
      }

      for (HighsInt i = 0; i < current_k; i++) {
        L[i * current_k_max + current_k] = l.value[i];
      }
      L[current_k * current_k_max + current_k] = sqrt(lambda);

      current_k++;
    } else {
      return QpSolverStatus::NOTPOSITIVDEFINITE;

      //     |LL' 0|
      // M = |0'  0| + bb' -aa'
      // a = (k * m, alpha), b = (k * m, beta)
      // b*b -a*a = mu
      // k(b-a) = 1
      // b + a = k*mu
      // Commented out unreachable code
      //      const double tolerance = 0.001;
      //
      //      double beta = max(tolerance, sqrt(m.norm2() / L[0] + fabs(mu)));
      //      double k = 1 / (beta + sqrt(beta * beta - mu));
      //      double alpha = k * mu - beta;
      //
      //      printf("k = %d, alpha = %lf, beta = %lf, k = %lf\n",
      //      (int)current_k, alpha,
      //             beta, k);
      //
      //      a.clear();
      //      a.resize(current_k + 1);
      //      for (HighsInt i = 0; i < current_k; i++) {
      //        a[i] = k * m.value[i];
      //      }
      //      a[current_k] = alpha;
      //
      //      std::vector<double> b(current_k + 1);
      //      for (HighsInt i = 0; i < current_k; i++) {
      //        b[i] = k * m.value[i];
      //      }
      //      b[current_k] = beta;
      //
      //      if (current_k_max <= current_k + 1) {
      //        resize(current_k_max * 2);
      //      }
      //
      //      // append b to the left of L
      //      for (HighsInt row = current_k; row > 0; row--) {
      //        // move row one position down
      //        for (HighsInt i = 0; i < current_k; i++) {
      //          L[row * current_k_max + i] = L[(row - 1) * current_k_max + i];
      //        }
      //      }
      //      for (HighsInt i = 0; i < current_k + 1; i++) {
      //        L[i] = b[i];
      //      }
      //
      //      // re-triangulize
      //      for (HighsInt i = 0; i < current_k + 1; i++) {
      //        eliminate(L, i, i + 1, current_k_max, current_k + 1);
      //      }
      //
      //      current_k = current_k + 1;
    }
    return QpSolverStatus::OK;
  }

  QpSolverStatus solveL(QpVector& rhs) {
    if (!uptodate) {
      // Recompute Cholesky
      QpSolverStatus status = recompute();
      if (status != QpSolverStatus::OK) return status;
    }

    if (current_k != rhs.dim) {
      assert(current_k == rhs.dim);
      return QpSolverStatus::kError;
    }

    for (HighsInt r = 0; r < rhs.dim; r++) {
      for (HighsInt j = 0; j < r; j++) {
        rhs.value[r] -= rhs.value[j] * L[j * current_k_max + r];
      }

      rhs.value[r] /= L[r * current_k_max + r];
    }
    return QpSolverStatus::OK;
  }

  // solve L' u = v
  void solveLT(QpVector& rhs) {
    for (HighsInt i = rhs.dim - 1; i >= 0; i--) {
      double sum = 0.0;
      for (HighsInt j = rhs.dim - 1; j > i; j--) {
        sum += rhs.value[j] * L[i * current_k_max + j];
      }
      rhs.value[i] = (rhs.value[i] - sum) / L[i * current_k_max + i];
    }
  }

  QpSolverStatus solve(QpVector& rhs) {
    QpSolverStatus status = QpSolverStatus::OK;
    if (!uptodate || (numberofreduces >= runtime.instance.num_var / 2 &&
                      !has_negative_eigenvalue)) {
      // Recompute Cholesky
      status = recompute();
      if (status != QpSolverStatus::OK) return status;
    }
    status = solveL(rhs);
    if (status != QpSolverStatus::OK) return status;
    solveLT(rhs);

    rhs.resparsify();
    return QpSolverStatus::OK;
  }

  void eliminate(std::vector<double>& m, HighsInt i, HighsInt j, HighsInt kmax,
                 HighsInt currentk) {
    // i = col, j = row
    if (m[j * kmax + i] == 0.0) {
      return;
    }
    double z = sqrt(m[i * kmax + i] * m[i * kmax + i] +
                    m[j * kmax + i] * m[j * kmax + i]);
    double cos_, sin_;
    if (z == 0) {
      cos_ = 1.0;
      sin_ = 0.0;
    } else {
      cos_ = m[i * kmax + i] / z;
      sin_ = -m[j * kmax + i] / z;
    }

    if (sin_ == 0.0) {
      if (cos_ > 0.0) {
        // nothing
      } else {
        for (HighsInt k = 0; k < current_k; k++) {
          // update entry i and j of column k
          double a_ik = m[i * kmax + k];
          // entry i
          m[i * kmax + k] = -a_ik;
          m[j * kmax + k] = -m[j * kmax + k];
        }
      }
    } else if (cos_ == 0.0) {
      if (sin_ > 0.0) {
        for (HighsInt k = 0; k < current_k; k++) {
          // update entry i and j of column k
          double a_ik = m[i * kmax + k];
          // entry i
          m[i * kmax + k] = -m[j * kmax + k];
          m[j * kmax + k] = a_ik;
        }
      } else {
        for (HighsInt k = 0; k < current_k; k++) {
          // update entry i and j of column k
          double a_ik = m[i * kmax + k];
          // entry i
          m[i * kmax + k] = m[j * kmax + k];
          m[j * kmax + k] = -a_ik;
        }
      }
    } else {
      // #pragma omp parallel for
      for (HighsInt k = 0; k < current_k; k++) {
        // update entry i and j of column k
        double a_ik = m[i * kmax + k];
        // entry i
        m[i * kmax + k] = cos_ * a_ik - sin_ * m[j * kmax + k];
        m[j * kmax + k] = sin_ * a_ik + cos_ * m[j * kmax + k];
      }
    }
    m[j * kmax + i] = 0.0;
  }

  void reduce(const QpVector& buffer_d, const HighsInt maxabsd, bool p_in_v) {
    if (current_k == 0) {
      return;
    }
    if (!uptodate) {
      return;
    }
    numberofreduces++;

    unsigned p = maxabsd;  // col we push to the right and remove

    // start situation: p=3, current_k = 5
    // |1 x  | |x    |       |1   | |xxxxx|
    // | 1x  | |xx   |  ===  | 1  | | xxxx|
    // |  x1 | |xxx  |       |xxxx| |  xxx|
    // |  x 1| |xxxx |       |  1 | |   xx|
    //         |xxxxx|       |   1| |    x|
    // next step: move row/col p to the bottom/right

    //> save row p
    std::vector<double> row_p(current_k, 0.0);
    for (HighsInt i = 0; i < current_k; i++) {
      row_p[i] = L[p * current_k_max + i];
    }

    //> move all rows > p up by one row
    for (HighsInt row = p; row < current_k - 1; row++) {
      for (HighsInt i = 0; i < current_k; i++) {
        L[row * current_k_max + i] = L[(row + 1) * current_k_max + i];
      }
    }

    //> load row p
    for (HighsInt i = 0; i < current_k; i++) {
      L[(current_k - 1) * current_k_max + i] = row_p[i];
    }

    //> now move col p to the right in each row
    for (HighsInt row = 0; row < current_k; row++) {
      double p_entry = L[row * current_k_max + p];
      for (HighsInt col = p; col < current_k - 1; col++) {
        L[row * current_k_max + col] = L[row * current_k_max + col + 1];
      }
      L[row * current_k_max + current_k - 1] = p_entry;
    }

    if (current_k == 1) {
      current_k--;
      return;
    }

    if (!p_in_v) {
      // situation now:
      // |1   x| |x    |       |1   | |xxxxx|
      // | 1  x| |xx   |  ===  | 1  | | xxxx|
      // |  1 x| |xxx x|       |  1 | |  xx |
      // |   1x| |xxxxx|       |   1| |   x |
      //         |xx  x|       |xxxx| |  xxx|
      // next: remove nonzero entries in last column except for diagonal element
      for (HighsInt r = (HighsInt)p - 1; r >= 0; r--) {  // to current_k-1
        eliminate(L, current_k - 1, r, current_k_max, current_k);
      }

      // situation now:
      // |1   x| |x   x|        |xxxx | |1   |
      // | 1  x| |xx  x|  ===   | xxx | | 1  |
      // |  1 x| |xxx x|        |  xx | |  1 |
      // |   1x| |xxxxx|        |   x | |   1|
      //         |    x|        |xxxxx| |xxxx|
      // next: multiply product
      // new last row: old last row (first current_k-1 elements) + r *
      // R_current_k_current_k

      for (HighsInt i = 0; i < buffer_d.num_nz; i++) {
        HighsInt idx = buffer_d.index[i];
        if (idx == maxabsd) {
          continue;
        }
        if (idx < maxabsd) {
          L[(current_k - 1) * current_k_max + idx] +=
              -buffer_d.value[idx] / buffer_d.value[maxabsd] *
              L[(current_k - 1) * current_k_max + current_k - 1];
        } else {
          L[(current_k - 1) * current_k_max + idx - 1] +=
              -buffer_d.value[idx] / buffer_d.value[maxabsd] *
              L[(current_k - 1) * current_k_max + current_k - 1];
        }
      }
      // situation now: as above, but no more product
    }
    // next: eliminate last row
    for (HighsInt i = 0; i < current_k - 1; i++) {
      eliminate(L, i, current_k - 1, current_k_max, current_k);
    }
    current_k--;
  }

  void report(const std::string& name = "") {
    printf("%s\n", name.c_str());
    for (HighsInt i = 0; i < current_k; i++) {
      for (HighsInt j = 0; j < current_k; j++) {
        printf("%lf ", L[i * current_k_max + j]);
      }
      printf("\n");
    }
  }

  double density() {
    if (current_k == 0) {
      return 0.0;
    }

    HighsInt num_nz = 0;
    for (HighsInt i = 0; i < current_k; i++) {
      for (HighsInt j = 0; j < current_k; j++) {
        if (fabs(L[i * current_k_max + j]) > 1e-7) {
          num_nz++;
        }
      }
    }
    return (double)num_nz / (current_k * (current_k + 1) / 2.0);
  }
};

#endif
