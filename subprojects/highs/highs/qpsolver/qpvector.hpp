/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_VECTOR_HPP__
#define __SRC_LIB_VECTOR_HPP__

#include <util/HighsType.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

struct QpVector {
  HighsInt num_nz;
  HighsInt dim;
  std::vector<HighsInt> index;
  std::vector<double> value;

  QpVector(HighsInt d) : dim(d) {
    index.resize(dim);
    value.resize(dim, 0.0);
    num_nz = 0;
  }

  QpVector(const QpVector& vec)
      : num_nz(vec.num_nz), dim(vec.dim), index(vec.index), value(vec.value) {}

  void reset() {
    for (HighsInt i = 0; i < num_nz; i++) {
      value[index[i]] = 0;
      index[i] = 0;
    }
    num_nz = 0;
  }

  QpVector& repopulate(const QpVector& other) {
    reset();
    for (HighsInt i = 0; i < other.num_nz; i++) {
      index[i] = other.index[i];
      value[index[i]] = other.value[index[i]];
    }
    num_nz = other.num_nz;
    return *this;
  }

  QpVector& operator=(const QpVector& other) {
    num_nz = other.num_nz;
    dim = other.dim;
    index = other.index;
    value = other.value;
    return *this;
  }

  static QpVector& unit(HighsInt dim, HighsInt u, QpVector& target) {
    target.reset();
    target.index[0] = u;
    target.value[u] = 1.0;
    target.num_nz = 1;
    return target;
  }

  static QpVector unit(HighsInt dim, HighsInt u) {
    QpVector vec(dim);
    vec.index[0] = u;
    vec.value[u] = 1.0;
    vec.num_nz = 1;
    return vec;
  }

  void report(const std::string& name = "") const {
    if (name != "") {
      printf("%s: ", name.c_str());
    }
    for (HighsInt i = 0; i < num_nz; i++) {
      printf("[%" HIGHSINT_FORMAT "] %lf ", index[i], value[index[i]]);
    }
    printf("\n");
  }

  double norm2() {
    double val = 0.0;

    for (HighsInt i = 0; i < num_nz; i++) {
      val += value[index[i]] * value[index[i]];
    }

    return val;
  }

  void sanitize(double threshold = 1e-14) {
    HighsInt new_idx = 0;

    for (HighsInt i = 0; i < num_nz; i++) {
      if (fabs(value[index[i]]) > threshold) {
        index[new_idx++] = index[i];
      } else {
        value[index[i]] = 0.0;
        index[i] = 0;
      }
    }
    num_nz = new_idx;
  }

  void resparsify() {
    num_nz = 0;
    for (HighsInt i = 0; i < dim; i++) {
      if (value[i] != 0.0) {
        index[num_nz++] = i;
      }
    }
  }

  QpVector& scale(double a) {
    for (HighsInt i = 0; i < num_nz; i++) {
      value[index[i]] *= a;
    }
    return *this;
  }

  QpVector& saxpy(double a, double b, const QpVector& x) {
    scale(a);
    saxpy(b, x);
    return *this;
  }

  QpVector& saxpy(double a, const QpVector& x) {
    sanitize(0.0);
    for (HighsInt i = 0; i < x.num_nz; i++) {
      if (value[x.index[i]] == 0.0) {
        index[num_nz++] = x.index[i];
      }
      value[x.index[i]] += a * x.value[x.index[i]];
    }
    resparsify();
    // sanitize(0.0);
    return *this;
  }

  // void saxpy(double a, HighsInt* idx, double* val, HighsInt nnz) {
  //    for (HighsInt i=0; i<nnz; i++) {
  //       value[idx[i]] += a * val[i];
  //    }
  //    resparsify();
  // }

  QpVector operator+(const QpVector& other) const {
    QpVector result(dim);

    for (HighsInt i = 0; i < dim; i++) {
      result.value[i] = value[i] + other.value[i];
      if (result.value[i] != 0.0) {
        result.index[result.num_nz++] = i;
      }
    }

    return result;
  }

  QpVector operator-(const QpVector& other) const {
    QpVector result(dim);

    for (HighsInt i = 0; i < dim; i++) {
      result.value[i] = value[i] - other.value[i];
      if (result.value[i] != 0.0) {
        result.index[result.num_nz++] = i;
      }
    }

    return result;
  }

  QpVector operator-() const {
    QpVector result(dim);

    for (HighsInt i = 0; i < num_nz; i++) {
      result.index[i] = index[i];
      result.value[index[i]] = -value[index[i]];
    }
    result.num_nz = num_nz;

    return result;
  }

  QpVector operator*(const double d) const {
    QpVector result(dim);

    for (HighsInt i = 0; i < num_nz; i++) {
      result.index[i] = index[i];
      result.value[index[i]] = d * value[index[i]];
    }
    result.num_nz = num_nz;

    return result;
  }

  double dot(const QpVector& other) const {
    double dot = 0.0;
    for (HighsInt i = 0; i < num_nz; i++) {
      dot += value[index[i]] * other.value[index[i]];
    }

    return dot;
  }

  double operator*(const QpVector& other) const { return dot(other); }

  double dot(const HighsInt* idx, const double* val, HighsInt nnz) const {
    double dot = 0.0;
    for (HighsInt i = 0; i < nnz; i++) {
      dot += value[idx[i]] * val[i];
    }

    return dot;
  }

  QpVector& operator+=(const QpVector& other) {
    // sanitize();
    for (HighsInt i = 0; i < other.num_nz; i++) {
      // if (value[other.index[i]] == 0.0) {
      //    index[num_nz++] = other.index[i];
      // }
      value[other.index[i]] += other.value[other.index[i]];
    }
    resparsify();
    return *this;
  }

  QpVector& operator*=(const double d) {
    for (HighsInt i = 0; i < num_nz; i++) {
      value[index[i]] *= d;
    }

    return *this;
  }
};

#endif
