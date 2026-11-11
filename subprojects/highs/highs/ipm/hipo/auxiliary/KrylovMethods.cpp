#include "KrylovMethods.h"

#include <cmath>
#include <iostream>

#include "VectorOperations.h"

namespace hipo {

static void applyRotation(double& x, double& y, double c, double s) {
  double t = c * x + s * y;
  y = -s * x + c * y;
  x = t;
}
static void getRotation(double x, double y, double& c, double& s) {
  if (y == 0.0) {
    c = 1.0;
    s = 0.0;
  } else if (std::abs(y) > std::abs(x)) {
    double t = x / y;
    s = 1.0 / sqrt(1.0 + t * t);
    c = t * s;
  } else {
    double t = y / x;
    c = 1.0 / sqrt(1.0 + t * t);
    s = t * c;
  }
}
static void update(std::vector<double>& x, Int k, std::vector<double>& H,
                   Int ldh, std::vector<double>& s,
                   std::vector<std::vector<double>>& V) {
  // Solve H * y = s
  std::vector<double> y = s;
  for (Int i = k; i >= 0; --i) {
    y[i] /= H[i + ldh * i];
    for (Int j = i - 1; j >= 0; --j) {
      y[j] -= H[j + ldh * i] * y[i];
    }
  }
  // x += V * y
  for (Int j = 0; j <= k; ++j) vectorAdd(x, 1.0, V[j], y[j]);
}

Int Gmres(const AbstractMatrix* M, const AbstractMatrix* P,
          const std::vector<double>& b, std::vector<double>& x, double tol,
          Int maxit) {
  // Attempt to solve M * x = b using GMRES, with the preconditioner P.
  // Return the number of iterations taken.
  //
  // Based on Netlib implementation and Kelley "Iterative Methods for Linear and
  // Nonlinear Equations":
  // https://www.netlib.org/templates/cpp/gmres.h
  // https://www.netlib.org/templates/templates.pdf
  // https://en.wikipedia.org/wiki/Generalized_minimal_residual_method#Example_code

  // sizes
  const Int n = x.size();
  const Int m = maxit;

  // vectors
  std::vector<std::vector<double>> V;

  // Hessenberg matrix
  std::vector<double> H((m + 1) * m);
  const Int ldh = m + 1;

  // Givens rotations
  std::vector<double> sn(m + 1);
  std::vector<double> cs(m + 1);

  // workspace
  std::vector<double> w(n);

  // preconditioned residual P^-1 * (b - M * x)
  std::vector<double> r = x;
  M->apply(r);
  vectorAdd(r, -1.0, b, 1.0);
  if (P) P->apply(r);
  double beta = norm2(r);

  // preconditioned rhs P^-1 * b
  w = b;
  if (P) P->apply(w);
  double normb = norm2(w);

  if (normb == 0.0) normb = 1.0;
  if (beta <= tol) return 0;

  V.push_back(r);
  vectorScale(V[0], 1.0 / beta);

  // s = beta * e1
  std::vector<double> s(m + 1);
  s[0] = beta;

  // main loop
  Int i;
  for (i = 0; i < maxit; ++i) {
    // Arnoldi
    w = V[i];
    M->apply(w);
    if (P) P->apply(w);
    V.push_back(w);
    for (Int k = 0; k <= i; ++k) {
      H[k + ldh * i] = dotProd(V.back(), V[k]);
      vectorAdd(V.back(), 1.0, V[k], -H[k + ldh * i]);
    }
    H[i + 1 + ldh * i] = norm2(V.back());

    // re-orthogonalize
    // if (norm2(V.back()) / norm2(w) < 1e-6) {
    // printf("Re-orthogonalize============\n");
    for (Int k = 0; k <= i; ++k) {
      double temp = dotProd(V.back(), V[k]);
      H[k + ldh * i] += temp;
      vectorAdd(V.back(), 1.0, V[k], -temp);
    }
    H[i + 1 + ldh * i] = norm2(V.back());
    // }

    // normalize
    vectorScale(V.back(), 1.0 / H[i + 1 + ldh * i]);

    // apply Givens rotations
    for (Int k = 0; k < i; ++k)
      applyRotation(H[k + ldh * i], H[k + 1 + ldh * i], cs[k], sn[k]);

    // get latest rotation and apply it also to s
    getRotation(H[i + ldh * i], H[i + 1 + ldh * i], cs[i], sn[i]);
    applyRotation(H[i + ldh * i], H[i + 1 + ldh * i], cs[i], sn[i]);
    applyRotation(s[i], s[i + 1], cs[i], sn[i]);

    // check termination
    if (std::abs(s[i + 1]) < tol) break;
  }

  // compute solution
  update(x, i, H, ldh, s, V);

  return i + 1;
}

Int Cg(const AbstractMatrix* M, const AbstractMatrix* P,
       const std::vector<double>& b, std::vector<double>& x, double tol,
       Int maxit) {
  // Attempt to solve M * x = b using CG, with the preconditioner P.
  // Return the number of iterations taken.

  Int n = b.size();

  std::vector<double> w(n);

  std::vector<double> r = x;
  M->apply(r);
  vectorAdd(r, -1.0, b, 1.0);

  std::vector<double> z = r;
  if (P) P->apply(z);

  std::vector<double> p(z);

  double rho_old = dotProd(r, z);
  double rho_new;
  double norm_b = norm2(b);

  Int iter = 0;
  while (iter < maxit) {
    w = p;
    M->apply(w);
    double alpha = rho_old / dotProd(p, w);
    vectorAdd(x, 1.0, p, alpha);
    vectorAdd(r, 1.0, w, -alpha);

    z = r;
    if (P) P->apply(z);
    rho_new = dotProd(r, z);

    if (norm2(r) < tol * norm_b) break;

    vectorAdd(p, rho_new / rho_old, z, 1.0);
    rho_old = rho_new;
    ++iter;

    if (isNanVector(x)) {
      break;
    }
  }

  return iter;
}

}  // namespace hipo