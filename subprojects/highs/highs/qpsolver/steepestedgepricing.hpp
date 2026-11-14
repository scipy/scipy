/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_PRICING_STEEPESTEDGEPRICING_HPP__
#define __SRC_LIB_PRICING_STEEPESTEDGEPRICING_HPP__

#include "qpsolver/basis.hpp"
#include "qpsolver/pricing.hpp"
#include "qpsolver/runtime.hpp"

//

class SteepestEdgePricing : public Pricing {
 private:
  Runtime& runtime;
  Basis& basis;
  ReducedCosts& redcosts;
  std::vector<double> weights;

  HighsInt chooseconstrainttodrop(const QpVector& lambda) {
    auto active_constraint_index = basis.getactive();
    auto constraintindexinbasisfactor = basis.getindexinfactor();

    HighsInt minidx = -1;
    double maxval = 0.0;
    for (size_t i = 0; i < active_constraint_index.size(); i++) {
      HighsInt iVar = active_constraint_index[i];
      HighsInt indexinbasis = constraintindexinbasisfactor[iVar];
      if (indexinbasis == -1) {
        printf("error\n");
      }
      assert(indexinbasis != -1);
      if (runtime.instance.isEquality(iVar)) continue;
      double val = lambda.value[indexinbasis] * lambda.value[indexinbasis] /
                   weights[indexinbasis];
      if (val > maxval && fabs(lambda.value[indexinbasis]) >
                              runtime.settings.lambda_zero_threshold) {
        if (basis.getstatus(iVar) == BasisStatus::kActiveAtLower &&
            -lambda.value[indexinbasis] > 0) {
          minidx = iVar;
          maxval = val;
        } else if (basis.getstatus(iVar) == BasisStatus::kActiveAtUpper &&
                   lambda.value[indexinbasis] > 0) {
          minidx = iVar;
          maxval = val;
        } else {
          // TODO
        }
      }
    }

    return minidx;
  }

 public:
  SteepestEdgePricing(Runtime& rt, Basis& bas, ReducedCosts& rc)
      : runtime(rt),
        basis(bas),
        redcosts(rc),
        weights(std::vector<double>(rt.instance.num_var, 1.0)) {
    compute_exact_weights();
  };

  HighsInt price(const QpVector& x, const QpVector& gradient) {
    HighsInt minidx = chooseconstrainttodrop(redcosts.getReducedCosts());
    return minidx;
  }

  double compute_exact_weight(HighsInt i) {
    QpVector y_i = basis.btran(QpVector::unit(runtime.instance.num_var, i));
    return y_i.dot(y_i);
  }

  void compute_exact_weights() {
    for (int i = 0; i < runtime.instance.num_var; i++) {
      weights[i] = compute_exact_weight(i);
    }
  }

  bool check_weight(HighsInt i) {
    double weight_is = weights[i];
    double weight_comp = compute_exact_weight(i);
    if (fabs(weight_comp - weight_is) > 1e-2) {
      // printf("weights[%d] = %lf, should be %lf\n", i, weight_is,
      // weight_comp);
      return false;
    }
    return true;
  }

  bool check_weights() {
    std::vector<int> correct_weights;
    std::vector<int> incorrect_weights;
    for (int i = 0; i < runtime.instance.num_var; i++) {
      bool correct = check_weight(i);
      if (correct) {
        correct_weights.push_back(i);
      } else {
        incorrect_weights.push_back(i);
      }
    }

    printf("correct weights: ");
    for (int& i : correct_weights) {
      printf("%d ", i);
    }
    printf("\n");

    printf("incorrect weights: ");
    for (int& i : incorrect_weights) {
      printf("%d ", i);
    }
    printf("\n");

    return incorrect_weights.size() == 0;
  }

  void recompute() { compute_exact_weights(); }

  void update_weights(const QpVector& aq, const QpVector& ep, HighsInt p,
                      HighsInt q) {
    HighsInt rowindex_p = basis.getindexinfactor()[p];
    // printf("Update weights, p = %d, rowindex = %d, q = %d\n", p, rowindex_p,
    // q);

    // if (!check_weights()) {
    //   printf("weight check failed\n");
    //   exit(1);
    // }

    QpVector delta = basis.ftran(aq);

    // double old_weight_p_updated = weights[rowindex_p];
    //  exact weight coming in needs to come in before update.
    double old_weight_p_computed = ep.dot(ep);

    // if (fabs(old_weight_p_computed - old_weight_p_updated) >= 1e-2) {
    //   printf("old weight[p] discrepancy: updated = %lf, computed=%lf\n",
    //   old_weight_p_updated, old_weight_p_computed);
    // }

    double weight_p = old_weight_p_computed;

    double t_p = aq.value[rowindex_p];
    for (HighsInt i = 0; i < runtime.instance.num_var; i++) {
      if (i != rowindex_p) {
        double t_i = aq.value[i];
        weights[i] = weights[i] - 2 * (t_i / t_p) * delta.value[i] +
                     ((t_i * t_i) / (t_p * t_p)) * weight_p;
        // printf("weights[%d] = %lf\n", i, weights[i]);
      }
    }
    // QpVector new_ep = basis.btran(QpVector::unit(runtime.instance.num_var,
    // rowindex_p)); double computed_weight = new_ep.dot(new_ep);
    double new_weight_p_updated = weight_p / (t_p * t_p);

    // if (fabs(updated_weight - computed_weight) > 1e-4) {
    //   printf("updated weight %lf vs computed weight %lf. aq[p] = %lf\n",
    //   updated_weight, computed_weight, t_p); printf("old weight = %lf, aq[p]
    //   = %lf, ^2 = %lf, new weight = %lf\n", weight_p, t_p, t_p*t_p,
    //   updated_weight);
    // }
    weights[rowindex_p] = new_weight_p_updated;
  }
};

#endif
