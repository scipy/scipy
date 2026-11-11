/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_PRICING_DANTZIGPRICING_HPP__
#define __SRC_LIB_PRICING_DANTZIGPRICING_HPP__

#include "basis.hpp"
#include "pricing.hpp"
#include "reducedcosts.hpp"
#include "runtime.hpp"

// 51561, 78965838.346823, 559, 213.280772, 0.000812, 801

class DantzigPricing : public Pricing {
 private:
  Runtime& runtime;
  Basis& basis;
  ReducedCosts& redcosts;

  HighsInt chooseconstrainttodrop(const QpVector& lambda) {
    auto active_constraint_index = basis.getactive();
    auto constraintindexinbasisfactor = basis.getindexinfactor();

    HighsInt minidx = -1;
    double maxabslambda = 0.0;
    for (size_t i = 0; i < active_constraint_index.size(); i++) {
      HighsInt iVar = active_constraint_index[i];
      HighsInt indexinbasis = constraintindexinbasisfactor[iVar];
      if (indexinbasis == -1) {
        printf("error\n");
      }
      assert(indexinbasis != -1);
      if (runtime.instance.isEquality(iVar)) continue;
      if (basis.getstatus(iVar) == BasisStatus::kActiveAtLower &&
          -lambda.value[indexinbasis] > maxabslambda) {
        minidx = iVar;
        maxabslambda = -lambda.value[indexinbasis];
      } else if (basis.getstatus(iVar) == BasisStatus::kActiveAtUpper &&
                 lambda.value[indexinbasis] > maxabslambda) {
        minidx = iVar;
        maxabslambda = lambda.value[indexinbasis];
      } else {
        // TODO
      }
    }

    if (maxabslambda <= runtime.settings.lambda_zero_threshold) {
      // printf("maxabslambda %lf\n", log(maxabslambda));
      return -1;
    }
    return minidx;
  }

 public:
  DantzigPricing(Runtime& rt, Basis& bas, ReducedCosts& rc)
      // clang-format off
      : runtime(rt), basis(bas), redcosts(rc) {};
  // clang-format on
  HighsInt price(const QpVector& x, const QpVector& gradient) {
    HighsInt minidx = chooseconstrainttodrop(redcosts.getReducedCosts());
    return minidx;
  }

  void recompute() {
    // do nothing
  }

  void update_weights(const QpVector& aq, const QpVector& ep, HighsInt p,
                      HighsInt q) {
    // does nothing
  }
};

#endif
