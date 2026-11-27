/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_REDUCEDCOSTS_HPP__
#define __SRC_LIB_REDUCEDCOSTS_HPP__

#include "qpsolver/basis.hpp"
#include "qpsolver/gradient.hpp"
#include "qpsolver/qpvector.hpp"
#include "qpsolver/runtime.hpp"

class ReducedCosts {
  Basis& basis;

  Gradient& gradient;

  QpVector reducedcosts;
  bool uptodate;

 public:
  ReducedCosts(Runtime& rt, Basis& bas, Gradient& grad)
      : basis(bas),
        gradient(grad),
        reducedcosts(QpVector(rt.instance.num_var)),
        uptodate(false) {}

  void recompute() {
    basis.ftran(gradient.getGradient(), reducedcosts);
    uptodate = true;
  }

  QpVector& getReducedCosts() {
    if (!uptodate) {
      recompute();
    }
    return reducedcosts;
  }

  void update() { uptodate = false; }
};

#endif
