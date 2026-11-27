/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_PRICING_HPP__
#define __SRC_LIB_PRICING_HPP__

#include "qpvector.hpp"

class Pricing {
 public:
  virtual HighsInt price(const QpVector& x, const QpVector& gradient) = 0;
  virtual void update_weights(const QpVector& aq, const QpVector& ep,
                              HighsInt p, HighsInt q) = 0;
  virtual void recompute() = 0;
  virtual ~Pricing() {}
};

#endif
