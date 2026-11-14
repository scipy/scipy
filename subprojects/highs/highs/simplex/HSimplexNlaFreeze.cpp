/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSimplexNlaFreeze.cpp
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#include <stdio.h>

#include "simplex/HSimplexNla.h"

void SimplexIterate::clear() {
  this->valid_ = false;
  this->basis_.clear();
  this->invert_.clear();
  this->dual_edge_weight_.clear();
}

void HSimplexNla::putInvert() {
  simplex_iterate_.valid_ = true;
  simplex_iterate_.invert_ = factor_.getInvert();
}
void HSimplexNla::getInvert() { factor_.setInvert(simplex_iterate_.invert_); }
