/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_DOMAIN_CHANGE_H_
#define HIGHS_DOMAIN_CHANGE_H_

#include "util/HighsType.h"

enum class HighsBoundType { kLower, kUpper };

struct HighsDomainChange {
  double boundval;
  HighsInt column;
  HighsBoundType boundtype;

  bool operator<(const HighsDomainChange& other) const {
    if (column < other.column) return true;
    if (other.column < column) return false;
    if ((HighsInt)boundtype < (HighsInt)other.boundtype) return true;
    if ((HighsInt)other.boundtype < (HighsInt)boundtype) return false;
    if (boundval < other.boundval) return true;
    return false;
  }

  bool operator==(const HighsDomainChange& other) const {
    return boundtype == other.boundtype && column == other.column &&
           boundval == other.boundval;
  }

  bool operator!=(const HighsDomainChange& other) const {
    return boundtype != other.boundtype || column != other.column ||
           boundval != other.boundval;
  }
};

struct HighsSubstitution {
  HighsInt substcol;
  HighsInt staycol;
  double scale;
  double offset;
};

#endif
