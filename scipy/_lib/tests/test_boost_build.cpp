// Testing availability of boost headers

#include <iostream>
#include <boost/math/distributions.hpp>

void test() {
  boost::math::binomial_distribution<double> d(10, 0.5);
}
